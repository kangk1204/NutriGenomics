from __future__ import annotations
import os
import sqlite3
import logging
import subprocess
import sys
import threading
from contextlib import asynccontextmanager
from pathlib import Path
from typing import Literal
from fastapi import FastAPI, HTTPException, Query, UploadFile, File
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, ConfigDict, Field
from . import __version__
from .jobs import JobStore, Runner

class Submission(BaseModel):
    model_config = ConfigDict(extra='forbid')
    algorithm: Literal['evidence', 'methylation', 'intervention', 'compound-target']
    action: Literal['validate', 'query', 'train', 'evaluate', 'predict', 'analyze']
    parameters: dict = Field(default_factory=dict)
    protocol_id: str | None = Field(default=None,pattern=r'^[a-z0-9-]{1,96}$')

class EvidencePair(BaseModel):
    model_config=ConfigDict(extra='forbid')
    source_id:str=Field(min_length=1,max_length=200)
    text:str=Field(min_length=1,max_length=50000)
    chemical_start:int=Field(ge=0)
    chemical_end:int=Field(gt=0)
    gene_start:int=Field(ge=0)
    gene_end:int=Field(gt=0)

SUPPORTED = {'evidence': {'validate', 'query'}, 'methylation': {'train', 'evaluate'},
             'intervention': {'analyze'}, 'compound-target': {'predict'}}

def safe_path(root: Path, value, must_exist=True):
    path = (root / str(value)).resolve()
    if not path.is_relative_to(root.resolve()) or (must_exist and not path.exists()):
        raise ValueError('Input must exist inside configured research root')
    return path

def command_for(root: Path, python: str, job, prepare: bool = True):
    params, algorithm, action = job['parameters'], job['algorithm'], job['action']
    output = Path(job['output'])
    if action not in SUPPORTED.get(algorithm, set()):
        raise ValueError('Unsupported algorithm/action combination')
    if algorithm == 'evidence':
        if set(params) - {'database', 'disease', 'tier', 'limit'}:
            raise ValueError('Unexpected evidence parameters')
        db = safe_path(root, params['database'])
        command = [python, '-m', 'nutriomics_evidence', action, '--database', str(db)]
        if action == 'query':
            # Export a bounded, reproducible evidence query to the job result.
            command = [python, '-m', 'nutriomics_evidence', 'export', '--database', str(db),
                       '--output', str(output / 'evidence.json'), '--limit', str(max(1,min(int(params.get('limit',100)),10000)))]
            for key in ('disease', 'tier'):
                if key in params:
                    command += ['--'+key, str(params[key])]
        return command
    if algorithm == 'methylation':
        if set(params) - {'root', 'cores', 'spaces', 'methods', 'seed'}:
            raise ValueError('Unexpected methylation parameters')
        study = safe_path(root, params['root'])
        if action == 'evaluate' and set(params)!={'root'}:
            raise ValueError('Evaluation scope is fixed by the existing training contract; provide only root')
        spaces = params.get('spaces', ['common', 'full'])
        methods = params.get('methods', ['nested5', 'loocv'])
        if not spaces or set(spaces) - {'common','full'} or not methods or set(methods) - {'nested5','loocv'}:
            raise ValueError('Unsupported feature space or evaluation method')
        # Prepared study files remain under the audited study root; preserve run isolation.
        if action == 'train' and prepare:
            data = study / 'data'
            if not data.is_dir():
                raise ValueError('Prepared data directory missing')
            import shutil
            shutil.copytree(data / 'prepared', output / 'data' / 'prepared')
            if (data / 'source_manifest.json').exists():
                shutil.copy2(data / 'source_manifest.json', output / 'data' / 'source_manifest.json')
            study = output
        if action == 'evaluate':
            required=[study/'results/training_config.json',study/'data/prepared/external_metadata.tsv']
            if any(not p.is_file() for p in required):
                raise ValueError('Existing frozen prediction/configuration bundle is required')
            if prepare:
                import hashlib
                import json
                import shutil
                inputs=required+[p for p in (study/'results').rglob('*predictions.tsv')]
                inputs += list((study/'results').glob('*/final_training_audit.json'))
                if (study/'data/source_manifest.json').is_file():inputs.append(study/'data/source_manifest.json')
                copied=[]
                for path in inputs:
                    original=safe_path(root,path)
                    destination=output/path.relative_to(study)
                    destination.parent.mkdir(parents=True,exist_ok=True)
                    shutil.copy2(original,destination)
                    original_hash=hashlib.sha256(original.read_bytes()).hexdigest()
                    if hashlib.sha256(destination.read_bytes()).hexdigest()!=original_hash:
                        raise ValueError('Copied prediction input checksum changed')
                    copied.append({'source':str(original),'copy':str(destination.relative_to(output)),
                                   'sha256':original_hash})
                (output/'evaluation_inputs.json').write_text(json.dumps(copied,indent=2),encoding='utf-8')
                study=output
        return [python, '-m', 'nutriomics_methylation.cli', action, '--root', str(study),
                '--cores', str(max(1,min(int(params.get('cores',4)),12))), '--spaces', *spaces,
                '--methods', *methods, '--seed', str(int(params.get('seed',20261001)))]
    if algorithm == 'compound-target':
        if set(params) != {'input', 'model_dir'}:
            raise ValueError('Prediction requires input and model_dir')
        return [python,'-m','nutriomics_dti.cli','predict','--input',str(safe_path(root,params['input'])),
                '--model-dir',str(safe_path(root,params['model_dir'])), '--out',str(output/'predictions.csv')]
    if set(params) != {'input_dir'}:
        raise ValueError('Intervention analysis requires audited ST001257 input_dir')
    return [python,'-m','nutriomics_atlas.cli','analyze-metabolomics','--input-dir',str(safe_path(root,params['input_dir'])),
            '--output-dir',str(output)]

def create_app(root: Path, database: Path | None=None, output: Path | None=None, command_builder=None,
               literature_database: Path | None=None, research_manifest: Path | None=None,
               food_database: Path | None=None, image_model: Path | None=None, web_directory: Path | None=None,
               evidence_model:Path|None=None,evidence_receipt:Path|None=None):
    root = Path(root).resolve()
    store = JobStore(database or root/'integration_jobs.sqlite', output or root/'job_results')
    def dispatch_command(job):
        import json
        path=Path(job['output'])/'job.json'
        path.write_text(json.dumps(job,ensure_ascii=False),encoding='utf-8')
        from .execution_sources import bound_command
        return bound_command(root,job['algorithm'],[sys.executable,'-m','nutriomics_final.dispatch','--root',str(root),'--job',str(path)])
    builder = command_builder or dispatch_command
    runner = Runner(store, builder)
    stop = threading.Event()
    queue_status = {'last_error': None}
    def worker():
        while not stop.is_set():
            try:
                for job in store.pending():
                    if job['state']=='queued' and not stop.is_set():
                        runner.run(job['id'])
                queue_status['last_error'] = None
            except sqlite3.OperationalError as exc:
                queue_status['last_error'] = str(exc)
                logging.exception('Metadata queue operation failed; queued work retained')
            stop.wait(0.2)
    @asynccontextmanager
    async def lifespan(app):
        store.recover()
        thread = threading.Thread(target=worker, daemon=True)
        app.state.worker_thread = thread
        thread.start()
        yield
        stop.set()
        with runner.lock:
            owned_processes = list(runner.processes)
        for job_id in owned_processes:
            runner.cancel(job_id)
        thread.join(timeout=3)
    app = FastAPI(title='NutriOmics research jobs',version=__version__,lifespan=lifespan)
    app.state.store, app.state.runner = store, runner
    from .food_registry import FoodRegistry
    food_registry=FoodRegistry(root,food_database or root/'expansion_20261003/derived/expansion.duckdb')
    from .image_candidates import FoodImageMatcher
    image_matcher=FoodImageMatcher(root,image_model or root/'expansion_20261003/models/clip',food_registry)
    @app.post('/food-images/candidates')
    async def food_image(image:UploadFile=File(...)):
        raw=await image.read(8*1024*1024+1)
        from starlette.concurrency import run_in_threadpool
        try:return await run_in_threadpool(image_matcher.predict,raw)
        except FileNotFoundError as exc:raise HTTPException(503,str(exc)) from exc
        except ValueError as exc:raise HTTPException(422,str(exc)) from exc
    @app.get('/foods')
    def foods(q:str=Query('',max_length=200),limit:int=Query(25,ge=1,le=100)):
        try:return food_registry.search(q,limit)
        except FileNotFoundError as exc:raise HTTPException(503,str(exc)) from exc
    @app.get('/foods/{food_id}/nutrients')
    def food_nutrients(food_id:str,grams:float=Query(...,gt=0,le=10000)):
        try:return food_registry.composition(food_id,grams)
        except KeyError as exc:raise HTTPException(404,str(exc)) from exc
        except FileNotFoundError as exc:raise HTTPException(503,str(exc)) from exc
        except ValueError as exc:raise HTTPException(422,str(exc)) from exc
    literature = literature_database or root/'nutriomics-evidence-engine/artifacts/literature_snapshot_v3.sqlite'
    evidence_predictor=None
    @app.post('/literature/predict-pair')
    async def predict_pair(value:EvidencePair):
        nonlocal evidence_predictor
        from nutriomics_evidence.neural_inference import EvidencePredictor
        from starlette.concurrency import run_in_threadpool
        try:
            if evidence_predictor is None:
                model_path=evidence_model or root/'expansion_20261003/models/evidence'
                receipt_path=evidence_receipt or root/'expansion_20261003/models/evidence_receipt.json'
                if not model_path.resolve().is_relative_to(root) or not receipt_path.resolve().is_relative_to(root):
                    raise HTTPException(503,'Evidence model must remain inside research root')
                evidence_predictor=EvidencePredictor(model_path,receipt_path)
            return await run_in_threadpool(evidence_predictor.predict,**value.model_dump())
        except FileNotFoundError as exc:raise HTTPException(503,str(exc)) from exc
        except ValueError as exc:raise HTTPException(422,str(exc)) from exc
    def literature_path():
        try:
            return safe_path(root,literature)
        except ValueError as exc:
            raise HTTPException(503,'Frozen literature snapshot is not configured') from exc
    @app.get('/literature/status')
    def literature_status():
        from nutriomics_evidence.literature_bridge import status
        return status(literature_path())
    @app.get('/literature/papers')
    def literature_papers(q:str='',limit:int=Query(25,ge=1,le=500),offset:int=Query(0,ge=0)):
        from nutriomics_evidence.literature_bridge import search
        return search(literature_path(),query=q,limit=limit,offset=offset)
    @app.get('/literature/papers/{pmid}')
    def literature_paper(pmid:str):
        from nutriomics_evidence.literature_bridge import search
        result=search(literature_path(),pmid=pmid)
        if not result['items']:
            raise HTTPException(404,'Unknown PMID')
        return result['items'][0]
    @app.get('/literature/evidence')
    def literature_evidence(pmid:str|None=None,reviewed:bool|None=None,q:str='',
                            limit:int=Query(25,ge=1,le=500),offset:int=Query(0,ge=0)):
        from nutriomics_evidence.literature_bridge import search
        return search(literature_path(),pmid=pmid,query=q,reviewed=reviewed,
                      limit=limit,offset=offset,assertions=True)
    @app.get('/health')
    def health():
        alive = getattr(app.state,'worker_thread',None)
        return {'status':'degraded' if queue_status['last_error'] or (alive and not alive.is_alive()) else 'ok',
                'version':__version__,'research_only':True,'single_worker':True,
                'queue_error':queue_status['last_error'],'worker_alive':bool(alive and alive.is_alive())}
    @app.post('/jobs',status_code=202)
    def submit(value:Submission):
        from .protocols import protocol_id
        fixed_protocol=protocol_id(value.algorithm,value.action)
        if value.protocol_id is not None and value.protocol_id!=fixed_protocol:
            raise HTTPException(422,'Protocol does not match this algorithm and action')
        if value.action not in SUPPORTED[value.algorithm]:
            raise HTTPException(422,'Unsupported algorithm/action combination')
        # Validate before queueing, without mutating study files.
        try:
            if command_builder is None:
                dry = {'algorithm':value.algorithm,'action':value.action,'parameters':value.parameters,
                       'output':str(store.output_root/'validation')}
                command_for(root,sys.executable,dry,prepare=False)
        except (ValueError,KeyError,TypeError) as exc:
            raise HTTPException(422,str(exc)) from exc
        return store.create(value.algorithm,value.action,value.parameters,__version__,protocol_id=fixed_protocol)
    @app.get('/jobs')
    def listing():
        return store.list()
    @app.get('/jobs/{job_id}')
    def status(job_id:str):
        try:return store.get(job_id)
        except KeyError as exc:raise HTTPException(404,'Unknown job') from exc
    @app.post('/jobs/{job_id}/cancel')
    def cancel(job_id:str):
        try:return runner.cancel(job_id)
        except KeyError as exc:raise HTTPException(404,'Unknown job') from exc
    @app.post('/jobs/{job_id}/retry',status_code=202)
    def retry(job_id:str):
        try:
            job=store.get(job_id)
            if job['state'] not in {'failed','cancelled','interrupted'}:
                raise HTTPException(409,'Only failed, cancelled, or interrupted jobs can be retried')
            return store.create(job['algorithm'],job['action'],job['parameters'],__version__,job_id,job.get('protocol_id'))
        except KeyError as exc:raise HTTPException(404,'Unknown job') from exc
    @app.get('/jobs/{job_id}/result')
    def result(job_id:str):
        try:return store.result(job_id)
        except KeyError as exc:raise HTTPException(404,'Unknown job') from exc
        except ValueError as exc:raise HTTPException(409,str(exc)) from exc
    @app.get('/jobs/{job_id}/files/{artifact_path:path}')
    def download(job_id:str, artifact_path:str):
        try:
            value = store.result(job_id)
            allowed = {item['path'] for item in value['files']}
            if artifact_path not in allowed:
                raise HTTPException(404,'Unknown result artifact')
            target = safe_path(Path(value['job']['output']), artifact_path)
            return FileResponse(target, filename=target.name)
        except KeyError as exc:raise HTTPException(404,'Unknown job') from exc
        except ValueError as exc:raise HTTPException(409,str(exc)) from exc
    from .research_registry import ResearchRegistry
    registry = ResearchRegistry(root, research_manifest or root/'nutriomics-final-report/docs/completion_20261003/bundle/manifest.json')
    @app.get('/research/status')
    def research_status():
        try:return registry.status()
        except FileNotFoundError as exc:raise HTTPException(503,str(exc)) from exc
        except ValueError as exc:raise HTTPException(409,str(exc)) from exc
    @app.get('/research/evidence/{key}')
    def research_evidence(key: str):
        try:return registry.evidence(key)
        except KeyError as exc:raise HTTPException(404,'Unknown evidence identifier') from exc
        except FileNotFoundError as exc:raise HTTPException(503,str(exc)) from exc
        except ValueError as exc:raise HTTPException(409,str(exc)) from exc
    web=web_directory or root/'nutriomics-final-report/clients/web/dist'
    if web.is_dir():
        if not web.resolve().is_relative_to(root):raise ValueError('Web build must remain inside research root')
        app.mount('/app',StaticFiles(directory=web,html=True),name='research_web')
    return app

def create_default_app():
    return create_app(Path(os.environ.get('NUTRIOMICS_RESEARCH_ROOT',Path.cwd())))
