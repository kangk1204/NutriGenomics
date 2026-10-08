"""Allowlisted algorithm subprocesses with original input and model provenance."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
from .api import command_for,safe_path
from . import __version__
from .execution_sources import PACKAGES,source_layout,bound_command,current_module_origins,verify_module_origins

def fingerprint(path):
    with path.open('rb') as handle:digest=hashlib.file_digest(handle,'sha256').hexdigest()
    return {'path':str(path),'sha256':digest,'bytes':path.stat().st_size}


def code_snapshot(root, algorithm):
    """Identify synchronized source bytes even when the server has no Git checkout."""
    snapshots=[]
    for source in source_layout(root,algorithm):
        name=source['package'];repo=source['directory']
        files=[]
        for path in sorted((repo/'src').rglob('*.py')):
            if '__pycache__' not in path.parts:
                if path.is_symlink() or not path.resolve().is_relative_to(repo):raise ValueError('Unauditable source file: '+str(path))
                files.append({**fingerprint(path),'relative_path':path.relative_to(repo).as_posix()})
        if not files:
            raise ValueError('Auditable source checkout missing: '+name)
        content=[(f['relative_path'],f['sha256']) for f in files]
        snapshots.append({'package':name,'source_tree_sha256':hashlib.sha256(json.dumps(content).encode()).hexdigest(),
                          'files':files,'git_commit':None,'version_policy':'Synchronized source byte identity; no fabricated Git HEAD'})
    return snapshots

def execute(root,job):
    root=Path(root).resolve();output=Path(job['output'])
    before=code_snapshot(root,job['algorithm'])
    from .protocols import input_identity,protocol_id
    fixed_protocol=protocol_id(job['algorithm'],job['action'])
    if job.get('protocol_id') not in {None,fixed_protocol}:raise ValueError('Unknown execution protocol')
    input_before,input_hash=input_identity(root,job['algorithm'],job['action'],job['parameters'])
    (output/'input_manifest.json').write_text(json.dumps(input_before,ensure_ascii=False,indent=2),encoding='utf-8')
    command=command_for(root,sys.executable,job)
    receipt=output/'execution_modules.json'
    subprocess.run(bound_command(root,job['algorithm'],command,receipt),check=True)
    module_origins=json.loads(receipt.read_text(encoding='utf-8'))['origins']
    module_origins.extend(current_module_origins(root,job['algorithm']))
    module_origins.append({'module':'nutriomics_final.dispatch',**fingerprint(Path(__file__).resolve())})
    verify_module_origins(before,module_origins)
    algorithm=job['algorithm'];params=job['parameters'];sources=[];models=[]
    if algorithm=='methylation':
        study=output
        if job['action']=='train':
            for action in ('evaluate','export'):
                followup=output/('execution_modules_'+action+'.json')
                subprocess.run(bound_command(root,job['algorithm'],[sys.executable,'-m','nutriomics_methylation.cli',action,'--root',str(study)],followup),check=True)
                origins=json.loads(followup.read_text(encoding='utf-8'))['origins'];verify_module_origins(before,origins);module_origins.extend(origins)
        if job['action']=='train':
            validation=study/'docs'/'validation'
            if validation.exists():shutil.copytree(validation,output/'report',dirs_exist_ok=True)
        else:
            # Evaluation writes new summaries of copied frozen predictions;
            # returning an older docs/validation directory would be stale.
            report=output/'report';report.mkdir(exist_ok=True)
            for path in (study/'results').rglob('*'):
                if path.is_file() and (path.name in {'metrics_all.json','performance.tsv','metrics.json','external_metrics.json','preHT_summary.json'} or path.name.endswith('roc_calibration.svg')):
                    target=report/path.relative_to(study/'results');target.parent.mkdir(parents=True,exist_ok=True)
                    shutil.copy2(path,target)
            (report/'REPORT_SCOPE.md').write_text('Evaluation of copied frozen predictions only. No model retraining or new cohort validation. Bootstrap intervals condition on realized predictions. Source study files remain unchanged.\n',encoding='utf-8')
            sources.extend({**r,'kind':'frozen original prediction/configuration input'} for r in json.loads((output/'evaluation_inputs.json').read_text(encoding='utf-8')))
        manifest=study/'data'/'source_manifest.json'
        if manifest.exists():sources.append({**fingerprint(manifest),'kind':'original source manifest'})
        models=[fingerprint(p) for p in (study/'results').glob('*/final_training_audit.json')]
    elif algorithm=='compound-target':
        sources.append({**fingerprint(safe_path(root,params['input'])),'kind':'user supplied molecular pairs'})
        model_dir=safe_path(root,params['model_dir'])
        training=json.loads((model_dir/'training.json').read_text(encoding='utf-8'))
        weight='baseline.joblib' if training['model']=='baseline' else 'cnn.pt'
        models=[fingerprint(safe_path(root,model_dir/name)) for name in ('training.json',weight,'training_entities.csv.gz')]
    elif algorithm=='evidence':
        import sqlite3
        database=safe_path(root,params['database'])
        sources.append({**fingerprint(database),'kind':'versioned evidence database'})
        with sqlite3.connect(database) as db:
            original=[json.loads(r[0]) for r in db.execute('SELECT metadata_json FROM source')]
        (output/'source_manifest.json').write_text(json.dumps(original,ensure_ascii=False,indent=2),encoding='utf-8')
        sources.extend({'id':r.get('id'),'source':r.get('source'),'sha256':r.get('sha256'),'version':r.get('version'),'kind':'original source'} for r in original)
    else:
        original=safe_path(root,params['input_dir'])
        sources=[{**fingerprint(p),'kind':'original metabolomics snapshot'} for p in sorted(original.iterdir()) if p.is_file()]
    model_version=hashlib.sha256(json.dumps(models,sort_keys=True).encode()).hexdigest() if models else 'no learned model in this operation'
    code=[fingerprint(p) for p in sorted(Path(__file__).parent.glob('*.py'))]
    after=code_snapshot(root,algorithm)
    if before!=after:
        raise RuntimeError('Source files changed during this job; result is not a frozen-code execution')
    input_after,after_hash=input_identity(root,algorithm,job['action'],params)
    if input_hash!=after_hash:raise RuntimeError('Protocol inputs changed during execution')
    provenance={'model_version':model_version,'model_manifests':models,'sources':sources,'code_files':code,'code_snapshots':before,
                'module_origins':module_origins,
                'code_version':hashlib.sha256(json.dumps([(s['package'],s['source_tree_sha256']) for s in before]).encode()).hexdigest(),
                'service_version':__version__,'research_only':True,'verified_treatment_effect':False,
                'protocol_id':fixed_protocol,'input_hash':input_hash,'output_manifest':'output_manifest.json'}
    (output/'provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2),encoding='utf-8')
    artifacts=[{**fingerprint(p),'path':p.relative_to(output).as_posix()} for p in sorted(output.rglob('*'))
               if p.is_file() and not p.is_symlink() and p.name not in {'output_manifest.json','run.log'}]
    (output/'output_manifest.json').write_text(json.dumps({'protocol_id':fixed_protocol,'input_hash':input_hash,
                'model_version':model_version,'files':artifacts},ensure_ascii=False,indent=2),encoding='utf-8')
    return provenance

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True);parser.add_argument('--job',type=Path,required=True)
    args=parser.parse_args();execute(args.root,json.loads(args.job.read_text(encoding='utf-8')))
