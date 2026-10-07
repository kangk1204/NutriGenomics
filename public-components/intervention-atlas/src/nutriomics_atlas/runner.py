"""Execute adapted R workflows against existing raw input, with fresh outputs."""
from __future__ import annotations
import json
import os
import shutil
import subprocess
from pathlib import Path
from .provenance import file_record,now,write_json
from .geo import prepare_design
from .counts import audit

PACKAGES={
 'GSE127530':['edgeR','limma','AnnotationDbi','fgsea','org.Hs.eg.db','reactome.db','jsonlite'],
 'GSE56960':['affy','limma','AnnotationDbi','hgu219.db','fgsea','reactome.db','jsonlite'],
 'GSE27385':['affy','hgu133plus2cdf','limma','AnnotationDbi','fgsea','org.Hs.eg.db','reactome.db','jsonlite']}

def preflight(rscript:Path,libraries:list[Path],accession:str,expected_version:str='4.5.3'):
    if accession not in PACKAGES:raise ValueError('unsupported accession')
    if not rscript.is_file() or any(not p.is_dir() for p in libraries):raise FileNotFoundError('Rscript or R library missing')
    environment={**os.environ,'R_LIBS_USER':os.pathsep.join(str(p.resolve()) for p in libraries),'NUTRIGENOMICS_R_LIB':str(libraries[0].resolve()),'R_THREADS':'1','OMP_NUM_THREADS':'2','OPENBLAS_NUM_THREADS':'2'}
    packages=','.join(json.dumps(p) for p in PACKAGES[accession])
    code=f'needed<-c({packages});missing<-needed[!vapply(needed,requireNamespace,logical(1),quietly=TRUE)];if(length(missing))stop(paste("Missing",paste(missing,collapse=",")));if(as.character(getRversion())!={json.dumps(expected_version)})stop("Unexpected R version");cat(jsonlite::toJSON(list(R=as.character(getRversion()),libraries=.libPaths(),packages=as.list(setNames(vapply(needed,function(x)as.character(packageVersion(x)),character(1)),needed))),auto_unbox=TRUE))'
    process=subprocess.run([str(rscript),'-e',code],env=environment,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,timeout=180)
    if process.returncode:raise RuntimeError('R preflight failed: '+process.stderr)
    report=json.loads(process.stdout);return environment,report


def link_raw(source:Path,target:Path):
    if not source.is_file():raise FileNotFoundError(source)
    target.parent.mkdir(parents=True,exist_ok=True)
    if target.exists():
        if target.resolve()!=source.resolve():raise ValueError('staged raw path conflicts')
        return
    target.symlink_to(source.resolve())


def run_geo(accession:str,data_dir:Path,expansion_root:Path|None,output:Path,rscript:Path,libraries:list[Path],expected_version:str='4.5.3',preflight_only:bool=False):
    project=Path(__file__).resolve().parents[2];script=project/'scripts'/('analyze_'+accession.lower()+'.R')
    output=output.resolve()
    if (output/'run_manifest.json').exists():raise FileExistsError('use a new output directory for each independent run')
    output.mkdir(parents=True,exist_ok=True)
    environment,r_environment=preflight(rscript,libraries,accession,expected_version)
    status={'accession':accession,'started_at':now(),'status':'preflight_passed','R_environment':r_environment,'rscript':str(rscript),'script':file_record(script),'helper':file_record(project/'scripts/bridge_omics_helpers.R'),'old_analysis_results_reused':False}
    if preflight_only:write_json(output/'preflight.json',status);return status
    if accession=='GSE127530':
        report=audit(data_dir,output)
        if report['sample_design']['geo_sample_crosswalk_status']!='complete' or report['status'] not in {'identical_payloads','verified_fixed_cleanup'}:raise ValueError('count payload/sample crosswalk audit failed')
        arguments=['--data-dir',str(data_dir.resolve()),'--output-dir',str(output)]
    elif accession=='GSE56960':
        family=data_dir/'raw/geo/GSE56960/GSE56960_family.soft.gz';prepare_design(accession,family,output/'design')
        arguments=['--data-dir',str(data_dir.resolve()),'--output-dir',str(output),'--sample-sheet',str(output/'design/sample_design.tsv'),'--bridge-contrasts','true','--analysis-version','gse56960-final-dose-response-v1']
    else:
        if expansion_root is None:raise ValueError('GSE27385 requires original raw expansion root')
        original=expansion_root/'raw/GSE27385';stage=output/'staged-inputs';raw=stage/'raw/GSE27385';raw.mkdir(parents=True,exist_ok=True)
        for name in ('GSE27385_family.soft.gz','GSE27385_series_matrix.txt.gz','GSE27385_RAW.tar'):link_raw(original/name,raw/name)
        prepare_design(accession,raw/'GSE27385_family.soft.gz',raw)
        arguments=['--root',str(stage),'--output-dir',str(output),'--raw-cel','true']
    status.update(status='running',command=[str(rscript),str(script),*arguments]);write_json(output/'run_manifest.json',status)
    with (output/'analysis.log').open('w',encoding='utf-8') as log:
        process=subprocess.run(status['command'],cwd=project,env=environment,stdout=log,stderr=subprocess.STDOUT)
    status.update(returncode=process.returncode,finished_at=now(),status='completed' if process.returncode==0 else 'failed')
    if process.returncode==0:
        status['outputs']=[file_record(path) for path in sorted(output.iterdir()) if path.is_file() and path.name!='run_manifest.json']
    write_json(output/'run_manifest.json',status)
    if process.returncode:raise RuntimeError('R analysis failed; inspect '+str(output/'analysis.log'))
    return status
