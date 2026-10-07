"""Complete CSV/schema/CRC audit without materializing additional data tables."""
from __future__ import annotations
import json
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
from .ctd import read_ctd
from .sources import atomic_json,now,sha256


def audit_file(path_string:str) -> dict:
    path=Path(path_string)
    original_hash=sha256(path)
    stream=read_ctd(path);header=next(stream);rows=sum(1 for _ in stream)
    if sha256(path)!=original_hash:raise ValueError('source changed during audit')
    return {'name':path.name,'sha256':original_hash,'bytes':path.stat().st_size,
            'rows':rows,'fields':header['fields'],'report_created':header['report_created'],
            'csv_schema_width':'passed','gzip_crc':'passed' if path.suffix=='.gz' else 'not_applicable',
            'validation':'complete_stream_to_eof'}


def audit_ctd(source_dir:Path,output:Path,previous_manifest:Path|None=None,workers:int=1) -> dict:
    if workers<1 or workers>16:raise ValueError('workers must be 1..16')
    files=sorted([*source_dir.glob('CTD_*.csv.gz'),*source_dir.glob('CTD_*.csv')])
    if not files:raise ValueError('no completed CTD CSV files')
    previous=json.loads(previous_manifest.read_text(encoding='utf-8'))['files'] if previous_manifest else {}
    result={'schema_version':1,'source_dir':str(source_dir.resolve()),'created_at':now(),
            'files':{},'errors':[],'excluded_incomplete_files':'*.crdownload are not read'}
    pending=[]
    for path in files:
        prior=previous.get(path.name)
        if prior and prior.get('rows') is not None and prior.get('gzip_crc') in {'passed','not_applicable'} and sha256(path)==prior['sha256']:
            result['files'][path.name]={key:prior[key] for key in ('sha256','bytes','rows','fields','report_created','gzip_crc')}
            result['files'][path.name].update({'csv_schema_width':'passed','validation':'prior_complete_ingestion_with_current_source_hash_verified'})
        else:pending.append(path)
    def store(path,data=None,error=None):
        if error is not None:result['errors'].append({'file':path.name,'error':str(error)})
        else:result['files'][path.name]=data
        atomic_json(output,result)
    if workers==1:
        for path in pending:
            try:store(path,audit_file(str(path)))
            except Exception as exc:store(path,error=exc)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures={pool.submit(audit_file,str(path)):path for path in pending}
            for future in as_completed(futures):
                try:store(futures[future],future.result())
                except Exception as exc:store(futures[future],error=exc)
    result['files']=dict(sorted(result['files'].items()));result['expected_files']=len(files)
    result['validated_files']=len(result['files']);result['complete']=not result['errors'] and len(files)==len(result['files'])
    atomic_json(output,result);return result
