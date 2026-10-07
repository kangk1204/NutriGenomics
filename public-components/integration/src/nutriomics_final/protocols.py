"""Allowlisted execution contracts and content-based input identity."""
import hashlib
import json
from pathlib import Path


def protocol_id(algorithm,action):
    return f'{algorithm}-{action}-20261003-v1'


def input_identity(root,algorithm,action,parameters):
    root=Path(root).resolve()
    files=[]
    def add(path):
        path=Path(path).resolve()
        if not path.is_relative_to(root) or not path.is_file():raise ValueError('Protocol input must be a file inside research root')
        with path.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
        files.append({'path':path.relative_to(root).as_posix(),'bytes':path.stat().st_size,'sha256':digest})
    def local(value):
        path=(root/str(value)).resolve()
        if not path.is_relative_to(root):raise ValueError('Protocol input escaped research root')
        return path
    if algorithm=='evidence':
        database=local(parameters['database'])
        wal=Path(str(database)+'-wal')
        if wal.exists() and wal.stat().st_size>0:raise ValueError('A checkpointed evidence database snapshot is required')
        add(database)
    elif algorithm=='compound-target':
        add(local(parameters['input']))
        for p in sorted(local(parameters['model_dir']).iterdir()):
            if p.is_file() and p.suffix in {'.json','.pt','.joblib','.gz'}:add(p)
    elif algorithm=='intervention':
        for p in sorted(local(parameters['input_dir']).iterdir()):
            if p.is_file():add(p)
    elif algorithm=='methylation':
        study=local(parameters['root'])
        paths=list((study/'data/prepared').glob('*'))
        if action=='evaluate':
            paths += list((study/'results').rglob('*predictions.tsv'))
            paths += [study/'results/training_config.json']
            paths += list((study/'results').glob('*/final_training_audit.json'))
        if (study/'data/source_manifest.json').is_file():paths.append(study/'data/source_manifest.json')
        for p in sorted(set(paths)):
            if p.is_file():add(p)
    if not files:raise ValueError('No auditable protocol inputs found')
    payload={'protocol_id':protocol_id(algorithm,action),'parameters':parameters,'files':files}
    return payload,hashlib.sha256(json.dumps(payload,sort_keys=True,ensure_ascii=False).encode()).hexdigest()
