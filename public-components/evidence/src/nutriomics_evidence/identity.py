"""Reviewed identity mappings; names and regulatory identifiers are not structures."""
from __future__ import annotations
import json
import re
from pathlib import Path
from .graph import canonical,connect,node,source,validate
from .sources import file_record,registry_hook

KEY=re.compile(r'^[A-Z]{14}-[A-Z]{10}-[A-Z]$')


def validate_mapping(row:dict) -> dict:
    required={'subject','namespace','identifier','relation','basis'}
    if not required.issubset(row):raise ValueError('missing identity mapping fields')
    if row['relation'] not in {'source_reported','exact_structure','related','candidate_name','regulatory_identity'}:raise ValueError('unsupported identity relation')
    if row['relation']=='exact_structure':
        if row.get('method')!='same_standard_inchikey' or not KEY.fullmatch(row.get('subject_inchikey','')) or row.get('subject_inchikey')!=row.get('object_inchikey'):
            raise ValueError('exact structure requires identical valid standard InChIKeys')
        if row.get('reviewed') is not True or row.get('mixture') is True:raise ValueError('exact structure requires reviewed isolated chemical form')
    if row['namespace']=='FDA-UNII' and row['relation']=='exact_structure' and row.get('method')!='same_standard_inchikey':raise ValueError('UNII alone is regulatory identity')
    return row


def ingest_identity(database:Path,path:Path,format:str,registry=None):
    payload=json.loads(path.read_text(encoding='utf-8'));record=file_record(path,format,'https://pubchem.ncbi.nlm.nih.gov/' if format=='pubchem' else payload.get('source_url',''),payload.get('version') if isinstance(payload,dict) else None,'Retain original source conditions; no name-only exact identity')
    db=connect(database);count=0
    try:
        with db:
            source(db,record)
            if format=='pubchem':
                props=payload.get('PropertyTable',{}).get('Properties',[])
                if not props:raise ValueError('PubChem PUG property JSON required')
                for p in props:
                    identifier='PubChem:'+str(p['CID']);node(db,identifier,'chemical',p.get('IUPACName',identifier),{'InChIKey':p.get('InChIKey'),'SMILES':p.get('ConnectivitySMILES',p.get('CanonicalSMILES')),'source_cid':p['CID']})
                    if p.get('InChIKey'):
                        if not KEY.fullmatch(p['InChIKey']):raise ValueError('invalid PubChem InChIKey')
                        db.execute('INSERT OR IGNORE INTO identifier_mapping VALUES (?,?,?,?,?,?)',(identifier,'InChIKey',p['InChIKey'],'source_reported',record['id'],'PubChem CID structure property'))
                    count+=1
            elif format=='reviewed-mappings':
                for row in payload['mappings']:
                    validate_mapping(row)
                    if not db.execute('SELECT 1 FROM node WHERE id=?',(row['subject'],)).fetchone():raise ValueError('mapping subject does not exist')
                    if row['relation']=='exact_structure':
                        known={r[0] for r in db.execute("SELECT identifier FROM identifier_mapping WHERE subject=? AND namespace='InChIKey' AND relation='source_reported'",(row['subject'],))}
                        target=row['namespace']+':'+row['identifier']
                        target_known={r[0] for r in db.execute("SELECT identifier FROM identifier_mapping WHERE subject=? AND namespace='InChIKey' AND relation='source_reported'",(target,))}
                        if row['subject_inchikey'] not in known or row['object_inchikey'] not in target_known:raise ValueError('exact structure requires separately ingested source structure properties for both nodes')
                    db.execute('INSERT OR IGNORE INTO identifier_mapping VALUES (?,?,?,?,?,?)',(row['subject'],row['namespace'],row['identifier'],row['relation'],record['id'],canonical(row)))
                    count+=1
            else:raise ValueError('unsupported identity format')
        registry_hook(registry,record);result=validate(db);result['mappings_ingested']=count;return result
    finally:db.close()
