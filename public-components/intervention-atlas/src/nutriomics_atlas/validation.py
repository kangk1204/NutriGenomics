from __future__ import annotations
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
from .metabolomics import bh
from .provenance import sha256

def validate_output(directory:Path) -> dict:
    metadata=json.loads((directory/'analysis_metadata.json').read_text(encoding='utf-8'));errors=[];tables={}
    expected=['metabolite_effects.tsv.gz'] if metadata['accession']=='ST001257' else ['gene_effects.tsv.gz','pathway_effects.tsv.gz']
    for name in expected:
        path=directory/name
        if not path.is_file():errors.append('missing result '+name);continue
        table=pd.read_csv(path,sep='\t');required={'feature_id','contrast','effect_estimate','p_value','q_value','n_people','n_samples'} if name!='metabolite_effects.tsv.gz' else {'feature_id','contrast','effect_estimate','p_value','q_value','n_people'}
        if not required.issubset(table):errors.append('missing result fields '+name);continue
        if table.duplicated(['feature_id','contrast']).any():errors.append('duplicate feature-contrast '+name)
        if ((table.p_value.dropna()<0)|(table.p_value.dropna()>1)).any() or ((table.q_value.dropna()<0)|(table.q_value.dropna()>1)).any():errors.append('invalid p/q range '+name)
        if (table.n_people>metadata['independent_people']).any():errors.append('person count exceeds independent cohort '+name)
        if 'q_value_study_family' not in table:errors.append('primary study-level testing family absent '+name)
        elif not np.allclose(table.q_value,bh(table.p_value),atol=1e-10,equal_nan=True):errors.append('BH family inconsistent '+name)
        if name!='pathway_effects.tsv.gz' and not {'confidence_low','confidence_high','standard_error'}.issubset(table):errors.append('gene/metabolite uncertainty absent '+name)
        finite=table.p_value.notna();tables[name]={'rows':len(table),'features':table.feature_id.nunique(),'contrasts':table.contrast.nunique(),'estimable_tests':int(finite.sum()),'q_below_005':int((table.q_value<.05).sum()),'sha256':sha256(path)}
    run=directory/'run_manifest.json'
    for record in metadata.get('result_files',[]):
        if not isinstance(record,dict):continue
        filename=re.split(r'[\\/]',record['path'])[-1]
        if not (directory/filename).is_file() or sha256(directory/filename)!=record['sha256']:errors.append('result provenance hash mismatch '+filename)
    if run.exists():
        status=json.loads(run.read_text(encoding='utf-8'))
        if status['status']!='completed' or status['returncode']!=0:errors.append('R run is not completed')
        for record in status.get('outputs',[]):
            filename=re.split(r'[\\/]',record['path'])[-1]
            if not (directory/filename).is_file() or sha256(directory/filename)!=record['sha256']:errors.append('output provenance hash mismatch '+filename)
    return {'accession':metadata['accession'],'independent_people':metadata['independent_people'],'tables':tables,'errors':errors,'patient_multiomics_fusion_confirmed':False}
