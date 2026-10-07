"""Read original GEO metadata; never adopt old statistical outputs."""
from __future__ import annotations
import gzip
import json
import re
from pathlib import Path
import pandas as pd
from .provenance import file_record,write_json

def parse_soft_samples(path:Path):
    samples={};current=None
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8') as handle:
        for line in handle:
            start=re.match(r'^\^SAMPLE\s*=\s*(GSM\d+)',line,re.I)
            if start:
                current=start[1].upper()
                if current in samples:raise ValueError('duplicate GEO sample accession')
                samples[current]={'characteristics':[],'metadata':{}};continue
            if current is None:continue
            field=re.match(r'^!Sample_([^=]+?)\s*=\s*(.*)$',line)
            if not field:continue
            key,value=field[1].strip(),field[2].strip()
            if key.casefold()=='title':
                if 'title' in samples[current]:raise ValueError('duplicate sample title')
                samples[current]['title']=value
            elif key.casefold().startswith('characteristics'):samples[current]['characteristics'].append(value)
            else:samples[current]['metadata'].setdefault(key,[]).append(value)
    if not samples or any('title' not in row for row in samples.values()):raise ValueError('GEO sample titles missing')
    return samples


def prepare_design(accession:str,family:Path,output:Path):
    samples=parse_soft_samples(family);output.mkdir(parents=True,exist_ok=True)
    rows=[]
    for gsm,row in samples.items():
        if accession=='GSE56960':
            match=re.fullmatch(r'Subject\s+([^,]+),\s*Meal\s+([A-Z])\s*\(([0-9]+)\s*kcal\),\s*(Fasting|Postprandial)\s*\(([0-9]+)h\)',row['title'],re.I)
            if not match:raise ValueError('GSE56960 title cannot map participant, meal, dose and time')
            rows.append({'sample_accession':gsm,'title':row['title'],'subject_candidate':match[1],'meal':match[2].upper(),'dose_kcal':int(match[3]),'timepoint':match[5]+'h','characteristics_json':json.dumps(row['characteristics'],ensure_ascii=False)})
        elif accession=='GSE27385':
            match=re.fullmatch(r'HPBMC_([0-9]+)([ABC])_(FISHOIL|FIBRATE|PLACEBO)',row['title'])
            if not match:raise ValueError('GSE27385 title cannot map person, period and arm')
            rows.append({'sample_accession':gsm,'title':row['title'],'subject_candidate':match[1],'period':match[2],'arm':match[3],'characteristics_json':json.dumps(row['characteristics'],ensure_ascii=False)})
        else:raise ValueError('unsupported design accession')
    frame=pd.DataFrame(rows)
    if accession=='GSE56960':
        if len(frame)!=168 or frame.subject_candidate.nunique()!=14 or set(frame.dose_kcal)!={500,1000,1500} or set(frame.timepoint)!={'0h','2h','4h','6h'}:raise ValueError('GSE56960 expected 14 people x 3 meals x 4 times')
        if frame.duplicated(['subject_candidate','meal','timepoint']).any() or not (frame.groupby('subject_candidate').size()==12).all():raise ValueError('incomplete repeated meal-time design')
    else:
        if len(frame)!=33 or frame.subject_candidate.nunique()!=11 or frame.duplicated(['subject_candidate','arm']).any() or not (frame.groupby('subject_candidate').size()==3).all():raise ValueError('GSE27385 expected 11 people x 3 crossover arms')
        records=[{'sample_accession':gsm,'title':[row['title']],'characteristics_ch1':row['characteristics'],**row['metadata']} for gsm,row in samples.items()]
        write_json(output/'sample_metadata.json',records)
    frame.to_csv(output/'sample_design.tsv',sep='\t',index=False)
    report={'accession':accession,'samples':len(frame),'independent_people':frame.subject_candidate.nunique(),'source':file_record(family),'metadata_rederived_from_original_soft':True,'old_analysis_results_reused':False}
    write_json(output/'design_audit.json',report);return report
