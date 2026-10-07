"""Official ST001257 source acquisition with immutable byte-level provenance."""
from __future__ import annotations
import os
from pathlib import Path
import requests
from .provenance import file_record,now,write_json

URLS={
 'untargeted.tsv':'https://www.metabolomicsworkbench.org/rest/study/analysis_id/AN002086/untarg_data',
 'mwtab.txt':'https://www.metabolomicsworkbench.org/data/study_textformat_view.php?ANALYSIS_ID=AN002086&MODE=d&STUDY_ID=ST001257',
 'factors.json':'https://www.metabolomicsworkbench.org/rest/study/study_id/ST001257/factors',
 'study-summary.json':'https://www.metabolomicsworkbench.org/rest/study/study_id/ST001257/summary',
 'analysis.json':'https://www.metabolomicsworkbench.org/rest/study/study_id/ST001257/analysis',
 'study-page.html':'https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST001257'}

def fetch_st001257(output:Path) -> dict:
    output.mkdir(parents=True,exist_ok=True)
    manifest={'accession':'ST001257','analysis_id':'AN002086','project':'PR000843','project_doi':'10.21228/M8QH5G','retrieved_at':now(),'files':[],'errors':[],
              'source_kind':'official_repository','raw_instrument_reprocessing':False,'rights':'Retain NMDR source attribution and terms; no independent redistribution determination'}
    for name,url in URLS.items():
        path=output/name
        try:
            if not path.exists():
                response=requests.get(url,timeout=(15,60),headers={'User-Agent':'Nutriomics research source audit/0.1'})
                response.raise_for_status()
                if not response.content:raise ValueError('empty repository response')
                if name.endswith('.json'):response.json()
                if name=='mwtab.txt' and b'#METABOLOMICS WORKBENCH' not in response.content[:2000]:raise ValueError('response is not mwTab')
                if name=='untargeted.tsv' and not response.content.startswith(b'Samples\tgroup\t'):raise ValueError('response is not the expected analysis data matrix')
                temporary=path.with_suffix(path.suffix+'.part');temporary.write_bytes(response.content);os.replace(temporary,path)
            manifest['files'].append({'name':name,**file_record(path,url)})
        except Exception as exc:manifest['errors'].append({'name':name,'error':str(exc)})
        write_json(output/'download_manifest.json',manifest)
    return manifest
