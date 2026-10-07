"""Separate source-metadata inventory from newly executed cohort analyses."""
from __future__ import annotations
import gzip
import json
import re
from pathlib import Path
from .geo import parse_soft_samples
from .provenance import file_record,now,write_json

COHORTS=[
 {'cohort_id':'S01','series':['GSE127530'],'reported_people':{'GSE127530':5}},
 {'cohort_id':'S02','series':['GSE56960'],'reported_people':{'GSE56960':14}},
 {'cohort_id':'S03','series':['GSE54325','GSE54643'],'containers':['GSE54690'],'reported_people':{'GSE54325':7,'GSE54643':10},'person_level_cross_assay_mapping':'unconfirmed; 7 expression people and 10 methylation people cannot be summed or automatically fused'},
 {'cohort_id':'S04','series':['GSE107205'],'reported_people':{'GSE107205':36}},
 {'cohort_id':'S05','series':['GSE219217'],'reported_people':{'GSE219217':18}},
 {'cohort_id':'S06','series':['GSE13466'],'reported_people':{'GSE13466':21}},
 {'cohort_id':'S07','series':['GSE27385'],'containers':['GSE27387'],'other_children':['GSE27384','GSE27386'],'reported_people':{'GSE27385':11}},
 {'cohort_id':'S08','series':['GSE98645'],'reported_people':{'GSE98645':7},'subset_scope':'7 RNA-seq people; acute/chronic visits are one participant cohort, not independent cohorts'}]


def series_fields(path:Path) -> dict:
    result={};opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',encoding='utf-8') as handle:
        for line in handle:
            match=re.match(r'^!Series_([^=]+?)\s*=\s*(.*)$',line)
            if match:result.setdefault(match[1].strip(),[]).append(match[2].strip())
    return result


def participant(sample:dict,accession:str):
    title=sample['title']
    patterns={
       'GSE127530':r'^(S\d+)-D\d+-(?:Fast|3hr|6hr)$',
       'GSE56960':r'^Subject\s+([^,]+),',
       'GSE54325':r'Volunteer\s*(\d+)',
       'GSE54643':r'^SN\s*(\d+)\s*V\d+$',
       'GSE13466':r'^PBMC (\d+) (?:PUFA|SFA) ',
       'GSE27385':r'^HPBMC_(\d+)[ABC]_(?:FISHOIL|FIBRATE|PLACEBO)$'}
    if accession in patterns:
        match=re.search(patterns[accession],title,re.I)
        if match:return match[1],'study-specific source title'
    for char in sample.get('characteristics',[]):
        key,separator,value=char.partition(':')
        if separator and key.strip().casefold() in {'subject id','subject','participant id','participant','individual','individual id','patient id','patient','person id','donor id'}:
            return value.strip(),'explicit source characteristic '+key.strip()
    return None,'unresolved; no numeric suffix inference'


def locate(accession:str,roots:list[Path]):
    candidates=[]
    for root in roots:
        candidates += [root/'raw/geo'/accession/(accession+'_family.soft.gz'),root/'raw'/accession/(accession+'_family.soft.gz'),root/'geo_metadata'/(accession+'.soft.txt'),root/'snapshots'/(accession+'.soft')]
    return next((path for path in candidates if path.is_file()),None)


def audit_inventory(roots:list[Path],analysis_dirs:list[Path],output:Path):
    executed={}
    for directory in analysis_dirs:
        manifest=directory/'run_manifest.json';metadata=directory/'analysis_metadata.json'
        if not manifest.is_file() or not metadata.is_file():raise ValueError('completed new run manifest and metadata required: '+str(directory))
        run=json.loads(manifest.read_text());info=json.loads(metadata.read_text())
        if run.get('status')!='completed' or run.get('returncode')!=0 or run.get('old_analysis_results_reused') is not False:raise ValueError('run is not newly completed: '+str(directory))
        executed[info['accession']]={'analysis_directory':str(directory.resolve()),'run_manifest':file_record(manifest),'analysis_version':info['analysis_version'],'people':info['independent_people']}
    report={'schema_version':1,'created_at':now(),'cohort_groups':[],'cohort_group_count':8,'patient_multiomics_fusion_confirmed':False,'unique_people_total':'not summed; platform subsets and containers are not independent cohorts'}
    for original in COHORTS:
        cohort={**original,'sources':[]}
        for accession in [*original['series'],*original.get('containers',[]),*original.get('other_children',[])]:
            path=locate(accession,roots);record={'accession':accession,'source_metadata_status':'not_located','fresh_analysis_status':'completed' if accession in executed else 'not_reanalyzed_in_final_atlas'}
            if accession in executed:record['fresh_analysis']=executed[accession]
            container=accession in original.get('containers',[])
            record['contributes_assay_subset_to_cohort_group']=not container and accession in original['series']
            record['count_as_additional_independent_cohort']=False
            if path:
                fields=series_fields(path);record.update(source=file_record(path),source_metadata_status='series_metadata_read',reported_GEO_sample_accessions=len(set(fields.get('sample_id',[]))),series_relationships=fields.get('relation',[]),is_superseries_container=container)
                try:
                    samples=parse_soft_samples(path)
                    subjects=[participant(sample,accession)[0] for sample in samples.values()]
                    resolved=[value for value in subjects if value is not None]
                    record.update(source_metadata_status='full_sample_metadata_read',sample_records=len(samples),resolved_subject_candidates=len(set(resolved)),unresolved_sample_subjects=sum(value is None for value in subjects),subject_count_basis='study-specific source titles or explicit subject characteristics; no patient prediction')
                    if container:record['resolved_subject_candidates']=None;record['subject_count_basis']='SuperSeries container; child records are not new participants'
                    if accession=='GSE54325':record['two_color_arrays']=len(samples);record['paired_biological_samples_if_channels_verified']=2*len(samples) if all(sample['metadata'].get('channel_count')==['2'] for sample in samples.values()) else None
                    if accession in original['reported_people']:record['reported_vs_resolved_agreement']=len(set(resolved))==original['reported_people'][accession] and all(value is not None for value in subjects)
                except ValueError as exc:record['sample_metadata_limit']=str(exc)
            cohort['sources'].append(record)
        report['cohort_groups'].append(cohort)
    write_json(output,report);return report
