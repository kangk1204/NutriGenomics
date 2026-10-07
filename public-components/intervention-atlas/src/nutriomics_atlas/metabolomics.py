"""Processed ST001257 abundance audit and person-paired response estimates."""
from __future__ import annotations
import csv
import json
import math
import re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from .provenance import file_record,now,sha256,write_json

FEATURE=re.compile(r'^([0-9]+(?:\.[0-9]+)?)_([0-9]+(?:\.[0-9]+)?)$')

def factors(value:str) -> dict:
    result={}
    for token in value.split('|'):
        key,separator,item=token.partition(':')
        if not separator or key.strip() in result:raise ValueError('invalid or duplicate sample factor')
        result[key.strip()]=item.strip()
    return result


def audit_matrix(directory:Path,output:Path|None=None) -> tuple[pd.DataFrame,pd.DataFrame,dict]:
    manifest=json.loads((directory/'download_manifest.json').read_text(encoding='utf-8'))
    records={row['name']:row for row in manifest['files']}
    for filename in ('untargeted.tsv','factors.json','analysis.json','mwtab.txt'):
        if filename not in records or sha256(directory/filename)!=records[filename]['sha256']:raise ValueError('source manifest missing or hash mismatch: '+filename)
    analysis=json.loads((directory/'analysis.json').read_text(encoding='utf-8'))
    if analysis.get('study_id')!='ST001257' or analysis.get('analysis_id')!='AN002086' or analysis.get('units')!='Abundance (Counts Log2 Transformed)':raise ValueError('wrong study, analysis or unverified abundance units')
    metadata=json.loads((directory/'factors.json').read_text(encoding='utf-8'))
    metadata=list(metadata.values()) if isinstance(metadata,dict) else metadata
    if any(row['study_id']!='ST001257' for row in metadata):raise ValueError('foreign study factor')
    mapping={row['local_sample_id']:row for row in metadata}
    if len(mapping)!=len(metadata) or len(set(row['mb_sample_id'] for row in metadata))!=len(metadata):raise ValueError('duplicate sample identity')
    path=directory/'untargeted.tsv'
    with path.open(encoding='utf-8-sig',newline='') as handle:columns=next(csv.reader(handle,delimiter='\t'))
    if columns[:2]!=['Samples','group'] or len(columns)!=len(set(columns)) or any(not FEATURE.fullmatch(name) for name in columns[2:]):raise ValueError('invalid or duplicate analytical feature headers')
    table=pd.read_csv(path,sep='\t',dtype={'Samples':str,'group':str})
    if table['Samples'].duplicated().any() or set(table['Samples'])!=set(mapping):raise ValueError('matrix and source sample identities are not one-to-one')
    values=table.iloc[:,2:].apply(pd.to_numeric,errors='raise')
    if np.isinf(values.to_numpy(dtype=float)).any():raise ValueError('infinite analytical values')
    samples=[]
    for _,row in table.iloc[:,:2].iterrows():
        original=mapping[row['Samples']];recorded=factors(original['factors']);matrix=factors(row['group'])
        if recorded!=matrix:raise ValueError('abundance group and official factors disagree')
        item={'sample_id':row['Samples'],'repository_sample_id':original['mb_sample_id'],'sample_type':recorded['Sample type'],'subject_id':None,'arm':None,'timepoint':None,'source_factors':original['factors']}
        if item['sample_type']=='Urine':
            subject=re.match(r'^S31-([0-9]+)\s',row['Samples']);condition=re.fullmatch(r'(Pre|Post) (Chicken|Pork)',recorded['Pre or Post and food'])
            if subject is None or condition is None:raise ValueError('unresolved participant or source arm/timepoint')
            item.update(subject_id='S31-'+subject[1],timepoint=condition[1].lower(),arm=condition[2].lower())
            suffix=re.search(r'(?:-|\s)(Ch|Po)$',row['Samples'],re.I)
            item['name_factor_arm_discordance']=bool(suffix and {'ch':'chicken','po':'pork'}[suffix[1].lower()]!=item['arm'])
            item['arm_mapping_basis']='official factor; name supplies participant only'
        elif item['sample_type']!='Food':raise ValueError('unexpected sample type')
        samples.append(item)
    sample_frame=pd.DataFrame(samples);urine=sample_frame[sample_frame.sample_type=='Urine']
    cells=urine.groupby(['subject_id','arm','timepoint']).size()
    if len(urine)!=76 or urine.subject_id.nunique()!=19 or len(cells)!=76 or (cells!=1).any() or set(urine.arm)!= {'chicken','pork'} or set(urine.timepoint)!={'pre','post'}:raise ValueError('expected complete 19 participant x 2 arm x 2 timepoint urine design')
    if len(sample_frame)!=90 or (sample_frame.sample_type=='Food').sum()!=14:raise ValueError('expected 76 urine and 14 food preparations')
    values.index=table.Samples
    report={'accession':'ST001257','analysis_id':'AN002086','created_at':now(),'independent_people':19,'urine_samples':76,'food_preparations':14,'analytical_features':values.shape[1],
            'finite_values':int(np.isfinite(values.to_numpy(dtype=float)).sum()),'missing_values':int(values.isna().sum().sum()),'zero_values':int((values==0).sum().sum()),
            'name_factor_discordance':urine[urine.name_factor_arm_discordance.fillna(False)].sample_id.tolist(),
            'units':analysis['units'],'second_log_transform_applied':False,'sequence_and_period_metadata_available':False,
            'identified_metabolites_confirmed':0,'features_are_neutral_mass_retention_time_pairs':True,
            'patient_multiomics_pairing_confirmed':False,'inputs':[records[name] for name in ('untargeted.tsv','factors.json','analysis.json','mwtab.txt')]}
    if output:
        output.mkdir(parents=True,exist_ok=True);sample_frame.to_csv(output/'sample_metadata.tsv',sep='\t',index=False)
        pd.DataFrame([{'feature_id':'AN002086:'+name,'source_feature':name,'neutral_mass':float(FEATURE.fullmatch(name)[1]),'retention_time_minutes':float(FEATURE.fullmatch(name)[2]),'identified_structure':False} for name in values.columns]).to_csv(output/'feature_crosswalk.tsv',sep='\t',index=False)
        write_json(output/'data_audit.json',report)
    return values,sample_frame,report


def bh(pvalues):
    pvalues=np.asarray(pvalues,dtype=float);out=np.full(len(pvalues),np.nan);valid=np.flatnonzero(np.isfinite(pvalues))
    if not len(valid):return out
    if ((pvalues[valid]<0)|(pvalues[valid]>1)).any():raise ValueError('p outside 0..1')
    order=valid[np.argsort(pvalues[valid],kind='stable')];rank=np.arange(1,len(order)+1)
    adjusted=pvalues[order]*len(order)/rank
    out[order]=np.minimum(1,np.minimum.accumulate(adjusted[::-1])[::-1]);return out


def paired_effect(before:np.ndarray,after:np.ndarray,min_pairs:int=16,zero_policy:str='exclude') -> dict:
    before=np.asarray(before,dtype=float);after=np.asarray(after,dtype=float)
    if before.shape!=after.shape or before.ndim!=2:raise ValueError('paired arrays must have the same person x feature shape')
    if zero_policy not in {'exclude','retain'}:raise ValueError('unknown zero policy')
    valid=np.isfinite(before)&np.isfinite(after)
    if zero_policy=='exclude':valid &= (before!=0)&(after!=0)
    delta=np.where(valid,after-before,np.nan);n=valid.sum(axis=0)
    effect=np.divide(np.nansum(delta,axis=0),n,out=np.full(before.shape[1],np.nan),where=n>0)
    residual=np.where(valid,delta-effect,0.)
    variance=np.divide(np.sum(residual**2,axis=0),n-1,out=np.full(before.shape[1],np.nan),where=n>1)
    se=np.sqrt(variance/np.maximum(n,1));estimable=(n>=min_pairs)&np.isfinite(se)&(se>0)
    statistic=np.divide(effect,se,out=np.full(before.shape[1],np.nan),where=estimable)
    p=2*stats.t.sf(np.abs(statistic),np.maximum(n-1,1));critical=stats.t.ppf(.975,np.maximum(n-1,1))
    return {'n_people':n,'effect_estimate':effect,'standard_error':se,'confidence_low':np.where(estimable,effect-critical*se,np.nan),'confidence_high':np.where(estimable,effect+critical*se,np.nan),'p_value':p,'status':np.where(estimable,'estimable','insufficient_pairs_or_constant_delta')}


def analyze(directory:Path,output:Path,min_pairs:int=16) -> dict:
    if not 3<=min_pairs<=19:raise ValueError('min_pairs must be 3..19')
    if (output/'analysis_metadata.json').exists():raise FileExistsError('use a new directory for each independent analysis')
    values,samples,audit=audit_matrix(directory,output);tables=[];sensitivities=[]
    for arm in ('chicken','pork'):
        subset=samples[(samples.sample_type=='Urine')&(samples.arm==arm)]
        index=subset.pivot(index='subject_id',columns='timepoint',values='sample_id').sort_index()
        before=values.loc[index['pre']].to_numpy(dtype=float);after=values.loc[index['post']].to_numpy(dtype=float)
        for policy,target in [('exclude',tables),('retain',sensitivities)]:
            result=paired_effect(before,after,min_pairs,policy);table=pd.DataFrame(result)
            table['feature_id']=['AN002086:'+name for name in values.columns];table['source_feature']=values.columns
            table['contrast']=arm+'_dash_post_vs_pre';table['intervention']='6-week DASH-style '+arm+' protein arm';table['tissue']='24-hour urine'
            table['effect_unit']='difference in depositor log2 abundance';table['zero_policy']=policy;table['q_value_contrast']=bh(table.p_value)
            target.append(table)
    primary=pd.concat(tables,ignore_index=True);sensitivity=pd.concat(sensitivities,ignore_index=True)
    for table,name in [(primary,'metabolite_effects.tsv.gz'),(sensitivity,'zero_retention_sensitivity.tsv.gz')]:
        table['q_value_study_family']=bh(table.p_value);table['q_value']=table['q_value_study_family'];table.to_csv(output/name,sep='\t',index=False,compression='gzip')
    row_index=samples.set_index('sample_id').loc[values.index]
    sample_qc=pd.DataFrame({'sample_id':values.index,'sample_type':row_index.sample_type.to_numpy(),'finite_features':np.isfinite(values).sum(axis=1),'zero_features':(values==0).sum(axis=1),'median_source_log2':values.median(axis=1)})
    sample_qc.to_csv(output/'sample_qc.tsv',sep='\t',index=False)
    metadata={**audit,'analysis_version':'st001257-paired-processed-v1','independent_unit':'participant; 19 paired people per arm, not 76 independent observations',
              'primary_method':'two-sided one-sample Student t test of within-person post-minus-pre differences; finite nonzero source values; minimum '+str(min_pairs)+' complete pairs',
              'zero_policy':'Zero meaning not documented in processed source; primary treats zeros as unknown, sensitivity retains them. No imputation or second log transformation.',
              'multiple_testing':'BH across all estimable feature x two within-arm contrasts; within-contrast q also retained; sensitivity separate',
              'identified_metabolite_claims':False,'blood_pressure_outcome_modeled':False,'causal_between_arm_effect_estimated':False,
              'limitations':['Processed abundance; raw instrument and identification QC not repeated.','Participant ID taken from sample names; official factor gives arm despite two name discrepancies.','Period/sequence metadata absent; within-arm before-after changes are exploratory and may include time/period effects.','Anonymous mass/RT features are not validated compound identities; no patient transcriptome/metabolome fusion.'],
              'result_rows':len(primary),'estimable_rows':int((primary.status=='estimable').sum()),'q_below_005':int((primary.q_value<.05).sum()),
              'result_files':[file_record(output/name) for name in ('metabolite_effects.tsv.gz','zero_retention_sensitivity.tsv.gz','sample_metadata.tsv','feature_crosswalk.tsv','sample_qc.tsv')]}
    write_json(output/'analysis_metadata.json',metadata);return metadata
