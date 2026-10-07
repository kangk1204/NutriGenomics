import gzip
import json
from pathlib import Path
import numpy as np
import pytest
from scipy import stats
from nutriomics_atlas.metabolomics import factors,bh,paired_effect,audit_matrix
from nutriomics_atlas.provenance import sha256
from nutriomics_atlas.geo import parse_soft_samples
from nutriomics_atlas.geo import prepare_design
from nutriomics_atlas.provenance import file_record
from nutriomics_atlas.counts import _parse_header

def test_paired_effect_uses_people_and_matches_student_t():
    before=np.array([[10],[11],[12],[13],[14]],dtype=float);delta=np.array([1,3,2,4,2.])
    result=paired_effect(before,before+delta[:,None],3)
    assert result['n_people'][0]==5
    assert result['effect_estimate'][0]==pytest.approx(delta.mean())
    assert result['standard_error'][0]==pytest.approx(delta.std(ddof=1)/np.sqrt(5))
    assert result['p_value'][0]==pytest.approx(stats.ttest_1samp(delta,0).pvalue)
    assert result['confidence_low'][0]<delta.mean()<result['confidence_high'][0]

def test_missing_and_zero_are_not_imputed():
    before=np.array([[10],[0],[np.nan],[11],[12],[14.]])
    after=np.array([[12],[3],[20],[13],[15],[15.]])
    primary=paired_effect(before,after,3);sensitivity=paired_effect(before,after,3,'retain')
    assert primary['n_people'][0]==4 and sensitivity['n_people'][0]==5
    assert primary['effect_estimate'][0]==pytest.approx(2.)
    with pytest.raises(ValueError):paired_effect(before,after[:3],3)

def test_insufficient_or_constant_delta_not_false_significance():
    result=paired_effect(np.ones((19,1)),np.ones((19,1))*2,16)
    assert np.isnan(result['p_value'][0]) and result['status'][0]=='insufficient_pairs_or_constant_delta'

def test_bh_joint_family_keeps_missing_unknown():
    q=bh([.01,.04,.03,np.nan])
    assert np.allclose(q[:3],[.03,.04,.04]) and np.isnan(q[3])
    with pytest.raises(ValueError):bh([-.1])

def test_factor_parser_rejects_ambiguous_keys():
    assert factors('Sample type:Urine | Food type:urine')['Sample type']=='Urine'
    with pytest.raises(ValueError):factors('Food type:one | Food type:two')

def test_soft_parser_preserves_original_titles_and_characteristics(tmp_path):
    path=tmp_path/'source.soft.gz';path.write_bytes(gzip.compress(b'^SAMPLE = GSM1\n!Sample_title = Test participant\n!Sample_characteristics_ch1 = participant: 001\n!Sample_platform_id = GPL1\n'))
    result=parse_soft_samples(path)
    assert result['GSM1']['title']=='Test participant' and result['GSM1']['characteristics']==['participant: 001']
    assert result['GSM1']['metadata']['platform_id']==['GPL1']

def test_soft_duplicate_accessions_rejected(tmp_path):
    path=tmp_path/'x.soft';path.write_text('^SAMPLE = GSM1\n!Sample_title = first\n^SAMPLE = GSM1\n!Sample_title = second\n')
    with pytest.raises(ValueError):parse_soft_samples(path)

def test_matrix_hash_mismatch_blocks_analysis(tmp_path):
    (tmp_path/'untargeted.tsv').write_text('Samples\tgroup\t137.0_1.5\n')
    (tmp_path/'download_manifest.json').write_text(json.dumps({'files':[{'name':'untargeted.tsv','sha256':'wrong'}]}))
    with pytest.raises(ValueError,match='hash mismatch'):audit_matrix(tmp_path)


def test_official_factor_overrides_name_suffix_and_log2_is_not_transformed(tmp_path):
    rows=[];records=[]
    for person in range(101,120):
        for arm,short in [('Chicken','Ch'),('Pork','Po')]:
            for time in ['Pre','Post']:
                suffix=({'Ch':'Po','Po':'Ch'}[short] if person==119 and time=='Post' else short)
                name=f'S31-{person} {time}-{suffix}';factor=f'Sample type:Urine | Food type:urine | Pre or Post and food:{time} {arm}'
                rows.append([name,factor,7.]);records.append({'study_id':'ST001257','local_sample_id':name,'mb_sample_id':'SA'+str(len(rows)),'factors':factor})
    for i in range(14):
        name='food'+str(i);factor='Sample type:Food | Food type:Food | Pre or Post and food:na';rows.append([name,factor,7.]);records.append({'study_id':'ST001257','local_sample_id':name,'mb_sample_id':'SA'+str(len(rows)),'factors':factor})
    import pandas as pd
    pd.DataFrame(rows,columns=['Samples','group','137.0_1.5']).to_csv(tmp_path/'untargeted.tsv',sep='\t',index=False)
    (tmp_path/'factors.json').write_text(json.dumps(records));(tmp_path/'analysis.json').write_text(json.dumps({'study_id':'ST001257','analysis_id':'AN002086','units':'Abundance (Counts Log2 Transformed)'}));(tmp_path/'mwtab.txt').write_text('#METABOLOMICS WORKBENCH')
    (tmp_path/'download_manifest.json').write_text(json.dumps({'files':[{'name':name,**file_record(tmp_path/name)} for name in ['untargeted.tsv','factors.json','analysis.json','mwtab.txt']]}))
    values,metadata,audit=audit_matrix(tmp_path)
    assert values.iloc[0,0]==7. and audit['second_log_transform_applied'] is False
    assert len(audit['name_factor_discordance'])==2
    match=metadata[metadata.sample_id=='S31-119 Post-Po'].iloc[0]
    assert match.arm=='chicken' and match.subject_id=='S31-119'


def test_count_header_independent_people_not_sample_columns():
    labels=[f'S{s}-D{d}-{t}_S{index}' for index,(s,d,t) in enumerate((s,d,t) for s in range(1,6) for d in range(1,4) for t in ['Fast','3hr','6hr'])]
    assert len(_parse_header('\t'.join(labels),Path('source.gz')))==45
    labels[-1]=labels[0]
    with pytest.raises(ValueError):_parse_header('\t'.join(labels),Path('source.gz'))


def test_geo_period_and_person_structure_from_original_metadata(tmp_path):
    path=tmp_path/'source.soft.gz';lines=[]
    for person in range(1,12):
        for index,arm in enumerate(['FISHOIL','FIBRATE','PLACEBO']):
            lines.extend([f'^SAMPLE = GSM{person*10+index}',f'!Sample_title = HPBMC_{person}{"ABC"[index]}_{arm}','!Sample_characteristics_ch1 = source: PBMC'])
    path.write_bytes(gzip.compress('\n'.join(lines).encode()))
    result=prepare_design('GSE27385',path,tmp_path/'out')
    assert result['independent_people']==11 and result['samples']==33 and result['old_analysis_results_reused'] is False
