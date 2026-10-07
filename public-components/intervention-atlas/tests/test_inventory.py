import json
from nutriomics_atlas.inventory import audit_inventory,participant

def test_superseries_and_unexecuted_assays_not_promoted(tmp_path):
    source=tmp_path/'source';(source/'geo_metadata').mkdir(parents=True)
    (source/'geo_metadata/GSE54690.soft.txt').write_text('^SERIES = GSE54690\n!Series_sample_id = GSM1\n!Series_relation = SuperSeries of: GSE54325\n!Series_relation = SuperSeries of: GSE54643\n')
    output=tmp_path/'inventory.json';result=audit_inventory([source],[],output)
    assert result['cohort_group_count']==8
    group=next(item for item in result['cohort_groups'] if item['cohort_id']=='S03')
    container=next(item for item in group['sources'] if item['accession']=='GSE54690')
    assert container['is_superseries_container'] and container['count_as_additional_independent_cohort'] is False
    assert container['fresh_analysis_status']=='not_reanalyzed_in_final_atlas'
    assert result['patient_multiomics_fusion_confirmed'] is False

def test_participant_numbers_require_study_specific_or_explicit_basis():
    assert participant({'title':'HPBMC_19B_FISHOIL'},'GSE27385')[0]=='19'
    assert participant({'title':'Unknown sample 19','characteristics':[]},'GSE219217')[0] is None
    assert participant({'title':'Unknown sample 19','characteristics':['subject id: P1']},'GSE219217')[0]=='P1'
