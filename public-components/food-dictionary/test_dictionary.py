import json
from pathlib import Path
import sqlite3
import unicodedata

import pytest

import dictionary as d
from evaluate import evaluate


@pytest.fixture
def store(tmp_path):
    s=d.Store(tmp_path/'fixture.sqlite',{})
    yield s
    s.con.close()


def molecule(key='AAAAAAAAAAAAAA-UHFFFAOYSA-N',inchi='InChI=1S/fixture'):
    return {'semantic_type':'molecule','full_inchikey':key,'original_inchi':inchi,'single_component':True,'stereo_unspecified':False,'form_scope':'fixture-free-form'}


def proof():
    return {'source_verified':True,'scope_review_passed':True}


def test_unicode_normalization_preserves_original_and_semantic_qualifiers():
    assert d.normalize(unicodedata.normalize('NFD','  캠페롤  '))=='캠페롤'
    assert d.normalize('  Citric   Acid ')=='citric acid'
    assert d.normalize('D-glucose')!=d.normalize('L-glucose')
    assert d.normalize('citric acid')!=d.normalize('citric acid monohydrate')
    assert d.normalize('Chicken, breast, raw')!=d.normalize('Chicken, thigh, cooked')
    assert d.normalize('cis-18:1')!=d.normalize('trans-18:1')


def test_ids_survive_label_changes_and_release_order():
    first=d.identifier('concept','USDA:FDC:food','321358')
    assert first==d.identifier('concept','USDA:FDC:food','321358')
    assert first!=d.identifier('concept','MFDS:food-material-code','321358')
    assert first!=d.identifier('concept','USDA:FDC:component','321358')


@pytest.mark.parametrize('field,value',[('full_inchikey','AAAAAAAAAAAAAA-ABCDEFGHIJ-N'),('original_inchi','InChI=1S/other-tautomer'),('form_scope','fixture-hydrate')])
def test_no_stereo_tautomer_salt_or_hydrate_collapse(field,value):
    left=molecule();right={**left,field:value}
    assert d.decide_equivalence(left,right,proof())[0]=='conflict'


@pytest.mark.parametrize('field,value',[('single_component',False),('stereo_unspecified',True)])
def test_unspecified_or_mixture_forms_abstain(field,value):
    left=molecule();right={**left,field:value}
    assert d.decide_equivalence(left,right,proof())[0]=='candidate'


def test_exact_name_and_complete_key_still_need_scope_evidence():
    assert d.decide_equivalence(molecule(),molecule(),{})[0]=='candidate'
    assert d.decide_equivalence(molecule(),molecule(),proof())[0]=='source_verified'


@pytest.mark.parametrize('field,value',[('organism','other species'),('part','leaf'),('preparation','pickled'),('authority_id','different native code')])
def test_food_organism_part_preparation_and_native_ids_stay_distinct(field,value):
    a={'semantic_type':'food_material','authority_namespace':'fixture','authority_id':'1','organism':'species A','part':'root','preparation':'raw','state':'fresh'}
    assert d.decide_equivalence(a,{**a,field:value},proof())[0]=='conflict'
    assert d.decide_equivalence(a,{**a,'state':None},proof())[0]=='conflict'
    assert d.decide_equivalence({**a,'state':None},{**a,'state':None},proof())[0]=='candidate'


def test_analytical_nutrient_and_molecule_not_equivalent():
    assert d.decide_equivalence(molecule(),{'semantic_type':'analytical_nutrient'},proof())[0]=='conflict'


def test_homonyms_keep_separate_ids_and_one_preferred_per_language(store):
    ev=store.evidence('fixture','local:fixture','1','fixture',{})
    cids=[store.concept('fixture',str(n),'food_material','fixture',{}) for n in (1,2)]
    for cid in cids:store.label(cid,'en','Shared food label','source_reported',ev,role='preferred')
    result=d.retrieve(store.con,'Shared food label',semantic_type='food_material')
    assert {r['concept_id'] for r in result}==set(cids)
    assert all(r['mapping_status']=='candidate' for r in result)
    store.label(cids[0],'en','Official fixture label','official_source',ev,role='preferred')
    assert store.con.execute("SELECT count(*) FROM label WHERE concept_id=? AND language='en' AND role='preferred'",(cids[0],)).fetchone()[0]==1
    assert store.con.execute('SELECT count(*) FROM label').fetchone()[0]==3


def test_duplicate_records_idempotent_but_conflicting_payload_rejected(store):
    ev=store.evidence('fixture','local:fixture','1','fixture',{})
    cid=store.concept('fixture','1','food_material','fixture',{})
    store.record(cid,'fixture','1','v1',{'raw':'unaltered'},ev,'private')
    store.record(cid,'fixture','1','v1',{'raw':'unaltered'},ev,'private')
    assert store.con.execute('SELECT count(*) FROM source_record').fetchone()[0]==1
    with pytest.raises(ValueError,match='Conflicting duplicate'):
        store.record(cid,'fixture','1','v1',{'raw':'changed'},ev,'private')


def test_changed_food_facets_do_not_silently_reuse_a_concept(store):
    store.concept('fixture','1','food_material','fixture',{'part':'leaf'})
    with pytest.raises(ValueError,match='Conflicting concept'):
        store.concept('fixture','1','food_material','fixture',{'part':'root'})


def test_exact_mappings_pinned_guard_idempotence_and_reversal(store):
    key=molecule()['full_inchikey']
    ev=store.evidence('fixture','local:fixture','1','fixture',{'full_inchikey':key,'form_review_passed':True,'original_inchi_exact_match':True})
    cids=[store.concept('fixture',str(n),'molecule','fixture',{'full_inchikey':key}) for n in (1,2,3)]
    decision={'reason':'Synthetic source fixture','full_inchikey':key,'scope_review_passed':True,'identity_evidence_verified':True}
    mid=store.mapping(cids[0],cids[1],'skos:exactMatch','source_verified',ev,decision)
    store.mapping(cids[0],cids[1],'skos:exactMatch','source_verified',ev,decision)
    store.mapping(cids[1],cids[2],'skos:exactMatch','source_verified',ev,decision)
    assert store.con.execute('SELECT count(*) FROM mapping').fetchone()[0]==2
    assert store.con.execute('SELECT count(*) FROM concept').fetchone()[0]==3
    assert store.con.execute('SELECT count(*) FROM mapping WHERE subject_id=? AND object_id=?',(cids[0],cids[2])).fetchone()[0]==0
    with pytest.raises(ValueError,match='Unverified exact'):
        store.mapping(cids[0],cids[2],'skos:exactMatch','candidate',ev,decision)
    store.revoke(mid,'Fixture incorrect scope','fixture-reviewer')
    assert store.con.execute('SELECT status FROM mapping WHERE mapping_id=?',(mid,)).fetchone()[0]=='revoked'
    assert store.con.execute("SELECT count(*) FROM change_event WHERE action='revoke_mapping'").fetchone()[0]==1
    assert store.con.execute('SELECT count(*) FROM concept').fetchone()[0]==3


def test_exact_mapping_cannot_be_fabricated_by_boolean_flags(store):
    ev=store.evidence('fixture','local:fixture','1','fixture',{})
    a=store.concept('fixture','a','molecule','fixture',{'full_inchikey':'key'})
    b=store.concept('fixture','b','molecule','fixture',{'full_inchikey':'key'})
    with pytest.raises(ValueError,match='pinned'):
        store.mapping(a,b,'skos:exactMatch','source_verified',ev,{'reason':'Unsupported','full_inchikey':'key','scope_review_passed':True,'identity_evidence_verified':True})


@pytest.mark.parametrize('predicate',['skos:closeMatch','skos:broadMatch','skos:narrowMatch','skos:relatedMatch'])
def test_nonexact_relations_remain_separate_candidates(store,predicate):
    ev=store.evidence('fixture','local:fixture','1','fixture',{})
    a=store.concept('fixture','a','food_material','fixture',{})
    b=store.concept('fixture','b','food_material','fixture',{})
    store.mapping(a,b,predicate,'candidate',ev,{'reason':'Fixture unreviewed relation'})
    assert store.con.execute('SELECT predicate_id,status FROM mapping').fetchone()==(predicate,'candidate')


def test_no_human_gold_means_no_accuracy_claim():
    result=evaluate([{'reviewer':None,'gold_equivalent':None}])
    assert result['reviewed_rows']==0
    assert result['precision'] is result['recall'] is result['false_merge_rate_among_accepted'] is None
    fixture=[{'reviewer':'fixture-person','gold_equivalent':True,'prediction_equivalent':True,'object_id':'a','retrieved_object_ids':['a']},
             {'reviewer':'fixture-person','gold_equivalent':False,'prediction_equivalent':True,'object_id':'b'},
             {'reviewer':'fixture-person','gold_equivalent':True,'prediction_equivalent':None,'object_id':'c','retrieved_object_ids':[]}]
    result=evaluate(fixture)
    assert result['precision']==result['recall']==result['false_merge_rate_among_accepted']==0.5
    assert result['retrieval_recall_at_k']==0.5


