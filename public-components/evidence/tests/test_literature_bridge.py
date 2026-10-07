import hashlib
import json
from pathlib import Path
import pytest
from nutriomics_evidence.literature_bridge import import_snapshot, search, status


def report(tmp_path, *, approved=False, invalid_review=False, stale=False):
    quote='DASH lowered blood pressure.'
    assertion={'assertion_id':'a1','fingerprint':'f1','subject':'local:DASH','predicate':'decreases',
        'object':'local:BP','effective':{'subject':'local:DASH','predicate':'decreases','object':'local:BP'},
        'quote':quote,'start':0,'end':len(quote),'source_span_valid':True,'payload':json.dumps({'explicit_relation':True}),
        'display_status':'approved' if approved else 'candidate',
        'reviews':[{'decision':'approved','fingerprint':'wrong' if invalid_review else 'f1'}] if approved else []}
    extraction={'source_hash':'h','source_text_hash':'h','source_kind':'pubmed','source_url':'https://pubmed.ncbi.nlm.nih.gov/1/',
        'source_license':None,'source_id':'s1','extraction_id':'e1','engine_version':'heuristic-v1',
        'dictionary_version':'dict1','result':{'rules_version':'rules1','summary':{'results':[quote]}}}
    p={'pmid':'1','doi':None,'pmcid':None,'title':'DASH experiment','pub_date':'2026','organism':'human',
        'oa_status':'unavailable','license':None,'retracted':False,'topics':[],'extraction':extraction,
        'extraction_current':not stale,'assertions':[assertion]}
    data={'schema_version':'literature-report-v1','generated_at':'2026-10-01T10:00:00+00:00',
        'counts':{'papers_in_report':1,'relation_candidates':1,'approved_relations':int(approved)},
        'quality':{'status':'not_evaluated','f1':None},'papers':[p]}
    path=tmp_path/'report.json';path.write_text(json.dumps(data),encoding='utf-8')
    return path,hashlib.sha256(path.read_bytes()).hexdigest()


def test_candidates_remain_extracted_with_original_provenance(tmp_path):
    path,sha=report(tmp_path);db=tmp_path/'snapshot.sqlite'
    original=path.read_bytes()
    result=import_snapshot(path,db,expected_sha256=sha)
    item=search(db,assertions=True)['items'][0]
    assert item['tier']=='extracted' and item['relation_is_causal'] is False and item['binding_label'] is False
    assert search(db,reviewed=True,assertions=True)['total']==0
    assert path.read_bytes()==original and result['CTD_graph_modified'] is False
    assert status(db)['scientific_accuracy']['f1'] is None


def test_matching_human_review_still_does_not_become_ctd_curated(tmp_path):
    path,sha=report(tmp_path,approved=True);db=tmp_path/'snapshot.sqlite'
    import_snapshot(path,db,expected_sha256=sha)
    assert search(db,assertions=True,reviewed=True)['items'][0]['tier']=='extracted'


def test_mismatched_approval_is_rejected_without_output(tmp_path):
    path,sha=report(tmp_path,approved=True,invalid_review=True);db=tmp_path/'snapshot.sqlite'
    with pytest.raises(ValueError,match='matching explicit review'):
        import_snapshot(path,db,expected_sha256=sha)
    assert not db.exists()


def test_changed_input_hash_is_rejected(tmp_path):
    path,sha=report(tmp_path);path.write_text('{}',encoding='utf-8')
    with pytest.raises(ValueError,match='SHA-256 mismatch'):
        import_snapshot(path,tmp_path/'snapshot.sqlite',expected_sha256=sha)


def test_effective_relation_cannot_replace_approved_original(tmp_path):
    path,_=report(tmp_path,approved=True)
    data=json.loads(path.read_text());data['papers'][0]['assertions'][0]['effective']['object']='local:OTHER'
    path.write_text(json.dumps(data));sha=hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError,match='bound to its matching review'):
        import_snapshot(path,tmp_path/'snapshot.sqlite',expected_sha256=sha)


def test_corrected_review_must_bind_exact_effective_triple(tmp_path):
    path,_=report(tmp_path,approved=True)
    data=json.loads(path.read_text());a=data['papers'][0]['assertions'][0]
    corrected={'subject':'local:DASH','predicate':'decreases','object':'local:REVIEWED'}
    a['effective']=corrected;a['reviews'][0].update({'decision':'corrected','corrected_payload':json.dumps(corrected)})
    path.write_text(json.dumps(data));sha=hashlib.sha256(path.read_bytes()).hexdigest()
    db=tmp_path/'snapshot.sqlite';import_snapshot(path,db,expected_sha256=sha)
    assert search(db,assertions=True,reviewed=True)['items'][0]['effective']['object']=='local:REVIEWED'


def test_single_report_byte_read_binds_declared_hash(monkeypatch,tmp_path):
    path,sha=report(tmp_path);read=Path.read_bytes;count=0
    def intercept(p):
        nonlocal count
        if p==path:
            count+=1
            if count>1:return b'{}'
        return read(p)
    monkeypatch.setattr(Path,'read_bytes',intercept)
    result=import_snapshot(path,tmp_path/'snapshot.sqlite',expected_sha256=sha)
    assert count==1 and result['papers']==1


def test_stale_extractions_not_imported_as_current_assertions(tmp_path):
    path,sha=report(tmp_path,stale=True);db=tmp_path/'snapshot.sqlite'
    r=import_snapshot(path,db,expected_sha256=sha)
    assert r['skipped_stale_assertions']==1 and search(db,assertions=True)['total']==0
    assert search(db)['items'][0]['summary'] is None


def test_pmid_union_deduplicates_collections_and_literal_search(tmp_path):
    path,sha=report(tmp_path);db=tmp_path/'snapshot.sqlite'
    old=tmp_path/'old.json';old.write_text(json.dumps({'unique_documents':2,'searches':{'t1':{'selected_pmids':['1','2']},'t2':{'selected_pmids':['1']}}}))
    r=import_snapshot(path,db,expected_sha256=sha,previous_manifest=old)
    assert r['overlap']['union_unique_pmids']==2 and r['overlap']['shared_pmids']==1
    assert search(db,query='%')['total']==0
    with pytest.raises(FileExistsError):
        import_snapshot(path,db,expected_sha256=sha)
