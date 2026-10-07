import json
import pytest
from nutriomics_evidence.extraction import rule_extract,parse_qwen,load_drugprot,write_relations,evaluate

@pytest.fixture
def doc():
    return {'pmid':'123','text':'X inhibits ACE.','entities':[{'id':'T1','type':'CHEMICAL','start':0,'end':1,'text':'X'},{'id':'T2','type':'GENE-Y','start':11,'end':14,'text':'ACE'}]}

def response():return {'relations':[{'relation':'INHIBITOR','arg1':'T1','arg2':'T2','evidence_text':'X inhibits ACE.','evidence_start':0,'evidence_end':15}]}

def test_rule_extract_retains_span_without_affinity(doc):
    result=rule_extract(doc);assert len(result)==1 and result[0]['relation']=='INHIBITOR';assert result[0]['binding_affinity'] is None

def test_rule_negation_abstains(doc):
    doc['text']='X does not inhibit ACE.';doc['entities'][1].update(start=19,end=22)
    assert rule_extract(doc)==[]

def test_qwen_supported_span_still_not_semantic_validation(doc):
    result=parse_qwen(response(),doc);assert len(result['accepted'])==1;assert result['semantic_entailment_validated'] is False

@pytest.mark.parametrize('field,value',[('relation','CURES_HYPERTENSION'),('arg2','T99'),('evidence_text','fabricated'),('evidence_start',True)])
def test_qwen_rejects_unsupported_relation_or_span(doc,field,value):
    payload=response();payload['relations'][0][field]=value;result=parse_qwen(payload,doc);assert result['accepted']==[] and result['rejected']

def test_qwen_extra_affinity_rejected(doc):
    payload=response();payload['relations'][0]['Kd']=0.01;assert parse_qwen(payload,doc)['accepted']==[]

def test_drugprot_roundtrip_offsets_and_micro_counts(tmp_path,doc):
    abstracts=tmp_path/'abstracts.tsv';entities=tmp_path/'entities.tsv';gold=tmp_path/'gold.tsv';pred=tmp_path/'pred.tsv'
    abstracts.write_text('123\tX\tinhibits ACE.\n');entities.write_text('123\tT1\tCHEMICAL\t0\t1\tX\n123\tT2\tGENE-Y\t11\t14\tACE\n')
    docs=load_drugprot(abstracts,entities);rows=rule_extract(docs[0]);write_relations(gold,rows);write_relations(pred,rows)
    result=evaluate(gold,pred,docs);assert result['micro']['tp']==1 and result['micro']['fp']==0 and result['gold_mentions_used']
    entities.write_text('123\tT1\tCHEMICAL\t1\t2\tX\n')
    with pytest.raises(ValueError):load_drugprot(abstracts,entities)
