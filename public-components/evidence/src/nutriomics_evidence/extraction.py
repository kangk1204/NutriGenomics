"""Conservative rule baseline and schema-checked Qwen evidence parsing.

DrugProt gold mentions are inputs, so this evaluates relation extraction, not NER.
The official evaluation library remains the submission-grade scoring authority.
"""
from __future__ import annotations
import csv
import json
import re
from collections import Counter
from pathlib import Path

RELATIONS={'ACTIVATOR','AGONIST','AGONIST-ACTIVATOR','AGONIST-INHIBITOR','ANTAGONIST',
           'DIRECT-REGULATOR','INDIRECT-DOWNREGULATOR','INDIRECT-UPREGULATOR','INHIBITOR',
           'PART-OF','PRODUCT-OF','SUBSTRATE','SUBSTRATE_PRODUCT-OF'}
RULES=(('ANTAGONIST',r'\bantagonist\b'),('AGONIST',r'\bagonist\b'),('INHIBITOR',r'\binhibit(?:s|ed|or|ion|ing)?\b'),('ACTIVATOR',r'\bactivat(?:es?|ed|or|ion|ing)\b'),('DIRECT-REGULATOR',r'\bbind(?:s|ing)?\b'))
NEGATION=re.compile(r'\b(?:no|not|neither|without|failed to|doesn.t|didn.t)\b',re.I)


def load_drugprot(abstracts:Path,entities:Path) -> list[dict]:
    # Official DrugProt TSV fields are literal text, not CSV-quoted fields.
    # In v1.2 PMID 23194825 starts its abstract with '"Ecstasy"'; CSV's
    # default quote handling deletes those source characters and shifts offsets.
    documents={}
    with abstracts.open(encoding='utf-8',newline='') as handle:
        for row in csv.reader(handle,delimiter='\t',quoting=csv.QUOTE_NONE):
            if len(row)!=3 or row[0] in documents:raise ValueError('invalid or duplicate DrugProt abstract')
            documents[row[0]]={'pmid':row[0],'text':row[1]+' '+row[2],'entities':[]}
    seen=set()
    with entities.open(encoding='utf-8',newline='') as handle:
        for row in csv.reader(handle,delimiter='\t',quoting=csv.QUOTE_NONE):
            if len(row)!=6 or row[0] not in documents:raise ValueError('invalid DrugProt entity row')
            pmid,identifier,kind,start,end,text=row;start,end=int(start),int(end)
            if (pmid,identifier) in seen:raise ValueError('duplicate entity ID')
            seen.add((pmid,identifier))
            if kind not in {'CHEMICAL','GENE','GENE-Y','GENE-N'} or not 0<=start<end<=len(documents[pmid]['text']) or documents[pmid]['text'][start:end]!=text:raise ValueError('entity offset/type/text mismatch')
            documents[pmid]['entities'].append({'id':identifier,'type':kind,'start':start,'end':end,'text':text})
    return list(documents.values())


def sentence_spans(text:str):
    start=0
    for match in re.finditer(r'[.!?](?:\s+|$)',text):
        yield start,match.end();start=match.end()
    if start<len(text):yield start,len(text)


def rule_extract(document:dict) -> list[dict]:
    result=[]
    for start,end in sentence_spans(document['text']):
        sentence=document['text'][start:end]
        if NEGATION.search(sentence):continue
        mentions=[e for e in document['entities'] if start<=e['start'] and e['end']<=end]
        chemicals=[e for e in mentions if e['type']=='CHEMICAL'];genes=[e for e in mentions if e['type'].startswith('GENE')]
        # Abstain on ambiguous multi-pair sentences rather than hallucinate attribution.
        if len(chemicals)!=1 or len(genes)!=1:continue
        for relation,pattern in RULES:
            if re.search(pattern,sentence,re.I):
                result.append({'pmid':document['pmid'],'relation':relation,'arg1':chemicals[0]['id'],'arg2':genes[0]['id'],'evidence_text':sentence,'evidence_start':start,'evidence_end':end,'method':'rule_baseline','tier':'extracted','binding_affinity':None})
                break
    return result


def qwen_prompt(document:dict) -> list[dict]:
    schema={'relations':[{'relation':'one allowed DrugProt type','arg1':'CHEMICAL entity ID','arg2':'GENE entity ID','evidence_text':'exact source substring','evidence_start':'integer char offset','evidence_end':'integer exclusive char offset'}]}
    instruction=('Extract only explicitly stated chemical-protein relations supported by the supplied text. '
                 'Use supplied entity IDs and original character offsets. Do not infer therapeutic benefit, '
                 'binding affinity, causality or species. Abstain on negated/speculative unsupported relations. '
                 'Return JSON only; no markdown. Allowed types: '+','.join(sorted(RELATIONS))+'. Schema: '+json.dumps(schema))
    return [{'role':'system','content':instruction},{'role':'user','content':json.dumps(document,ensure_ascii=False)}]


def parse_qwen(response:str|dict,document:dict) -> dict:
    try:payload=json.loads(response) if isinstance(response,str) else response
    except json.JSONDecodeError as exc:raise ValueError('Qwen output is not JSON') from exc
    if not isinstance(payload,dict) or set(payload)!= {'relations'} or not isinstance(payload['relations'],list):raise ValueError('Qwen JSON must contain only a relations list')
    entities={e['id']:e for e in document['entities']};accepted=[];rejected=[];seen=set()
    for index,row in enumerate(payload['relations']):
        reasons=[]
        if not isinstance(row,dict):rejected.append({'index':index,'reasons':['relation_not_object']});continue
        if set(row)!= {'relation','arg1','arg2','evidence_text','evidence_start','evidence_end'}:reasons.append('unsupported_or_missing_fields')
        if row.get('relation') not in RELATIONS:reasons.append('unsupported_relation')
        a,b=entities.get(row.get('arg1')),entities.get(row.get('arg2'))
        if a is None or b is None or a['type']!='CHEMICAL' or not b['type'].startswith('GENE'):reasons.append('unknown_or_wrong_entity_types')
        start,end=row.get('evidence_start'),row.get('evidence_end')
        span_valid=isinstance(start,int) and not isinstance(start,bool) and isinstance(end,int) and not isinstance(end,bool) and 0<=start<end<=len(document['text'])
        if not span_valid or document['text'][start:end]!=row.get('evidence_text'):reasons.append('invalid_evidence_span')
        elif a and b:
            if not all(start<=e['start'] and e['end']<=end for e in [a,b]):reasons.append('evidence_does_not_contain_both_mentions')
            if NEGATION.search(row['evidence_text']):reasons.append('negation_requires_manual_review')
        key=(row.get('relation'),row.get('arg1'),row.get('arg2'))
        if key in seen:reasons.append('duplicate_relation')
        if reasons:rejected.append({'index':index,'reasons':reasons})
        else:
            seen.add(key);accepted.append({**row,'pmid':document['pmid'],'method':'qwen_json','tier':'extracted','binding_affinity':None,'semantic_support_status':'span_validated_requires_expert_review'})
    return {'pmid':document['pmid'],'accepted':accepted,'rejected':rejected,'semantic_entailment_validated':False}


def write_relations(path:Path,rows:list[dict]) -> None:
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',encoding='utf-8',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        for row in sorted(rows,key=lambda x:(x['pmid'],x['relation'],x['arg1'],x['arg2'])):
            if row['relation'] not in RELATIONS:raise ValueError('unsupported relation for DrugProt export')
            writer.writerow([row['pmid'],row['relation'],'Arg1:'+row['arg1'],'Arg2:'+row['arg2']])


def read_relations(path:Path) -> set[tuple]:
    result=set()
    with path.open(encoding='utf-8',newline='') as handle:
        for row in csv.reader(handle,delimiter='\t'):
            if len(row)!=4 or row[1] not in RELATIONS or not row[2].startswith('Arg1:') or not row[3].startswith('Arg2:'):raise ValueError('invalid DrugProt relation TSV')
            result.add((row[0],row[1],row[2][5:],row[3][5:]))
    return result


def evaluate(gold:Path,predictions:Path,documents:list[dict]) -> dict:
    expected,observed=read_relations(gold),read_relations(predictions)
    lookup={(d['pmid'],e['id']):e for d in documents for e in d['entities']}
    for pmid,relation,a,b in expected|observed:
        if (pmid,a) not in lookup or (pmid,b) not in lookup or lookup[(pmid,a)]['type']!='CHEMICAL' or not lookup[(pmid,b)]['type'].startswith('GENE'):raise ValueError('prediction/gold relation references invalid corpus entities')
    def metrics(a,b):
        tp=len(a&b);fp=len(b-a);fn=len(a-b);p=tp/(tp+fp) if tp+fp else 0.;r=tp/(tp+fn) if tp+fn else 0.
        return {'tp':tp,'fp':fp,'fn':fn,'precision':p,'recall':r,'f1':2*p*r/(p+r) if p+r else 0.}
    return {'micro':metrics(expected,observed),'by_relation':{relation:metrics({r for r in expected if r[1]==relation},{r for r in observed if r[1]==relation}) for relation in sorted(RELATIONS)},'gold_mentions_used':True,'corpus_documents':len(documents),'scope':'exact relation tuples; compare official evaluator before submission','official_evaluator':'https://github.com/tonifuc3m/drugprot-evaluation-library'}
