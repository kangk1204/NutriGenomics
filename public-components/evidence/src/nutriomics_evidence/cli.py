from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
from .sources import atomic_json


def main(argv=None):
    parser=argparse.ArgumentParser(description='Source-backed nutrition evidence research CLI')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('ingest-ctd');p.add_argument('--source-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--files',nargs='+');p.add_argument('--source-registry',type=Path)
    p=sub.add_parser('audit-ctd');p.add_argument('--source-dir',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--previous-manifest',type=Path);p.add_argument('--workers',type=int,default=1)
    p=sub.add_parser('build-graph');p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--database',type=Path,required=True);p.add_argument('--include-inferred',action='store_true');p.add_argument('--overwrite',action='store_true')
    p=sub.add_parser('ingest-food');p.add_argument('--format',choices=['fdc','kfind','rda'],required=True);p.add_argument('--input',nargs='+',type=Path,required=True);p.add_argument('--database',type=Path,required=True);p.add_argument('--selection',type=Path);p.add_argument('--mapping',type=Path);p.add_argument('--source-registry',type=Path)
    p=sub.add_parser('ingest-identities');p.add_argument('--format',choices=['pubchem','reviewed-mappings'],required=True);p.add_argument('--input',type=Path,required=True);p.add_argument('--database',type=Path,required=True);p.add_argument('--source-registry',type=Path)
    for command in ['validate','query','export','meal']:
        p=sub.add_parser(command);p.add_argument('--database',type=Path,required=True)
        if command in ['query','export']:
            p.add_argument('--subject');p.add_argument('--disease');p.add_argument('--tier');p.add_argument('--limit',type=int,default=20 if command=='query' else 100000);p.add_argument('--offset',type=int,default=0)
        if command=='export':p.add_argument('--output',type=Path,required=True)
        if command=='meal':p.add_argument('--input',type=Path,required=True)
    p=sub.add_parser('extract-drugprot');p.add_argument('--abstracts',type=Path,required=True);p.add_argument('--entities',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--evidence-output',type=Path);p.add_argument('--prompt-output',type=Path)
    p=sub.add_parser('parse-qwen');p.add_argument('--documents',type=Path,required=True);p.add_argument('--responses',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--relations-output',type=Path)
    p=sub.add_parser('evaluate-drugprot');p.add_argument('--abstracts',type=Path,required=True);p.add_argument('--entities',type=Path,required=True);p.add_argument('--gold',type=Path,required=True);p.add_argument('--predictions',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    if args.command=='audit-ctd':
        from .audit import audit_ctd
        result=audit_ctd(args.source_dir,args.output,args.previous_manifest,args.workers)
        result={key:result[key] for key in ('expected_files','validated_files','complete','errors')}
    elif args.command=='ingest-ctd':
        from .ctd import ingest_ctd,DEFAULT_FILES
        result=ingest_ctd(args.source_dir,args.output_dir,args.files or DEFAULT_FILES,args.source_registry)
        result={'files':{name:{'rows':r['rows'],'sha256':r['sha256'],'gzip_crc':r['gzip_crc']} for name,r in result['files'].items()}}
    elif args.command=='build-graph':
        from .graph import build_graph
        result=build_graph(args.input_dir,args.database,args.overwrite,args.include_inferred)
    elif args.command=='ingest-food':
        from .food import ingest_food
        result=ingest_food(args.database,args.format,args.input,args.selection,args.mapping,args.source_registry)
    elif args.command=='ingest-identities':
        from .identity import ingest_identity
        result=ingest_identity(args.database,args.input,args.format,args.source_registry)
    elif args.command in ['validate','query','export','meal']:
        from .graph import connect,validate,query
        db=connect(args.database,readonly=True)
        try:
            if args.command=='validate':result=validate(db)
            elif args.command=='meal':
                from .food import meal
                result=meal(db,json.loads(args.input.read_text(encoding='utf-8'))['foods'])
            else:
                rows=query(db,args.subject,args.disease,args.tier,args.limit,args.offset)
                if args.command=='export':atomic_json(args.output,{'schema_version':1,'evidence':rows,'redistribution_verified':False});result={'rows':len(rows),'output':str(args.output)}
                else:result=rows
        finally:db.close()
    elif args.command=='extract-drugprot':
        from .extraction import load_drugprot,rule_extract,write_relations,qwen_prompt
        documents=load_drugprot(args.abstracts,args.entities);rows=[r for d in documents for r in rule_extract(d)];write_relations(args.output,rows)
        if args.evidence_output:atomic_json(args.evidence_output,{'documents':documents,'relations':rows})
        if args.prompt_output:
            args.prompt_output.parent.mkdir(parents=True,exist_ok=True)
            with args.prompt_output.open('w',encoding='utf-8') as handle:
                for d in documents:handle.write(json.dumps({'pmid':d['pmid'],'messages':qwen_prompt(d)},ensure_ascii=False)+'\n')
        result={'documents':len(documents),'rule_relations':len(rows),'gold_mentions_used':True,'Qwen_inference_performed':False}
    elif args.command=='parse-qwen':
        from .extraction import parse_qwen,write_relations
        payload=json.loads(args.documents.read_text(encoding='utf-8'));documents=payload['documents'] if isinstance(payload,dict) else payload
        responses=[json.loads(line) for line in args.responses.read_text(encoding='utf-8').splitlines() if line.strip()];lookup={d['pmid']:d for d in documents};reports=[];seen=set()
        for response in responses:
            if response['pmid'] in seen or response['pmid'] not in lookup:raise ValueError('duplicate or unknown response PMID')
            seen.add(response['pmid']);reports.append(parse_qwen(response['response'],lookup[response['pmid']]))
        result={'documents_responded':len(reports),'documents_missing':sorted(set(lookup)-seen),'accepted':sum(len(r['accepted']) for r in reports),'rejected':sum(len(r['rejected']) for r in reports),'reports':reports};atomic_json(args.output,result)
        if args.relations_output:write_relations(args.relations_output,[r for report in reports for r in report['accepted']])
        result={k:v for k,v in result.items() if k!='reports'}
    else:
        from .extraction import load_drugprot,evaluate
        result=evaluate(args.gold,args.predictions,load_drugprot(args.abstracts,args.entities));atomic_json(args.output,result)
    print(json.dumps(result,indent=2,ensure_ascii=False))
    if isinstance(result,dict) and result.get('errors'):raise SystemExit(1)


if __name__=='__main__':main()
