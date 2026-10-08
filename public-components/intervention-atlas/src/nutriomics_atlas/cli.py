from __future__ import annotations
import argparse
import json
from pathlib import Path

def main(argv=None):
    parser=argparse.ArgumentParser(description='Source-verified intervention omics atlas')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('fetch-metabolomics');p.add_argument('--output-dir',type=Path,required=True)
    for name in ('audit-metabolomics','analyze-metabolomics'):
        p=sub.add_parser(name);p.add_argument('--input-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
        if name=='analyze-metabolomics':p.add_argument('--min-pairs',type=int,default=16)
    p=sub.add_parser('prepare-geo');p.add_argument('--accession',choices=['GSE56960','GSE27385'],required=True);p.add_argument('--family-soft',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    p=sub.add_parser('audit-counts');p.add_argument('--data-dir',type=Path,required=True);p.add_argument('--output-dir',type=Path,required=True)
    p=sub.add_parser('run-geo');p.add_argument('--accession',choices=['GSE127530','GSE56960','GSE27385'],required=True);p.add_argument('--data-dir',type=Path,required=True);p.add_argument('--expansion-root',type=Path);p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--rscript',type=Path,required=True);p.add_argument('--r-library',nargs='+',type=Path,required=True);p.add_argument('--expected-r-version',default='4.5.3');p.add_argument('--preflight-only',action='store_true')
    p=sub.add_parser('validate');p.add_argument('--input-dir',type=Path,required=True)
    p=sub.add_parser('audit-inventory');p.add_argument('--source-roots',nargs='+',type=Path,required=True);p.add_argument('--analysis-dirs',nargs='*',type=Path,default=[]);p.add_argument('--output',type=Path,required=True)
    args=parser.parse_args(argv)
    if args.command=='audit-inventory':
        from .inventory import audit_inventory
        result=audit_inventory(args.source_roots,args.analysis_dirs,args.output)
    elif args.command=='fetch-metabolomics':
        from .download import fetch_st001257
        result=fetch_st001257(args.output_dir)
    elif args.command in {'audit-metabolomics','analyze-metabolomics'}:
        from .metabolomics import audit_matrix,analyze
        result=audit_matrix(args.input_dir,args.output_dir)[2] if args.command=='audit-metabolomics' else analyze(args.input_dir,args.output_dir,args.min_pairs)
    elif args.command=='prepare-geo':
        from .geo import prepare_design
        result=prepare_design(args.accession,args.family_soft,args.output_dir)
    elif args.command=='audit-counts':
        from .counts import audit
        result=audit(args.data_dir,args.output_dir)
    elif args.command=='run-geo':
        from .runner import run_geo
        result=run_geo(args.accession,args.data_dir,args.expansion_root,args.output_dir,args.rscript,args.r_library,args.expected_r_version,args.preflight_only)
    else:
        from .validation import validate_output
        result=validate_output(args.input_dir)
    print(json.dumps(result,ensure_ascii=False,indent=2))
    if result.get('errors'):raise SystemExit(1)
    if args.command=='validate' and result.get('full_validation_passed') is not True:raise SystemExit(2)

if __name__=='__main__':main()
