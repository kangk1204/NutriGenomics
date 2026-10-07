import argparse
import json
from pathlib import Path
from .api import create_app

def main():
    parser=argparse.ArgumentParser(description='Research job service and evidence artifact inspection')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('serve');p.add_argument('--root',type=Path,required=True);p.add_argument('--host',default='127.0.0.1');p.add_argument('--port',type=int,default=8765)
    p=sub.add_parser('inspect');p.add_argument('--root',type=Path,required=True)
    args=parser.parse_args()
    if args.command=='serve':
        import uvicorn
        uvicorn.run(create_app(args.root),host=args.host,port=args.port,workers=1)
    else:
        application=create_app(args.root)
        print(json.dumps(application.state.store.list(),ensure_ascii=False,indent=2))

if __name__=='__main__':main()
