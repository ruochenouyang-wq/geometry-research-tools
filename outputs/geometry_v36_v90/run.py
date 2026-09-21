"""Legacy mathematical CLI plus compact request/response tools and JSONL sessions."""
import argparse,json,sys
from pathlib import Path
import legacy_cli
import research_registry
import protocol
import token_meter as meter

verify=research_registry.verify
load_library=legacy_cli.load_library

OPERATIONS={
    'register':['goal'],'revise':['goal_ref','changes'],'solve':['goal_ref','options?'],
    'submit':['goal','options?','token_limit?'],'inspect':['ref','pointer?','offset?','limit?'],
    'batch':['goal_refs','options?','exceptions_only?','deduplicate?'],'add_lemma':['lemma'],
    'lemma_index':[],'poll':['receipt','known_evidence'],'budget':['receipt','limit?','encoding?'],
    'checkpoint':['goals','receipts'],'restore':['ref'],'context':['receipts','pending'],
    'optimize':['spec','tolerance?','max_rounds?','max_range_leaves?','token_limit?'],
    'lemma':['spec','token_limit?'],
    'extrema':['spec','tolerance?','max_nodes?','token_limit?'],
    'positive':['spec','options?','token_limit?'],
    'uniform':['spec','rule?','token_limit?'],
    'refine':['spec','tolerance?','budget?','max_nodes?','token_limit?'],
    'research':['goal','rule?','coefficient?','token_limit?'],
}


def call(service,request):
    if not isinstance(request,dict):raise ValueError('Request object required')
    args=dict(request);op=args.pop('op')
    if op=='capabilities':return {'profile':protocol.PROFILE_TEXT,'operations':OPERATIONS}
    if op not in OPERATIONS:raise ValueError('Unknown operation')
    mappings={'goal_ref':'ref','goal_refs':'refs','goal':'raw','lemma':'raw'}
    args={mappings.get(k,k):v for k,v in args.items()}
    return getattr(service,op)(**args)


def main():
    if len(sys.argv)>1 and sys.argv[1]=='verify':
        parser=argparse.ArgumentParser();parser.add_argument('command');parser.add_argument('paths',type=Path,nargs='+');args=parser.parse_args()
        checked=0;failed=[]
        for path in args.paths:
            files=sorted(path.glob('*.json')) if path.is_dir() else [path]
            if not files:failed.append(str(path))
            for file in files:
                checked+=1
                try:
                    raw=json.loads(file.read_text());protocol.checked(raw)
                except (ValueError,TypeError,KeyError,IndexError,OSError,ArithmeticError,AttributeError,RecursionError):failed.append(str(file))
        print(meter.wire({'checked':checked,'failed':failed}))
        if failed:raise SystemExit(1)
        return
    if len(sys.argv)>1 and sys.argv[1] not in ('call','serve'):
        return legacy_cli.main()
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['call','serve']);parser.add_argument('request',nargs='?',type=Path);parser.add_argument('--store',type=Path,default=Path('.geometry-store'));args=parser.parse_args()
    service=protocol.Service(args.store)
    if args.command=='call':
        if args.request is None:raise ValueError('Request JSON file required')
        print(meter.wire(call(service,json.loads(args.request.read_text()))));return
    for line in sys.stdin:
        if not line.strip():continue
        try:result=call(service,json.loads(line))
        except (ValueError,TypeError,KeyError,IndexError,OSError,RuntimeError,SyntaxError,ArithmeticError,AttributeError) as exc:result={'error':type(exc).__name__,'detail':str(exc),'mathematical_verdict':None}
        print(meter.wire(result),flush=True)


if __name__=='__main__':
    try:main()
    except (ValueError,TypeError,KeyError,IndexError,OSError,RuntimeError,SyntaxError,ArithmeticError,AttributeError) as exc:
        print(meter.wire({'error':type(exc).__name__,'detail':str(exc),'mathematical_verdict':None}));raise SystemExit(2)
