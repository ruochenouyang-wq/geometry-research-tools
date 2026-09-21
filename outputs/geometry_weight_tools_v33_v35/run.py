"""Structured geometry goals and certified matrix-envelope optimization."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
import compiler_v23 as compiler
import reductions
import family_v29 as family
import matrix_v30 as matrix
import transport_v31 as transport
import decision_v28 as decision
import planner_v32 as baseline
import planner_v35 as planner
import dual_v33 as dual
import optimize_v34 as optimizer


def verify(c):
    if not isinstance(c,dict):return False
    method=c.get('method')
    if method=='expression_compiler_v23':return compiler.verify(c)
    if 'rule' in c:return reductions.verify(c)
    if method in ('poisson_scalar_v29','poisson_matrix_v30'):return family.verify(c)
    if method=='lemma_transport_v31':return transport.verify(c)
    if method=='threshold_decision_v28':return decision.verify(c)
    if method=='proof_search_v32':return baseline.verify(c)
    if method=='matrix_dual_lower_v33':return dual.verify_lower(c)
    if method=='global_matrix_envelope_v33':return dual.verify(c)
    if method=='penalty_budget_v35':return planner.verify_budget(c)
    return False


def load_library(raw):
    if isinstance(raw,list):
        if len(raw)>32 or any(not family.verify(t) for t in raw):raise ValueError('Invalid lemma list')
        return raw
    if not isinstance(raw,dict):raise ValueError('Invalid library input')
    cert=raw.get('certificate',raw)
    if not verify(cert):raise ValueError('Invalid library certificate')
    if cert['method']=='global_matrix_envelope_v33':return [cert['template']]
    if cert['method'] in ('poisson_scalar_v29','poisson_matrix_v30'):return [cert]
    raise ValueError('Expected a family lemma or a checked envelope optimization result')


def main():
    parser=argparse.ArgumentParser();sub=parser.add_subparsers(dest='command',required=True)
    s=sub.add_parser('solve');s.add_argument('goal',type=Path);s.add_argument('--library',type=Path);s.add_argument('--no-synthesis',action='store_true');s.add_argument('--no-penalty-weights',action='store_true');s.add_argument('--max-leaves',type=int,default=64);s.add_argument('--output',type=Path,default=Path('result.json'))
    o=sub.add_parser('optimize');o.add_argument('spec',type=Path);o.add_argument('--tolerance',default='1/100000');o.add_argument('--max-rounds',type=int,default=8);o.add_argument('--max-range-leaves',type=int,default=128);o.add_argument('--output',type=Path,default=Path('weights.json'))
    l=sub.add_parser('lemma');l.add_argument('spec',type=Path);l.add_argument('--output',type=Path,default=Path('lemma.json'))
    c=sub.add_parser('compile');c.add_argument('goal',type=Path);c.add_argument('--output',type=Path,default=Path('compiled.json'))
    v=sub.add_parser('verify');v.add_argument('paths',type=Path,nargs='+')
    args=parser.parse_args()
    if args.command=='verify':
        checked=0;failed=[]
        for path in args.paths:
            paths=sorted(path.glob('*.json')) if path.is_dir() else [path]
            if not paths:failed.append(str(path))
            for file in paths:
                checked+=1
                try:
                    raw=json.loads(file.read_text());valid=verify(raw.get('certificate',raw) if isinstance(raw,dict) else raw)
                    if isinstance(raw,dict) and 'diagnostics' in raw:valid=valid and all(planner.verify_budget(c) for c in raw['diagnostics'])
                except (ValueError,TypeError,KeyError,OSError):valid=False
                if not valid:failed.append(str(file))
        print(json.dumps({'checked':checked,'failed':failed}))
        if failed:raise SystemExit(1)
        return
    if args.command=='solve':
        raw=json.loads(args.goal.read_text());library=load_library(json.loads(args.library.read_text())) if args.library else []
        result=planner.solve(raw,library=library,synthesize=not args.no_synthesis,use_penalty_weights=not args.no_penalty_weights,max_leaves=args.max_leaves)
    elif args.command=='compile':result={'certificate':compiler.compile_goal(json.loads(args.goal.read_text()))}
    elif args.command=='optimize':
        spec=json.loads(args.spec.read_text());result=optimizer.solve(spec['directions'],spec.get('cost'),F(args.tolerance),args.max_rounds,args.max_range_leaves)
    else:
        spec=json.loads(args.spec.read_text());result={'certificate':matrix.synthesize(spec['directions'],spec.get('weight'))}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
    c=result['certificate'];print(json.dumps({'output':str(args.output),'method':c['method'],'status':c.get('status'),'statistics':result.get('statistics')},ensure_ascii=False))

if __name__=='__main__':
    try:main()
    except (ValueError,TypeError,KeyError,SyntaxError,ZeroDivisionError) as exc:
        print(json.dumps({'status':'unsupported_or_invalid_input','reason':str(exc)},ensure_ascii=False));raise SystemExit(2)
