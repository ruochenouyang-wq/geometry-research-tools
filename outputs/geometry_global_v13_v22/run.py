"""Command-line entry point. Python 3.9+, standard library only."""
import argparse
import json
from pathlib import Path
from fractions import Fraction as F
import global_core as core
import equality_v19 as equality
import robust_v20 as robust
import constrained_v21 as constrained


def verify(c):
    if not isinstance(c,dict):return False
    method=c.get('method')
    if method=='global_regularized_ground_v11':return core.g.verify(c,True)
    if method=='global_function_bounds_v13_v22':return core.verify(c,True)
    if method=='affine_equality_global_v19':return equality.verify(c,True)
    if method=='robust_spectral_dual_v20':return robust.verify(c,True)
    if method in ('constrained_rayleigh_dual_v21','moment_infeasible_v21'):return constrained.verify(c,True)
    return False


def demo(version):
    if version==19:
        p=core.problem([0],[[0,1],[0,0,1],[0,0,0,1]],[[-1,1]]*3,[0]*3,[[int(i==j) for j in range(3)] for i in range(3)])
        return equality.search(p,[[1,0,1]],[0])
    if version==20:
        p=core.problem([0],[[0,1]],[[-1,1]],[0],[[1]])
        return robust.search(p,[[0,1]],[[-1,2]])
    if version==21:return constrained.search([0],[0,1],F(1,2))
    if version==16:
        p=core.all_real_problem([0,8],[[0,1],[0,0,1]],[0,0],[[F(1,5),0],[0,1]])
        return core.search(p,F(1,10000),version=16,domain='all_real')
    if version in (17,18):
        p=core.problem([0],[[0,1],[0,0,1],[0,0,0,1]],[[-1,1]]*3,[0]*3,[[2*int(i==j) for j in range(3)] for i in range(3)])
        return core.search(p,F(1,10),version=version)
    p=core.problem([0],[[0,1]],[[-6,6]],[F(1,50)],[[F(1,5)]])
    return core.search(p,version=version)


def main():
    parser=argparse.ArgumentParser(description='Certified geometry research tools V13--V22')
    sub=parser.add_subparsers(dest='command',required=True)
    d=sub.add_parser('demo');d.add_argument('--version',type=int,choices=range(13,23),default=22);d.add_argument('--output',type=Path,default=Path('result.json'))
    v=sub.add_parser('verify');v.add_argument('paths',type=Path,nargs='+')
    args=parser.parse_args()
    if args.command=='demo':
        result=demo(args.version);args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2)+'\n')
        c=result['certificate'];print(json.dumps({'version':args.version,'status':c['status'],'gap':c.get('global_gap'),'output':str(args.output)},ensure_ascii=False));return
    failed=[];count=0
    for path in args.paths:
        paths=sorted(path.glob('*.json')) if path.is_dir() else [path]
        if not paths:failed.append(str(path))
        for item in paths:
            try:
                raw=json.loads(item.read_text());valid=verify(raw.get('certificate',raw) if isinstance(raw,dict) else raw)
            except (ValueError,KeyError,TypeError,OSError):valid=False
            count+=1
            if not valid:failed.append(str(item))
    print(json.dumps({'checked':count,'failed':failed},ensure_ascii=False))
    if failed:raise SystemExit(1)

if __name__=='__main__':main()
