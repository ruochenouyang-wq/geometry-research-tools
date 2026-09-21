#!/usr/bin/env python3
"""CLI for V11/V12 global research tools."""
import argparse
from fractions import Fraction as F
from pathlib import Path
import json
import sys
import global_search as g
import coercive
import moment_v12


def verify_document(doc):
    if not isinstance(doc,dict):return False
    cert=doc.get('certificate',doc)
    if not isinstance(cert,dict):return False
    method=cert.get('method')
    if method==g.METHOD:return g.verify(cert)
    if method==coercive.METHOD:return coercive.verify(cert)
    if method==moment_v12.METHOD:return moment_v12.verify(cert)
    return False


def main(argv=None):
    p=argparse.ArgumentParser(description='连续参数与无限维函数的全局最优差距认证')
    sub=p.add_subparsers(dest='command',required=True)
    s=sub.add_parser('optimize');s.add_argument('--problem',type=Path,required=True)
    s.add_argument('--domain',choices=['box','all-real'],default='all-real')
    s.add_argument('--all-real-method',choices=['elimination','growth'],default='elimination')
    s.add_argument('--bound',choices=['concavity','lipschitz'],default='concavity')
    s.add_argument('--epsilon',type=F,default=F(1,1000));s.add_argument('--modes',type=int,default=10)
    s.add_argument('--max-leaves',type=int,default=256);s.add_argument('--out',type=Path,required=True)
    s=sub.add_parser('verify');s.add_argument('path',type=Path)
    a=p.parse_args(argv)
    try:
        if a.command=='verify':
            doc=json.loads(a.path.read_text());ok=verify_document(doc)
            cert=doc.get('certificate',doc) if isinstance(doc,dict) else {}
            print(json.dumps({'valid':ok,'status':cert.get('status') if ok else None},ensure_ascii=False));return 0 if ok else 1
        raw=json.loads(a.problem.read_text());cost=raw['penalty']
        if a.domain=='box':
            problem=g.problem(raw['q0'],raw['directions'],raw['box'],cost['linear'],cost['hessian'],cost.get('constant',0))
            result=g.search(problem,a.epsilon,a.modes,a.max_leaves,a.bound)
        else:
            fn=moment_v12.search if a.all_real_method=='elimination' else coercive.search
            result=fn(raw['q0'],raw['directions'],cost['linear'],cost['hessian'],cost.get('constant',0),
                      a.epsilon,a.modes,a.max_leaves,a.bound)
        if not verify_document(result):raise ArithmeticError('全局证书复核未通过')
        g.v3.base.write_json(a.out,result);c=result['certificate']
        print(json.dumps({'status':c['status'],'parameter_domain':a.domain,
                          'candidate_parameters_approx':[float(F(x)) for x in c['candidate_parameters']],
                          'global_lower_approx':float(F(c['global_lower'])),
                          'candidate_upper_approx':float(F(c['candidate_upper'])),
                          'global_gap_approx':float(F(c['global_gap'])),
                          'exact_values_in':str(a.out.resolve())},ensure_ascii=False,indent=2))
        return 0 if c['status']=='epsilon_global' else 2
    except (OSError,ValueError,KeyError,TypeError,ArithmeticError) as exc:
        print(str(exc),file=sys.stderr);return 1


if __name__=='__main__':raise SystemExit(main())
