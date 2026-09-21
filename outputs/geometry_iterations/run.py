#!/usr/bin/env python3
"""Unified command-line interface for the eight research tool iterations."""
import argparse
import json
from pathlib import Path
from fractions import Fraction as F
import sys
import v03_tail as v3
import v04_search as v4
import v05_range as v5
import v06_parity as v6
import v07_precision as v7
import v08_cluster as v8
import v09_family as v9
import v10_adaptive as v10
from check_document import verify_document


def main(argv=None):
    parser=argparse.ArgumentParser(description='几何谱研究工具：有理证书、误差估计与预算搜索')
    commands=parser.add_subparsers(dest='command',required=True)
    p=commands.add_parser('verify',help='独立复核保存的证书和精度声明');p.add_argument('path',type=Path)
    for command in ['spectrum','range','refine','cluster','family','function']:
        p=commands.add_parser(command);p.add_argument('--out',type=Path,required=True)
        if command!='family':p.add_argument('--q',required=True,help='从常数项起的有理系数，例如 1,-2,3')
        if command in ('spectrum','function'):p.add_argument('--tolerance',type=F,default=F(1,10**10))
        if command in ('spectrum','function'):p.add_argument('--max-modes',type=int,default=32)
        if command in ('spectrum','family'):p.add_argument('--index',type=int,default=1)
        if command in ('spectrum','refine','cluster','family','function'):p.add_argument('--modes',type=int,default=12)
        if command=='spectrum':
            p.add_argument('--version',type=int,choices=[3,4,5,6,10],default=10)
            p.add_argument('--bits',type=int,default=112)
        if command=='range':
            p.add_argument('--epsilon',type=F,default=F(1,10**6));p.add_argument('--max-leaves',type=int,default=64)
        if command in ('refine','function'):p.add_argument('--bits',type=int,default=160)
        if command=='refine':p.add_argument('--steps',type=int,default=1)
        if command=='cluster':p.add_argument('--size',type=int,default=2)
        if command=='family':
            p.add_argument('--q0',required=True);p.add_argument('--direction',required=True)
            p.add_argument('--interval',default='-1,1');p.add_argument('--threshold',type=F,required=True)
            p.add_argument('--max-cells',type=int,default=64);p.add_argument('--uniform',action='store_true')
    a=parser.parse_args(argv)
    try:
        if a.command=='verify':
            data=json.loads(a.path.read_text());valid=verify_document(data)
            print(json.dumps({'valid':valid,'independent_dense_checks':True},ensure_ascii=False));return 0 if valid else 1
        if a.command!='family':q=v3.base.potential(a.q.split(','))
        if a.command=='spectrum':
            if a.version==3:
                c=v3.certify(q,a.index,a.modes,a.bits)
                result={'algorithm_version':3,'certificate':c,'tolerance':str(a.tolerance),
                        'status':'target_met' if F(c['exact_width'])<=a.tolerance else 'target_not_met'}
            else:
                method={4:v4.search,5:v5.search,6:v6.search,10:v10.adaptive}[a.version]
                result=method(q,a.index,a.tolerance,a.max_modes)
        elif a.command=='range':result={'method':'polynomial_range_document_v5','q_coefficients':[str(x) for x in q],
                                       'range_proof':v5.make_range(q,a.epsilon,a.max_leaves)}
        elif a.command=='refine':result=v7.refine(q,a.modes,a.bits,a.steps)
        elif a.command=='cluster':result=v8.cluster(q,a.size,a.modes)
        elif a.command=='function':result=v10.function(q,a.tolerance,a.modes,a.max_modes,a.bits)
        else:result=v9.family(v3.base.potential(a.q0.split(',')),v3.base.potential(a.direction.split(',')),
                             a.interval.split(','),a.index,a.threshold,a.modes,a.max_cells,not a.uniform)
        cert=result if 'method' in result else result['certificate']
        if cert is not None and not verify_document(result):raise ArithmeticError('统一复核未通过')
        v3.base.write_json(a.out,result)
        status=result.get('status',cert.get('status','certified') if cert else 'target_not_met')
        if a.command=='range' and result['range_proof'].get('range_target_met') is False:status='target_not_met'
        brief={'status':status,'output':str(a.out.resolve())}
        if cert is not None:
            for key in ['lower','upper','exact_width','projector_hilbert_schmidt_squared_upper','nontrivial_bound','uniform_lower']:
                if key in cert:brief[key]=cert[key]
        print(json.dumps(brief,ensure_ascii=False,indent=2))
        return 2 if status in ('target_not_met','unresolved') else 0
    except (OSError,ValueError,TypeError,ArithmeticError,KeyError) as exc:
        print(str(exc),file=sys.stderr);return 1


if __name__=='__main__':raise SystemExit(main())
