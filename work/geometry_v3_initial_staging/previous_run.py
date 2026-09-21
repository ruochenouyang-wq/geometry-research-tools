#!/usr/bin/env python3
"""Run a certified ground-state computation, with an explicit error target."""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys
from time import perf_counter
import error_bounds as e


def search(q, method='auto', tolerance=F(1,10**10), max_modes=24):
    q = e.base.potential([str(x) for x in q])
    tolerance = F(tolerance)
    if not 0 < tolerance <= 1:
        raise ValueError('Tolerance must be in (0,1]')
    e.base.check_sizes(1,max_modes)
    if max_modes < 4: raise ValueError('Adaptive search needs max_modes >= 4')
    if method not in ('auto','temple','symmetry','exponential','banded','baseline'):
        raise ValueError('Unknown method')
    start, attempts, best, best_size, latest = perf_counter(), [], None, None, None
    def consider(name,size,fn):
        nonlocal best,best_size,latest
        t = perf_counter()
        latest = None
        try:
            cert = fn()
            if not e.verify(cert): raise ArithmeticError('Certificate rejected')
            width = F(cert['exact_width'])
            latest = cert
            attempts.append({'method':name,'size':size,'seconds':perf_counter()-t,
                             'exact_width':str(width),'status':'certified'})
            if best is None or width < F(best['exact_width']): best,best_size = cert,size
            return width <= tolerance
        except ValueError as exc:
            attempts.append({'method':name,'size':size,'seconds':perf_counter()-t,
                             'status':'not_certified','reason':str(exc)})
            return False
    def result(met):
        return {'status':'target_met' if met else 'target_not_met','tolerance':str(tolerance),
                'requested_method':method,'size':best_size,'seconds':perf_counter()-start,
                'attempts':attempts,'certificate':best}
    if method == 'auto':
        # Cheap structural proposal, useful when q is close to the exactly
        # representable exponential family. It is always certified afresh.
        slope = -q[1]/2 if len(q)>1 else F(0)
        if consider('exponential_structural',1,lambda:e.barta_certificate(q,[0,slope])):
            return result(True)
    if method == 'exponential':
        for degree in [1,2,4,6,8,10,12,14,16]:
            if consider('exponential',degree,lambda:e.certify_exponential(q,degree)):
                return result(True)
        return result(False)
    sizes = list(range(4,max_modes+1,2))
    if sizes[-1] != max_modes: sizes.append(max_modes)
    for n in sizes:
        if method in ('auto','temple','symmetry'):
            symmetry = method in ('auto','symmetry')
            label = 'symmetry_temple' if symmetry else 'temple'
            if consider(label,n,lambda:e.certify_temple(q,n,symmetry=symmetry)): return result(True)
            if method in ('temple','symmetry'): continue
            # Both the missing gap and an amplified residual are meaningful
            # reasons to try a certificate without a neighbouring-eigenvalue
            # assumption. Fallback is charged to total search time.
            near_cluster = latest is None or (F(latest['gap_lower'])-F(latest['rayleigh_quotient'])
                              < F(1,1000)*max(F(1),abs(F(latest['rayleigh_quotient']))))
            if near_cluster and consider('banded',n,lambda:e.certify_banded(q,1,n,48)):
                return result(True)
        elif method == 'banded':
            if consider('banded',n,lambda:e.certify_banded(q,1,n,48)): return result(True)
        else:
            if consider('baseline',n,lambda:e.base.certify(q,1,n,48)): return result(True)
    return result(False)


def main(argv=None):
    parser = argparse.ArgumentParser(description='精确误差界：残差、指数试探函数与带状 Schur 方法')
    sub = parser.add_subparsers(dest='command',required=True)
    calc = sub.add_parser('bound')
    calc.add_argument('--q',default='0,1')
    calc.add_argument('--method',choices=['auto','temple','symmetry','exponential','banded','baseline'],default='auto')
    calc.add_argument('--tolerance',default='1e-10')
    calc.add_argument('--max-modes',type=int,default=24)
    calc.add_argument('--out',required=True)
    check = sub.add_parser('verify')
    check.add_argument('file')
    args = parser.parse_args(argv)
    try:
        if args.command == 'verify':
            path = Path(args.file)
            if path.stat().st_size > 2_000_000: raise ValueError('File too large')
            data = json.loads(path.read_text())
            if not isinstance(data,dict): raise ValueError('A certificate object is required')
            passed = e.verify(data.get('certificate',data))
            if 'certificate' in data:
                # A valid enclosure must not masquerade as a reached target.
                if data.get('status') not in ('target_met','target_not_met'): passed = False
                if passed:
                    tolerance = F(data['tolerance'])
                    met = F(data['certificate']['exact_width']) <= tolerance
                    passed = 0 < tolerance <= 1 and met == (data['status']=='target_met')
            print('证书复核通过。' if passed else '证书复核未通过。')
            return 0 if passed else 1
        result = search(e.base.potential(args.q.split(',')),args.method,F(args.tolerance),args.max_modes)
        e.base.write_json(args.out,result)
        print('达到目标精度。' if result['status']=='target_met' else '未达到目标精度；保留已有证据。')
        if result['certificate']:
            print('严格区间：['+', '.join(result['certificate']['decimal_enclosure'])+']')
        print('耗时：%.6f 秒；结果：%s' % (result['seconds'],args.out))
        return 0 if result['status']=='target_met' else 2
    except (ValueError,TypeError,OSError,KeyError) as exc:
        print('错误：'+str(exc),file=sys.stderr)
        return 2


if __name__ == '__main__':
    sys.exit(main())
