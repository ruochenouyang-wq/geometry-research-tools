"""A: complete negative spectrum counts, with optional signed-sum enclosures.

The constant-potential cases use exactly the same original-moment Kernel as
nonconstant potentials.  No finite Ritz spectrum is accepted as a PDE proof.
"""
from copy import deepcopy
from fractions import Fraction as F
import importlib.util
import json
from math import isfinite
from pathlib import Path
import sys
from time import perf_counter

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location(
    'parallel_a_shared_baseline', ROOT.parent/'shared_baseline.py')
shared = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(shared)
import spectral_gap as sg
import power_forms as pf

direct = sg.direct
FORMAT = 'complete_axis_negative_spectrum_v1'
VERSION = 'parallel-A-multispectrum-v1'
QUANTITY = 'sum_of_strictly_negative_eigenvalues_with_real_multiplicity'


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def normalize_task(task):
    if type(task) is not dict or set(task)-{'kind','function','mean_zero','tolerance','budget','id'}:
        raise ValueError('Malformed negative-spectrum task')
    if task.get('kind') != 'negative_spectrum':
        raise ValueError('Expected negative_spectrum')
    q = shared.baseline.normalize_function(task['function'])
    q = direct.normalize(q)
    if q['profile'] == 'abs_power' and F(q['exponent']) <= -F(1,2):
        raise ValueError('Complete original q^2 moments require alpha > -1/2')
    mean_zero = task.get('mean_zero', True)
    if type(mean_zero) is not bool:
        raise ValueError('mean_zero must be boolean')
    tol = shared.baseline.rational(task.get('tolerance','1/1000000'))
    if not F(1,10**30) <= tol <= 1:
        raise ValueError('Sum-width tolerance must lie in 1e-30..1')
    budget = task.get('budget',{})
    if type(budget) is not dict or set(budget)-{'wall_seconds'}:
        raise ValueError('Only wall_seconds is supported')
    wall = budget.get('wall_seconds',10)
    if type(wall) not in (int,float) or not isfinite(wall) or not 0 <= wall <= 86400:
        raise ValueError('Invalid wall_seconds')
    return {'kind':'negative_spectrum','function':q,'mean_zero':mean_zero,
            'tolerance':str(tol),'budget':{'wall_seconds':float(wall)}}


def form_for(q, modes):
    try:
        return direct.form_bound(q,48), None
    except ValueError:
        fc = pf.select_form(q,modes)
        return deepcopy(fc['form']),fc


def checked_form(q, form, source):
    if source is None:
        expected = direct.form_bound(q,48)
    else:
        if not pf.verify_form(source,expected_function=q):
            raise ValueError('Custom original-function form failed verification')
        expected = source['form']
    if canonical(form) != canonical(expected) or F(form['energy_factor']) <= 0:
        raise ValueError('Invalid positive-energy form bound')
    return expected


def kernel(q,m,mean_zero,n,near,form,source):
    k = sg._kernel(q,m,mean_zero,n,near,48,source,
                   pf.verify_form if source is not None else None)
    if canonical(k.form) != canonical(form):
        raise ValueError('Kernel form mismatch')
    return k


def angular_cutoff(form):
    # Include m=0 explicitly, even for a nonnegative potential.
    for cutoff in range(1,14):
        if direct.tail_bound(form,cutoff) >= 0:
            return cutoff
    raise ValueError('Angular exclusion exceeds the bounded m<=12 budget')


def endpoint(k,j,lower,upper):
    direct.integer(j,1,k.n,'eigenvalue_index')
    low,high = direct.f.rational(lower),direct.f.rational(upper)
    if low > high or high > 0:
        raise ValueError('Negative-eigenvalue endpoints must satisfy lower<=upper<=0')
    lc = direct.base.inertia(k.matrix(low,'lower'))
    uc = direct.base.inertia(k.matrix(high,'upper'))
    if lc[0] > j-1 or uc[0]+uc[1] < j:
        raise ValueError('Ordered eigenvalue endpoint inertia failed')
    return {'index':j,'lower':str(low),'upper':str(high),'width':str(high-low),
            'lower_comparison_inertia':list(lc),'ritz_upper_inertia':list(uc)}


def radial_evidence(k, intervals=None):
    # Strict N(<0): zero directions are deliberately not added here.
    ri = direct.base.inertia(k.matrix(F(0),'upper'))
    li = direct.base.inertia(k.matrix(F(0),'lower'))
    if ri[0] > li[0]:
        raise ArithmeticError('Contradictory strict count enclosure')
    complete = ri[0] == li[0]
    rebuilt = None
    if intervals is not None:
        if not complete or type(intervals) is not list or len(intervals) != ri[0]:
            raise ValueError('Sum intervals require every strictly negative eigenvalue')
        rebuilt=[]
        for j, interval in enumerate(intervals,1):
            if type(interval) is not dict or interval.get('index') != j:
                raise ValueError('Missing, duplicated or reordered eigenvalue index')
            item=endpoint(k,j,interval['lower'],interval['upper'])
            if canonical(item) != canonical(interval):
                raise ValueError('Endpoint evidence does not replay')
            rebuilt.append(item)
    elif complete and ri[0] == 0:
        rebuilt=[]
    return {'azimuth_m':k.m,'multiplicity':1 if k.m==0 else 2,
            'projection':'remove_l0' if k.m==0 and k.mean_zero else 'none',
            'modes':k.n,'near_tail':k.near_tail,'sqrt_bits':48,
            'radial_tail_start':k.tail_start,'radial_tail_lower':str(k.beta),
            'strict_radial_tail_positive':k.beta>0,
            'kernel_evidence':k.evidence(),
            'ritz_at_zero_inertia':list(ri),'lower_at_zero_inertia':list(li),
            'count_lower':ri[0],'count_upper':li[0],
            'count_certified':complete,'eigenvalue_intervals':rebuilt}


def assemble(task,form,source,sectors):
    q,mean_zero=task['function'],task['mean_zero']
    checked_form(q,form,source)
    if type(sectors) is not list or not 1<=len(sectors)<=13:
        raise ValueError('Require a nonempty consecutive azimuth sequence')
    for m,row in enumerate(sectors):
        if row['azimuth_m'] != m or row['multiplicity'] != (1 if m==0 else 2):
            raise ValueError('Sector order or multiplicity mismatch')
    angular=direct.tail_bound(form,len(sectors))
    if angular<0:
        raise ValueError('Uncomputed angular sectors may contain negative spectrum')
    low=sum(row['multiplicity']*row['count_lower'] for row in sectors)
    high=sum(row['multiplicity']*row['count_upper'] for row in sectors)
    exact=low==high
    have_sum=exact and all(row['eigenvalue_intervals'] is not None for row in sectors)
    sum_low=sum_high=width=None
    if have_sum:
        sum_low=sum((row['multiplicity']*sum((F(e['lower']) for e in row['eigenvalue_intervals']),F(0))
                     for row in sectors),F(0))
        sum_high=sum((row['multiplicity']*sum((F(e['upper']) for e in row['eigenvalue_intervals']),F(0))
                      for row in sectors),F(0))
        width=sum_high-sum_low
        if width<0 or sum_high>0:
            raise ArithmeticError('Invalid signed negative-spectrum sum')
    met=have_sum and width<=F(task['tolerance'])
    return {'format':FORMAT,'kind':'negative_spectrum','function':q,
            'mean_zero':mean_zero,'geometry':'unit_S2','measure':'d_sigma/(4*pi)',
            'scope':'all_real_mean_zero_H1' if mean_zero else 'all_real_H1',
            'quantity':QUANTITY,'negative_threshold':'0','threshold_convention':'strictly_less',
            'tolerance':task['tolerance'],'form_bound':deepcopy(form),
            'form_certificate':deepcopy(source),'sectors':deepcopy(sectors),
            'angular_tail_m_start':len(sectors),'angular_tail_lower':str(angular),
            'angular_tail_nonnegative':True,'real_multiplicity_rule':'m0:1;m_positive:2',
            'count_lower':low,'count_upper':high,'count_certified':exact,
            'negative_count':low if exact else None,
            'sum_lower':str(sum_low) if have_sum else None,
            'sum_upper':str(sum_high) if have_sum else None,
            'sum_width':str(width) if have_sum else None,
            'abs_negative_sum_lower':str(-sum_high) if have_sum else None,
            'abs_negative_sum_upper':str(-sum_low) if have_sum else None,
            'sum_certified':have_sum,'sum_target_met':bool(met),
            'target_met':bool(exact and met),
            'full_infinite_space_covered':True,'original_function_moments_exact':True,
            'constant_shortcut_used':False,'formal_proof_assistant_checked':False}


def assessment(cert):
    keys=('count_lower','count_upper','count_certified','negative_count',
          'sum_lower','sum_upper','sum_width','abs_negative_sum_lower','abs_negative_sum_upper',
          'sum_certified','sum_target_met','target_met')
    return {'certificate_valid':True,**{key:cert[key] for key in keys},
            'status':'target_met' if cert['target_met'] else
                     ('count_certified_sum_open' if cert['count_certified'] else 'count_open'),
            'tolerance':cert['tolerance'],'fallback':False}


def verify(cert,task):
    if type(task) is dict and task.get('kind') != 'negative_spectrum':
        try:
            result=shared.verify(cert,task)
            return {**result,'fallback':True,'route':'baseline_delegation'}
        except Exception as exc:
            return {'certificate_valid':False,'target_met':False,'fallback':True,
                    'error':type(exc).__name__+': '+str(exc)}
    try:
        t=normalize_task(task)
        if not sg._json_native(cert) or type(cert) is not dict or cert.get('format')!=FORMAT:
            raise ValueError('Invalid certificate format')
        if cert['function']!=t['function'] or cert['mean_zero'] is not t['mean_zero']:
            raise ValueError('Original function or exact space binding mismatch')
        if cert['tolerance']!=t['tolerance']:
            raise ValueError('Sum tolerance binding mismatch')
        form,source=cert['form_bound'],cert['form_certificate']
        checked_form(t['function'],form,source)
        rows=cert['sectors']
        if type(rows) is not list or not 1<=len(rows)<=13:
            raise ValueError('Invalid sector count')
        rebuilt=[]
        for m,row in enumerate(rows):
            if type(row) is not dict or type(row.get('azimuth_m')) is not int or row['azimuth_m']!=m:
                raise ValueError('Require every consecutive m exactly once')
            k=kernel(t['function'],m,t['mean_zero'],row['modes'],row['near_tail'],form,source)
            rebuilt.append(radial_evidence(k,row['eigenvalue_intervals']))
        fresh=assemble(t,form,source,rebuilt)
        if canonical(fresh)!=canonical(cert):
            raise ValueError('Certificate evidence does not replay')
        return assessment(fresh)
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError,OverflowError,RecursionError) as exc:
        return {'certificate_valid':False,'target_met':False,'count_certified':False,
                'sum_target_met':False,'status':'verification_failed','fallback':False,
                'error':type(exc).__name__+': '+str(exc)}


class BudgetReached(Exception):
    pass


def check_deadline(deadline):
    if perf_counter()>=deadline:
        raise BudgetReached('Cooperative computation budget exhausted')


def eigenvalue_interval(k,j,bits,deadline):
    # Counts at zero already prove there are j strictly negative eigenvalues.
    analytic=direct.tail_bound(k.form,k.start)
    left,right=analytic-1,F(0)
    for _ in range(bits):
        check_deadline(deadline)
        mid=(left+right)/2
        neg,zero,_=k.count(mid,'upper',independent=True)
        if neg+zero>=j:right=mid
        else:left=mid
    upper=right
    distance=F(1)
    for _ in range(64):
        check_deadline(deadline)
        left=analytic-distance
        if k.count(left,'lower',independent=True)[0]<=j-1:break
        distance*=2
    else:raise ValueError('Lower bracket exhausted')
    right=upper
    for _ in range(bits):
        check_deadline(deadline)
        mid=(left+right)/2
        if k.count(mid,'lower',independent=True)[0]<=j-1:left=mid
        else:right=mid
    return endpoint(k,j,left,upper)


def negative_spectrum(task,modes=(8,12,16),near_tail=16,bits=28):
    """Bounded count-first search; all successful and failed stages are retained."""
    started=perf_counter();t=normalize_task(task)
    if type(modes) not in (list,tuple) or not modes:
        raise ValueError('Require a finite mode schedule')
    for n in modes:direct.integer(n,1,32,'modes')
    direct.integer(near_tail,0,32,'near_tail');direct.integer(bits,8,48,'bits')
    # Reserve time for final evidence replay and the caller's own replay.
    limit=t['budget']['wall_seconds'];deadline=started+min(limit*.7,7.0)
    attempts=[];best=None
    def prefer(candidate):
        nonlocal best
        def score(c):
            gap=c['count_upper']-c['count_lower']
            return (gap,0 if c['sum_certified'] else 1,
                    F(c['sum_width']) if c['sum_certified'] else F(10**30))
        if best is None or score(candidate)<score(best):best=candidate
    for n in modes:
        begin=perf_counter();row={'modes':n,'near_tail':near_tail}
        try:
            check_deadline(deadline)
            form,source=form_for(t['function'],n)
            cutoff=angular_cutoff(form);kernels=[];sectors=[]
            for m in range(cutoff):
                check_deadline(deadline)
                k=kernel(t['function'],m,t['mean_zero'],n,near_tail,form,source)
                kernels.append(k);sectors.append(radial_evidence(k))
            cert=assemble(t,form,source,sectors);prefer(cert)
            row.update(count_lower=cert['count_lower'],count_upper=cert['count_upper'],
                       count_certified=cert['count_certified'],angular_sectors=cutoff)
            if cert['count_certified']:
                refined=[]
                for k,sector in zip(kernels,sectors):
                    intervals=[eigenvalue_interval(k,j,bits,deadline)
                               for j in range(1,sector['count_lower']+1)]
                    refined.append(radial_evidence(k,intervals))
                cert=assemble(t,form,source,refined);prefer(cert)
                row.update(sum_width=cert['sum_width'],sum_target_met=cert['sum_target_met'])
            row['status']='target_met' if cert['target_met'] else 'certified_open'
        except (ValueError,TypeError,ArithmeticError,KeyError,BudgetReached) as exc:
            row.update(status='budget_exhausted' if isinstance(exc,BudgetReached) else 'stage_failed',
                       error=type(exc).__name__+': '+str(exc))
        row['elapsed_seconds']=perf_counter()-begin;attempts.append(row)
        if best is not None and best['target_met']:break
        if perf_counter()>=deadline:break
    cert=json.loads(json.dumps(best,allow_nan=False)) if best is not None else None
    checked=verify(cert,t) if cert is not None else {
        'certificate_valid':False,'target_met':False,'count_certified':False,
        'sum_target_met':False,'status':'no_verified_certificate','fallback':False}
    elapsed=perf_counter()-started
    result={'version':VERSION,'route':'complete_moment_multispectrum','fallback':False,
            'certificate':cert,**checked,'attempts':attempts,'total_elapsed_seconds':elapsed,
            'within_budget':elapsed<=limit,'budget_scope':'cooperative; caller enforces hard deadline',
            'actual_model_calls':0,'actual_model_tokens':None}
    if elapsed>limit:result.update(target_met=False,status='budget_exceeded')
    return result


def solve(task):
    if type(task) is dict and task.get('kind') in ('spectrum','approximation'):
        result=shared.solve(task)
        return {**result,'variant':VERSION,'fallback':True,'route':'baseline_delegation'}
    try:
        return negative_spectrum(task)
    except (ValueError,TypeError,KeyError,ArithmeticError) as exc:
        return {'version':VERSION,'certificate':None,'certificate_valid':False,
                'target_met':False,'count_certified':False,'sum_target_met':False,
                'fallback':False,'status':'invalid_or_unsupported_task',
                'error':type(exc).__name__+': '+str(exc),'attempts':[]}


if __name__=='__main__':
    for line in sys.stdin:
        if line.strip():print(json.dumps(solve(json.loads(line)),allow_nan=False),flush=True)
