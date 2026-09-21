"""V5: certified extremal range of a polynomial on [-1,1]."""
from fractions import Fraction as F
from math import comb
from time import perf_counter
import error_bounds as e


def quadratic(q):
    if len(q)>3:raise ValueError('Quadratic range requires degree <= 2')
    points=[F(-1),F(1)]
    if len(q)==3 and q[2]:
        t=-q[1]/(2*q[2])
        if -1<t<1:points.append(t)
    values=[e.evaluate(q,t) for t in points]
    return {'kind':'quadratic_range_v5','lower':str(min(values)),'upper':str(max(values))}


def interval_coefficients(q,a,b):
    n=len(q)-1
    power=[sum((q[j]*comb(j,k)*a**(j-k)*(b-a)**k for j in range(k,n+1)),F(0)) for k in range(n+1)]
    return [sum((power[j]*F(comb(i,j),comb(n,j)) for j in range(i+1)),F(0)) for i in range(n+1)]


def bounds_from_partition(q,intervals):
    if not isinstance(intervals,list) or not 1<=len(intervals)<=256:raise ValueError('Invalid partition')
    lower,upper=None,None;last=F(-1)
    for raw in intervals:
        if not isinstance(raw,list) or len(raw)!=2:raise ValueError('Invalid interval')
        a,b=map(e.base.rational,raw)
        if a!=last or not a<b<=1:raise ValueError('Partition must cover [-1,1] in order with no gaps')
        coeff=interval_coefficients(q,a,b)
        lower=min(coeff) if lower is None else min(lower,min(coeff))
        upper=max(coeff) if upper is None else max(upper,max(coeff))
        last=b
    if last!=1:raise ValueError('Partition misses the right endpoint')
    return lower,upper


def verify_bounds(q,proof):
    if proof.get('kind')=='quadratic_range_v5':
        if proof!=quadratic(q):raise ValueError('Incorrect quadratic extrema')
        return F(proof['lower']),F(proof['upper'])
    if proof.get('kind')!='adaptive_bernstein_v5':raise ValueError('Unknown polynomial range proof')
    lo,hi=bounds_from_partition(q,proof['intervals'])
    epsilon=e.base.rational(proof['requested_tolerance'])
    if not F(1,10**12)<=epsilon<=1:raise ValueError('Invalid range tolerance')
    witness=[]
    for key in ['minimum_witness','maximum_witness']:
        t=e.base.rational(proof[key]['at'])
        if not -1<=t<=1:raise ValueError('Witness is outside the domain')
        val=e.evaluate(q,t)
        if proof[key]!={'at':str(t),'value':str(val)}:raise ValueError('Incorrect witness value')
        witness.append(val)
    met=witness[0]-lo<=epsilon and hi-witness[1]<=epsilon
    expected={'kind':'adaptive_bernstein_v5','intervals':proof['intervals'],
              'lower':str(lo),'upper':str(hi),'requested_tolerance':str(epsilon),
              'minimum_witness':proof['minimum_witness'],'maximum_witness':proof['maximum_witness'],
              'range_target_met':met}
    if proof!=expected:raise ValueError('Incorrect range enclosure or target claim')
    return lo,hi


def make_range(q,epsilon=F(1,10**6),max_leaves=64):
    q=e.base.potential([str(x) for x in q]);epsilon=F(epsilon)
    if not F(1,10**12)<=epsilon<=1:raise ValueError('Range tolerance must be in [1e-12,1]')
    if type(max_leaves) is not int or not 1<=max_leaves<=256:raise ValueError('Invalid leaf budget')
    if len(q)<=3:return quadratic(q)
    leaves=[(F(-1),F(1),e.bernstein_coefficients(q))]
    samples={t:e.evaluate(q,t) for t in [F(-1),F(0),F(1)]}
    while True:
        tmin=min(samples,key=samples.get);tmax=max(samples,key=samples.get)
        low=min(min(c) for _,_,c in leaves);high=max(max(c) for _,_,c in leaves)
        met=samples[tmin]-low<=epsilon and high-samples[tmax]<=epsilon
        if met or len(leaves)>=max_leaves:break
        index=max(range(len(leaves)),key=lambda i:max(samples[tmin]-min(leaves[i][2]),
                                                    max(leaves[i][2])-samples[tmax]))
        a,b,c=leaves.pop(index);mid=(a+b)/2
        left,right=e.split_bernstein(c)
        leaves.extend([(a,mid,left),(mid,b,right)])
        for t in [(a+mid)/2,mid,(mid+b)/2]:samples[t]=e.evaluate(q,t)
    intervals=[[str(a),str(b)] for a,b,_ in sorted(leaves)]
    proof={'kind':'adaptive_bernstein_v5','intervals':intervals,'lower':str(low),'upper':str(high),
           'requested_tolerance':str(epsilon),'minimum_witness':{'at':str(tmin),'value':str(samples[tmin])},
           'maximum_witness':{'at':str(tmax),'value':str(samples[tmax])},'range_target_met':met}
    verify_bounds(q,proof)
    return proof


def search(q,k=1,tolerance=F(1,10**10),max_modes=32,epsilon=F(1,10**6),max_leaves=64):
    import v04_search as v4
    start=perf_counter();proof=make_range(q,epsilon,max_leaves)
    result=v4.search(q,k,tolerance,max_modes,range_proof=proof)
    result['algorithm_version']=5;result['seconds']=perf_counter()-start
    return result
