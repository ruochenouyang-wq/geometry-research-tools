"""Extend a box certificate to ALL real parameters using growth at infinity."""
from fractions import Fraction as F
from math import ceil
from time import perf_counter
import global_search as g

METHOD='all_real_regularized_ground_v11'


def exterior_data(p,q0_range,direction_ranges):
    d=len(p['directions']);cost=p['penalty'];hessian=[g.rationals(row) for row in cost['hessian']]
    h=min(hessian[i][i]-sum((abs(hessian[i][j]) for j in range(d) if j!=i),F(0)) for i in range(d))
    if h<=0:raise ValueError('This exterior proof requires a positive Gershgorin Hessian lower bound')
    qlo,_=g.v3.range_bounds(g.rationals(p['q0']),q0_range)
    ranges=[g.v3.range_bounds(g.rationals(direction),proof) for direction,proof in zip(p['directions'],direction_ranges)]
    if len(ranges)!=d:raise ValueError('Direction proofs missing')
    linear=g.rationals(cost['linear']);b=[abs(linear[i]+(lo+hi)/2)+(hi-lo)/2 for i,(lo,hi) in enumerate(ranges)]
    return h,b,qlo+F(cost['constant'])


def exterior_lower(p,h,b,offset):
    radii=[]
    for i,pair in enumerate(p['box']):
        lo,hi=g.rationals(pair)
        if lo!=-hi or hi<b[i]/h:raise ValueError('Exterior box must be symmetric and beyond quadratic turning points')
        radii.append(hi)
    bounds=[offset+h*r*r/2-b[j]*r-sum((b[i]**2/(2*h) for i in range(len(b)) if i!=j),F(0))
            for j,r in enumerate(radii)]
    return min(bounds)


def certificate(boxed,q0_range,independent=True):
    if not g.verify(boxed,independent):raise ValueError('Invalid bounded-domain proof')
    p=boxed['problem'];h,b,offset=exterior_data(p,q0_range,boxed['direction_range_proofs'])
    outside=exterior_lower(p,h,b,offset);upper=F(boxed['candidate_upper'])
    if outside<=upper:raise ValueError('A better solution outside the box has not been excluded')
    return {'method':METHOD,'model':g.v3.base.MODEL,'scope':g.v3.base.SCOPE,
            'optimization_scope':'all real parameter vectors and all nonzero axisymmetric form-domain functions',
            'boxed_certificate':boxed,'q0_range_proof':q0_range,'hessian_lower':str(h),
            'absolute_linear_bounds':[str(x) for x in b],'exterior_lower':str(outside),
            'global_lower':boxed['global_lower'],'candidate_upper':boxed['candidate_upper'],
            'global_gap':boxed['global_gap'],'tolerance':boxed['tolerance'],'status':boxed['status'],
            'candidate_parameters':boxed['candidate_parameters'],
            'possible_optimizer_boxes':boxed['possible_optimizer_boxes'],'formal_assistant_checked':False}


def verify(cert,independent=True):
    try:return cert==certificate(cert['boxed_certificate'],cert['q0_range_proof'],independent)
    except (KeyError,TypeError,ValueError,IndexError,ZeroDivisionError,RecursionError):return False


def search(q0,directions,linear,hessian,constant=0,tolerance=F(1,1000),modes=12,max_leaves=256,policy='concavity'):
    start=perf_counter();d=len(directions)
    p=g.problem(q0,directions,[[-1,1]]*d,linear,hessian,constant)
    q0_range=g.v5.make_range(g.rationals(p['q0']))
    ranges=[g.v5.make_range(g.rationals(direction)) for direction in p['directions']]
    h,b,offset=exterior_data(p,q0_range,ranges)
    trial=g.v3.old.polynomial_trial(g.rationals(p['q0']),modes)
    mu,_,_=g.v3.old.residual_statistics(g.rationals(p['q0']),trial)
    incumbent=mu+F(p['penalty']['constant'])
    radius=max(1,ceil(max(b)/h))
    while True:
        p['box']=[[str(-radius),str(radius)] for _ in range(d)]
        if exterior_lower(p,h,b,offset)>incumbent:break
        radius+=1
        if radius>1000:raise ValueError('Required exterior radius exceeds this prototype budget')
    p=g.normalized(p)
    result=g.search(p,tolerance,modes,max_leaves,policy,seeds=[[0]*d])
    cert=certificate(result['certificate'],q0_range,independent=True)
    return {'algorithm_version':11,'certificate':cert,'stop_reason':result['stop_reason'],
            'modes':modes,'point_evaluations':result['point_evaluations'],'history':result['history'],
            'seconds':perf_counter()-start}
