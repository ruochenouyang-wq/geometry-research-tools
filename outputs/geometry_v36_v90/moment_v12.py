"""V12: eliminate convex parameter variables and certify a nonlinear function optimum."""
from fractions import Fraction as F
from time import perf_counter
import global_search as g

METHOD='moment_eliminated_global_v12'


def confinement(p,direction_proofs):
    d=len(p['directions']);cost=p['penalty'];h=[g.rationals(row) for row in cost['hessian']]
    if g.v3.base.inertia(h)!=[0,0,d]:raise ValueError('Parameter elimination requires a positive definite Hessian')
    inverse_columns=[g.v7.solve(h,[F(i==j) for i in range(d)]) for j in range(d)]
    inv=[[inverse_columns[j][i] for j in range(d)] for i in range(d)]
    linear=g.rationals(cost['linear'])
    if len(direction_proofs)!=d:raise ValueError('Direction proofs required')
    ranges=[g.v3.range_bounds(g.rationals(q),proof) for q,proof in zip(p['directions'],direction_proofs)]
    box=[]
    for row in inv:
        endpoints=[[-a*(linear[j]+lo),-a*(linear[j]+hi)] for j,(a,(lo,hi)) in enumerate(zip(row,ranges))]
        box.append([sum((min(pair) for pair in endpoints),F(0)),sum((max(pair) for pair in endpoints),F(0))])
    return inv,ranges,box


def moments(p,c):
    # Difference of exact Rayleigh quotients cancels the same kinetic energy.
    mu0,_,_=g.v3.old.residual_statistics(g.rationals(p['q0']),c)
    kinetic,_,_=g.v3.old.residual_statistics([F(0)],c)
    m=[g.v3.old.residual_statistics(g.rationals(direction),c)[0]-kinetic for direction in p['directions']]
    return mu0,m


def eliminated_value(p,c,inverse):
    mu,m=moments(p,c);linear=g.rationals(p['penalty']['linear']);z=[a+b for a,b in zip(linear,m)]
    parameters=[-sum((a*b for a,b in zip(row,z)),F(0)) for row in inverse]
    energy=mu+F(p['penalty']['constant'])+sum((a*b for a,b in zip(z,parameters)),F(0))/2
    return energy,parameters,m


def certificate(boxed,independent=True):
    if not g.verify(boxed,independent):raise ValueError('Invalid box certificate')
    p=boxed['problem'];inverse,ranges,confined=confinement(p,boxed['direction_range_proofs'])
    box=[g.rationals(pair) for pair in p['box']]
    if any(lo>a or b>hi for (lo,hi),(a,b) in zip(box,confined)):
        raise ValueError('Search box does not contain every conditionally optimal parameter')
    candidates=[]
    for idx,record in enumerate(boxed['points']):
        coeff=g.v7.coefficients(record['trial_legendre']);energy,parameters,m=eliminated_value(p,coeff,inverse)
        if not g.inside(parameters,box):raise ArithmeticError('Moment confinement inconsistency')
        candidates.append((energy,idx,parameters,m))
    upper,index,parameters,m=min(candidates,key=lambda x:x[0]);lower=F(boxed['global_lower'])
    if upper<lower:raise ArithmeticError('Inconsistent global bracket')
    gap=upper-lower;tolerance=F(boxed['tolerance'])
    return {'method':METHOD,'model':g.v3.base.MODEL,'scope':g.v3.base.SCOPE,
            'optimization_scope':'all real parameter vectors and all nonzero axisymmetric form-domain functions',
            'boxed_certificate':boxed,'inverse_penalty_hessian':[[str(x) for x in row] for row in inverse],
            'moment_confined_box':[[str(x) for x in pair] for pair in confined],
            'candidate_trial_index':index,'candidate_parameters':[str(x) for x in parameters],
            'candidate_moments':[str(x) for x in m],
            'global_lower':str(lower),'candidate_upper':str(upper),'global_gap':str(gap),
            'tolerance':str(tolerance),'status':'epsilon_global' if gap<=tolerance else 'global_gap_open',
            'formal_assistant_checked':False}


def verify(cert,independent=True):
    try:return cert==certificate(cert['boxed_certificate'],independent)
    except (KeyError,ValueError,TypeError,IndexError,ZeroDivisionError,RecursionError):return False


def search(q0,directions,linear,hessian,constant=0,tolerance=F(1,1000),modes=12,max_leaves=256,policy='concavity'):
    start=perf_counter();d=len(directions)
    p=g.problem(q0,directions,[[-1,1]]*d,linear,hessian,constant)
    proofs=[g.v5.make_range(g.rationals(q)) for q in p['directions']]
    inverse,ranges,box=confinement(p,proofs)
    # If a moment coordinate is fixed, padding keeps a nondegenerate search box
    # without excluding the exact confined set. This is safe but not optimal.
    box=[[lo if lo<hi else lo-1,hi if lo<hi else hi+1] for lo,hi in box]
    p['box']=[[str(x) for x in pair] for pair in box];p=g.normalized(p)
    r=g.search(p,tolerance,modes,max_leaves,policy)
    cert=certificate(r['certificate'])
    reason='epsilon_global_after_elimination' if cert['status']=='epsilon_global' else r['stop_reason']
    return {'algorithm_version':12,'certificate':cert,'stop_reason':reason,
            'modes':modes,'point_evaluations':r['point_evaluations'],'history':r['history'],
            'seconds':perf_counter()-start}
