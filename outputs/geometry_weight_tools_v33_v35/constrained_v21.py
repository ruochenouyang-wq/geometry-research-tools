"""V21: constrained Rayleigh minimization using certified spectral duality.

min R_q(u), with <u,p u>/<u,u> >= target, over the full form domain.
No strong-duality or convergence assumption is used in certificate acceptance.
"""
from fractions import Fraction as F
from time import perf_counter
import global_core as core
import point_oracle


def expectation(p,c):
    return core.v3.old.residual_statistics(p,c)[0]-core.v3.old.residual_statistics([F(0)],c)[0]


def certificate(q,p,target,beta,spectral,trial,tolerance,independent=True):
    q=core.v3.base.potential(list(map(str,q)));p=core.v3.base.potential(list(map(str,p)))
    target=F(target);beta=F(beta);tolerance=F(tolerance)
    if beta<0 or not F(1,10**12)<=tolerance<=100:raise ValueError('Invalid multiplier or tolerance')
    tilted=core.v3.old.add(q,core.v3.old.scale(p,-beta))
    if spectral['eigenvalue_index']!=1 or core.v3.base.potential(spectral['q_coefficients'])!=tilted or not core.v3.verify(spectral,independent):raise ValueError('Incorrect tilted spectrum')
    lower=F(spectral['lower'])+beta*target;upper=None;actual=None
    if trial is not None:
        c=core.g.v7.coefficients(list(map(str,trial)));actual=expectation(p,c)
        if actual<target:raise ValueError('Trial violates the moment constraint')
        upper=core.v3.old.residual_statistics(q,c)[0];trial=list(map(str,c))
        if upper<lower:raise ValueError('Inconsistent constrained bounds')
    return {'method':'constrained_rayleigh_dual_v21','q':list(map(str,q)),'moment_polynomial':list(map(str,p)),
            'target':str(target),'multiplier':str(beta),'spectral_certificate':spectral,'trial_legendre':trial,
            'candidate_moment':None if actual is None else str(actual),'global_lower':str(lower),
            'candidate_upper':None if upper is None else str(upper),'global_gap':None if upper is None else str(upper-lower),
            'tolerance':str(tolerance),'status':'feasibility_not_found' if upper is None else ('epsilon_global' if upper-lower<=tolerance else 'global_gap_open'),
            'optimization_scope':'all nonzero axisymmetric form-domain functions satisfying the expectation inequality',
            'formal_assistant_checked':False}


def infeasible(q,p,target,proof):
    q=core.v3.base.potential(list(map(str,q)));p=core.v3.base.potential(list(map(str,p)))
    _,upper=core.v3.range_bounds(p,proof);target=F(target)
    if target<=upper:raise ValueError('Infeasibility not proved')
    return {'method':'moment_infeasible_v21','q':list(map(str,q)),'moment_polynomial':list(map(str,p)),
            'target':str(target),'range_proof':proof,'status':'certified_infeasible'}


def verify(c,independent=True):
    try:
        if c['method']=='moment_infeasible_v21':return c==infeasible(c['q'],c['moment_polynomial'],c['target'],c['range_proof'])
        return c==certificate(c['q'],c['moment_polynomial'],c['target'],c['multiplier'],c['spectral_certificate'],c['trial_legendre'],c['tolerance'],independent)
    except (ValueError,KeyError,TypeError,IndexError,ZeroDivisionError):return False


def blend_feasible(p,target,left,right,steps=32):
    n=max(len(left),len(right));left=list(left)+[F(0)]*(n-len(left));right=list(right)+[F(0)]*(n-len(right))
    if sum(a*b*F(2,2*j+1) for j,(a,b) in enumerate(zip(left,right)))<0:right=[-v for v in right]
    lo=F(0);hi=F(1);best=right
    if expectation(p,left)>=target or expectation(p,right)<target:raise ValueError('Feasibility bracket required')
    for _ in range(steps):
        mid=(lo+hi)/2;c=[(1-mid)*a+mid*b for a,b in zip(left,right)]
        if not any(c):lo=mid;continue
        if expectation(p,c)>=target:hi=mid;best=c
        else:lo=mid
    return best


def search(q,p,target,tolerance=F(1,1000),max_steps=32,max_modes=24):
    start=perf_counter();q=core.v3.base.potential(list(map(str,q)));p=core.v3.base.potential(list(map(str,p)));target=F(target);tolerance=F(tolerance)
    if type(max_steps) is not int or not 1<=max_steps<=64:raise ValueError('Invalid step budget')
    if not F(1,10**12)<=tolerance<=100:raise ValueError('Invalid tolerance')
    proof=core.g.v5.make_range(p)
    if target>core.v3.range_bounds(p,proof)[1]:
        return {'algorithm_version':21,'certificate':infeasible(q,p,target,proof),'steps':0,'seconds':perf_counter()-start}
    model=core.problem(q,[[-v for v in p]],[[0,1]])
    oracle=point_oracle.Oracle(model,min(F(1),tolerance/32),max_modes=max_modes)
    lo=F(0);hi=None;left=None;right=None;beta=F(0);bestlow=None;bestup=None;history=[]
    for step in range(max_steps):
        idx=oracle.point([beta]);r=oracle.records[idx];c=list(map(F,r['trial_legendre']));m=expectation(p,c)
        low=F(r['spectral_certificate']['lower'])+beta*target
        if bestlow is None or low>bestlow[0]:bestlow=(low,beta,r['spectral_certificate'])
        if m>=target:
            hi=beta;right=c
            candidates=[c]
            if left is not None:candidates.append(blend_feasible(p,target,left,right))
            for candidate in candidates:
                up=core.v3.old.residual_statistics(q,candidate)[0]
                if bestup is None or up<bestup[0]:bestup=(up,candidate)
        else:lo=beta;left=c
        history.append({'step':step,'beta':str(beta),'lower':str(bestlow[0]),'upper':None if bestup is None else str(bestup[0])})
        if bestup is not None and bestup[0]-bestlow[0]<=tolerance:break
        proposal=(lo+hi)/2 if hi is not None else (2*lo if lo else F(1))
        proposal=F(round(proposal*2**24),2**24)
        # A bounded search may remain unresolved; it never claims infeasibility.
        if proposal==beta or proposal>128:break
        beta=proposal
    _,beta,sc=bestlow
    cert=certificate(q,p,target,beta,sc,None if bestup is None else bestup[1],tolerance)
    return {'algorithm_version':21,'certificate':cert,'steps':len(history),'history':history,'statistics':oracle.statistics(),'seconds':perf_counter()-start}
