"""V20: continuous uncertainty reduces to vertices; rational primal-dual cuts.

min_a Q(a) - min_theta lambda_1(q0 + a.p + theta.r), all a in R^d.
The uncertainty set is a box. Q has a positive definite Hessian.
"""
from fractions import Fraction as F
from itertools import product
from time import perf_counter
import global_core as core
import moment_v12 as moment
import point_oracle
import simplex_v15 as simplex


def data(p,uncertainty_directions,uncertainty_box):
    p=core.normalized(p);d=len(p['directions'])
    if not 1<=d<=2:raise ValueError('Robust prototype supports one or two design parameters')
    if not 1<=len(uncertainty_directions)<=2 or len(uncertainty_box)!=len(uncertainty_directions):raise ValueError('One or two uncertainty parameters required')
    dirs=[core.v3.base.potential(list(map(str,q))) for q in uncertainty_directions]
    box=[list(map(F,row)) for row in uncertainty_box]
    if any(len(row)!=2 or row[0]>=row[1] for row in box):raise ValueError('Invalid uncertainty box')
    proofs=[core.g.v5.make_range(list(map(F,q))) for q in p['directions']]
    inv,_,_=moment.confinement(p,proofs)
    scenarios=[]
    for theta in product(*box):
        q=list(map(F,p['q0']))
        for t,r in zip(theta,dirs):q=core.v3.old.add(q,core.v3.old.scale(r,t))
        scenarios.append(core.problem(q,p['directions'],p['box'],p['penalty']['linear'],p['penalty']['hessian'],p['penalty']['constant']))
    return inv,dirs,box,scenarios


def cut(p,scenario,c):
    alpha,m=moment.moments(p,c)
    return {'scenario':scenario,'trial_legendre':list(map(str,c)),'constant':str(-alpha),'linear':list(map(str,[-x for x in m]))}


def dual(p,inv,cuts,weights):
    if len(cuts)!=len(weights) or not cuts or min(weights)<0 or sum(weights,F(0))!=1:raise ValueError('Invalid dual weights')
    ell=list(map(F,p['penalty']['linear']));d=len(ell)
    z=[ell[i]+sum(w*F(c['linear'][i]) for w,c in zip(weights,cuts)) for i in range(d)]
    a=[-sum(x*y for x,y in zip(row,z)) for row in inv]
    low=F(p['penalty']['constant'])+sum(w*F(c['constant']) for w,c in zip(weights,cuts))+sum(x*y for x,y in zip(z,a))/2
    return low,a


def best_weights(p,inv,cuts):
    n=len(cuts);d=len(inv);ell=list(map(F,p['penalty']['linear']));slopes=[list(map(F,c['linear'])) for c in cuts]
    # Minimize the negative concave dual quadratic over the probability simplex.
    c=-F(p['penalty']['constant'])+sum(ell[i]*inv[i][j]*ell[j] for i in range(d) for j in range(d))/2
    l=[-F(cut['constant'])+sum(s[i]*inv[i][j]*ell[j] for i in range(d) for j in range(d)) for cut,s in zip(cuts,slopes)]
    h=[[sum(s[i]*inv[i][j]*t[j] for i in range(d) for j in range(d)) for t in slopes] for s in slopes]
    vertices=[[F(i==j) for i in range(n)] for j in range(n)]
    return simplex.minimum(c,l,h,vertices)[1]


def certificate(p,dirs,box,cuts,weights,candidate,scenario_records,tolerance,independent=True):
    inv,dirs,box,scenarios=data(p,dirs,box);p=core.normalized(p);tolerance=F(tolerance)
    if not F(1,10**12)<=tolerance<=100 or not 1<=len(cuts)<=8:raise ValueError('Invalid robust certificate budget')
    for entry in cuts:
        s=entry['scenario']
        if type(s) is not int or not 0<=s<len(scenarios):raise ValueError('Invalid scenario index')
        c=core.g.v7.coefficients(entry['trial_legendre'])
        if entry!=cut(scenarios[s],s,c):raise ValueError('Incorrect affine Rayleigh cut')
    weights=list(map(F,weights));lower,_=dual(p,inv,cuts,weights);candidate=list(map(F,candidate))
    if len(scenario_records)!=len(scenarios):raise ValueError('All uncertainty vertices required')
    lows=[]
    for sp,record in zip(scenarios,scenario_records):
        sc=record['spectral_certificate'];q=core.at(sp,candidate)
        if record['at']!=list(map(str,candidate)) or sc['eigenvalue_index']!=1 or core.v3.base.potential(sc['q_coefficients'])!=q or not core.v3.verify(sc,independent):raise ValueError('Invalid robust upper certificate')
        coeff=core.g.v7.coefficients(record['trial_legendre']);mu=core.v3.old.residual_statistics(q,coeff)[0]
        if record!={'at':list(map(str,candidate)),'spectral_certificate':sc,'trial_legendre':list(map(str,coeff)),'rayleigh_quotient':str(mu)}:raise ValueError('Invalid robust trial')
        lows.append(F(sc['lower']))
    upper=core.penalty(p,candidate)-min(lows)
    if upper<lower:raise ValueError('Inconsistent primal-dual bounds')
    return {'method':'robust_spectral_dual_v20','problem':p,'uncertainty_directions':[list(map(str,q)) for q in dirs],
            'uncertainty_box':[list(map(str,row)) for row in box],'cuts':cuts,'weights':list(map(str,weights)),
            'candidate_parameters':list(map(str,candidate)),'scenario_records':scenario_records,
            'global_lower':str(lower),'candidate_upper':str(upper),'global_gap':str(upper-lower),'tolerance':str(tolerance),
            'status':'epsilon_global' if upper-lower<=tolerance else 'global_gap_open',
            'optimization_scope':'all real design parameters; every real uncertainty vector in its box; ground spectrum on the axisymmetric sphere model',
            'formal_assistant_checked':False}


def verify(c,independent=True):
    try:return c==certificate(c['problem'],c['uncertainty_directions'],c['uncertainty_box'],c['cuts'],c['weights'],c['candidate_parameters'],c['scenario_records'],c['tolerance'],independent)
    except (ValueError,KeyError,TypeError,IndexError,ZeroDivisionError):return False


def search(p,dirs,box,tolerance=F(1,1000),max_steps=20,max_modes=24):
    start=perf_counter();inv,dirs,box,scenarios=data(p,dirs,box);tolerance=F(tolerance)
    if type(max_steps) is not int or not 1<=max_steps<=100:raise ValueError('Invalid step budget')
    if not F(1,10**12)<=tolerance<=100:raise ValueError('Invalid tolerance')
    oracles=[point_oracle.Oracle(sp,min(F(1),tolerance/16),max_modes=max_modes,orbits=True) for sp in scenarios]
    a=[F(0)]*len(inv);cuts=[];bestlower=None;bestupper=None;history=[]
    for step in range(max_steps):
        # Quantization is only for proposed evaluation points, never dual bounds.
        a=[F(round(x*2**20),2**20) for x in a]
        records=[o.records[o.point(a)] for o in oracles]
        upper=core.penalty(p,a)-min(F(r['spectral_certificate']['lower']) for r in records)
        if bestupper is None or upper<bestupper[0]:bestupper=(upper,a,records)
        for s,(sp,r) in enumerate(zip(scenarios,records)):
            new=cut(sp,s,list(map(F,r['trial_legendre'])))
            if new not in cuts:cuts.append(new)
        if len(cuts)>8:raise ValueError('Active dual cut budget exceeded')
        w=best_weights(p,inv,cuts);lower,a=dual(p,inv,cuts,w)
        if bestlower is None or lower>bestlower[0]:bestlower=(lower,cuts[:],w[:])
        history.append({'step':step,'lower':str(bestlower[0]),'upper':str(bestupper[0]),'active_cuts':len(cuts)})
        if bestupper[0]-bestlower[0]<=tolerance:break
        cuts=[c for c,weight in zip(cuts,w) if weight>0]
    _,savedcuts,weights=bestlower;_,candidate,records=bestupper
    cert=certificate(p,dirs,box,savedcuts,weights,candidate,records,tolerance)
    return {'algorithm_version':20,'certificate':cert,'steps':len(history),'history':history,
            'statistics':{'spectral_searches':sum(o.solves for o in oracles),'points':sum(len(o.records) for o in oracles)},
            'seconds':perf_counter()-start}
