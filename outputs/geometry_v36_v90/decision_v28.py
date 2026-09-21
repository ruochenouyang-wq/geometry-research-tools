"""V28: threshold-directed whole-domain proof or verified counterexample search."""
from fractions import Fraction as F
from itertools import product
from time import perf_counter
import copy
import compiler_v23 as compiler
import global_core as core
import point_oracle
import algebra
import reductions
import moment_v12 as moment


def as_core(p):
    p=compiler.normalized(p)
    if p['function_space']!='axisymmetric':raise ValueError('Apply the full-sphere reduction first')
    q=p['penalty'];d=len(p['directions'])
    if not d:raise ValueError('Fixed operator uses the dedicated point rule')
    if p['domain']['kind']=='all_real':return core.all_real_problem(p['q0'],p['directions'],q['linear'],q['hessian'],q['constant'])
    return core.problem(p['q0'],p['directions'],p['domain']['box'],q['linear'],q['hessian'],q['constant'])


def fixed(p,tolerance=F(1,10**6),max_modes=24):
    q0,dirs,c,l,h=reductions.data(p)
    if p['function_space']!='axisymmetric':raise ValueError('Apply the full-sphere reduction first')
    if dirs:raise ValueError('Fixed rule requires zero parameters')
    fake=core.problem(q0,[[0]],[[-1,1]])
    oracle=point_oracle.Oracle(fake,tolerance,max_modes=max_modes)
    record=oracle.records[oracle.point([0])]
    return {'method':'fixed_spectrum','problem':p,'record':record,'lower_bound':str(F(record['spectral_certificate']['lower'])+c),
            'upper_bound':str(F(record['rayleigh_quotient'])+c),'candidate_parameters':[]}


def analytic(p):
    q0,dirs,c,l,h=reductions.data(p)
    if p['domain']['kind']=='all_real':
        cp=as_core(p);box=[list(map(F,r)) for r in cp['box']]
    else:box=[list(map(F,r)) for r in p['domain']['box']]
    proofs=[];lows=[]
    for x in product(*box):
        q=q0[:]
        for v,direction in zip(x,dirs):q=core.v3.old.add(q,core.v3.old.scale(direction,v))
        proof=core.g.v5.make_range(q);proofs.append(proof);lows.append(core.v3.range_bounds(q,proof)[0])
    low=algebra.quadratic_min(c,l,h,box)[0]+min(lows)
    return {'method':'range_quadratic_bound','problem':p,'corner_ranges':proofs,'lower_bound':str(low)}


def verify_fixed(cert):
    try:
        p=compiler.normalized(cert['problem']);q0,dirs,c,_,_=reductions.data(p);r=cert['record'];sc=r['spectral_certificate']
        if p['function_space']!='axisymmetric' or dirs or sc['eigenvalue_index']!=1 or list(map(F,sc['q_coefficients']))!=q0 or not core.v3.verify(sc,True):return False
        coef=core.g.v7.coefficients(r['trial_legendre']);mu=core.v3.old.residual_statistics(q0,coef)[0]
        return cert=={'method':'fixed_spectrum','problem':p,'record':r,'lower_bound':str(F(sc['lower'])+c),'upper_bound':str(mu+c),'candidate_parameters':[]} and r=={'at':['0'],'spectral_certificate':sc,'trial_legendre':list(map(str,coef)),'rayleigh_quotient':str(mu)}
    except (ValueError,KeyError,TypeError,IndexError,ZeroDivisionError):return False


def result(p,threshold,inner,kind):
    p=compiler.normalized(p);threshold=F(threshold)
    if kind=='fixed':
        if inner['problem']!=p or not verify_fixed(inner):raise ValueError('Invalid fixed spectrum')
        low=F(inner['lower_bound']);up=F(inner['upper_bound']);x=[]
    elif kind=='analytic':
        if inner!=analytic(p):raise ValueError('Invalid analytic floor')
        low=F(inner['lower_bound']);up=None;x=None
    elif kind=='partition':
        cp=as_core(p)
        if inner['problem']!=cp or inner['domain']!=p['domain']['kind'] or not core.verify(inner,True):raise ValueError('Invalid full-domain decision proof')
        low=F(inner['global_lower']);up=F(inner['candidate_upper']);x=inner['candidate_parameters']
    else:raise ValueError('Unknown proof kind')
    status='proved' if low>=threshold else ('refuted' if up is not None and up<threshold else 'unresolved')
    return {'method':'threshold_decision_v28','problem':p,'threshold':str(threshold),'proof_kind':kind,'inner':inner,
            'lower_bound':str(low),'upper_bound':None if up is None else str(up),'candidate_parameters':x,'status':status,'formal_assistant_checked':False}


def verify(cert):
    try:return cert==result(cert['problem'],cert['threshold'],cert['inner'],cert['proof_kind'])
    except (ValueError,KeyError,TypeError,IndexError,ZeroDivisionError):return False


def search(p,threshold,max_leaves=128,max_modes=24,tolerance=F(1,10**6),try_analytic=True):
    start=perf_counter();p=compiler.normalized(p);threshold=F(threshold)
    if p['function_space']!='axisymmetric':raise ValueError('Apply the full-sphere reduction first')
    if type(max_leaves) is not int or not 1<=max_leaves<=512:raise ValueError('Decision leaf budget')
    if not p['directions']:
        cert=result(p,threshold,fixed(p,tolerance,max_modes),'fixed')
        return {'certificate':cert,'seconds':perf_counter()-start,'statistics':{'points':1,'parameter_splits':0}}
    cp=as_core(p);d=len(p['directions']);root=[list(map(F,r)) for r in cp['box']]
    if try_analytic:
        quick=result(p,threshold,analytic(p),'analytic')
        if quick['status']=='proved':return {'certificate':quick,'seconds':perf_counter()-start,'statistics':{'points':0,'parameter_splits':0}}
    proofs=[core.g.v5.make_range(list(map(F,q))) for q in cp['directions']]
    inv=moment.confinement(cp,proofs)[0] if p['domain']['kind']=='all_real' else None
    oracle=point_oracle.Oracle(cp,F(tolerance),max_modes=max_modes,orbits=True);active={};nodes={}
    kind='center_simplex' if d<=2 else 'refined_bernstein'
    def add(path,box):
        ids=[oracle.point(x) for x in core.positions(box)];node={'kind':kind,'indices':ids}
        if d>2:
            poly=core.bernstein.lower_polynomial(cp,box,[F(oracle.records[i]['spectral_certificate']['lower']) for i in ids[:-1]])
            low,tree=core.bernstein.refine(poly,d,8);node['refinement']=tree
        else:low=core.cell_bound(cp,box,oracle.records,ids,kind)
        node['lower']=str(low);nodes[path]=node;active[path]=(low,box)
    add('',root);steps=0
    while True:
        path=min(active,key=lambda k:active[k][0]);low,box=active[path]
        upper=min(moment.eliminated_value(cp,list(map(F,r['trial_legendre'])),inv)[0] for r in oracle.records) if inv is not None else min(F(r['rayleigh_quotient'])+core.penalty(cp,list(map(F,r['at']))) for r in oracle.records)
        if low>=threshold or upper<threshold or len(active)>=max_leaves or len(path)>=20:break
        axis=max(range(d),key=lambda j:box[j][1]-box[j][0]);mid=sum(box[axis])/2
        left=copy.deepcopy(box);right=copy.deepcopy(box);left[axis][1]=mid;right[axis][0]=mid
        nodes[path]={'axis':axis,'split':str(mid)};del active[path];add(path+'0',left);add(path+'1',right);steps+=1
    def tree(path):
        node=nodes[path]
        return {**node,'left':tree(path+'0'),'right':tree(path+'1')} if 'axis' in node else node
    boxed=core.certificate(cp,p['domain']['kind'],proofs,oracle.records,tree(''),F(1,10**12))
    cert=result(p,threshold,boxed,'partition')
    return {'certificate':cert,'seconds':perf_counter()-start,'statistics':{**oracle.statistics(),'parameter_splits':steps}}
