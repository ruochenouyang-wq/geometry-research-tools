"""Certified global minimization of a regularized infinite-dimensional ground energy.

min_{a in box, u != 0} <u,H(q0 + sum a_i p_i)u>/<u,u> + Q(a).
The box has one or two real parameters. All acceptance arithmetic is rational.
"""
from fractions import Fraction as F
from itertools import product
from time import perf_counter
import v03_tail as v3
import v04_search as v4
import v05_range as v5
import v07_precision as v7
import global_quadratic as quadratic

METHOD='global_regularized_ground_v11'


def rationals(values):return [v3.base.rational(str(x)) for x in values]


def problem(q0,directions,box,linear=None,hessian=None,constant=0):
    q0=v3.base.potential([str(x) for x in q0])
    if not isinstance(directions,list) or len(directions) not in (1,2):raise ValueError('One or two affine directions required')
    d=len(directions);directions=[v3.base.potential([str(x) for x in p]) for p in directions]
    if not isinstance(box,list) or len(box)!=d or any(len(x)!=2 for x in box):raise ValueError('Box dimension mismatch')
    box=[rationals(x) for x in box]
    if any(lo>=hi or max(abs(lo),abs(hi))>1000 for lo,hi in box):raise ValueError('Nonempty box with endpoints bounded by 1000 required')
    linear=rationals([0]*d if linear is None else linear)
    hessian=[rationals(row) for row in ([[0]*d for _ in range(d)] if hessian is None else hessian)]
    if len(linear)!=d or len(hessian)!=d or any(len(row)!=d for row in hessian):raise ValueError('Penalty dimension mismatch')
    if any(hessian[i][j]!=hessian[j][i] for i in range(d) for j in range(d)):raise ValueError('Penalty Hessian must be symmetric')
    constant=v3.base.rational(str(constant))
    result={'q0':[str(x) for x in q0],'directions':[[str(x) for x in p] for p in directions],
            'box':[[str(x) for x in pair] for pair in box],
            'penalty':{'constant':str(constant),'linear':[str(x) for x in linear],
                       'hessian':[[str(x) for x in row] for row in hessian]}}
    for corner in product(*box):at(result,corner)
    return result


def normalized(raw):
    p=raw['penalty']
    result=problem(raw['q0'],raw['directions'],raw['box'],p['linear'],p['hessian'],p['constant'])
    if result!=raw:raise ValueError('Problem must use canonical rational fields')
    return result


def at(p,x):
    if len(x)!=len(p['directions']):raise ValueError('Parameter dimension mismatch')
    q=rationals(p['q0'])
    for a,direction in zip(x,p['directions']):q=v3.old.add(q,v3.old.scale(rationals(direction),F(a)))
    return v3.base.potential([str(z) for z in q])


def penalty(p,x):
    cost=p['penalty']
    return quadratic.value(F(cost['constant']),rationals(cost['linear']),[rationals(row) for row in cost['hessian']],x)


def inside(x,box):return len(x)==len(box) and all(lo<=z<=hi for z,(lo,hi) in zip(x,box))


class Oracle:
    def __init__(self,p,modes,tolerance):
        self.p=p;self.modes=modes;self.tolerance=tolerance;self.records=[];self.cache={};self.values=[]

    def point(self,x):
        key=tuple(x)
        if key in self.cache:return self.cache[key]
        q=at(self.p,x)
        cert=v4.guided(q,1,self.modes,self.tolerance,range_proof=v5.make_range(q))['certificate']
        even=not any(q[1::2]);c=v3.old.polynomial_trial(q,self.modes,even=even)
        mu,_,_=v3.old.residual_statistics(q,c)
        record={'at':[str(z) for z in x],'spectral_certificate':cert,
                'trial_legendre':[str(z) for z in c],'rayleigh_quotient':str(mu)}
        index=len(self.records);self.cache[key]=index;self.records.append(record)
        self.values.append((F(cert['lower'])+penalty(self.p,x),mu+penalty(self.p,x)))
        return index


def relaxation(p,box,records,indices,policy,ranges):
    d=len(box);cost=p['penalty']
    c,l,h=quadratic.on_unit_box(F(cost['constant']),rationals(cost['linear']),[rationals(row) for row in cost['hessian']],box)
    if policy=='concavity':
        corners=list(product(*box))
        if len(indices)!=2**d:raise ValueError('All corners are required')
        data={}
        for bits,x,idx in zip(product((0,1),repeat=d),corners,indices):
            if rationals(records[idx]['at'])!=list(x):raise ValueError('Corner point mismatch')
            data[bits]=F(records[idx]['spectral_certificate']['lower'])
        origin=(0,)*d;base=data[origin];c+=base
        for i in range(d):
            bits=tuple(int(j==i) for j in range(d));l[i]+=data[bits]-base
        if d==2:
            mixed=data[(1,1)]-data[(1,0)]-data[(0,1)]+base
            h[0][1]+=mixed;h[1][0]+=mixed
    elif policy=='lipschitz':
        mid=[(lo+hi)/2 for lo,hi in box]
        if len(indices)!=1 or rationals(records[indices[0]]['at'])!=mid:raise ValueError('Center point mismatch')
        c+=F(records[indices[0]]['spectral_certificate']['lower'])
        for i,(lo,hi) in enumerate(ranges):
            # Remove exact constant shifts before the norm perturbation bound.
            shift=(lo+hi)/2;radius=(hi-lo)/2;width=box[i][1]-box[i][0]
            c-=width*(shift+radius)/2;l[i]+=width*shift
    else:raise ValueError('Unknown global lower-bound policy')
    return quadratic.minimum(c,l,h,[(F(0),F(1))]*d)


def certificate(p,policy,range_proofs,records,tree,tolerance,independent=True):
    p=normalized(p);root=[rationals(pair) for pair in p['box']];d=len(root);tolerance=F(tolerance)
    if not F(1,10**12)<=tolerance<=100:raise ValueError('Global tolerance must be in [1e-12,100]')
    if len(range_proofs)!=d:raise ValueError('Direction range proofs required')
    ranges=[v3.range_bounds(rationals(direction),proof) for direction,proof in zip(p['directions'],range_proofs)]
    if not isinstance(records,list) or not 1<=len(records)<=8192:raise ValueError('Point record budget exceeded')
    upper=[];seen=set()
    for record in records:
        x=rationals(record['at']);q=at(p,x);key=tuple(x)
        if not inside(x,root) or key in seen:raise ValueError('Point outside box or duplicate')
        seen.add(key);sc=record['spectral_certificate']
        if sc['eigenvalue_index']!=1 or v3.base.potential(sc['q_coefficients'])!=q or not v3.verify(sc,independent):
            raise ValueError('Ground-state certificate mismatch')
        coeff=v7.coefficients(record['trial_legendre']);mu,_,_=v3.old.residual_statistics(q,coeff)
        if record!={'at':[str(z) for z in x],'spectral_certificate':sc,
                    'trial_legendre':[str(z) for z in coeff],'rayleigh_quotient':str(mu)}:raise ValueError('Incorrect trial value')
        upper.append(mu+penalty(p,x))
    best=min(range(len(upper)),key=upper.__getitem__);u=upper[best];leaves=[]

    def check_node(node,box,depth=0):
        if depth>20 or not isinstance(node,dict):raise ValueError('Invalid partition tree')
        if 'axis' in node:
            if set(node)!= {'axis','split','left','right'}:raise ValueError('Incomplete split node')
            axis=node['axis'];split=v3.base.rational(node['split'])
            if type(axis) is not int or not 0<=axis<d or not box[axis][0]<split<box[axis][1]:raise ValueError('Invalid split')
            left=[pair[:] for pair in box];right=[pair[:] for pair in box]
            left[axis][1]=split;right[axis][0]=split
            check_node(node['left'],left,depth+1);check_node(node['right'],right,depth+1)
        else:
            if set(node)!= {'oracle_indices','lower_bound'}:raise ValueError('Invalid leaf fields')
            ids=node['oracle_indices']
            if not isinstance(ids,list) or any(type(i) is not int or not 0<=i<len(records) for i in ids):raise ValueError('Invalid point reference')
            lower,_=relaxation(p,box,records,ids,policy,ranges)
            if node['lower_bound']!=str(lower):raise ValueError('Invalid whole-cell lower bound')
            leaves.append((lower,box))
            if len(leaves)>2048:raise ValueError('Leaf budget exceeded')
    check_node(tree,root)
    lower=min(v for v,box in leaves)
    if lower>u:raise ArithmeticError('Global bounds are inconsistent')
    return {'method':METHOD,'model':v3.base.MODEL,'scope':v3.base.SCOPE,'problem':p,
            'optimization_scope':'all real parameters in the declared box and all nonzero axisymmetric form-domain functions',
            'bound_policy':policy,'direction_range_proofs':range_proofs,'points':records,'partition_tree':tree,
            'candidate_index':best,'candidate_parameters':records[best]['at'],
            'global_lower':str(lower),'candidate_upper':str(u),'global_gap':str(u-lower),
            'tolerance':str(tolerance),'status':'epsilon_global' if u-lower<=tolerance else 'global_gap_open',
            'leaf_count':len(leaves),'excluded_leaf_count':sum(v>u for v,box in leaves),
            'possible_optimizer_boxes':[[[str(z) for z in pair] for pair in box] for v,box in leaves if v<=u],
            'formal_assistant_checked':False}


def verify(cert,independent=True):
    try:
        return cert==certificate(cert['problem'],cert['bound_policy'],cert['direction_range_proofs'],
                                 cert['points'],cert['partition_tree'],cert['tolerance'],independent)
    except (KeyError,ValueError,TypeError,ZeroDivisionError,IndexError,RecursionError):return False


def search(p,tolerance=F(1,1000),modes=12,max_leaves=256,policy='concavity',seeds=None):
    start=perf_counter();p=normalized(p);root=[rationals(pair) for pair in p['box']];d=len(root)
    v3.base.check_sizes(1,modes);tolerance=F(tolerance)
    if not F(1,10**12)<=tolerance<=100:raise ValueError('Global tolerance must be in [1e-12,100]')
    if type(max_leaves) is not int or not 1<=max_leaves<=2048:raise ValueError('Leaf budget must be 1..2048')
    if policy not in ('concavity','lipschitz'):raise ValueError('Unknown bound policy')
    proofs=[v5.make_range(rationals(p)) for p in p['directions']]
    ranges=[v3.range_bounds(rationals(direction),proof) for direction,proof in zip(p['directions'],proofs)]
    oracle=Oracle(p,modes,min(F(1,10**8),tolerance/100));active={};nodes={};history=[]
    for seed in seeds or []:
        x=rationals(seed)
        if not inside(x,root):raise ValueError('Seed outside the search box')
        oracle.point(x)

    def add_leaf(path,box):
        locations=list(product(*box)) if policy=='concavity' else [[(lo+hi)/2 for lo,hi in box]]
        indices=[oracle.point(x) for x in locations]
        bound,_=relaxation(p,box,oracle.records,indices,policy,ranges)
        nodes[path]={'oracle_indices':indices,'lower_bound':str(bound)}
        active[path]=(bound,box)
        # Feasible center samples also improve the incumbent, never the cell proof.
        oracle.point([(lo+hi)/2 for lo,hi in box])
    add_leaf('',root);reason='leaf_budget'
    while True:
        path=min(active,key=lambda x:active[x][0]);lower,box=active[path]
        upper=min(v[1] for v in oracle.values)
        history.append({'leaves':len(active),'evaluations':len(oracle.records),'lower':str(lower),'upper':str(upper),'gap':str(upper-lower)})
        if upper-lower<=tolerance:reason='epsilon_global';break
        if len(active)>=max_leaves:break
        if len(path)>=18:reason='parameter_resolution_budget';break
        axis=max(range(d),key=lambda i:box[i][1]-box[i][0]);split=sum(box[axis],F(0))/2
        left=[pair[:] for pair in box];right=[pair[:] for pair in box]
        left[axis][1]=split;right[axis][0]=split
        nodes[path]={'axis':axis,'split':str(split)};del active[path]
        add_leaf(path+'0',left);add_leaf(path+'1',right)

    def build(path):
        node=nodes[path]
        if 'axis' not in node:return node
        return {**node,'left':build(path+'0'),'right':build(path+'1')}
    cert=certificate(p,policy,proofs,oracle.records,build(''),tolerance,independent=True)
    return {'algorithm_version':11,'certificate':cert,'stop_reason':reason,'modes':modes,
            'point_evaluations':len(oracle.records),'history':history,'seconds':perf_counter()-start}
