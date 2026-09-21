"""V13--V18 / V22: rational certificates for global function/parameter minima."""
from fractions import Fraction as F
from itertools import product
from time import perf_counter
import copy
import global_search as g
import global_quadratic as quadratic
import moment_v12 as moment
import point_oracle
import simplex_v15 as simplex
import bernstein_v17 as bernstein
v3=g.v3
rationals=g.rationals
at=g.at
penalty=g.penalty
inside=g.inside

def problem(q0,directions,box,linear=None,hessian=None,constant=0):
    q0=v3.base.potential([str(x) for x in q0])
    if not isinstance(directions,list) or len(directions) not in (1,2,3,4):raise ValueError('One to four affine directions required')
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


def all_real_problem(q0,directions,linear,hessian,constant=0):
    p=problem(q0,directions,[[-1,1]]*len(directions),linear,hessian,constant)
    proofs=[g.v5.make_range(rationals(q)) for q in p['directions']]
    _,_,box=moment.confinement(p,proofs)
    p['box']=[[str(lo if lo<hi else lo-1),str(hi if lo<hi else hi+1)] for lo,hi in box]
    return normalized(p)


def moment_exclusion(p,box):
    if p['directions']!=[['0','1'],['0','0','1']]:return None
    ell=rationals(p['penalty']['linear']);k=[rationals(row) for row in p['penalty']['hessian']]
    # f(a)=m_2-m_1^2, m=-ell-Ka. The maximum of f below zero excludes all realizable moments.
    c=-ell[1]-ell[0]**2
    l=[-k[1][i]-2*ell[0]*k[0][i] for i in range(2)]
    h=[[-2*k[0][i]*k[0][j] for j in range(2)] for i in range(2)]
    negative_min,_=quadratic.minimum(-c,[-z for z in l],[[-z for z in row] for row in h],box)
    upper=-negative_min
    return {'kind':'moment_exclusion','variance_upper':str(upper)} if upper<0 else None


def positions(box):return [list(x) for x in product(*box)]+[[sum(pair)/2 for pair in box]]


def cell_bound(p,box,records,ids,kind,refinement=None):
    d=len(box);loc=positions(box)
    if len(ids)!=len(loc) or any(type(i) is not int or not 0<=i<len(records) for i in ids):raise ValueError('Invalid point references')
    if any(rationals(records[idx]['at'])!=x for idx,x in zip(ids,loc)):raise ValueError('Cell coordinates mismatch')
    values=[F(records[i]['spectral_certificate']['lower']) for i in ids]
    if kind=='quadratic':
        return g.relaxation(p,box,records,ids[:-1],'concavity',[])[0]
    if kind=='center_simplex':
        if d not in (1,2):raise ValueError('Center triangulation supports one or two parameters')
        cost=p['penalty'];c,l,h=quadratic.on_unit_box(F(cost['constant']),rationals(cost['linear']),[rationals(r) for r in cost['hessian']],box)
        vertices=list(product((F(0),F(1)),repeat=d))+[(F(1,2),)*d]
        faces=[(0,2),(2,1)] if d==1 else [(0,1,4),(1,3,4),(3,2,4),(2,0,4)]
        lows=[]
        for face in faces:
            v=[vertices[i] for i in face];a,b=simplex.affine_interpolant(v,[values[i] for i in face])
            lows.append(simplex.minimum(c+a,[x+y for x,y in zip(l,b)],h,v)[0])
        previous=g.relaxation(p,box,records,ids[:-1],'concavity',[])[0]
        return max(previous,min(lows))
    poly=bernstein.lower_polynomial(p,box,values[:-1])
    if kind=='bernstein':return min(bernstein.coefficients(poly,d).values())
    if kind=='refined_bernstein':return bernstein.verify_refinement(poly,d,refinement)
    raise ValueError('Unknown cell proof')


def verify_points(p,records,root,independent):
    if not isinstance(records,list) or not 1<=len(records)<=16384:raise ValueError('Point budget exceeded')
    upper=[];seen=set()
    for record in records:
        x=rationals(record['at']);q=at(p,x);sc=record['spectral_certificate']
        if not inside(x,root) or tuple(x) in seen:raise ValueError('Duplicate or outside point')
        seen.add(tuple(x))
        if sc['eigenvalue_index']!=1 or v3.base.potential(sc['q_coefficients'])!=q or not v3.verify(sc,independent):raise ValueError('Invalid ground spectrum')
        coeff=g.v7.coefficients(record['trial_legendre']);mu=v3.old.residual_statistics(q,coeff)[0]
        if record!={'at':list(map(str,x)),'spectral_certificate':sc,'trial_legendre':list(map(str,coeff)),'rayleigh_quotient':str(mu)}:raise ValueError('Invalid trial')
        upper.append((mu+penalty(p,x),x))
    return upper


def certificate(p,domain,proofs,records,tree,tolerance,independent=True):
    p=normalized(p);d=len(p['directions']);root=[rationals(r) for r in p['box']];tolerance=F(tolerance)
    if not F(1,10**12)<=tolerance<=100:raise ValueError('Invalid tolerance')
    if domain not in ('box','all_real'):raise ValueError('Invalid domain')
    if len(proofs)!=d:raise ValueError('Missing direction proofs')
    for q,proof in zip(p['directions'],proofs):v3.range_bounds(rationals(q),proof)
    inv=None
    if domain=='all_real':
        inv,_,confined=moment.confinement(p,proofs)
        if any(a>lo or b<hi for (a,b),(lo,hi) in zip(root,confined)):raise ValueError('Incomplete all-real confinement')
    candidates=verify_points(p,records,root,independent)
    if inv is not None:
        candidates=[moment.eliminated_value(p,g.v7.coefficients(r['trial_legendre']),inv)[:2] for r in records]
    best=min(range(len(candidates)),key=lambda i:candidates[i][0]);upper,x=candidates[best];leaves=[];excluded=[]
    def visit(node,box,depth):
        if depth>24:raise ValueError('Partition depth exceeded')
        if 'axis' in node:
            if set(node)!= {'axis','split','left','right'}:raise ValueError('Incomplete split')
            axis=node['axis'];mid=F(node['split'])
            if type(axis) is not int or not 0<=axis<d or not box[axis][0]<mid<box[axis][1]:raise ValueError('Invalid split')
            left=copy.deepcopy(box);right=copy.deepcopy(box);left[axis][1]=mid;right[axis][0]=mid
            visit(node['left'],left,depth+1);visit(node['right'],right,depth+1)
        elif node.get('kind')=='moment_exclusion':
            if domain!='all_real' or node!=moment_exclusion(p,box):raise ValueError('Invalid moment exclusion')
            excluded.append(box)
        else:
            fields={'kind','indices','lower'}|({'refinement'} if node.get('kind')=='refined_bernstein' else set())
            if set(node)!=fields:raise ValueError('Invalid cell fields')
            low=cell_bound(p,box,records,node['indices'],node['kind'],node.get('refinement'))
            if node['lower']!=str(low):raise ValueError('Invalid cell bound')
            leaves.append((low,box))
        if len(leaves)+len(excluded)>2048:raise ValueError('Leaf budget exceeded')
    visit(tree,root,0)
    if not leaves:raise ValueError('No realizable parameter region remains')
    lower=min(z[0] for z in leaves)
    if upper<lower:raise ValueError('Inconsistent global bracket')
    return {'method':'global_function_bounds_v13_v22','model':v3.base.MODEL,'scope':v3.base.SCOPE,
            'problem':p,'domain':domain,'direction_proofs':proofs,'points':records,'tree':tree,
            'global_lower':str(lower),'candidate_upper':str(upper),'global_gap':str(upper-lower),
            'candidate_index':best,'candidate_parameters':list(map(str,x)),
            'tolerance':str(tolerance),'status':'epsilon_global' if upper-lower<=tolerance else 'global_gap_open',
            'leaf_count':len(leaves)+len(excluded),'moment_excluded_leaves':len(excluded),
            'formal_assistant_checked':False}


def verify(cert,independent=True):
    try:return cert==certificate(cert['problem'],cert['domain'],cert['direction_proofs'],cert['points'],cert['tree'],cert['tolerance'],independent)
    except (ValueError,TypeError,KeyError,IndexError,ZeroDivisionError,RecursionError):return False


def search(p,tolerance=F(1,1000),version=22,domain='box',max_leaves=256,max_modes=24,
           start_modes=4,adaptive=True,orbits=None,kind=None,moment_cuts=None,refinement_budget=8,max_actions=2048):
    start=perf_counter();p=normalized(p);d=len(p['directions']);root=[rationals(r) for r in p['box']];tolerance=F(tolerance)
    if type(version) is not int or not 13<=version<=22:raise ValueError('Invalid version')
    if version in (20,21):raise ValueError('Use the robust_v20 or constrained_v21 entry point for its distinct objective')
    if d>2 and version<17:raise ValueError('Three and four parameters require V17 or later')
    if not F(1,10**12)<=tolerance<=100:raise ValueError('Invalid tolerance')
    if type(max_leaves) is not int or not 1<=max_leaves<=2048:raise ValueError('Invalid leaf budget')
    if type(max_actions) is not int or not 1<=max_actions<=8192:raise ValueError('Invalid action budget')
    if type(refinement_budget) is not int or not 1<=refinement_budget<=128:raise ValueError('Invalid refinement budget')
    if kind is None:kind=('quadratic' if version<15 else 'center_simplex') if d<=2 else ('bernstein' if version<18 else 'refined_bernstein')
    if orbits is None:orbits=version>=14
    if moment_cuts is None:moment_cuts=version>=16 and domain=='all_real'
    if moment_cuts and domain!='all_real':raise ValueError('Moment cuts require all-real parameter elimination')
    proofs=[g.v5.make_range(rationals(q)) for q in p['directions']]
    inv=None
    if domain=='all_real':
        inv,_,confined=moment.confinement(p,proofs)
        if any(a>lo or b<hi for (a,b),(lo,hi) in zip(root,confined)):raise ValueError('Search box fails confinement')
    elif domain!='box':raise ValueError('Invalid domain')
    scheduled=version==22
    oracle=point_oracle.Oracle(p,min(F(1),tolerance/(2 if scheduled else 16)),start_modes,max_modes,adaptive,orbits)
    active={};nodes={};budgets={};history=[];actions={'parameter_splits':0,'point_refinements':0,'polynomial_refinements':0}
    def upper():
        if inv is not None:return min(moment.eliminated_value(p,rationals(r['trial_legendre']),inv)[0] for r in oracle.records)
        return min(F(r['rayleigh_quotient'])+penalty(p,rationals(r['at'])) for r in oracle.records)
    def add(path,box):
        exclusion=moment_exclusion(p,box) if moment_cuts else None
        if exclusion is not None:nodes[path]=exclusion;return
        ids=[oracle.point(x) for x in positions(box)]
        node={'kind':kind,'indices':ids};refinement=None
        if kind=='refined_bernstein':
            values=[F(oracle.records[i]['spectral_certificate']['lower']) for i in ids[:-1]]
            poly=bernstein.lower_polynomial(p,box,values)
            low,refinement=bernstein.refine(poly,d,budgets.get(path,1 if scheduled else refinement_budget));node['refinement']=refinement
        else:low=cell_bound(p,box,oracle.records,ids,kind,refinement)
        # The final verifier independently expands every retained Bernstein leaf.
        # Repeating that expansion at each search proposal adds no proof strength.
        node['lower']=str(low)
        nodes[path]=node;active[path]=(low,box)
    add('',root);reason='action_budget';failed_refinements=set()
    for step in range(max_actions):
        if not active:raise ArithmeticError('All cells were excluded')
        path=min(active,key=lambda k:active[k][0]);low,box=active[path];up=upper()
        history.append({'action':step,'leaves':sum('axis' not in n for n in nodes.values()),'points':len(oracle.records),'lower':str(low),'upper':str(up),'gap':str(up-low)})
        if up-low<=tolerance:reason='epsilon_global';break
        if scheduled:
            ids=nodes[path]['indices'];i=max(ids,key=lambda j:point_oracle.width(oracle.records[j]))
            if point_oracle.width(oracle.records[i])>tolerance/8 and i not in failed_refinements:
                improved=oracle.refine(i,tolerance/32);actions['point_refinements']+=1
                if not improved:failed_refinements.add(i)
                for key,(val,b) in list(active.items()):
                    if i in nodes[key]['indices']:add(key,b)
                continue
            if kind=='refined_bernstein' and budgets.get(path,1)<refinement_budget:
                # Only refine a polynomial whose computable value range exceeds the target.
                values=[F(oracle.records[i]['spectral_certificate']['lower']) for i in ids[:-1]]
                poly=bernstein.lower_polynomial(p,box,values)
                sample=min(bernstein.evaluate(poly,x) for x in product((F(0),F(1,2),F(1)),repeat=d))
                if sample-low>tolerance/4:
                    # Batch the algebraic work once instead of regenerating the
                    # same 1, 2, 4, ... leaves at every scheduler action.
                    budgets[path]=refinement_budget;add(path,box);actions['polynomial_refinements']+=1;continue
        if sum('axis' not in n for n in nodes.values())>=max_leaves:reason='leaf_budget';break
        if len(path)>=20:reason='parameter_resolution_budget';break
        axis=max(range(d),key=lambda j:box[j][1]-box[j][0]);mid=sum(box[axis])/2
        left=copy.deepcopy(box);right=copy.deepcopy(box);left[axis][1]=mid;right[axis][0]=mid
        nodes[path]={'axis':axis,'split':str(mid)};del active[path]
        add(path+'0',left);add(path+'1',right);actions['parameter_splits']+=1
    def tree(path):
        n=nodes[path]
        return {**n,'left':tree(path+'0'),'right':tree(path+'1')} if 'axis' in n else n
    cert=certificate(p,domain,proofs,oracle.records,tree(''),tolerance)
    return {'algorithm_version':version,'certificate':cert,'stop_reason':reason,'statistics':oracle.statistics(),
            'actions':actions,'history':history,'seconds':perf_counter()-start}
