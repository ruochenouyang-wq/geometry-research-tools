"""V17 tensor Bernstein lower bounds; V18 certified adaptive subdivision."""
from fractions import Fraction as F
from itertools import product
from math import comb


def add_term(p,key,value):
    p[key]=p.get(key,F(0))+F(value)
    if not p[key]:del p[key]


def lower_polynomial(p,box,corner_values):
    import global_search as g
    from global_quadratic import on_unit_box
    d=len(box);cost=p['penalty']
    c,l,h=on_unit_box(F(cost['constant']),list(map(F,cost['linear'])),[list(map(F,row)) for row in cost['hessian']],box)
    poly={};zero=(0,)*d;add_term(poly,zero,c)
    for i in range(d):
        e=list(zero);e[i]=1;add_term(poly,tuple(e),l[i])
        for j in range(d):
            e=list(zero);e[i]+=1;e[j]+=1;add_term(poly,tuple(e),h[i][j]/2)
    for bits,value in zip(product((0,1),repeat=d),corner_values):
        terms={zero:F(value)}
        for i,bit in enumerate(bits):
            new={}
            for e,v in terms.items():
                if not bit:add_term(new,e,v)
                k=list(e);k[i]+=1;add_term(new,tuple(k),v if bit else -v)
            terms=new
        for e,v in terms.items():add_term(poly,e,v)
    return poly


def evaluate(poly,x):
    return sum((v*product_value([x[i]**n for i,n in enumerate(e)]) for e,v in poly.items()),F(0))


def product_value(values):
    result=F(1)
    for v in values:result*=v
    return result


def affine_change(poly,box):
    d=len(box);out={}
    for exponent,coef in poly.items():
        for powers in product(*(range(n+1) for n in exponent)):
            v=coef*product_value([F(comb(n,k))*lo**(n-k)*(hi-lo)**k for n,k,(lo,hi) in zip(exponent,powers,box)])
            add_term(out,powers,v)
    return out


def coefficients(poly,d):
    return {k:sum((v*product_value([F(comb(a,b),comb(2,b)) for a,b in zip(k,e)])
                   for e,v in poly.items() if all(a>=b for a,b in zip(k,e))),F(0))
            for k in product(range(3),repeat=d)}


def split(coeff,d,axis):
    left={};right={}
    for rest in product(range(3),repeat=d-1):
        keys=[]
        for k in range(3):
            idx=list(rest);idx.insert(axis,k);keys.append(tuple(idx))
        a,b,c=[coeff[k] for k in keys];mid=(a+2*b+c)/4
        for k,v,w in zip(keys,[a,(a+b)/2,mid],[mid,(b+c)/2,c]):left[k]=v;right[k]=w
    return left,right


def refine(poly,d,max_leaves=1):
    if type(max_leaves) is not int or not 1<=max_leaves<=128:raise ValueError('Invalid Bernstein budget')
    root=[(F(0),F(1))]*d;nodes={};active={'':(coefficients(poly,d),root)}
    while len(active)<max_leaves:
        path=min(active,key=lambda k:min(active[k][0].values()))
        coeff,box=active.pop(path)
        # The axis with the largest Bernstein first-difference is a proposal only.
        scores=[]
        for axis in range(d):
            scores.append(max(abs(coeff[tuple(k[j]+int(j==axis) for j in range(d))]-v)
                              for k,v in coeff.items() if k[axis]<2))
        axis=max(range(d),key=scores.__getitem__)
        if not scores[axis]:active[path]=(coeff,box);break
        a,b=split(coeff,d,axis);mid=sum(box[axis])/2
        left=list(box);right=list(box);left[axis]=(box[axis][0],mid);right[axis]=(mid,box[axis][1])
        nodes[path]={'axis':axis};active[path+'0']=(a,left);active[path+'1']=(b,right)
    def tree(path):
        if path in active:return {'lower':str(min(active[path][0].values()))}
        return {**nodes[path],'left':tree(path+'0'),'right':tree(path+'1')}
    return min(min(v[0].values()) for v in active.values()),tree('')


def verify_refinement(poly,d,tree):
    leaves=[]
    def visit(node,box,depth):
        if depth>32:raise ValueError('Refinement depth exceeded')
        if set(node)=={'lower'}:
            # Independent direct affine monomial expansion, not de Casteljau.
            low=min(coefficients(affine_change(poly,box),d).values())
            if str(low)!=node['lower']:raise ValueError('Incorrect Bernstein coefficient bound')
            leaves.append(low)
            if len(leaves)>128:raise ValueError('Refinement budget exceeded')
        elif set(node)=={'axis','left','right'}:
            axis=node['axis']
            if type(axis) is not int or not 0<=axis<d:raise ValueError('Invalid refinement split')
            mid=sum(box[axis])/2;left=list(box);right=list(box)
            left[axis]=(box[axis][0],mid);right[axis]=(mid,box[axis][1])
            visit(node['left'],left,depth+1);visit(node['right'],right,depth+1)
        else:raise ValueError('Invalid refinement tree')
    visit(tree,[(F(0),F(1))]*d,0)
    return min(leaves)
