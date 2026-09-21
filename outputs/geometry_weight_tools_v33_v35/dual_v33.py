"""Exact primal/dual certificates for minimum-cost Poisson matrix envelopes.

P = inf trace(C G), G >= (1-t*t) v(t)v(t)^T for ALL -1<=t<=1.
Any finite PSD atoms Y_k with sum Y_k <= C give a GLOBAL lower bound.
Neither a floating optimization log nor its stopping flag is trusted.
"""
from fractions import Fraction as F
import algebra as a
import error_bounds as e
import family_v29 as family


def matrix(raw,d,positive=False):
    out=[list(map(F,row)) for row in raw]
    if len(out)!=d or any(len(row)!=d for row in out):raise ValueError('Matrix size')
    if any(out[i][j]!=out[j][i] for i in range(d) for j in range(d)):raise ValueError('Symmetric matrix required')
    signs=a.inertia(out)
    if signs[0] or (positive and signs[1]):raise ValueError('Positive (semi)definite matrix required')
    return out


def inner(x,y):return sum((x[i][j]*y[j][i] for i in range(len(x)) for j in range(len(x))),F(0))
def encode(x):return [list(map(str,row)) for row in x]


def setup(directions,cost=None):
    directions=[e.base.potential(list(map(str,p))) for p in directions];d=len(directions)
    if not 1<=d<=4:raise ValueError('One to four directions required')
    v=[e.derivative(family.poisson(p)[0]) for p in directions]
    n=max(map(len,v));rows=[p+[F(0)]*(n-len(p)) for p in v]
    if len(a.independent_rows(rows))!=d:raise ValueError('Directions must be independent modulo constants')
    c=matrix(cost if cost is not None else [[int(i==j) for j in range(d)] for i in range(d)],d,True)
    return directions,v,c


def at(derivatives,t):
    t=F(t)
    if not -1<=t<=1:raise ValueError('Dual sample outside sphere coordinate')
    values=[e.evaluate(p,t) for p in derivatives]
    return [[(1-t*t)*x*y for y in values] for x in values]


def lower_certificate(directions,cost,atoms):
    dirs,v,c=setup(directions,cost);d=len(dirs)
    if not isinstance(atoms,list) or not 1<=len(atoms)<=128:raise ValueError('One to 128 dual atoms required')
    total=[[F(0)]*d for _ in range(d)];low=F(0);normalized=[]
    for atom in atoms:
        t=F(atom['t']);y=matrix(atom['matrix'],d);aa=at(v,t)
        low+=inner(y,aa)
        for i in range(d):
            for j in range(d):total[i][j]+=y[i][j]
        normalized.append({'t':str(t),'matrix':encode(y)})
    slack=[[c[i][j]-total[i][j] for j in range(d)] for i in range(d)]
    matrix(slack,d)
    return {'method':'matrix_dual_lower_v33','directions':[list(map(str,p)) for p in dirs],
            'cost':encode(c),'atoms':normalized,'cost_slack':encode(slack),'lower_bound':str(low),
            'scope':'minimum trace(C G) over all matrix Poisson envelopes on the full interval [-1,1]',
            'formal_assistant_checked':False}


def verify_lower(cert):
    try:return cert==lower_certificate(cert['directions'],cert['cost'],cert['atoms'])
    except (ValueError,TypeError,KeyError,IndexError,ZeroDivisionError):return False


def certificate(template,dual,tolerance):
    tolerance=F(tolerance)
    if not F(1,10**10)<=tolerance<=1:raise ValueError('Optimization tolerance budget')
    if not family.verify(template) or not verify_lower(dual):raise ValueError('Invalid primal or dual proof')
    if template['directions']!=dual['directions']:raise ValueError('Primal and dual describe different directions')
    g=[list(map(F,row)) for row in template['quadratic_bound']];c=[list(map(F,row)) for row in dual['cost']]
    low=F(dual['lower_bound']);upper=inner(c,g)
    if low>upper:raise ValueError('Inconsistent primal and dual bounds')
    return {'method':'global_matrix_envelope_v33','template':template,'dual':dual,'lower_bound':str(low),
            'upper_bound':str(upper),'gap':str(upper-low),'tolerance':str(tolerance),
            'status':'global_gap_closed' if upper-low<=tolerance else 'global_gap_open',
            'objective':'trace(C G); G dominates (1-t^2) v(t)v(t)^T for all t in [-1,1]',
            'limitation':'Optimality is within the fixed Poisson exponential ansatz, not the sharp spectral inequality.',
            'formal_assistant_checked':False}


def verify(cert):
    try:return cert==certificate(cert['template'],cert['dual'],cert['tolerance'])
    except (ValueError,TypeError,KeyError,IndexError,ZeroDivisionError):return False


def repair_dual(directions,cost,atoms,accuracy=F(1,10**10)):
    """Scale proposed rational PSD atoms until their sum is <= C, EXACTLY."""
    dirs,_,c=setup(directions,cost);d=len(dirs);accuracy=F(accuracy)
    if not F(1,10**14)<=accuracy<=1:raise ValueError('Dual repair tolerance')
    ys=[matrix(atom['matrix'],d) for atom in atoms]
    if not 1<=len(ys)<=128:raise ValueError('Dual atom budget')
    total=[[sum(y[i][j] for y in ys) for j in range(d)] for i in range(d)]
    def feasible(scale):return a.inertia([[scale*c[i][j]-total[i][j] for j in range(d)] for i in range(d)])[0]==0
    lo=F(0);hi=F(1)
    for _ in range(128):
        if feasible(hi):break
        hi*=2
    else:raise ValueError('Dual scaling overflow budget')
    while hi-lo>accuracy:
        mid=(lo+hi)/2
        if feasible(mid):hi=mid
        else:lo=mid
    hi=max(hi,F(1))
    scaled=[{'t':str(F(atom['t'])),'matrix':encode([[z/hi for z in row] for row in y])} for atom,y in zip(atoms,ys)]
    return lower_certificate(dirs,c,scaled)
