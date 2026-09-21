"""V19: exact affine equality elimination followed by all-real global search."""
from fractions import Fraction as F
from time import perf_counter
import global_core as core


def affine_space(matrix,rhs,n):
    if not isinstance(matrix,list) or not 1<=len(matrix)<=8 or len(rhs)!=len(matrix) or any(len(row)!=n for row in matrix):raise ValueError('Invalid equality dimensions')
    a=[list(map(F,row))+[F(b)] for row,b in zip(matrix,rhs)];pivots=[];r=0
    for col in range(n):
        pivot=next((k for k in range(r,len(a)) if a[k][col]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];scale=a[r][col];a[r]=[x/scale for x in a[r]]
        for k in range(len(a)):
            if k!=r:
                scale=a[k][col];a[k]=[x-scale*y for x,y in zip(a[k],a[r])]
        pivots.append(col);r+=1
    if any(not any(row[:n]) and row[n] for row in a):raise ValueError('Inconsistent equalities')
    free=[j for j in range(n) if j not in pivots]
    if not free:raise ValueError('This prototype requires at least one free parameter')
    origin=[F(0)]*n;basis=[[F(i==j) for j in free] for i in range(n)]
    for k,i in enumerate(pivots):
        origin[i]=a[k][n]
        for z,j in enumerate(free):basis[i][z]=-a[k][j]
    return origin,basis


def transform(p,matrix,rhs):
    p=core.normalized(p);n=len(p['directions']);origin,basis=affine_space(matrix,rhs,n);d=len(basis[0])
    k=[list(map(F,row)) for row in p['penalty']['hessian']];ell=list(map(F,p['penalty']['linear']))
    q0=core.at(p,origin);directions=[]
    for j in range(d):
        q=[F(0)]
        for i in range(n):q=core.v3.old.add(q,core.v3.old.scale(list(map(F,p['directions'][i])),basis[i][j]))
        directions.append(q)
    newell=[sum(basis[i][j]*(ell[i]+sum(k[i][z]*origin[z] for z in range(n))) for i in range(n)) for j in range(d)]
    newk=[[sum(basis[i][j]*k[i][z]*basis[z][h] for i in range(n) for z in range(n)) for h in range(d)] for j in range(d)]
    reduced=core.all_real_problem(q0,directions,newell,newk,core.penalty(p,origin))
    return reduced,origin,basis


def certificate(original,matrix,rhs,inner,independent=True):
    original=core.normalized(original);matrix=[list(map(str,map(F,row))) for row in matrix];rhs=list(map(str,map(F,rhs)))
    reduced,origin,basis=transform(original,matrix,rhs)
    if inner['problem']!=reduced or inner['domain']!='all_real' or not core.verify(inner,independent):raise ValueError('Invalid reduced problem proof')
    z=list(map(F,inner['candidate_parameters']));x=[a+sum(b*t for b,t in zip(row,z)) for a,row in zip(origin,basis)]
    if any(sum(F(c)*t for c,t in zip(row,x))!=F(b) for row,b in zip(matrix,rhs)):raise ValueError('Candidate violates equality')
    return {'method':'affine_equality_global_v19','original_problem':original,'equality_matrix':matrix,'equality_rhs':rhs,
            'origin':list(map(str,origin)),'basis':[list(map(str,row)) for row in basis],
            'inner_certificate':inner,'candidate_parameters':list(map(str,x)),
            **{key:inner[key] for key in ('global_lower','candidate_upper','global_gap','tolerance','status')},
            'optimization_scope':'all real parameters satisfying the declared equalities and all nonzero axisymmetric form-domain functions',
            'formal_assistant_checked':False}


def verify(c,independent=True):
    try:return c==certificate(c['original_problem'],c['equality_matrix'],c['equality_rhs'],c['inner_certificate'],independent)
    except (ValueError,KeyError,TypeError,IndexError,ZeroDivisionError):return False


def search(p,matrix,rhs,tolerance=F(1,1000),**options):
    start=perf_counter();reduced,_,_=transform(p,matrix,rhs)
    result=core.search(reduced,tolerance,version=19,domain='all_real',**options)
    cert=certificate(p,matrix,rhs,result['certificate'])
    return {'algorithm_version':19,'certificate':cert,'statistics':result['statistics'],'stop_reason':result['stop_reason'],'seconds':perf_counter()-start}
