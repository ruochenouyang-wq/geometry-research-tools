"""V7: rational Newton refinement of trial eigenfunctions."""
from fractions import Fraction as F
from time import perf_counter
import v03_tail as v3
import v04_search as v4
import v05_range as v5
e=v3.old


def solve(a,b):
    rows=[list(map(F,row))+[F(value)] for row,value in zip(a,b)]
    n=len(b)
    if len(rows)!=n or any(len(row)!=n+1 for row in rows):raise ValueError('Square system required')
    for k in range(n):
        pivot=max(range(k,n),key=lambda i:abs(rows[i][k]))
        rows[k],rows[pivot]=rows[pivot],rows[k]
        if not rows[k][k]:raise ValueError('Singular Newton system')
        for i in range(k+1,n):
            if not rows[i][k]:continue
            factor=rows[i][k]/rows[k][k]
            rows[i][k]=F(0)
            for j in range(k+1,n+1):rows[i][j]-=factor*rows[k][j]
    x=[F(0)]*n
    for i in reversed(range(n)):
        x[i]=(rows[i][-1]-sum((rows[i][j]*x[j] for j in range(i+1,n)),F(0)))/rows[i][i]
    return x


def coefficients(raw):
    if not isinstance(raw,list) or not 1<=len(raw)<=32:raise ValueError('At most 32 trial coefficients')
    c=[v3.base.rational(x) for x in raw]
    if any(abs(x)>10**6 or x.denominator>2**256 for x in c):raise ValueError('Trial coefficient exceeds budget')
    return c


def newton_step(q,c,bits=96):
    if type(bits) is not int or not 48<=bits<=224:raise ValueError('Coefficient bits must be in 48..224')
    n=len(c);data=v3.base.assemble(q,n);a,m=data['A'],data['M']
    pivot=max(range(n),key=lambda i:abs(c[i]))
    if not c[pivot]:raise ValueError('Zero function')
    c=[x/c[pivot] for x in c]
    mu,_,_=e.residual_statistics(q,c)
    residual=[sum((a[i][j]*c[j] for j in range(n)),F(0))-mu*m[i]*c[i] for i in range(n)]
    if not any(residual):return c
    free=[i for i in range(n) if i!=pivot]
    jac=[[a[i][j]-(mu*m[i] if i==j else 0) for j in free]+[-m[i]*c[i]] for i in range(n)]
    delta=solve(jac,[-x for x in residual])
    denominator=2**bits
    for j,idx in enumerate(free):c[idx]=F(round((c[idx]+delta[j])*denominator),denominator)
    return c


def gap_proof(q,c):
    mu,_,_=e.residual_statistics(q,c)
    even=not any(q[i] for i in range(1,len(q),2)) and not any(c[i] for i in range(1,len(c),2))
    range_proof=v5.make_range(q);lo,_=v3.range_bounds(q,range_proof)
    gamma=(6 if even else 2)+lo
    if gamma>mu:
        return {'kind':'range_minmax_v7','sector':'even' if even else 'all','range_proof':range_proof}
    if even:
        return {'kind':'even_schur_v7','proof':e.even_separation(q,max(4,len(c)),40)}
    g=v4.guided(q,2,max(4,len(c)),F(1,10**8),range_proof=range_proof)['certificate']
    return {'kind':'full_schur_v7','certificate':g}


def separation(q,c,proof):
    if not isinstance(proof,dict):raise ValueError('Gap proof required')
    kind=proof.get('kind')
    even=kind=='even_schur_v7' or (kind=='range_minmax_v7' and proof.get('sector')=='even')
    if even and (any(q[i] for i in range(1,len(q),2)) or any(c[i] for i in range(1,len(c),2))):
        raise ValueError('Parity condition failed')
    if kind=='range_minmax_v7':
        if proof.get('sector') not in ('all','even'):raise ValueError('Unknown sector')
        lo,_=v3.range_bounds(q,proof['range_proof']);return (6 if even else 2)+lo
    if kind=='even_schur_v7':return e.gap_value(q,proof['proof'])
    if kind=='full_schur_v7':
        cert=proof['certificate']
        if not v3.verify(cert) or cert['eigenvalue_index']!=2 or v3.base.potential(cert['q_coefficients'])!=q:
            raise ValueError('Incorrect full-space separation')
        return F(cert['lower'])
    raise ValueError('Unknown gap proof')


def certificate(q,c,proof):
    q=v3.base.potential([str(x) for x in q]);c=coefficients([str(x) for x in c])
    gamma=separation(q,c,proof);mu,inside,outside=e.residual_statistics(q,c)
    if gamma<=mu:raise ValueError('No certified spectral separation')
    rho2=inside+outside
    return {'method':'rational_trial_refinement_v7','model':v3.base.MODEL,'scope':v3.base.SCOPE,
            'eigenvalue_index':1,'q_coefficients':[str(x) for x in q],
            'trial_legendre_coefficients':[str(x) for x in c],'gap_proof':proof,
            'rayleigh_quotient':str(mu),'gap_lower':str(gamma),
            'inside_residual_squared':str(inside),'outside_residual_squared':str(outside),
            'excited_weight_upper':str(min(F(1),rho2/(gamma-mu)**2)),
            'formal_assistant_checked':False,**e.enclosure_fields(mu-rho2/(gamma-mu),mu)}


def verify(cert):
    try:
        if not isinstance(cert,dict):return False
        return cert==certificate(cert['q_coefficients'],cert['trial_legendre_coefficients'],cert['gap_proof'])
    except (KeyError,ValueError,TypeError,ZeroDivisionError):return False


def refine(q,modes=14,bits=96,steps=1):
    if type(steps) is not int or not 0<=steps<=3:raise ValueError('At most 3 Newton steps')
    if type(bits) is not int or not 48<=bits<=224:raise ValueError('Coefficient bits must be in 48..224')
    start=perf_counter();q=v3.base.potential([str(x) for x in q])
    even=not any(q[i] for i in range(1,len(q),2))
    c=e.polynomial_trial(q,modes,even=even)
    proof=gap_proof(q,c);before=certificate(q,c,proof)
    history=[{'step':0,'exact_width':before['exact_width'],'inside_residual_squared':before['inside_residual_squared']}]
    best=before
    for step in range(1,steps+1):
        c=newton_step(q,c,bits)
        cert=certificate(q,c,proof)
        history.append({'step':step,'exact_width':cert['exact_width'],'inside_residual_squared':cert['inside_residual_squared']})
        if F(cert['exact_width'])<F(best['exact_width']):best=cert
    if not verify(best):raise ArithmeticError('Refined function certificate rejected')
    return {'algorithm_version':7,'modes':modes,'coefficient_bits':bits,'steps':steps,
            'before_certificate':before,'certificate':best,'history':history,'seconds':perf_counter()-start}
