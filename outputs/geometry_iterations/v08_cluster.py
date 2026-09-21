"""V8: exact residual and projector bounds for the FIRST m eigenfunctions."""
from fractions import Fraction as F
from math import sqrt, copysign
from time import perf_counter
import v07_precision as v7
v3,v4,v5=v7.v3,v7.v4,v7.v5


def trial_basis(q,n,m):
    """Untrusted cyclic Jacobi proposal, returning columns in Legendre coordinates."""
    data=v3.base.assemble(q,n);mass=list(map(float,data['M']))
    a=[[float(data['A'][i][j])/sqrt(mass[i]*mass[j]) for j in range(n)] for i in range(n)]
    vectors=[[float(i==j) for j in range(n)] for i in range(n)]
    for _ in range(60):
        if max((abs(a[i][j]) for i in range(n) for j in range(i+1,n)),default=0)<1e-14:break
        for p in range(n):
            for r in range(p+1,n):
                apq=a[p][r]
                if abs(apq)<1e-16:continue
                tau=(a[r][r]-a[p][p])/(2*apq)
                t=copysign(1.0,tau)/(abs(tau)+sqrt(1+tau*tau));c=1/sqrt(1+t*t);s=t*c
                app,arr=a[p][p],a[r][r]
                for i in range(n):
                    if i not in (p,r):
                        aip,air=a[i][p],a[i][r]
                        a[i][p]=a[p][i]=c*aip-s*air
                        a[i][r]=a[r][i]=s*aip+c*air
                    vip,vir=vectors[i][p],vectors[i][r]
                    vectors[i][p],vectors[i][r]=c*vip-s*vir,s*vip+c*vir
                a[p][p],a[r][r],a[p][r],a[r][p]=app-t*apq,arr+t*apq,0.0,0.0
    indices=sorted(range(n),key=lambda j:a[j][j])[:m];columns=[]
    for j in indices:
        col=[vectors[i][j]/sqrt(mass[i]) for i in range(n)];norm=max(map(abs,col))
        columns.append([F(round(x/norm*2**44),2**44) for x in col])
    return columns


def statistics(q,columns):
    if not isinstance(columns,list) or not 1<=len(columns)<=4:raise ValueError('Cluster size must be 1..4')
    columns=[v7.coefficients([str(x) for x in col]) for col in columns]
    n=len(columns[0]);m=len(columns)
    if any(len(col)!=n for col in columns) or m>n:raise ValueError('Inconsistent basis')
    mass=[F(2,2*i+1) for i in range(n+len(q)-1)]
    gram=[[sum((mass[i]*u[i]*w[i] for i in range(n)),F(0)) for w in columns] for u in columns]
    if v3.base.inertia(gram)!=[0,0,m]:raise ValueError('Trial basis must have full rank')
    actions=[]
    for col in columns:
        h=[F(0)]*len(mass)
        for j,c in enumerate(col):
            h[j]+=j*(j+1)*c
            for i,value in v3.base.q_times_p(q,j).items():h[i]+=c*value
        actions.append(h)
    ritz=[[sum((mass[i]*u[i]*h[i] for i in range(n)),F(0)) for h in actions] for u in columns]
    projected=[v7.solve(gram,[ritz[i][j] for i in range(m)]) for j in range(m)]
    residuals=[[h[i]-(sum((columns[l][i]*projected[j][l] for l in range(m)),F(0)) if i<n else 0)
                for i in range(len(mass))] for j,h in enumerate(actions)]
    rg=[[sum((mass[i]*u[i]*w[i] for i in range(len(mass))),F(0)) for w in residuals] for u in residuals]
    residual_squared=sum((v7.solve(gram,[rg[i][j] for i in range(m)])[j] for j in range(m)),F(0))
    if residual_squared<0:raise ArithmeticError('Negative residual norm')
    return gram,ritz,residual_squared


def ritz_upper(q,n,gram,ritz):
    proof=v5.make_range(q);lo,hi=v3.range_bounds(q,proof);hi+=(n-1)*n+1;lo-=1
    m=len(gram)
    for _ in range(40):
        mid=(lo+hi)/2
        counts=v3.base.inertia([[ritz[i][j]-mid*gram[i][j] for j in range(m)] for i in range(m)])
        if counts[0]+counts[1]==m:hi=mid
        else:lo=mid
    return hi


def certificate(q,columns,alpha,gap_certificate):
    q=v3.base.potential([str(x) for x in q]);alpha=v3.base.rational(str(alpha))
    gram,ritz,rho2=statistics(q,columns);m=len(columns)
    gap=gap_certificate
    if not v3.verify(gap) or gap['eigenvalue_index']!=m+1 or v3.base.potential(gap['q_coefficients'])!=q:
        raise ValueError('External gap must bound the (m+1)-th eigenvalue of this operator')
    gamma=F(gap['lower'])
    counts=v3.base.inertia([[ritz[i][j]-alpha*gram[i][j] for j in range(m)] for i in range(m)])
    if counts[0]+counts[1]!=m or gamma<=alpha:raise ValueError('No certified external cluster gap')
    bound=min(F(2*m),2*rho2/(gamma-alpha)**2)
    return {'method':'first_cluster_projector_v8','model':v3.base.MODEL,'scope':v3.base.SCOPE,
            'cluster_indices':list(range(1,m+1)),
            'q_coefficients':[str(x) for x in q],
            'trial_columns':[[str(F(x)) for x in col] for col in columns],
            'ritz_upper':str(alpha),'external_lower':str(gamma),'gap_certificate':gap,
            'full_residual_frobenius_squared':str(rho2),
            'projector_hilbert_schmidt_squared_upper':str(bound),
            'nontrivial_bound':bound<2*m,'formal_assistant_checked':False}


def verify(cert):
    try:
        return cert==certificate(cert['q_coefficients'],cert['trial_columns'],cert['ritz_upper'],cert['gap_certificate'])
    except (KeyError,TypeError,ValueError,ZeroDivisionError,IndexError):return False


def cluster(q,m=2,modes=16):
    if type(m) is not int or not 1<=m<=4:raise ValueError('Cluster size must be 1..4')
    v3.base.check_sizes(m+1,modes)
    start=perf_counter();q=v3.base.potential([str(x) for x in q])
    cols=trial_basis(q,modes,m);gram,ritz,_=statistics(q,cols)
    alpha=ritz_upper(q,modes,gram,ritz)
    gap=v4.guided(q,m+1,modes,F(1,10**8),range_proof=v5.make_range(q))['certificate']
    cert=certificate(q,cols,alpha,gap)
    if not verify(cert):raise ArithmeticError('Cluster certificate rejected')
    return {'algorithm_version':8,'modes':modes,'cluster_size':m,'certificate':cert,'seconds':perf_counter()-start}
