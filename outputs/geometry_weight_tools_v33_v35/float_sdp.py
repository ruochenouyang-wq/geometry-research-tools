"""Untrusted small dense log-barrier candidate search (standard library only).

The final acceptance happens in dual_v33 and family_v29 using rational arithmetic.
"""
import math


def inv_logdet(a):
    n=len(a);l=[[0.0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            value=a[i][j]-sum(l[i][k]*l[j][k] for k in range(j))
            if i==j:
                if not value>0 or not math.isfinite(value):raise ValueError('Floating candidate not positive definite')
                l[i][j]=math.sqrt(value)
            else:l[i][j]=value/l[j][j]
    inv=[]
    for k in range(n):
        y=[]
        for i in range(n):y.append(((1.0 if i==k else 0.0)-sum(l[i][j]*y[j] for j in range(i)))/l[i][i])
        x=[0.0]*n
        for i in reversed(range(n)):x[i]=(y[i]-sum(l[j][i]*x[j] for j in range(i+1,n)))/l[i][i]
        inv.append(x)
    return [[(inv[i][j]+inv[j][i])/2 for j in range(n)] for i in range(n)],2*sum(math.log(l[i][i]) for i in range(n))


def linear(a,b):
    rows=[row[:]+[v] for row,v in zip(a,b)];n=len(b)
    for j in range(n):
        k=max(range(j,n),key=lambda i:abs(rows[i][j]))
        if abs(rows[k][j])<1e-28:raise ValueError('Singular floating Newton matrix')
        rows[j],rows[k]=rows[k],rows[j];v=rows[j][j]
        rows[j]=[z/v for z in rows[j]]
        for i in range(j+1,n):
            v=rows[i][j];rows[i]=[x-v*y for x,y in zip(rows[i],rows[j])]
    x=[0.0]*n
    for i in reversed(range(n)):x[i]=rows[i][-1]-sum(rows[i][j]*x[j] for j in range(i+1,n))
    return x


def solve(samples,cost,final_mu,max_steps=240):
    d=len(cost);indices=[(i,j) for i in range(d) for j in range(i,d)];n=len(indices)
    def unpack(x):
        g=[[0.0]*d for _ in range(d)]
        for value,(i,j) in zip(x,indices):g[i][j]=g[j][i]=value
        return g
    objective=[cost[i][j]*(1 if i==j else 2) for i,j in indices]
    alpha=max(1.0,2*max(sum(aa[i][i] for i in range(d)) for aa in samples));x=[alpha if i==j else 0.0 for i,j in indices]
    def evaluate(x,mu,derivatives=False):
        g=unpack(x);value=sum(u*v for u,v in zip(objective,x));gradient=objective[:];hessian=[[0.0]*n for _ in range(n)]
        for aa in samples:
            inv,logdet=inv_logdet([[g[i][j]-aa[i][j] for j in range(d)] for i in range(d)])
            value-=mu*logdet
            if derivatives:
                for k,(i,j) in enumerate(indices):gradient[k]-=mu*inv[i][j]*(1 if i==j else 2)
                for k,(i,j) in enumerate(indices):
                    for l,(u,v) in enumerate(indices):
                        if i==j and u==v:z=inv[i][u]**2
                        elif i==j:z=2*inv[i][u]*inv[i][v]
                        elif u==v:z=2*inv[i][u]*inv[j][u]
                        else:z=2*(inv[i][u]*inv[j][v]+inv[i][v]*inv[j][u])
                        hessian[k][l]+=mu*z
        return value,gradient,hessian
    mu=max(0.1,final_mu);steps=0;stages=0;converged=False
    while steps<max_steps:
        stages+=1
        for _ in range(40):
            if steps>=max_steps:break
            value,gradient,hessian=evaluate(x,mu,True);direction=linear(hessian,[-z for z in gradient]);slope=sum(u*v for u,v in zip(gradient,direction));steps+=1
            if -slope/2<max(1e-18,mu*1e-10):break
            rate=1.0;accepted=False
            for _ in range(60):
                candidate=[z+rate*v for z,v in zip(x,direction)]
                try:new_value=evaluate(candidate,mu)[0]
                except ValueError:new_value=float('inf')
                if new_value<=value+0.01*rate*slope:
                    x=candidate;accepted=True;break
                rate*=0.5
            if not accepted:break
        if mu<=final_mu:
            converged=True;break
        mu=max(final_mu,mu/5)
    return unpack(x),mu,{'newton_steps':steps,'barrier_stages':stages,'reached_final_barrier':converged}
