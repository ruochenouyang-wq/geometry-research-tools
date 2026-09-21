"""Small exact linear/quadratic algorithms used by the proof transformations."""
from fractions import Fraction as F
from itertools import product
from v07_precision import solve
from spectral_certifier import inertia
from global_quadratic import value


def transpose(a):return [list(row) for row in zip(*a)]
def matmul(a,b):return [[sum((x*y for x,y in zip(row,col)),F(0)) for col in zip(*b)] for row in a]
def matvec(a,x):return [sum((b*c for b,c in zip(row,x)),F(0)) for row in a]
def inverse(a):return transpose([solve(a,[F(i==j) for i in range(len(a))]) for j in range(len(a))])


def rectangular(a,b,n):
    """Return one exact solution and pivot columns; reject inconsistent systems."""
    if len(a)!=len(b) or any(len(row)!=n for row in a):raise ValueError('Linear dimensions')
    rows=[list(map(F,row))+[F(y)] for row,y in zip(a,b)];pivots=[];r=0
    for col in range(n):
        k=next((k for k in range(r,len(rows)) if rows[k][col]),None)
        if k is None:continue
        rows[r],rows[k]=rows[k],rows[r];v=rows[r][col];rows[r]=[z/v for z in rows[r]]
        for i in range(len(rows)):
            if i!=r:
                v=rows[i][col];rows[i]=[x-v*y for x,y in zip(rows[i],rows[r])]
        pivots.append(col);r+=1
    if any(not any(row[:n]) and row[n] for row in rows):raise ValueError('Inconsistent linear system')
    x=[F(0)]*n
    for i,j in enumerate(pivots):x[j]=rows[i][n]
    return x,pivots


def independent_rows(a):
    if not a:return []
    return rectangular(transpose(a),[0]*len(a[0]),len(a))[1]


def quadratic_min(c,l,h,box=None):
    d=len(l);c=F(c);l=list(map(F,l));h=[list(map(F,row)) for row in h]
    if not 0<=d<=4 or len(h)!=d or any(len(r)!=d for r in h):raise ValueError('Quadratic size')
    if any(h[i][j]!=h[j][i] for i in range(d) for j in range(d)):raise ValueError('Symmetric Hessian required')
    if not d:return c,[]
    if box is None:
        if inertia(h)[0]:raise ValueError('Quadratic unbounded below')
        try:x,_=rectangular(h,[-z for z in l],d)
        except ValueError:raise ValueError('Linear term in null direction: unbounded below')
        return value(c,l,h,x),x
    if len(box)!=d or any(len(r)!=2 or r[0]>r[1] for r in box):raise ValueError('Invalid box')
    candidates=[]
    for state in product((-1,0,1),repeat=d):
        free=[i for i,s in enumerate(state) if not s];fixed=[i for i,s in enumerate(state) if s]
        x=[F(box[i][int(s==1)]) if s else F(0) for i,s in enumerate(state)]
        try:
            z=solve([[h[i][j] for j in free] for i in free],[-l[i]-sum(h[i][j]*x[j] for j in fixed) for i in free]) if free else []
        except ValueError:continue
        for i,v in zip(free,z):x[i]=v
        if all(lo<=v<=hi for v,(lo,hi) in zip(x,box)):candidates.append((value(c,l,h,x),x))
    return min(candidates,key=lambda row:row[0])
