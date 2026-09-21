"""Exact quadratic minimum over the convex hull of at most eight vertices.

Enumerate every face. A singular stationary set has a direction of constant
value meeting its boundary; hence a minimum is retained on a smaller face.
"""
from fractions import Fraction as F
from itertools import combinations
from v07_precision import solve
from global_quadratic import value


def minimum(c,l,h,vertices):
    vertices=[list(map(F,v)) for v in vertices];n=len(vertices);d=len(l)
    if not 1<=n<=8 or any(len(v)!=d for v in vertices):raise ValueError('Invalid vertices')
    if len(h)!=d or any(len(row)!=d for row in h) or any(h[i][j]!=h[j][i] for i in range(d) for j in range(d)):
        raise ValueError('Symmetric quadratic required')
    candidates=[]
    for size in range(1,n+1):
        for face in combinations(range(n),size):
            a=vertices[face[0]];basis=[[vertices[k][i]-a[i] for i in range(d)] for k in face[1:]]
            gradient=[F(l[i])+sum(F(h[i][j])*a[j] for j in range(d)) for i in range(d)]
            rhs=[-sum(b[i]*gradient[i] for i in range(d)) for b in basis]
            matrix=[[sum(b[i]*F(h[i][j])*e[j] for i in range(d) for j in range(d)) for e in basis] for b in basis]
            try:z=solve(matrix,rhs) if basis else []
            except (ValueError,ZeroDivisionError):continue
            w=[1-sum(z,F(0))]+z
            if min(w)<0:continue
            x=[sum(w[j]*vertices[k][i] for j,k in enumerate(face)) for i in range(d)]
            weights=[F(0)]*n
            for k,b in zip(face,w):weights[k]=b
            candidates.append((value(c,l,h,x),x,weights))
    return min(candidates,key=lambda item:item[0])


def affine_interpolant(vertices,values):
    d=len(vertices[0])
    solution=solve([[F(1)]+list(v) for v in vertices],list(values))
    return solution[0],solution[1:]
