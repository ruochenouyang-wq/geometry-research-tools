"""Exact minimization of a quadratic on a rational box of dimension <= 2."""
from fractions import Fraction as F
from itertools import product


def value(constant,linear,hessian,x):
    return F(constant)+sum((F(a)*F(b) for a,b in zip(linear,x)),F(0))+sum(
        (F(hessian[i][j])*F(x[i])*F(x[j])/2 for i in range(len(x)) for j in range(len(x))),F(0))


def minimum(constant,linear,hessian,box):
    """Enumerate relative-interior stationary points of all faces.

    Singular stationary sets reach a lower-dimensional face with equal value.
    Hence their omission here cannot omit every global minimizer.
    """
    d=len(box)
    if d not in (1,2) or len(linear)!=d or len(hessian)!=d or any(len(row)!=d for row in hessian):
        raise ValueError('Quadratic dimension must be 1 or 2')
    if any(hessian[i][j]!=hessian[j][i] for i in range(d) for j in range(d)):
        raise ValueError('Symmetric Hessian required')
    if any(lo>hi for lo,hi in box):raise ValueError('Reversed box')
    candidates=[]
    for state in product((-1,0,1),repeat=d):
        free=[i for i,s in enumerate(state) if s==0]
        x=[box[i][s==1] if s else F(0) for i,s in enumerate(state)]
        rhs=[-linear[i]-sum((hessian[i][j]*x[j] for j in range(d) if j not in free),F(0)) for i in free]
        if len(free)==1:
            i=free[0]
            if not hessian[i][i]:continue
            x[i]=rhs[0]/hessian[i][i]
        elif len(free)==2:
            a,b=hessian[0];c,e=hessian[1];det=a*e-b*c
            if not det:continue
            x=[(e*rhs[0]-b*rhs[1])/det,(a*rhs[1]-c*rhs[0])/det]
        if all(lo<=z<=hi for z,(lo,hi) in zip(x,box)):
            candidates.append((value(constant,linear,hessian,x),tuple(x)))
    return min(candidates)


def on_unit_box(constant,linear,hessian,box):
    origin=[lo for lo,hi in box];width=[hi-lo for lo,hi in box];d=len(box)
    c=value(constant,linear,hessian,origin)
    l=[width[i]*(linear[i]+sum((hessian[i][j]*origin[j] for j in range(d)),F(0))) for i in range(d)]
    h=[[width[i]*width[j]*hessian[i][j] for j in range(d)] for i in range(d)]
    return c,l,h
