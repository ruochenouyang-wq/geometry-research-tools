"""V30: rank-one matrix domination reduces all-amplitude proof synthesis to 1D.

For fixed SPD W, the best delta in (1-t^2)r'(t)r'(t)^T <= delta W is
max_t (1-t^2)r'(t)^T W^-1 r'(t). Search and witness bound this constant.
"""
from fractions import Fraction as F
import family_v29 as family


def synthesize(directions,weight=None,tolerance=F(1,10**6),max_leaves=64):
    d=len(directions)
    if weight is None:weight=[[F(i==j) for j in range(d)] for i in range(d)]
    dirs,w,_,_,g=family.assemble(directions,weight)
    return family.build(dirs,w,family.maximum(g,tolerance,max_leaves),'poisson_matrix_v30')

verify=family.verify
