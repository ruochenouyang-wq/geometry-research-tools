"""V6: exact parity-block counts for every global eigenvalue index."""
from fractions import Fraction as F
from time import perf_counter
import v03_tail as v3
import v04_search as v4
import v05_range as v5


class ParityKernel(v3.Kernel):
    def __init__(self,q,modes,policy='sandwich',range_proof=None):
        super().__init__(q,modes,policy,range_proof)
        if any(self.q[i] for i in range(1,len(self.q),2)):
            raise ValueError('Parity decomposition requires exact reflection symmetry')
        self.blocks=[list(range(0,modes,2)),list(range(1,modes,2))]

    def count(self,x,kind,dense=False):
        self.exact_calls+=1
        a=self.matrix(x,kind)
        if dense:return v3.base.inertia(a)
        if any(a[i][j] for i in self.blocks[0] for j in self.blocks[1]):
            raise ArithmeticError('Claimed invariant blocks are coupled')
        counts=[0,0,0]
        for indices in self.blocks:
            current=v3.old.banded_inertia([[a[i][j] for j in indices] for i in indices])
            counts=[x+y for x,y in zip(counts,current)]
        return counts


def search(q,k=1,tolerance=F(1,10**10),max_modes=32):
    start=perf_counter()
    q=v3.base.potential([str(x) for x in q])
    if any(q[i] for i in range(1,len(q),2)):
        raise ValueError('Odd coefficients must be exactly zero')
    proof=v5.make_range(q)
    result=v4.search(q,k,tolerance,max_modes,range_proof=proof,kernel_factory=ParityKernel)
    result['algorithm_version']=6;result['seconds']=perf_counter()-start
    return result
