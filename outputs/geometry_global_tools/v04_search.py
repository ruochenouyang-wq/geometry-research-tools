"""V4: inexpensive endpoint proposals followed by exact certification."""
from fractions import Fraction as F
from math import ceil,log2,nextafter,inf,isfinite
from time import perf_counter
import v03_tail as v3


class FloatProposal:
    """Never used to accept a theorem: every result is only a proposal."""
    def __init__(self,kernel):
        self.n=kernel.n
        self.a=[[float(x) for x in row] for row in kernel.a]
        self.mass=list(map(float,kernel.mass))
        self.qlo,self.qhi,self.beta=map(float,(kernel.qlo,kernel.qhi,kernel.beta))
        self.policy,self.width=kernel.policy,kernel.degree
        self.blocks=getattr(kernel,'blocks',None)
        self.terms=[]
        for j,mass,column in kernel.couplings:
            entries=[(i,l,float(c*d/mass)) for i,c in column for l,d in column]
            low=self.beta if kernel.policy=='scalar_ritz' else j*(j+1)+self.qlo
            self.terms.append((low,j*(j+1)+self.qhi,entries))

    def count(self,x,kind):
        a=[row[:] for row in self.a]
        for i in range(self.n):a[i][i]-=x*self.mass[i]
        if kind!='upper_ritz':
            for low,high,entries in self.terms:
                denominator=(high if kind=='upper_tail' else low)-x
                if denominator<=0:return self.n
                for i,j,value in entries:a[i][j]-=value/denominator
        negative=0
        blocks=self.blocks if self.blocks is not None else [list(range(self.n))]
        width=(self.width+1)//2 if self.blocks is not None else self.width
        for indices in blocks:
            block=[[a[i][j] for j in indices] for i in indices]
            for k in range(len(block)):
                d=block[k][k]
                # Proposal only. Exact verification never replaces zero.
                if abs(d)<1e-30:d=-1e-30
                negative+=d<0
                stop=min(len(block),k+width+1)
                for i in range(k+1,stop):
                    for j in range(i,stop):
                        block[i][j]-=block[i][k]*block[k][j]/d
                        block[j][i]=block[i][j]
        return negative

    def root(self,k,kind):
        distance=1.0
        lo=self.qlo-distance
        while self.count(lo,kind)>=k:
            distance*=2
            lo=self.qlo-distance
            if not isfinite(lo):raise ValueError('Float proposal overflow')
        hi=self.qhi+(k-1)*k+1
        if kind!='upper_ritz':hi=min(hi,nextafter(self.beta,-inf))
        if kind=='upper_tail' and self.count(hi,kind)<k:
            return self.root(k,'upper_ritz'),'upper_ritz'
        for _ in range(64):
            mid=(lo+hi)/2
            if mid in (lo,hi):break
            if self.count(mid,kind)>=k:hi=mid
            else:lo=mid
        return (lo+hi)/2


def proposals(kernel,k):
    f=FloatProposal(kernel)
    low=f.root(k,'lower')
    kind='upper_tail' if kernel.policy=='sandwich' else 'upper_ritz'
    up=f.root(k,kind)
    if isinstance(up,tuple):up,kind=up
    return low,up,kind


def guided(q,k=1,modes=8,tolerance=F(1,10**10),policy='sandwich',range_proof=None,kernel_factory=v3.Kernel):
    v3.sizes(k,modes)
    tolerance=F(tolerance)
    if not F(1,10**40)<=tolerance<=1:raise ValueError('Tolerance must be in [1e-40,1]')
    kernel=kernel_factory(q,modes,policy,range_proof)
    start=perf_counter()
    try:
        cl,cu,kind=proposals(kernel,k)
        if not all(isfinite(x) for x in (cl,cu)):raise ValueError('Invalid float proposal')
        cl,cu=F(cl),F(cu)
        guard=max(tolerance/16,(1+max(abs(cl),abs(cu)))/2**48)
        for _ in range(12):
            lo,hi=cl-guard,cu+guard
            if kind=='upper_tail' and hi>=kernel.beta:
                kind='upper_ritz';cu=F(FloatProposal(kernel).root(k,kind));hi=cu+guard
            if lo<=hi and kernel.lower_holds(k,lo) and kernel.upper_holds(k,hi,kind):break
            guard*=2
        else:raise ValueError('Proposals failed exact checks')
        path='guided'
        if hi-lo>tolerance:
            inner_low=min(cl+guard,kernel.beta)
            inner_up=cu-guard
            low_fails=not kernel.lower_holds(k,inner_low)
            up_fails=not kernel.upper_holds(k,inner_up,kind)
            if low_fails and up_fails and inner_up-inner_low>tolerance:
                path='tail_bound_limited'
            else:
                right=inner_low
                while right-lo>tolerance/16:
                    mid=(lo+right)/2
                    if kernel.lower_holds(k,mid):lo=mid
                    else:right=mid
                left=inner_up
                while hi-left>tolerance/16:
                    mid=(left+hi)/2
                    if kernel.upper_holds(k,mid,kind):hi=mid
                    else:left=mid
                path='exact_local_refinement'
        cert=v3.certificate(kernel,k,lo,hi,kind)
        if not v3.verify(cert):raise ArithmeticError('Guided certificate rejected')
        return {'certificate':cert,'target_met':hi-lo<=tolerance,'path':path,
                'exact_count_calls':kernel.exact_calls,'seconds':perf_counter()-start}
    except (ValueError,OverflowError,ZeroDivisionError):
        bits=min(160,max(48,ceil(-log2(float(tolerance)))+16))
        cert=v3.certify(q,k,modes,bits,policy,range_proof)
        return {'certificate':cert,'target_met':F(cert['exact_width'])<=tolerance,
                'path':'rational_fallback','exact_count_calls':None,'seconds':perf_counter()-start}


def search(q,k=1,tolerance=F(1,10**10),max_modes=32,policy='sandwich',range_proof=None,kernel_factory=v3.Kernel):
    v3.sizes(k,max_modes)
    if max_modes<max(4,k+2):raise ValueError('Budget needs at least max(4,k+2) modes')
    start=perf_counter();attempts=[];best=None
    sizes=list(range(max(4,k+2),max_modes+1,2))
    if sizes[-1]!=max_modes:sizes.append(max_modes)
    for n in sizes:
        result=guided(q,k,n,tolerance,policy,range_proof,kernel_factory)
        attempts.append({'modes':n,**{key:result[key] for key in ['path','exact_count_calls','seconds']},
                         'exact_width':result['certificate']['exact_width']})
        if best is None or F(result['certificate']['exact_width'])<F(best['certificate']['exact_width']):best=result
        if result['target_met']:break
    return {'algorithm_version':4,'status':'target_met' if best['target_met'] else 'target_not_met',
            'tolerance':str(F(tolerance)),'certificate':best['certificate'],'attempts':attempts,
            'seconds':perf_counter()-start}
