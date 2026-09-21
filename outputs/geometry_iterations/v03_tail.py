"""V3: mode-dependent two-sided Schur bounds, no internal spectral gap.

All certificate arithmetic is rational. See V03.md for the full proof.
"""
from fractions import Fraction as F
import error_bounds as old
base = old.base

POLICIES = ('scalar_ritz','graded_ritz','sandwich')


def sizes(k,n):
    if type(k) is not int or type(n) is not int or not 1 <= k <= n <= 64:
        raise ValueError('Require 1 <= k <= modes <= 64')


def coefficient_range(q):
    r = sum((abs(x) for x in q[1:]),F(0))
    return {'kind':'coefficient','lower':str(q[0]-r),'upper':str(q[0]+r)}


def range_bounds(q,proof):
    if not isinstance(proof,dict): raise ValueError('Range proof required')
    if proof.get('kind') == 'coefficient':
        if proof != coefficient_range(q): raise ValueError('Incorrect coefficient bound')
        return F(proof['lower']),F(proof['upper'])
    if proof.get('kind') in ('quadratic_range_v5','adaptive_bernstein_v5'):
        from v05_range import verify_bounds
        return verify_bounds(q,proof)
    raise ValueError('Unsupported range proof')


class Kernel:
    def __init__(self,q,modes,policy='sandwich',range_proof=None):
        sizes(1,modes)
        if policy not in POLICIES: raise ValueError('Unknown tail policy')
        self.q = base.potential([str(x) for x in q])
        self.n,self.policy = modes,policy
        self.range_proof = coefficient_range(self.q) if range_proof is None else range_proof
        self.qlo,self.qhi = range_bounds(self.q,self.range_proof)
        data = base.assemble(self.q,modes)
        self.a,self.mass = data['A'],data['M']
        self.beta = modes*(modes+1)+self.qlo
        self.degree = len(self.q)-1
        self.couplings = []
        for j in range(modes,modes+self.degree):
            column = [(i,self.mass[i]*c) for i,c in base.q_times_p(self.q,j).items() if i < modes and c]
            self.couplings.append((j,F(2,2*j+1),column))
        self.exact_calls = 0

    def matrix(self,x,kind):
        x = F(x)
        if kind not in ('lower','upper_ritz','upper_tail'): raise ValueError('Unknown matrix kind')
        if kind != 'upper_ritz' and x >= self.beta:
            raise ValueError('Tail positivity not established')
        a = [row[:] for row in self.a]
        for i in range(self.n): a[i][i] -= x*self.mass[i]
        if kind == 'upper_ritz': return a
        if kind == 'upper_tail' and self.policy != 'sandwich':
            raise ValueError('Tail upper bound disabled by this policy')
        for j,mass,column in self.couplings:
            energy = j*(j+1)+(self.qhi if kind=='upper_tail' else self.qlo)
            if kind=='lower' and self.policy=='scalar_ritz': energy = self.beta
            denom = mass*(energy-x)
            for i,c in column:
                for l,d in column: a[i][l] -= c*d/denom
        return a

    def count(self,x,kind,dense=False):
        self.exact_calls += 1
        a = self.matrix(x,kind)
        return base.inertia(a) if dense else old.banded_inertia(a)

    def lower_holds(self,k,x):
        return x < self.beta and self.count(x,'lower')[0] < k

    def upper_holds(self,k,x,kind):
        if kind=='upper_tail' and x >= self.beta: return False
        c = self.count(x,kind)
        return c[0]+c[1] >= k


def certificate(kernel,k,lo,hi,upper_kind,dense=False):
    sizes(k,kernel.n)
    if upper_kind not in ('upper_ritz','upper_tail'):
        raise ValueError('An upper endpoint requires an upper comparison matrix')
    lo,hi = F(lo),F(hi)
    if lo > hi: raise ValueError('Reversed enclosure')
    low_count = kernel.count(lo,'lower',dense)
    up_count = kernel.count(hi,upper_kind,dense)
    if low_count[0] >= k or up_count[0]+up_count[1] < k:
        raise ValueError('Endpoint failed its exact spectral count')
    return {'method':'modewise_schur_v3','model':base.MODEL,'scope':base.SCOPE,
            'q_coefficients':[str(x) for x in kernel.q],'modes':kernel.n,'eigenvalue_index':k,
            'tail_policy':kernel.policy,'range_proof':kernel.range_proof,'upper_kind':upper_kind,
            'tail_lower':str(kernel.beta),'lower_inertia':low_count,'upper_inertia':up_count,
            'formal_assistant_checked':False,**old.enclosure_fields(lo,hi)}


def verify(cert,independent=False):
    try:
        if not isinstance(cert,dict): return False
        kernel = Kernel(cert['q_coefficients'],cert['modes'],cert['tail_policy'],cert['range_proof'])
        expected = certificate(kernel,cert['eigenvalue_index'],base.rational(cert['lower']),
                               base.rational(cert['upper']),cert['upper_kind'],dense=independent)
        return expected == cert
    except (KeyError,ValueError,TypeError,ZeroDivisionError):
        return False


def bisect_upper(kernel,k,lo,hi,kind,bits):
    if not kernel.upper_holds(k,hi,kind): raise ValueError('Upper bracket not certified')
    for _ in range(bits):
        mid = (lo+hi)/2
        if kernel.upper_holds(k,mid,kind): hi = mid
        else: lo = mid
    return hi


def certify(q,k=1,modes=8,bits=48,policy='sandwich',range_proof=None):
    sizes(k,modes)
    if type(bits) is not int or not 8 <= bits <= 160: raise ValueError('Bits must be in 8..160')
    kernel = Kernel(q,modes,policy,range_proof)
    low,high = kernel.qlo-1,kernel.qhi+(k-1)*k+1
    upper_kind = 'upper_ritz'
    if policy=='sandwich' and high < kernel.beta: upper_kind = 'upper_tail'
    up = bisect_upper(kernel,k,low,high,upper_kind,bits)
    if policy=='sandwich' and upper_kind=='upper_ritz' and up < kernel.beta:
        upper_kind = 'upper_tail'
        up = bisect_upper(kernel,k,low,up,upper_kind,bits)
    distance = F(1)
    lo = kernel.qlo-distance
    while not kernel.lower_holds(k,lo):
        distance *= 2
        lo = kernel.qlo-distance
    hi = up
    for _ in range(bits):
        mid = (lo+hi)/2
        if kernel.lower_holds(k,mid): lo = mid
        else: hi = mid
    cert = certificate(kernel,k,lo,up,upper_kind)
    if not verify(cert): raise ArithmeticError('V3 certificate rejected')
    return cert
