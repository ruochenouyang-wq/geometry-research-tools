"""Exact original-function moments and infinite-tail spectral comparisons.

Finite radial matrices integrate q itself, not a polynomial approximation of q.
Q q P has a finite Gram representation from q^2 moments and Parseval. All
spectral endpoints use rational inertia; the complete m tail is retained.
"""
from copy import deepcopy
from functools import lru_cache
from math import comb, factorial
from bridge import F, same, base, arithmetic, wide, previous_functions as f
from bridge import previous_transfer as old

SECTOR = 'direct_axis_moment_sector_v1'
FULL = 'direct_axis_moment_full_v1'


def integer(x, lo, hi, name):
    if type(x) is not int or not lo <= x <= hi:
        raise ValueError(f'{name} must be an integer in {lo}..{hi}')
    return x


def normalize(function):
    result = f.normalize_function(function)
    if result['kind'] != 'axis_profile':
        raise ValueError('Direct moments require an axis profile')
    values = [F(result['amplitude']), F(result['offset'])]
    if result['profile'] == 'abs_power': values.append(F(result['exponent']))
    if any(abs(v) > 1000 or v.denominator > 10**6 for v in values):
        raise ValueError('Profile rational input budget exceeded')
    return result


def moment(function, degree, power=1):
    """Exact integral over dt on [-1,1]; power is 0, 1, or 2."""
    q = normalize(function)
    integer(degree, 0, 256, 'moment_degree'); integer(power, 0, 2, 'power')
    basic = F(0) if degree % 2 else F(2, degree+1)
    if power == 0: return basic
    a, b = F(q['amplitude']), F(q['offset'])
    if q['profile'] == 'step':
        first = second = F(1, degree+1)
    else:
        alpha = F(q['exponent'])
        first = F(0) if degree % 2 else 2/(alpha+degree+1)
        second = F(0) if degree % 2 else 2/(2*alpha+degree+1)
    if power == 1: return a*first+b*basic
    return a*a*second+2*a*b*first+b*b*basic


def multiply(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): out[i+j] += a*b
    return out


@lru_cache(maxsize=512)
def derivative_legendre(l, m):
    integer(l, 0, 96, 'degree_l'); integer(m, 0, l, 'azimuth_m')
    # Rodrigues' closed polynomial formula, with common (-1)^m omitted.
    p = [F(0)]*(l-m+1)
    for k in range((l-m)//2+1):
        j = l-2*k
        c = F((-1)**k*factorial(2*l-2*k),
              2**l*factorial(k)*factorial(l-k)*factorial(j))
        p[j-m] = c*F(factorial(j), factorial(j-m))
    return tuple(p)


@lru_cache(maxsize=4096)
def radial_product(m, l, k):
    p = multiply(derivative_legendre(l,m), derivative_legendre(k,m))
    weight = [F(0)]*(2*m+1)
    for j in range(m+1): weight[2*j] = F((-1)**j*comb(m,j))
    return tuple(multiply(p,weight))


def integral_product(q, m, l, k, power=1):
    return sum((v*moment(q,i,power) for i,v in enumerate(radial_product(m,l,k)) if v), F(0))


def form_bound(function, sqrt_bits=40):
    q=normalize(function); integer(sqrt_bits,8,160,'sqrt_bits')
    a,b=F(q['amplitude']),F(q['offset'])
    if q['profile']=='step':
        lower=min(b,a+b)
        return {'method':'pointwise_step','energy_factor':'1','mass_offset':str(lower),
                'pointwise_lower':str(lower),'measure':old.MEASURE}
    center=moment(q,0)/2
    variance=moment(q,0,2)/2-center*center
    eta=f.sqrt_upper(variance,sqrt_bits)
    if eta>=1: raise ValueError('Centered L2 form bound requires eta < 1')
    return {'method':'centered_probability_L2','energy_factor':str(1-eta),
            'mass_offset':str(center-eta),'center':str(center),
            'variance':str(variance),'eta_upper':str(eta),'sqrt_bits':sqrt_bits,
            'measure':old.MEASURE,'embedding':deepcopy(old.EMBEDDING)}


def tail_bound(form, l):
    return F(form['energy_factor'])*l*(l+1)+F(form['mass_offset'])


def matrix_assembly(function, m=0, mean_zero=True, modes=8, near_tail=0):
    q=normalize(function)
    integer(m,0,12,'azimuth_m'); integer(modes,1,32,'modes')
    integer(near_tail,0,32,'near_tail')
    if type(mean_zero) is not bool: raise ValueError('mean_zero must be boolean')
    start=1 if m==0 and mean_zero else m
    degrees=list(range(start,start+modes)); mass=[wide.basis_mass(m,l) for l in degrees]
    v=[[integral_product(q,m,l,k) for k in degrees] for l in degrees]
    t=[[integral_product(q,m,l,k,2) for k in degrees] for l in degrees]
    a=[row[:] for row in v]
    for i,l in enumerate(degrees): a[i][i]+=l*(l+1)*mass[i]
    c=[[t[i][j]-sum((v[i][k]*v[j][k]/mass[k] for k in range(modes)),F(0))
        for j in range(modes)] for i in range(modes)]
    removed_constant=None
    if m==0 and mean_zero:
        column=[integral_product(q,m,l,0) for l in degrees]
        removed_constant={'mass':F(2),'column':column}
        for i in range(modes):
            for j in range(modes): c[i][j]-=column[i]*column[j]/2
    near=[]
    for k in range(start+modes,start+modes+near_tail):
        column=[integral_product(q,m,l,k) for l in degrees]
        mk=wide.basis_mass(m,k)
        near.append({'degree':k,'mass':mk,'column':column})
        for i in range(modes):
            for j in range(modes): c[i][j]-=column[i]*column[j]/mk
    if base.inertia(c)[0]: raise ArithmeticError('Exact residual Gram must be positive semidefinite')
    return {'degrees':degrees,'mass':mass,'A':a,'V':v,'T':t,'C':c,
            'removed_constant':removed_constant,'near_tail':near}


def encode(value):
    if isinstance(value,F): return str(value)
    if isinstance(value,dict): return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [encode(v) for v in value]
    return value


class Kernel:
    def __init__(self,function,m=0,mean_zero=True,modes=8,sqrt_bits=40,near_tail=0):
        self.function=normalize(function)
        self.m,self.mean_zero,self.n,self.near_tail=m,mean_zero,modes,near_tail
        self.sqrt_bits=sqrt_bits
        self.data=matrix_assembly(self.function,m,mean_zero,modes,near_tail)
        self.form=form_bound(self.function,sqrt_bits)
        self.start=self.data['degrees'][0]
        self.tail_start=self.start+modes
        self.beta=tail_bound(self.form,self.tail_start)
        self.remainder_beta=tail_bound(self.form,self.tail_start+near_tail)

    def evidence(self): return encode(self.data)

    def matrix(self,x,kind='lower'):
        x=f.rational(x)
        if kind not in ('lower','upper'): raise ValueError('Unknown comparison kind')
        if kind=='lower' and x>=self.beta: raise ValueError('Strictly positive radial tail required')
        result=[row[:] for row in self.data['A']]
        for i in range(self.n): result[i][i]-=x*self.data['mass'][i]
        if kind=='upper': return result
        for i in range(self.n):
            for j in range(self.n):
                result[i][j]-=self.data['C'][i][j]/(self.remainder_beta-x)
        for near in self.data['near_tail']:
            denominator=near['mass']*(tail_bound(self.form,near['degree'])-x)
            for i in range(self.n):
                for j in range(self.n): result[i][j]-=near['column'][i]*near['column'][j]/denominator
        return result

    def count(self,x,kind='lower',independent=False):
        matrix=self.matrix(x,kind)
        return base.inertia(matrix) if independent else arithmetic.banded_inertia(matrix)

    def lower_holds(self,x): return x<self.beta and self.count(x)[0]==0
    def upper_holds(self,x):
        neg,zero,_=self.count(x,'upper')
        return neg+zero>=1


def sector_certificate(kernel,lower,upper,independent=True):
    low,high=f.rational(lower),f.rational(upper)
    if low>high: raise ValueError('Inverted spectral interval')
    lc=kernel.count(low,'lower',independent)
    uc=kernel.count(high,'upper',independent)
    if lc[0] or uc[0]+uc[1]<1: raise ValueError('Inertia fails to establish spectral endpoints')
    return {'format':SECTOR,'function':kernel.function,'geometry':'unit_S2',
            'axis_isometry':'coordinate_permutation_preserves_energy_mass_and_mean',
            'scope':'single_azimuth_sector_only','azimuth_m':kernel.m,'mean_zero':kernel.mean_zero,
            'projection':'remove_l0' if kernel.m==0 and kernel.mean_zero else 'none',
            'eigenvalue_index':1,'modes':kernel.n,'sqrt_bits':kernel.sqrt_bits,
            'near_tail':kernel.near_tail,'kernel_evidence':kernel.evidence(),
            'form_bound':kernel.form,'radial_tail_start':kernel.tail_start,
            'radial_tail_lower':str(kernel.beta),'remainder_tail_lower':str(kernel.remainder_beta),
            'lower':str(low),'upper':str(high),'exact_width':str(high-low),
            'lower_inertia':list(lc),'upper_inertia':list(uc),
            'original_function_moments_exact':True,'function_approximation_error':'0',
            'formal_proof_assistant_checked':False}


def certify_sector(function,m=0,mean_zero=True,modes=8,bits=32,sqrt_bits=40,near_tail=0):
    integer(bits,8,96,'bits')
    k=Kernel(function,m,mean_zero,modes,sqrt_bits,near_tail)
    analytic_low=tail_bound(k.form,k.start)
    high=min(k.data['A'][i][i]/k.data['mass'][i] for i in range(k.n))
    low=analytic_low-1
    if not k.upper_holds(high): raise ArithmeticError('Ritz trial upper bracket failed')
    for _ in range(bits):
        mid=(low+high)/2
        if k.upper_holds(mid): high=mid
        else: low=mid
    upper=high
    distance=F(1);low=analytic_low-distance
    for _ in range(64):
        if k.lower_holds(low): break
        distance*=2;low=analytic_low-distance
    else: raise ValueError('Lower bracket budget exhausted')
    high=min(upper,k.beta)
    for _ in range(bits):
        mid=(low+high)/2
        if k.lower_holds(mid): low=mid
        else: high=mid
    return sector_certificate(k,low,upper)


def verify_sector(cert,expected_function=None,expected_m=None,expected_mean_zero=None):
    try:
        if not isinstance(cert,dict) or cert.get('format')!=SECTOR: return False
        if expected_function is not None and normalize(expected_function)!=cert['function']: return False
        if expected_m is not None and (type(expected_m) is not int or expected_m!=cert['azimuth_m']): return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool or expected_mean_zero is not cert['mean_zero']): return False
        k=Kernel(cert['function'],cert['azimuth_m'],cert['mean_zero'],cert['modes'],cert['sqrt_bits'],cert['near_tail'])
        return same(cert,sector_certificate(k,cert['lower'],cert['upper']))
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError): return False


def full_certificate(function,mean_zero,sectors,tolerance='1/100000000'):
    q=normalize(function)
    if type(mean_zero) is not bool or not isinstance(sectors,list) or not 1<=len(sectors)<=13:
        raise ValueError('Require consecutive sectors and boolean mean_zero')
    tolerance=f.rational(tolerance)
    if not F(1,10**30)<=tolerance<=1: raise ValueError('Require tolerance in 1e-30..1')
    for m,cert in enumerate(sectors):
        if not verify_sector(cert,q,m,mean_zero): raise ValueError('Invalid complete sector sequence')
    form=form_bound(q,sectors[0]['sqrt_bits'])
    angular=tail_bound(form,len(sectors))
    lower=min([angular]+[F(s['lower']) for s in sectors])
    upper=min(F(s['upper']) for s in sectors)
    if lower>upper: raise ArithmeticError('Global spectral contradiction')
    return {'format':FULL,'function':q,'geometry':'unit_S2','measure':old.MEASURE,
            'scope':old.sphere_scope(mean_zero),'mean_zero':mean_zero,'eigenvalue_index':1,
            'sectors':deepcopy(sectors),'angular_form_bound':form,
            'omitted_azimuth_m_start':len(sectors),'angular_tail_lower':str(angular),
            'angular_tail_cannot_improve_best_upper':angular>=upper,
            'lower':str(lower),'upper':str(upper),'exact_width':str(upper-lower),
            'tolerance':str(tolerance),'status':'target_met' if upper-lower<=tolerance else 'certified_open',
            'original_function_moments_exact':True,'function_approximation_error':'0',
            'full_infinite_space_covered':True,'formal_proof_assistant_checked':False}


def full_ground(function,mean_zero=True,modes=8,bits=32,max_m=4,
                tolerance='1/100000000',sqrt_bits=40,near_tail=0):
    q=normalize(function);integer(max_m,0,12,'max_m')
    form=form_bound(q,sqrt_bits);sectors=[]
    for m in range(max_m+1):
        sectors.append(certify_sector(q,m,mean_zero,modes,bits,sqrt_bits,near_tail))
        if tail_bound(form,m+1)>=min(F(c['upper']) for c in sectors): break
    return full_certificate(q,mean_zero,sectors,tolerance)


def verify_full(cert,expected_function=None,expected_mean_zero=None,expected_tolerance=None):
    try:
        if not isinstance(cert,dict) or cert.get('format')!=FULL: return False
        if expected_function is not None and normalize(expected_function)!=cert['function']: return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool or expected_mean_zero is not cert['mean_zero']): return False
        if expected_tolerance is not None and f.rational(expected_tolerance)!=F(cert['tolerance']): return False
        return same(cert,full_certificate(cert['function'],cert['mean_zero'],cert['sectors'],cert['tolerance']))
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError): return False
