"""Exact L2 certificates for finite piecewise polynomials of rational powers.

For alpha=-m/d, use t=R**d. Both source moments and scale transports are
rational. A compact finite block description avoids repeating huge scaled
coefficients. It describes polynomials, not copies of the singular source.
"""
from fractions import Fraction as F
from math import comb, isqrt, lcm
from time import perf_counter
from collections import OrderedDict
import re
import json

FORMAT = 'general_rational_power_piecewise_L2_v1'
MEASURE = 'd_sigma/(4*pi)_on_unit_S2'
MAX_DEGREE = 32
MAX_LEVELS = 16384


def rational(value):
    if type(value) in (int,F):
        return F(value)
    if type(value) is str and len(value)<=1000000 and re.fullmatch(r'-?\d+(?:/[1-9]\d*)?',value):
        parts = value.split('/')
        return F(_parse_integer(parts[0]),_parse_integer(parts[1]) if len(parts)==2 else 1)
    raise ValueError('Use exact integers or rational strings, not floats or booleans')


def _parse_integer(text):
    """Parse bounded decimal chunks without changing Python's global limits."""
    sign = -1 if text.startswith('-') else 1
    text = text.lstrip('-')
    value = 0
    for start in range(0,len(text),9):
        chunk = text[start:start+9]
        value = value*10**len(chunk)+int(chunk)
    return sign*value


def _integer_text(value):
    if value.bit_length()<12000:
        return str(value)
    sign = '-' if value<0 else ''
    value = abs(value)
    chunks = []
    while value:
        value,remainder = divmod(value,10**9)
        chunks.append(remainder)
    return sign+str(chunks[-1])+''.join(f'{x:09d}' for x in reversed(chunks[:-1]))


def rational_text(value):
    x = F(value)
    n,d = _integer_text(x.numerator),_integer_text(x.denominator)
    return n if x.denominator==1 else n+'/'+d


def integer(value,low,high,name):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(name+' exceeds its integer budget')
    return value


def normalize(function):
    allowed = {'kind','profile','axis','exponent','amplitude','offset'}
    if type(function) is not dict or set(function)-allowed:
        raise ValueError('Malformed or unknown original-function fields')
    if function.get('kind') != 'axis_profile' or function.get('profile') != 'abs_power':
        raise ValueError('Require an axis_profile abs_power function')
    alpha = rational(function['exponent'])
    if not -F(1,2) < alpha < 0:
        raise ValueError('Require -1/2 < alpha < 0 for this L2 approximation family')
    if alpha.denominator > 4096:
        raise ValueError('Exponent denominator exceeds the exact-arithmetic budget4096')
    amplitude,offset = rational(function.get('amplitude',1)),rational(function.get('offset',0))
    if any(abs(x)>10**6 or x.denominator>10**6 for x in (amplitude,offset)):
        raise ValueError('Coefficient input budget exceeded')
    return {'kind':'axis_profile','profile':'abs_power',
            'axis':integer(function.get('axis',2),0,2,'axis'),
            'exponent':str(alpha),'amplitude':str(amplitude),'offset':str(offset)}


def _ratio(value):
    r = rational(value)
    if not 0 < r < 1 or r.denominator > 10**6:
        raise ValueError('Root ratio must lie in (0,1), denominator at most1e6')
    return r


def sqrt_upper(value,bits):
    x = rational(value)
    integer(bits,0,256,'sqrt_bits')
    if x < 0:
        raise ValueError('Negative error square')
    den = 1 << bits
    n = isqrt(x.numerator*den*den//x.denominator)
    if F(n*n,den*den) < x:
        n += 1
    return F(n,den)


def power_moment(alpha,k,left_root,right_root,power=1):
    """Integral t**(k+power*alpha), endpoints are rational d-th powers."""
    alpha = rational(alpha)
    if not -F(1,2) < alpha < 0:
        raise ValueError('Invalid singular exponent')
    integer(k,0,2*MAX_DEGREE,'moment_degree')
    integer(power,0,2,'source_power')
    left,right = rational(left_root),rational(right_root)
    if not 0 <= left < right <= 1:
        raise ValueError('Invalid root endpoints')
    d,m = alpha.denominator,-alpha.numerator
    exponent = d*(k+1)-power*m
    if exponent <= 0:
        raise ValueError('Source moment is not integrable')
    return F(d,exponent)*(right**exponent-left**exponent)


def _local_moments(alpha,ratio,degree):
    d = alpha.denominator
    left,width = ratio**d,1-ratio**d
    raw = [power_moment(alpha,k,ratio,F(1)) for k in range(degree+1)]
    local = [sum((F(comb(k,j))*(-left)**(k-j)*raw[j]
                  for j in range(k+1)),F(0))/width**k for k in range(degree+1)]
    return left,width,raw,local


def _legendre(degree):
    return [F((-1)**(degree+k)*comb(degree,k)*comb(degree+k,k))
            for k in range(degree+1)]


def _integer_coefficients(coefficients):
    denominator = lcm(*(c.denominator for c in coefficients))
    return [c.numerator*(denominator//c.denominator) for c in coefficients],denominator


def _polynomial_norm(coefficients,width):
    # One large Fraction construction, not a quadratic number of large gcds.
    nums,den = _integer_coefficients(coefficients)
    common = lcm(*range(1,2*len(nums)))
    total = sum(x*y*(common//(i+j+1))
                for i,x in enumerate(nums) for j,y in enumerate(nums))
    return width*F(total,den*den*common)


def _dot(a,b):
    aa,ad = _integer_coefficients(a)
    bb,bd = _integer_coefficients(b)
    return F(sum(x*y for x,y in zip(aa,bb)),ad*bd)


class ShapeCache:
    """Optional immutable exact planning shapes, keyed by exponent AND ratio."""
    def __init__(self,capacity=64):
        self.capacity = integer(capacity,1,1024,'cache_capacity')
        self.values = OrderedDict()

    def get(self,alpha,ratio,degree):
        key = (F(alpha),F(ratio),degree)
        value = self.values.get(key)
        if value is not None:
            self.values.move_to_end(key)
        return value

    def put(self,alpha,ratio,degree,polynomial,error):
        key = (F(alpha),F(ratio),degree)
        self.values[key] = (tuple(polynomial),F(error))
        self.values.move_to_end(key)
        while len(self.values)>self.capacity:
            self.values.popitem(last=False)


def shape_series(alpha,ratio,max_degree,cache=None):
    """Yield nested exact local L2 projections, constructing only needed degrees."""
    alpha,ratio = rational(alpha),_ratio(ratio)
    integer(max_degree,0,MAX_DEGREE,'max_degree')
    if cache is not None and type(cache) is not ShapeCache:
        raise ValueError('Expected an explicit ShapeCache')
    norm = power_moment(alpha,0,ratio,F(1),2)
    polynomial,projection = [],F(0)
    moments = None
    for degree in range(max_degree+1):
        cached = cache.get(alpha,ratio,degree) if cache is not None else None
        if cached is not None:
            polynomial,error = list(cached[0]),cached[1]
            projection = norm-error
            yield tuple(polynomial),error
            continue
        if moments is None:
            _,width,_,moments = _local_moments(alpha,ratio,max_degree)
        basis = _legendre(degree)
        inner = _dot(basis,moments[:degree+1])
        mass = width/F(2*degree+1)
        coefficient = inner/mass
        projection += coefficient*inner
        polynomial += [F(0)]
        for k,c in enumerate(basis):
            polynomial[k] += coefficient*c
        error = norm-projection
        if error < 0:
            raise ArithmeticError('Negative exact projection error')
        if cache is not None:
            cache.put(alpha,ratio,degree,polynomial,error)
        yield tuple(polynomial),error


def _core_coefficient(alpha):
    return 1/(1+alpha)


def _core_error(alpha,constant):
    return 1/(1+2*alpha)-2*constant/(1+alpha)+constant*constant


def _first_level(coefficient,ratio_power,target,max_levels):
    """First j in [0,max_levels] with coefficient*ratio_power**j<=target."""
    if coefficient <= target:
        return 0
    if coefficient*ratio_power**max_levels > target:
        return max_levels+1
    lo,hi = 1,max_levels
    while lo < hi:
        mid = (lo+hi)//2
        if coefficient*ratio_power**mid <= target:
            hi = mid
        else:
            lo = mid+1
    return lo


def _block_weight(ratio_power,start,count):
    return ratio_power**start*(1-ratio_power**count)/(1-ratio_power)


def _summary(q,ratio,levels,core_constant,templates,blocks,tolerance,bits):
    alpha,A,b = (F(q[x]) for x in ('exponent','amplitude','offset'))
    d,m = alpha.denominator,-alpha.numerator
    gamma = d-2*m
    weight = ratio**gamma
    core = A*A*_core_error(alpha,core_constant)*weight**levels
    annuli = A*A*sum((rational(templates[str(block['degree'])]['error_squared'])
                      *_block_weight(weight,block['start'],block['count']) for block in blocks),F(0))
    error = core+annuli
    upper = sqrt_upper(error,bits)
    return {'format':FORMAT,'function':q,'scope':'all_unit_S2_positions_in_probability_L2',
            'measure':MEASURE,'norm':'L2_probability','tolerance':str(tolerance),
            'error_squared':rational_text(error),'error_upper':rational_text(upper),
            'sqrt_bits':bits,'sqrt_upper_squared':rational_text(upper*upper),
            'core_error_squared':rational_text(core),'annuli_error_squared':rational_text(annuli),
            'source_norm_squared':rational_text(A*A/(1+2*alpha)+2*A*b/(1+alpha)+b*b),
            'root_ratio':str(ratio),'root_degree':d,'singular_numerator':m,
            'squared_error_scale_exponent':gamma,'levels':levels,
            'templates':templates,'blocks':blocks,
            'core_normal_coefficient':str(core_constant),
            'model':{'kind':'finite_even_piecewise_polynomial_scaled_templates',
                     'coordinate':'t=abs(axis_coordinate)',
                     'boundary_rule':'t_j=root_ratio**(root_degree*j)',
                     'annulus_rule':'j in block; t in [t_(j+1),t_j]; s=(t/t_j-root_ratio**root_degree)/(1-root_ratio**root_degree)',
                     'polynomial_rule':'offset+amplitude*root_ratio**(-singular_numerator*j)*template_degree(s)',
                     'core_rule':'0<=t<=t_levels; offset+amplitude*root_ratio**(-singular_numerator*levels)*core_normal_coefficient',
                     'seam_rule':'outer cell at internal annulus seams; core includes t_levels and zero',
                     'all_output_basis_functions_are_piecewise_polynomials':True},
            'expanded_polynomial_coefficient_slots':1+sum(block['count']*(block['degree']+1) for block in blocks),
            'stored_template_coefficient_slots':1+sum(len(t['polynomial_in_s']) for t in templates.values()),
            'status':'target_met' if upper<=tolerance else 'certified_open',
            'spectral_transfer_claimed':False,'global_polynomial':False,
            'formal_proof_assistant_checked':False}


def verify(cert,function=None,tolerance=None):
    """True means the full error bound is valid, even if status is open.

    Re-integrates the supplied finite polynomials; no planner/cache/projection
    generator is used. The caller checks error_upper<=the requested tolerance.
    """
    try:
        if type(cert) is not dict or cert.get('format') != FORMAT:
            return False
        q = normalize(cert['function'])
        if function is not None and q != normalize(function):
            return False
        tol = rational(cert['tolerance'])
        if not 0 < tol <= 1 or (tolerance is not None and rational(tolerance) != tol):
            return False
        alpha,ratio = F(q['exponent']),_ratio(cert['root_ratio'])
        levels = integer(cert['levels'],1,MAX_LEVELS,'levels')
        bits = integer(cert['sqrt_bits'],0,256,'sqrt_bits')
        core = rational(cert['core_normal_coefficient'])
        if _core_error(alpha,core) < 0:
            return False
        blocks,templates = cert['blocks'],cert['templates']
        if type(blocks) is not list or not 1<=len(blocks)<=MAX_DEGREE+1 or type(templates) is not dict:
            return False
        cursor,used = 0,set()
        for block in blocks:
            if type(block) is not dict or set(block)!={'start','count','degree'}:
                return False
            start = integer(block['start'],0,levels-1,'block_start')
            count = integer(block['count'],1,levels,'block_count')
            degree = integer(block['degree'],0,MAX_DEGREE,'degree')
            if start != cursor or start+count > levels:
                return False
            cursor += count;used.add(str(degree))
        if cursor != levels or set(templates)!=used:
            return False
        degree = max(int(x) for x in used)
        _,width,_,moments = _local_moments(alpha,ratio,degree)
        source_norm = power_moment(alpha,0,ratio,F(1),2)
        rebuilt_templates = {}
        for name in sorted(templates,key=int):
            entry = templates[name]
            if type(entry) is not dict or set(entry)!={'polynomial_in_s','error_squared'}:
                return False
            raw = entry['polynomial_in_s']
            if type(raw) is not list or len(raw) != int(name)+1:
                return False
            coefficients = [rational(x) for x in raw]
            approximation = _polynomial_norm(coefficients,width)
            cross = _dot(coefficients,moments[:len(coefficients)])
            error = source_norm-2*cross+approximation
            if error < 0 or error != rational(entry['error_squared']):
                return False
            rebuilt_templates[name] = {'polynomial_in_s':list(map(rational_text,coefficients)),
                                       'error_squared':rational_text(error)}
        rebuilt = _summary(q,ratio,levels,core,rebuilt_templates,blocks,tol,bits)
        return json.dumps(cert,sort_keys=True,allow_nan=False)==json.dumps(rebuilt,sort_keys=True,allow_nan=False)
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError,OverflowError):
        return False


def solve(function,tolerance='1/100000000',max_degree=24,max_levels=2048,
          root_ratio=None,sqrt_bits=None,shape_cache=None):
    """Return one complete polynomial-model certificate, never the source itself."""
    q = normalize(function)
    tol = rational(tolerance)
    if not F(1,10**60) <= tol <= 1:
        raise ValueError('Require tolerance in 1e-60..1')
    integer(max_degree,0,MAX_DEGREE,'max_degree')
    integer(max_levels,1,MAX_LEVELS,'max_levels')
    if shape_cache is not None and type(shape_cache) is not ShapeCache:
        raise ValueError('Expected an explicit ShapeCache')
    alpha,A = F(q['exponent']),F(q['amplitude'])
    d,m = alpha.denominator,-alpha.numerator
    ratio = _ratio(root_ratio if root_ratio is not None else F(2*d-1,2*d))
    weight = ratio**(d-2*m)
    core_constant = _core_coefficient(alpha)
    normalized_target = tol*tol/(A*A) if A else F(1)
    levels = min(max_levels,max(1,_first_level(_core_error(alpha,core_constant),weight,
                                               normalized_target/4,max_levels)))
    per_cell_target = normalized_target/(2*levels)
    series = shape_series(alpha,ratio,max_degree,shape_cache)
    entries,errors = {},[]
    for degree,(polynomial,error) in enumerate(series):
        entries[str(degree)] = {'polynomial_in_s':list(map(rational_text,polynomial)),
                                'error_squared':rational_text(error)}
        errors.append(error)
        if error <= per_cell_target or not A:
            break
    # Inner annuli need fewer polynomial coefficients than outer annuli.
    starts = [_first_level(error,weight,per_cell_target,levels) for error in errors]
    starts[-1] = 0  # Budget-open case still covers every outer annulus.
    blocks = []
    for degree in range(len(errors)-1,-1,-1):
        start = min(levels,starts[degree])
        stop = min(levels,starts[degree-1]) if degree else levels
        if stop > start:
            blocks.append({'start':start,'count':stop-start,'degree':degree})
    used = {str(block['degree']) for block in blocks}
    templates = {key:value for key,value in entries.items() if key in used}
    if sqrt_bits is None:
        bits = 48
        while bits < 256 and F(1,2**bits)>tol/8:
            bits += 1
    else:
        bits = integer(sqrt_bits,0,256,'sqrt_bits')
    cert = _summary(q,ratio,levels,core_constant,templates,blocks,tol,bits)
    if not verify(cert,q,tol):
        raise ArithmeticError('Independent integration rejected the generated approximation')
    return cert


def evaluate(cert,coordinate):
    """Evaluate the finite polynomial representative at an exact coordinate."""
    if not verify(cert):
        raise ValueError('Invalid whole-domain polynomial certificate')
    t = abs(rational(coordinate))
    if t>1:
        raise ValueError('Axis coordinate must lie in [-1,1]')
    q,r,L = cert['function'],F(cert['root_ratio']),cert['levels']
    d,m = cert['root_degree'],cert['singular_numerator']
    A,b = F(q['amplitude']),F(q['offset'])
    if t <= r**(d*L):
        return b+A*r**(-m*L)*F(cert['core_normal_coefficient'])
    lo,hi = 0,L-1
    while lo<hi:
        mid = (lo+hi)//2
        if t >= r**(d*(mid+1)):
            hi = mid
        else:
            lo = mid+1
    j = lo
    block = next(block for block in cert['blocks'] if block['start']<=j<block['start']+block['count'])
    s = (t/r**(d*j)-r**d)/(1-r**d)
    value = F(0)
    for c in reversed(cert['templates'][str(block['degree'])]['polynomial_in_s']):
        value = value*s+rational(c)
    return b+A*r**(-m*j)*value
