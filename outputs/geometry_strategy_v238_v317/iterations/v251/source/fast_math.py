"""Exact optimizations with immutable caches and frozen public verification."""
from collections import Counter
from copy import deepcopy
from functools import lru_cache
from math import comb
import json
from time import perf_counter

from backend import F, canonical, digest, pieces, direct, enriched

_COUNTERS = Counter()
_CACHES = []


def _count(name, amount=1):
    _COUNTERS[name] += amount


def _cached(maxsize):
    def wrap(function):
        result = lru_cache(maxsize=maxsize)(function)
        _CACHES.append(result)
        return result
    return wrap


def reset_counters(clear_caches=False):
    if clear_caches:
        for cache in _CACHES:
            cache.cache_clear()
    _COUNTERS.clear()


def counters():
    return dict(_COUNTERS)


def _cell(left_root, right_root, degree, amplitude, offset, core=False):
    """Integrate the returned polynomial in its bounded local coordinate."""
    _count('cell_integrations')
    a, b = left_root**4, right_root**4
    width = b-a
    moments = [offset*(b**(k+1)-a**(k+1))/F(k+1)
               + amplitude*pieces._power_moment(k, left_root, right_root)
               for k in range(degree+1)]
    # s=(t-a)/width, expanded exactly once per degree instead of per basis.
    local_moments = [sum((F(comb(k,j))*(-a)**(k-j)*moments[j]
                          for j in range(k+1)), F(0))/width**k
                     for k in range(degree+1)]
    source_squared = (offset**2*width
                      + 2*offset*amplitude*pieces._power_moment(0,left_root,right_root)
                      + 2*amplitude**2*(right_root**2-left_root**2))
    entries, polynomial, projection_norm = [], [F(0)], F(0)
    for l, basis in enumerate(pieces._shifted_legendre(degree)):
        inner = sum((c*local_moments[k] for k,c in enumerate(basis)), F(0))
        mass = width/F(2*l+1)
        coefficient = inner/mass
        contribution = coefficient**2*mass
        projection_norm += contribution
        polynomial = pieces._add(polynomial, pieces._scale(basis,coefficient))
        entries.append({'degree':l,'shifted_legendre_in_s':list(map(str,basis)),
                        'source_inner_product_dt':str(inner),'basis_mass_dt':str(mass),
                        'coefficient':str(coefficient),'norm_squared_dt':str(contribution)})
    squared = pieces._mul(polynomial,polynomial)
    direct_squared = width*sum((c/F(k+1) for k,c in enumerate(squared)), F(0))
    cross = sum((c*local_moments[k] for k,c in enumerate(polynomial)), F(0))
    residual = source_squared-2*cross+direct_squared
    if cross != projection_norm or direct_squared != projection_norm or residual < 0:
        raise ArithmeticError('Exact local-coordinate projection identity failed')
    return {'left':str(a),'right':str(b),'left_fourth_root':str(left_root),
            'right_fourth_root':str(right_root),'core':core,'degree':degree,
            'coordinate':'s=(abs(t)-left)/(right-left)',
            'polynomial_in_s':list(map(str,polynomial)),
            'source_monomial_moments_dt':list(map(str,moments)),
            'projection':entries,'source_norm_squared_dt':str(source_squared),
            'approximation_norm_squared_dt':str(direct_squared),'cross_integral_dt':str(cross),
            'error_squared_dt':str(residual)}


@_cached(32)
def _shape_text(ratio, degree):
    _count('shape_integrations')
    return canonical(_cell(ratio, F(1), degree, F(1), F(0)))


def piecewise_model(function, levels=4, degree=3, root_ratio='1/2', sqrt_bits=40):
    """The legacy wire certificate, with a cached immutable unit shape."""
    function = pieces.normalize_function(function)
    pieces._integer(levels, 1, pieces.MAX_LEVELS, 'levels')
    pieces._integer(degree, 0, pieces.MAX_DEGREE, 'degree')
    pieces._integer(sqrt_bits, 0, 256, 'sqrt_bits')
    ratio = pieces.rational(root_ratio)
    if not 0 < ratio < 1 or ratio.denominator > 4096:
        raise ValueError('root_ratio must lie in (0,1) with denominator<=4096')
    if function['profile'] == 'step':
        return pieces.piecewise_model(function, levels, degree, root_ratio, sqrt_bits)
    amplitude, offset = F(function['amplitude']), F(function['offset'])
    shape = json.loads(_shape_text(ratio, degree))
    cells = [_cell(F(0), ratio**levels, 0, amplitude, offset, core=True)]
    for j in range(levels-1, -1, -1):
        _count('scaled_cells')
        cells.append(pieces._scaled_cell(shape, ratio**j, amplitude, offset))
    error = sum((F(cell['error_squared_dt']) for cell in cells), F(0))
    geometric = amplitude**2*(F(2, 9)*ratio**(2*levels)
        + F(shape['error_squared_dt'])*(1-ratio**(2*levels))/(1-ratio**2))
    source_norm = offset**2+F(8, 3)*offset*amplitude+2*amplitude**2
    if error != geometric or sum((F(c['source_norm_squared_dt']) for c in cells), F(0)) != source_norm:
        raise ArithmeticError('Complete scaled model identities failed')
    upper = pieces.sqrt_upper(error, sqrt_bits)
    return {'format': pieces.FORMAT, 'function': function, 'scope': pieces.SCOPE,
            'norm': 'L2_probability', 'error_squared': str(error), 'error_upper': str(upper),
            'sqrt_upper_squared': str(upper**2), 'source_norm_squared': str(source_norm),
            'levels': levels, 'degree': degree, 'root_ratio': str(ratio), 'sqrt_bits': sqrt_bits,
            'representation': 'even_reflection_of_positive_half_cells', 'cells': cells,
            'shape_proof': shape,
            'shape_rule': 't=R^4*x scales centered profile coefficients by 1/R and squared residual integral by R^2',
            'source_value_at_axis_zero': str(offset),
            'axis_zero_convention': 'step(0)=0; singular_power(0)=0 as an a.e. representative',
            'endpoint_convention': 'at nonzero seams choose the cell starting at abs(t); final endpoint included',
            'expanded_polynomial_coefficient_slots': 1+levels*(degree+1),
            'shape_coefficient_count': degree+1,
            'bound_statement': 'probability_L2_norm(original_function-piecewise_model)<=error_upper',
            'spectral_transfer_claimed': False, 'global_polynomial': False,
            'formal_assistant_checked': False}


verify_piecewise = pieces.verify_piecewise


def evaluate_many(certificate, coordinates):
    """Replay an owned certificate once, then evaluate all exact coordinates."""
    owned = deepcopy(certificate)
    _count('piecewise_batch_verifications')
    if not pieces.verify_piecewise(owned):
        raise ValueError('Invalid piecewise model')
    points = [pieces.rational(t) for t in coordinates]
    if any(not -1 <= t <= 1 for t in points):
        raise ValueError('Axis coordinate must be in [-1,1]')
    q = owned['function']
    if q['profile'] == 'step':
        a,b = F(q['amplitude']),F(q['offset'])
        return [b+(a if t > 0 else 0) for t in points]
    cells = [(F(c['left']),F(c['right']),tuple(map(F,c['polynomial_in_s'])))
             for c in owned['cells']]
    values = []
    for t in points:
        t = abs(t)
        for a,b,coefficients in cells:
            if a <= t < b or t == b == 1:
                s, value = (t-a)/(b-a),F(0)
                for coefficient in reversed(coefficients):
                    value = value*s+coefficient
                values.append(value)
                break
        else:
            raise ArithmeticError('A verified cover must contain this coordinate')
    return values


def evaluate(certificate, t):
    return evaluate_many(certificate, [t])[0]


@_cached(8192)
def _raw_moment(key, degree, power):
    _count('raw_moment_computations')
    profile,a,b,alpha = key
    basic = F(0) if degree % 2 else F(2,degree+1)
    if power == 0: return basic
    if profile == 'step':
        first = second = F(1,degree+1)
    else:
        first = F(0) if degree % 2 else 2/(alpha+degree+1)
        second = F(0) if degree % 2 else 2/(2*alpha+degree+1)
    return a*first+b*basic if power == 1 else a*a*second+2*a*b*first+b*b*basic


class MomentTable:
    """Normalize one owned source and deduplicate exact (degree,power) moments."""
    def __init__(self,function):
        _count('moment_source_normalizations')
        self.function = direct.normalize(function)
        q = self.function
        self.key = (q['profile'],F(q['amplitude']),F(q['offset']),F(q.get('exponent',0)))

    def moment(self,degree,power=1):
        direct.integer(degree,0,256,'moment_degree')
        direct.integer(power,0,2,'power')
        _count('raw_moment_requests')
        return _raw_moment(self.key,degree,power)

    def integral(self,m,l,k,power=1):
        _count('radial_integral_requests')
        return sum((v*self.moment(i,power)
                    for i,v in enumerate(direct.radial_product(m,l,k)) if v),F(0))


def moment(function,degree,power=1):
    return MomentTable(function).moment(degree,power)


def matrix_assembly(function,m=0,mean_zero=True,modes=8,near_tail=0):
    table = MomentTable(function)
    direct.integer(m,0,12,'azimuth_m'); direct.integer(modes,1,32,'modes')
    direct.integer(near_tail,0,32,'near_tail')
    if type(mean_zero) is not bool: raise ValueError('mean_zero must be boolean')
    start = 1 if m == 0 and mean_zero else m
    degrees = list(range(start,start+modes))
    mass = [direct.wide.basis_mass(m,l) for l in degrees]
    v = [[table.integral(m,l,k) for k in degrees] for l in degrees]
    t = [[table.integral(m,l,k,2) for k in degrees] for l in degrees]
    a = [row[:] for row in v]
    for i,l in enumerate(degrees): a[i][i] += l*(l+1)*mass[i]
    c = [[t[i][j]-sum((v[i][k]*v[j][k]/mass[k] for k in range(modes)),F(0))
          for j in range(modes)] for i in range(modes)]
    removed_constant = None
    if m == 0 and mean_zero:
        column = [table.integral(m,l,0) for l in degrees]
        removed_constant = {'mass':F(2),'column':column}
        for i in range(modes):
            for j in range(modes): c[i][j] -= column[i]*column[j]/2
    near = []
    for k in range(start+modes,start+modes+near_tail):
        column = [table.integral(m,l,k) for l in degrees]
        mk = direct.wide.basis_mass(m,k)
        near.append({'degree':k,'mass':mk,'column':column})
        for i in range(modes):
            for j in range(modes): c[i][j] -= column[i]*column[j]/mk
    if direct.base.inertia(c)[0]: raise ArithmeticError('Exact residual Gram must be positive semidefinite')
    return {'degrees':degrees,'mass':mass,'A':a,'V':v,'T':t,'C':c,
            'removed_constant':removed_constant,'near_tail':near}
