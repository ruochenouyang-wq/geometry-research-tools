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
    _count('cell_integrations')
    return pieces._cell(left_root, right_root, degree, amplitude, offset, core)


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
