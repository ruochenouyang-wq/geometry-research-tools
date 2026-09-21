"""Automatic choices for the frozen, exactly verified probability-L2 model.

All planning values are exact rationals. Only a replayed old-format certificate
can satisfy the requested tolerance; this module makes no spectral assertion.
"""
from fractions import Fraction as F
import backend

pieces = backend.pieces
FORMAT = 'adaptive_piecewise_search_v1'


def _limits(max_degree, max_levels):
    pieces._integer(max_degree, 0, pieces.MAX_DEGREE, 'max_degree')
    pieces._integer(max_levels, 1, pieces.MAX_LEVELS, 'max_levels')


def _ratio(value):
    ratio = pieces.rational(value)
    if not 0 < ratio < 1 or ratio.denominator > 4096:
        raise ValueError('root_ratio must lie in (0,1), denominator <=4096')
    return ratio


def _components(shape, levels, ratio, amplitude):
    core = amplitude**2*F(2, 9)*ratio**(2*levels)
    annuli = amplitude**2*shape*(1-ratio**(2*levels))/(1-ratio**2)
    return {'shape_error_squared': str(shape), 'core_error_squared': str(core),
            'annulus_error_squared': str(annuli), 'error_squared': str(core+annuli),
            'dominant_component': 'core' if core >= annuli else 'annuli'}


def error_terms(degree, levels, root_ratio='1/2', amplitude='1'):
    """Separate the exact core, annulus and outward-rounding-free errors."""
    _limits(degree, levels)
    ratio, amplitude = _ratio(root_ratio), pieces.rational(amplitude)
    shape = F(pieces._cell(ratio, F(1), degree, F(1), F(0))['error_squared_dt'])
    return _components(shape, levels, ratio, amplitude)


def minimum_levels(shape, tolerance, ratio, amplitude, max_levels):
    """Find the least admissible level count by exact monotone bisection."""
    floor = amplitude**2*shape/(1-ratio**2)
    initial = amplitude**2*F(2, 9)
    if floor > initial:
        raise ArithmeticError('Projection error must decrease under subdivision')
    if floor > tolerance**2 or (floor == tolerance**2 and initial > floor):
        return max_levels, 0, False, 'fixed_degree_floor'
    checks = 0
    def meets(levels):
        nonlocal checks
        checks += 1
        return floor+(initial-floor)*ratio**(2*levels) <= tolerance**2
    if not meets(max_levels):
        return max_levels, checks, False, 'core_level_budget'
    lo, hi = 1, max_levels
    while lo < hi:
        mid = (lo+hi)//2
        if meets(mid):
            hi = mid
        else:
            lo = mid+1
    return lo, checks, True, 'reachable'


def sqrt_precision(error_squared, tolerance, minimum_bits, maximum_bits):
    """Increase outward-square-root precision only when it blocks a valid plan."""
    checks = 0
    def upper(bits):
        nonlocal checks
        checks += 1
        return pieces.sqrt_upper(error_squared, bits)
    first = upper(minimum_bits)
    if first <= tolerance or error_squared > tolerance**2:
        return minimum_bits, first, checks, 'sufficient' if first <= tolerance else 'approximation_error'
    last = upper(maximum_bits)
    if last > tolerance:
        return maximum_bits, last, checks, 'sqrt_bit_budget'
    lo, hi = minimum_bits+1, maximum_bits
    while lo < hi:
        mid = (lo+hi)//2
        if upper(mid) <= tolerance:
            hi = mid
        else:
            lo = mid+1
    return lo, upper(lo), checks, 'increased_to_meet_target'


def _finish(function, tolerance, options, attempts, diagnostics, cost):
    certificate = pieces.piecewise_model(function, **options)
    cost['certificate_builds'] = cost.get('certificate_builds', 0)+1
    if not pieces.verify_piecewise(certificate, expected_function=function):
        raise ArithmeticError('Frozen verifier rejected the selected approximation')
    cost['certificate_replays'] = cost.get('certificate_replays', 0)+1
    met = F(certificate['error_upper']) <= tolerance
    return {'format': FORMAT, 'function': function, 'tolerance': str(tolerance),
            'status': 'met' if met else 'open', 'certificate': certificate,
            'attempts': attempts, 'diagnostics': diagnostics, 'cost': cost,
            'spectral_transfer_claimed': False}


def solve(function, tolerance='1/100000000', max_degree=24, max_levels=64, **budget):
    function = pieces.normalize_function(function)
    tolerance = pieces.rational(tolerance)
    if tolerance <= 0:
        raise ValueError('tolerance must be positive')
    _limits(max_degree, max_levels)
    if set(budget)-{'root_ratio', 'root_ratios', 'sqrt_bits', 'max_sqrt_bits'}:
        raise ValueError('Unknown approximation budget')
    if 'root_ratio' in budget and 'root_ratios' in budget:
        raise ValueError('Specify root_ratio or root_ratios, not both')
    choices = budget.get('root_ratios', [budget['root_ratio']] if 'root_ratio' in budget
                         else ['1/2', '2/3', '3/4'])
    if type(choices) not in (list, tuple) or not 1 <= len(choices) <= 8:
        raise ValueError('root_ratios must contain between one and eight choices')
    ratios = list(dict.fromkeys(_ratio(choice) for choice in choices))
    bits = pieces._integer(budget.get('sqrt_bits', 40), 0, 256, 'sqrt_bits')
    maximum_bits = pieces._integer(budget.get('max_sqrt_bits', 256), bits, 256, 'max_sqrt_bits')
    cost = {'shape_evaluations': 0, 'level_comparisons': 0}
    if function['profile'] == 'step':
        options = {'degree': 0, 'levels': 1, 'root_ratio': str(ratios[0]), 'sqrt_bits': bits}
        return _finish(function, tolerance, options, [], {'selection': 'exact_step'}, cost)
    amplitude = F(function['amplitude'])
    attempts, best = [], None
    for ratio in ratios:
        for degree in range(max_degree+1):
            shape = F(pieces._cell(ratio, F(1), degree, F(1), F(0))['error_squared_dt'])
            cost['shape_evaluations'] += 1
            levels, checks, reachable, reason = minimum_levels(shape, tolerance, ratio, amplitude, max_levels)
            cost['level_comparisons'] += checks
            terms = _components(shape, levels, ratio, amplitude)
            terms.update({'predicted_target_reachable': reachable, 'limiting_component': reason,
                          'infinite_levels_error_floor_squared': str(amplitude**2*shape/(1-ratio**2))})
            options = {'degree': degree, 'levels': levels, 'root_ratio': str(ratio), 'sqrt_bits': bits}
            attempt = dict(options, **terms)
            attempts.append(attempt)
            if best is None or F(attempt['error_squared']) < F(best['error_squared']):
                best = attempt
            if reachable:
                break
    chosen_bits, rounded, checks, reason = sqrt_precision(F(best['error_squared']), tolerance, bits, maximum_bits)
    best['sqrt_bits'] = chosen_bits
    cost['sqrt_comparisons'] = checks
    rounding = {'minimum_bits': bits, 'selected_bits': chosen_bits, 'reason': reason,
                'rounded_error_upper': str(rounded),
                'squared_rounding_excess': str(rounded**2-F(best['error_squared']))}
    options = {name: best[name] for name in ('degree', 'levels', 'root_ratio', 'sqrt_bits')}
    return _finish(function, tolerance, options, attempts,
                   {'error_decomposition': best, 'selection': 'compare_ratio_first_attainable_degrees',
                    'root_ratios_considered': list(map(str, ratios)), 'rounding': rounding}, cost)
