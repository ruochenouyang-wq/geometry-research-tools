"""Replayable L^(3/2) clipped-potential bounds for negative axis powers.

Only the proposal of two rational clipping bases uses floating point.  Every
accepted bound, cube-root enclosure and choice between those bases is exact.
"""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import math
import sys

sys.dont_write_bytecode = True
OLD = Path(__file__).resolve().parents[1] / 'geometry_strategy_v238_v317'
if str(OLD) not in sys.path:
    sys.path.insert(0, str(OLD))
from backend import direct

FORMAT = 'negative_power_clipped_holder_form_v1'
SELECTED_FORMAT = 'selected_negative_power_holder_form_v1'
GRID = 1000


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def _native(value, depth=0):
    if depth > 32:
        return False
    if type(value) in (str, int, bool, type(None)):
        return True
    if type(value) is list:
        return all(_native(x, depth+1) for x in value)
    if type(value) is dict:
        return all(type(k) is str and _native(v, depth+1) for k, v in value.items())
    return False


def normalize(function):
    q = direct.normalize(function)
    if q['profile'] != 'abs_power' or not -F(1, 2) < F(q['exponent']) < 0:
        raise ValueError('Complete residual route requires -1/2 < alpha < 0')
    if F(q['amplitude']) >= 0:
        raise ValueError('Require a negative amplitude')
    if F(q['exponent']).denominator > 64:
        raise ValueError('Exact clipping exponent denominator exceeds 64')
    return q


def cube_root_upper(value, bits=48):
    """Least nonnegative dyadic with denominator 2**bits whose cube >= value."""
    q = direct.f.rational(value)
    direct.integer(bits, 8, 160, 'root_bits')
    if q < 0:
        raise ValueError('Cube-root enclosure requires a nonnegative value')
    scaled_numerator = q.numerator << (3*bits)
    denominator = q.denominator
    integer_target = scaled_numerator//denominator
    low, high = 0, 1 << ((integer_target.bit_length()+2)//3)
    while low < high:
        mid = (low+high)//2
        if mid**3*denominator >= scaled_numerator:
            high = mid
        else:
            low = mid+1
    result = F(low, 1 << bits)
    if result**3 < q or (low > 0 and F(low-1, 1 << bits)**3 >= q):
        raise ArithmeticError('Outward cube-root invariant failed')
    return result


def _base(raw):
    value = direct.f.rational(raw)
    if not F(1) < value <= 8 or value.denominator > 1000000:
        raise ValueError('Clipping base must lie in (1,8], denominator <= 1000000')
    return value


def _raw_certificate(function, base, root_bits=48):
    q, base = normalize(function), _base(base)
    direct.integer(root_bits, 8, 160, 'root_bits')
    p, A, b = -F(q['exponent']), -F(q['amplitude']), F(q['offset'])
    r, s = p.numerator, p.denominator
    delta, cap = base**(-s), A*base**r
    D = 1/(1-F(3, 2)*p)-F(3, 2)/(1-p/2)+F(1, 2)/(1+p/2)
    if D <= 0:
        raise ArithmeticError('Integrated clipped-residual upper bound is not positive')
    eta_cube = A**3*base**(3*r-2*s)*D**2
    eta = cube_root_upper(eta_cube, root_bits)
    return {
        'format': FORMAT, 'function': q, 'geometry': 'unit_S2',
        'axis_coordinate': q['axis'], 'measure': 'd_sigma/(4*pi)',
        'scope': 'all_real_H1_on_unit_S2', 'root_bits': root_bits,
        'proof': {
            'base': str(base), 'delta': str(delta), 'cap': str(cap),
            'residual': '(abs(amplitude)*abs(axis_coordinate)^(-p)-cap)_+',
            'holder_exponent': '3/2', 'dual_function_exponent': '6',
            'sqrt_upper': 'sqrt(1-t)<=1-t/2 for 0<=t<=1',
            'integration_region': '0<=abs(axis_coordinate)<=delta_only',
            'integral_coefficient_D': str(D),
            'residual_norm_cube_upper': str(eta_cube),
            'eta_upper': str(eta), 'eta_upper_cube': str(eta**3),
            'previous_dyadic_cube': str((eta-F(1, 1 << root_bits))**3),
            'embedding': 'norm_L6_squared <= 2*Dirichlet_energy + L2_mass',
            'embedding_source': 'https://arxiv.org/pdf/1210.1853; equation (1), d=2,p=6',
        },
        'form': {'energy_factor': str(1-2*eta), 'mass_offset': str(b-cap-eta)},
        'positive_energy_factor': eta < F(1, 2),
    }


def certificate(function, base, root_bits=48):
    cert = _raw_certificate(function, base, root_bits)
    if not cert['positive_energy_factor']:
        raise ValueError('Clipped Holder bound has no positive energy factor')
    return cert


def _select_from_bases(function, tail_start, bases, root_bits=48):
    q = normalize(function)
    direct.integer(tail_start, 1, 96, 'tail_start')
    if not isinstance(bases, (list, tuple)) or len(bases) != 2:
        raise ValueError('Exactly two adjacent rational grid proposals are required')
    bases = [_base(x) for x in bases]
    if bases[1]-bases[0] != F(1, GRID) or any((b*GRID).denominator != 1 for b in bases):
        raise ValueError('Proposals must be adjacent on the fixed 1/1000 grid')
    attempts, accepted = [], []
    for base in bases:
        cert = _raw_certificate(q, base, root_bits)
        barrier = direct.tail_bound(cert['form'], tail_start)
        attempts.append({'base': str(base), 'eta_upper': cert['proof']['eta_upper'],
                         'energy_factor': cert['form']['energy_factor'],
                         'mass_offset': cert['form']['mass_offset'],
                         'initial_tail_lower': str(barrier),
                         'admissible': cert['positive_energy_factor']})
        if cert['positive_energy_factor']:
            accepted.append((barrier, -base, cert))
    if not accepted:
        error = ValueError('Both declared Holder form proposals failed positive-energy acceptance')
        error.attempts = attempts
        raise error
    selected = max(accepted, key=lambda x: (x[0], x[1]))[2]
    return {
        'format': SELECTED_FORMAT, 'function': q, 'geometry': 'unit_S2',
        'axis_coordinate': q['axis'], 'measure': 'd_sigma/(4*pi)',
        'scope': 'all_real_H1_on_unit_S2', 'root_bits': root_bits,
        'form': deepcopy(selected['form']), 'selected_form': selected,
        'selection': {
            'policy': 'maximum_exact_initial_tail_lower_of_two_adjacent_rational_bases',
            'tail_start': tail_start, 'grid_denominator': GRID,
            'candidate_bases': list(map(str, bases)), 'attempts': attempts,
            'selected_base': selected['proof']['base'],
            'continuous_stationary_point_is_proposal_only': True,
            'optimal_over_all_bases_claimed': False,
        },
    }


def select_form(function, tail_start, root_bits=48):
    """Propose around the stationary point; exact proof/selection then decides."""
    q = normalize(function)
    direct.integer(tail_start, 1, 96, 'tail_start')
    p = -F(q['exponent'])
    r, s = p.numerator, p.denominator
    D = 1/(1-F(3, 2)*p)-F(3, 2)/(1-p/2)+F(1, 2)/(1+p/2)
    d, T = F(2*s-3*r, 3), tail_start*(tail_start+1)
    # B**(2*s) = ((2*T+1)*d/r)**3 * D**2. This floating point
    # stationary value is never used as a proof premise or an accepted bound.
    logbase = (3*math.log(float((2*T+1)*d/r))+2*math.log(float(D)))/(2*s)
    proposal = math.exp(logbase)
    k = min(8*GRID-1, max(GRID+1, math.floor(proposal*GRID)))
    return _select_from_bases(q, tail_start, [F(k, GRID), F(k+1, GRID)], root_bits)


def verify_form(cert, expected_function=None):
    """No trust in recorded cubes, integrals, selected base or scalar form."""
    try:
        if not _native(cert) or type(cert) is not dict:
            return False
        if expected_function is not None and normalize(expected_function) != cert['function']:
            return False
        if cert.get('format') == FORMAT:
            rebuilt = certificate(cert['function'], cert['proof']['base'], cert['root_bits'])
        elif cert.get('format') == SELECTED_FORMAT:
            selection = cert['selection']
            rebuilt = _select_from_bases(cert['function'], selection['tail_start'],
                                          selection['candidate_bases'], cert['root_bits'])
        else:
            return False
        return canonical(rebuilt) == canonical(cert)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, OverflowError):
        return False
