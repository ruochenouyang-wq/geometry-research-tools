"""Second radial eigenvalue lower bounds from exact Schur-complement inertia.

The finite comparison retains every omitted mode via its exact residual Gram.
It is the same infinite-dimensional comparison used by the frozen ground-state
checker; the new observation is that at most one negative direction suffices
to bound the second eigenvalue of a SINGLE real Fourier radial component.
"""
from copy import deepcopy
from pathlib import Path
import sys

sys.dont_write_bytecode = True
_OLD = Path(__file__).resolve().parent.parent / 'geometry_strategy_v238_v317'
if str(_OLD) not in sys.path:
    sys.path.insert(0, str(_OLD))
from backend import F, canonical, direct

FORMAT = 'complete_radial_second_inertia_v1'
GROUND_FORMAT = 'complete_radial_ground_form_v1'


def _json_native(value, depth=0):
    if depth > 64:
        return False
    if type(value) in (str, int, bool, type(None)):
        return True
    if type(value) is list:
        return all(_json_native(v, depth+1) for v in value)
    if type(value) is dict:
        return all(type(k) is str and _json_native(v, depth+1) for k, v in value.items())
    return False


def _kernel(function, m, mean_zero, modes, near_tail, sqrt_bits,
            form_certificate=None, form_verifier=None):
    q = direct.normalize(function)
    if q['profile'] == 'abs_power' and F(q['exponent']) <= -F(1, 2):
        raise ValueError('Exact complete residual requires exponent > -1/2')
    direct.integer(sqrt_bits, 8, 160, 'sqrt_bits')
    if form_certificate is None:
        return direct.Kernel(q, m, mean_zero, modes, sqrt_bits, near_tail)
    if not _json_native(form_certificate) or not callable(form_verifier):
        raise ValueError('A custom form requires its independent verifier')
    fc = deepcopy(form_certificate)
    if form_verifier(fc, expected_function=q) is not True:
        raise ValueError('Custom form does not replay for the exact function')
    if fc.get('function') != q:
        raise ValueError('Custom form must explicitly bind the normalized function')
    form = deepcopy(fc['form'])
    if F(form['energy_factor']) <= 0:
        raise ValueError('Positive form energy factor is required')
    F(form['mass_offset'])
    # The frozen constructor insists on its own centered form. Construct the
    # same validated exact moment kernel with the separately checked form.
    k = direct.Kernel.__new__(direct.Kernel)
    k.function = q
    k.m, k.mean_zero, k.n, k.near_tail = m, mean_zero, modes, near_tail
    k.sqrt_bits = sqrt_bits
    k.data = direct.matrix_assembly(q, m, mean_zero, modes, near_tail)
    k.form = form
    k.start = k.data['degrees'][0]
    k.tail_start = k.start + modes
    k.beta = direct.tail_bound(form, k.tail_start)
    k.remainder_beta = direct.tail_bound(form, k.tail_start + near_tail)
    return k


def _certificate(k, lower, form_certificate=None):
    low = direct.f.rational(lower)
    comparison = k.matrix(low, 'lower')  # also checks strict positive tail
    counts = direct.base.inertia(comparison)
    if counts[0] > 1:
        raise ValueError('More than one negative comparison direction')
    return {
        'format': FORMAT, 'function': k.function, 'geometry': 'unit_S2',
        'axis_coordinate': k.function['axis'],
        'axis_isometry': 'coordinate_permutation_preserves_energy_mass_and_mean',
        'scope': 'one_real_radial_fourier_component', 'azimuth_m': k.m,
        'mean_zero': k.mean_zero, 'projection': 'remove_l0' if k.m == 0 and k.mean_zero else 'none',
        'eigenvalue_index': 2, 'modes': k.n, 'near_tail': k.near_tail,
        'sqrt_bits': k.sqrt_bits, 'lower': str(low), 'beta': str(low),
        'form_bound': deepcopy(k.form), 'form_certificate': deepcopy(form_certificate),
        'radial_tail_start': k.tail_start, 'radial_tail_lower': str(k.beta),
        'remainder_tail_lower': str(k.remainder_beta),
        'kernel_evidence': k.evidence(), 'comparison_matrix': direct.encode(comparison),
        'comparison_inertia': list(counts), 'max_negative_directions': 1,
        'full_radial_space_covered': True,
        'sine_cosine_equivalence': k.m > 0,
        'original_function_moments_exact': True,
        'function_approximation_error': '0', 'formal_proof_assistant_checked': False,
    }


def certified_gap(function, m=1, mean_zero=True, modes=8, near_tail=8,
                  bits=24, sqrt_bits=40, form_certificate=None, form_verifier=None):
    """Return beta <= lambda_2 of one whole real radial m component.

    No numerical eigenvalue is trusted. Every bisection decision uses exact
    rational inertia, and the returned endpoint is independently recounted.
    The upper search endpoint is a tail barrier, NOT an eigenvalue upper bound.
    """
    direct.integer(bits, 8, 96, 'bits')
    k = _kernel(function, m, mean_zero, modes, near_tail, sqrt_bits,
                form_certificate, form_verifier)
    analytic = direct.tail_bound(k.form, k.start)
    distance = F(1)
    for _ in range(64):
        low = analytic - distance
        if low < k.beta and k.count(low, 'lower', independent=True)[0] <= 1:
            break
        distance *= 2
    else:
        raise ValueError('Second-eigenvalue lower bracket budget exhausted')
    high = k.beta
    for _ in range(bits):
        mid = (low + high) / 2
        if k.count(mid, 'lower', independent=True)[0] <= 1:
            low = mid
        else:
            high = mid
    return _certificate(k, low, form_certificate)


def verify_gap(cert, expected_function=None, expected_m=None,
               expected_mean_zero=None, form_verifier=None):
    """Rebuild all original moments, tail terms, comparison, and exact inertia."""
    try:
        if not _json_native(cert) or type(cert) is not dict or cert.get('format') != FORMAT:
            return False
        if expected_function is not None and direct.normalize(expected_function) != cert['function']:
            return False
        if expected_m is not None and (type(expected_m) is not int or expected_m != cert['azimuth_m']):
            return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool or expected_mean_zero is not cert['mean_zero']):
            return False
        k = _kernel(cert['function'], cert['azimuth_m'], cert['mean_zero'],
                    cert['modes'], cert['near_tail'], cert['sqrt_bits'],
                    cert['form_certificate'], form_verifier)
        rebuilt = _certificate(k, cert['lower'], cert['form_certificate'])
        return canonical(cert) == canonical(rebuilt)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, OverflowError):
        return False


def _ground_certificate(k, lower, upper, form_certificate=None):
    cert = direct.sector_certificate(k, lower, upper, independent=True)
    cert['format'] = GROUND_FORMAT
    cert['form_certificate'] = deepcopy(form_certificate)
    return cert


def certified_ground(function, m=0, mean_zero=True, modes=12, near_tail=12,
                     bits=24, sqrt_bits=40, form_certificate=None, form_verifier=None):
    """A full radial ground interval, also supporting independently checked forms."""
    direct.integer(bits, 8, 96, 'bits')
    k = _kernel(function, m, mean_zero, modes, near_tail, sqrt_bits,
                form_certificate, form_verifier)
    analytic = direct.tail_bound(k.form, k.start)
    high = min(k.data['A'][i][i]/k.data['mass'][i] for i in range(k.n))
    low = analytic - 1
    if not k.upper_holds(high):
        raise ArithmeticError('Ritz upper bracket failed')
    for _ in range(bits):
        mid = (low+high)/2
        if k.upper_holds(mid):
            high = mid
        else:
            low = mid
    upper = high
    distance = F(1)
    for _ in range(64):
        low = analytic-distance
        if k.lower_holds(low):
            break
        distance *= 2
    else:
        raise ValueError('Ground lower bracket budget exhausted')
    high = min(upper, k.beta)
    for _ in range(bits):
        mid = (low+high)/2
        if k.lower_holds(mid):
            low = mid
        else:
            high = mid
    return _ground_certificate(k, low, upper, form_certificate)


def verify_ground(cert, expected_function=None, expected_m=None,
                  expected_mean_zero=None, form_verifier=None):
    """Replay custom-form ground evidence without invoking the old form policy."""
    try:
        if not _json_native(cert) or type(cert) is not dict or cert.get('format') != GROUND_FORMAT:
            return False
        if expected_function is not None and direct.normalize(expected_function) != cert['function']:
            return False
        if expected_m is not None and (type(expected_m) is not int or expected_m != cert['azimuth_m']):
            return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool or expected_mean_zero is not cert['mean_zero']):
            return False
        k = _kernel(cert['function'], cert['azimuth_m'], cert['mean_zero'],
                    cert['modes'], cert['near_tail'], cert['sqrt_bits'],
                    cert['form_certificate'], form_verifier)
        rebuilt = _ground_certificate(k, cert['lower'], cert['upper'], cert['form_certificate'])
        return canonical(cert) == canonical(rebuilt)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, OverflowError):
        return False
