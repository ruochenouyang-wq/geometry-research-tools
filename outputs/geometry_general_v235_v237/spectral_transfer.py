"""V235--V237: certified function error transported to full-sphere spectra.

The function model proves its own remainder. A source is a verified complete
sphere certificate, never a finite Ritz value or single angular sector.
All norms use probability area and all transformations keep the same space.
"""
from math import isqrt
from copy import deepcopy
from dependencies import F, exact, same, oldmath as a, constraints, wide
import function_models as functions

SOURCE = 'full_sphere_polynomial_source_v235'
LINF = 'linf_function_spectrum_v235'
L2 = 'l2_function_spectrum_v236'
FIXED = 'fixed_energy_function_spectrum_v237'
COERCIVITY = 'fixed_reference_coercivity_v237'
MEASURE = 'd_sigma/(4*pi)_on_unit_S2'
EMBEDDING = {'domain': 'all_real_H1_functions_on_unit_S2', 'measure': MEASURE,
             'inequality': 'norm_L4_squared <= Dirichlet_energy + L2_mass',
             'energy_coefficient': '1', 'mass_coefficient': '1',
             'source': 'https://arxiv.org/pdf/1210.1853; equation (1), d=2,p=4'}


def integer(x, low, high, name):
    if type(x) is not int or not low <= x <= high:
        raise ValueError(f'{name} must be an integer in {low}..{high}')
    return x


def poly(raw):
    p = a.polynomial(raw, maximum_degree=24)
    if any(abs(v) > 10**12 or v.denominator > 2**4096 for v in p.values()):
        raise ValueError('Polynomial coefficient budget exceeded')
    return a.encode(p)


def sphere_scope(mean_zero):
    if type(mean_zero) is not bool:
        raise ValueError('mean_zero must be a boolean')
    return 'all_real_mean_zero_H1_on_unit_S2' if mean_zero else 'all_real_H1_on_unit_S2'


def sqrt_upper(value, bits=40):
    value = exact(value)
    integer(bits, 8, 160, 'sqrt_bits')
    if value < 0: raise ValueError('Negative squared norm')
    scale = 2**bits
    ceil_square = (value.numerator*scale**2 + value.denominator-1)//value.denominator
    n = isqrt(ceil_square)
    if n*n < ceil_square: n += 1
    return F(n, scale)


def axis_polynomial(raw):
    p = a.polynomial(poly(raw), maximum_degree=24)
    active = {j for e in p for j in range(3) if e[j]}
    if len(active) > 1: return None
    axis = next(iter(active)) if active else 2
    degree = max((sum(e) for e in p), default=0)
    q = [F(0)]*(degree+1)
    for e, c in p.items(): q[e[axis]] = c
    return axis, list(map(str, q))


def _range_record(p, kind, proof):
    p = poly(p)
    if kind == 'axis_isometry':
        ap = axis_polynomial(p)
        if ap is None: raise ValueError('Polynomial is not one-axis')
        axis, q = ap
        lo, hi = wide.range_bounds(q, proof)
        return {'kind': kind, 'polynomial': p, 'axis': axis, 'proof': proof,
                'lower': str(lo), 'upper': str(hi),
                'identity': 'coordinate_permutation_isometry_of_full_S2'}
    if kind == 'coefficient_bound':
        pp = a.polynomial(p, maximum_degree=24)
        c = pp.get((0,0,0), F(0))
        radius = sum((abs(v) for e,v in pp.items() if e != (0,0,0)), F(0))
        if proof is not None: raise ValueError('Unexpected coefficient-bound payload')
        return {'kind': kind, 'polynomial': p, 'proof': None,
                'lower': str(c-radius), 'upper': str(c+radius),
                'identity': 'abs_each_coordinate_at_most_one_on_unit_S2'}
    raise ValueError('Unknown pointwise range evidence')


def pointwise_range(p):
    ap = axis_polynomial(p)
    if ap is not None:
        return _range_record(p, 'axis_isometry', wide.potential_range(ap[1]))
    return _range_record(p, 'coefficient_bound', None)


def verify_range(c, expected_polynomial=None):
    try:
        if not isinstance(c,dict): return False
        if expected_polynomial is not None and poly(expected_polynomial) != c['polynomial']: return False
        return same(c, _range_record(c['polynomial'], c['kind'], c['proof']))
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError): return False


def _source_record(p, mean_zero, backend, proof):
    p = poly(p)
    scope = sphere_scope(mean_zero)
    if backend == 'axis_isometry_wide':
        ap = axis_polynomial(p)
        if ap is None: raise ValueError('Invalid axis reduction')
        axis, q = ap
        if not wide.verify_full(proof, expected_q=q, expected_mean_zero=mean_zero):
            raise ValueError('Full-sphere source proof failed')
        r = _range_record(p, 'axis_isometry', proof['range_proof'])
        conversion = {'axis': axis, 'identity': 'coordinate_permutation_preserves_L2_energy_and_mean'}
    elif backend == 'cartesian_constraints':
        expected_constraints = [{'0,0,0': '1'}] if mean_zero else []
        if proof.get('format') != constraints.SPECTRAL or not constraints.verify(
                proof, expected_q=p, expected_constraints=expected_constraints):
            raise ValueError('Source must cover the requested complete sphere space')
        r = _range_record(p, 'coefficient_bound', None)
        conversion = {'constraints': expected_constraints,
                      'identity': 'complete_real_cartesian_harmonic_basis'}
    else: raise ValueError('Unknown source backend')
    lo, hi = exact(proof['lower']), exact(proof['upper'])
    if lo > hi: raise ValueError('Inverted source interval')
    return {'format': SOURCE, 'geometry': 'unit_S2', 'measure': MEASURE,
            'scope': scope, 'mean_zero': mean_zero, 'eigenvalue_index': 1,
            'polynomial': p, 'backend': backend, 'conversion': conversion,
            'source_proof': proof, 'pointwise_range': r,
            'lower': str(lo), 'upper': str(hi), 'exact_width': str(hi-lo),
            'full_infinite_space_covered': True}


def polynomial_source(p, mean_zero=True, modes=8, max_modes=16, max_m=8,
                      bits=40, L=3, tolerance='1/10000000000'):
    """Use frozen engines; original polynomial and full-space scope are rebound."""
    p = poly(p); sphere_scope(mean_zero)
    ap = axis_polynomial(p)
    if ap is not None:
        proof = wide.full_ground(ap[1], mean_zero=mean_zero, modes=modes,
                                 max_modes=max_modes, max_m=max_m, bits=bits,
                                 tolerance=exact(tolerance))
        return _source_record(p, mean_zero, 'axis_isometry_wide', proof)
    proof = constraints.certify(p, [{'0,0,0':'1'}] if mean_zero else [], L=L, bits=bits)
    return _source_record(p, mean_zero, 'cartesian_constraints', proof)


def verify_source(c, expected_polynomial=None, expected_mean_zero=None):
    try:
        if not isinstance(c,dict): return False
        if c.get('format') != SOURCE: return False
        if expected_polynomial is not None and poly(expected_polynomial) != c['polynomial']: return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool or
                                             c['mean_zero'] is not expected_mean_zero): return False
        return same(c, _source_record(c['polynomial'], c['mean_zero'], c['backend'], c['source_proof']))
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError): return False


def fixed_reference(p, reference, sqrt_bits=40, range_evidence=None):
    """V237: prove coercivity relative to a fixed bounded polynomial reference."""
    p, reference = poly(p), poly(reference)
    integer(sqrt_bits, 8, 160, 'sqrt_bits')
    r = pointwise_range(reference) if range_evidence is None else range_evidence
    if not verify_range(r, reference): raise ValueError('Reference range proof failed')
    difference = a.add(a.polynomial(p, maximum_degree=24),
                       a.scale(a.polynomial(reference, maximum_degree=24), F(-1)))
    theta2 = a.inner(difference, difference)
    theta = sqrt_upper(theta2, sqrt_bits)
    if theta >= 1: raise ValueError('Fixed-reference distance must certify theta < 1')
    shift = max(F(0), 1-exact(r['lower']))
    return {'format': COERCIVITY, 'polynomial': p, 'reference': reference,
            'measure': MEASURE, 'range_evidence': r, 'sqrt_bits': sqrt_bits,
            'distance_squared': str(theta2), 'distance_upper': str(theta),
            'shift': str(shift), 'coercivity_factor': str(1-theta),
            'embedding': deepcopy(EMBEDDING),
            'inequality': 'a_p + shift*M >= (1-theta)*(E+M)'}


def verify_fixed_reference(c):
    try:
        if not isinstance(c,dict): return False
        return same(c, fixed_reference(c['polynomial'], c['reference'], c['sqrt_bits'], c['range_evidence']))
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError): return False


def _transfer(kind, model, source, target_width, coercivity=None):
    if not functions.verify_model(model): raise ValueError('Function remainder certificate failed')
    if not verify_source(source, expected_polynomial=model['polynomial']):
        raise ValueError('Source polynomial or full-space proof does not match the function model')
    target = exact(target_width)
    if not F(1,10**30) <= target <= 100:
        raise ValueError('Target width outside 1e-30..100')
    eps = exact(model['error_upper'])
    if eps < 0: raise ValueError('Negative error bound')
    low, high = exact(source['lower']), exact(source['upper'])
    if kind == LINF:
        if model['norm'] != 'Linf' or coercivity is not None:
            raise ValueError('V235 requires a uniform remainder')
        lo, hi = low-eps, high+eps
        error = {'uniform_error': str(eps), 'rule': 'minmax_additive_potential_error'}
    elif kind in (L2, FIXED):
        if model['norm'] != 'L2_probability': raise ValueError('Probability-L2 remainder required')
        if kind == L2:
            if coercivity is not None: raise ValueError('Unexpected reference evidence')
            shift = max(F(0), 1-exact(source['pointwise_range']['lower']))
            rho = eps
            positive_floor = F(1)
        else:
            if not verify_fixed_reference(coercivity) or coercivity['polynomial'] != source['polynomial']:
                raise ValueError('Wrong or invalid fixed-reference evidence')
            shift = exact(coercivity['shift'])
            rho = eps/exact(coercivity['coercivity_factor'])
            positive_floor = exact(coercivity['coercivity_factor'])
        if rho >= 1: raise ValueError('Relative form error must be strictly below one')
        shifted_low = max(low+shift, positive_floor)
        if high+shift < positive_floor: raise ArithmeticError('Source contradicts certified coercivity')
        lo, hi = (1-rho)*shifted_low-shift, (1+rho)*(high+shift)-shift
        error = {'probability_L2_error': str(eps), 'relative_form_error': str(rho),
                 'shift': str(shift), 'embedding': deepcopy(EMBEDDING),
                 'coercivity_shifted_floor': str(positive_floor),
                 'source_shifted_lower_used': str(shifted_low),
                 'rule': 'minmax_relative_closed_form_comparison'}
    else: raise ValueError('Unknown transfer format')
    return {'format': kind, 'geometry': 'unit_S2', 'measure': MEASURE,
            'function': model['function'], 'approximant': source['polynomial'],
            'mean_zero': source['mean_zero'], 'scope': source['scope'], 'eigenvalue_index': 1,
            'model': model, 'source': source, 'coercivity': coercivity,
            'error_evidence': error, 'lower': str(lo), 'upper': str(hi),
            'exact_width': str(hi-lo), 'requested_width': str(target),
            'status': 'target_met' if hi-lo <= target else 'certified_open',
            'error_budget': {'source_width': str(high-low),
                             'lower_widening': str(low-lo), 'upper_widening': str(hi-high),
                             'representation_widening': str(hi-lo-(high-low))},
            'full_infinite_space_covered': True, 'formal_proof_assistant_checked': False}


def transfer_linf(model, source, target_width='1/100000000'):
    return _transfer(LINF, model, source, target_width)


def transfer_l2(model, source, target_width='1/100000000'):
    return _transfer(L2, model, source, target_width)


def transfer_fixed(model, source, reference, target_width='1/100000000', sqrt_bits=40):
    c = fixed_reference(model['polynomial'], reference, sqrt_bits)
    return _transfer(FIXED, model, source, target_width, c)


def verify(c, expected_function=None, expected_mean_zero=None, expected_width=None):
    try:
        if not isinstance(c,dict): return False
        if expected_function is not None and functions.normalize_function(expected_function) != c['function']: return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool or c['mean_zero'] is not expected_mean_zero): return False
        if expected_width is not None and exact(expected_width) != exact(c['requested_width']): return False
        return same(c, _transfer(c['format'], c['model'], c['source'], c['requested_width'], c['coercivity']))
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError): return False


def solve(function, method='linf', order=12, degree=8, mean_zero=True,
          target_width='1/100000000', reference=None, sqrt_bits=40, source_options=None):
    """End-to-end original-function evaluation; budget failures stay explicit."""
    if method not in ('linf','l2','fixed'): raise ValueError('Unknown transfer method')
    options = {} if source_options is None else source_options
    allowed = {'modes','max_modes','max_m','bits','L','tolerance'}
    if not isinstance(options,dict) or set(options)-allowed: raise ValueError('Unknown source options')
    if method != 'fixed' and reference is not None: raise ValueError('Reference applies only to fixed method')
    model = functions.analytic_model(function, order=order) if method == 'linf' else functions.l2_model(function, degree=degree, sqrt_bits=sqrt_bits)
    source = polynomial_source(model['polynomial'], mean_zero=mean_zero, **options)
    if method == 'linf': return transfer_linf(model, source, target_width)
    if method == 'l2': return transfer_l2(model, source, target_width)
    if reference is None: raise ValueError('A fixed reference polynomial must be supplied')
    return transfer_fixed(model, source, reference, target_width, sqrt_bits)
