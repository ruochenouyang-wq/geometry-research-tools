"""Exact probability-L2 piecewise models for a step and |t|**(-1/4).

The model is a function approximation, not a polynomial spectral certificate.
Only rational arithmetic and finite polynomial integration enter its replay.
"""
from fractions import Fraction as F
from math import isqrt
import json
import re

FORMAT = 'sphere_structured_piecewise_L2_v1'
SCOPE = 'unit_S2_probability_measure'
MAX_LEVELS = 64
MAX_DEGREE = 32


def rational(value):
    if type(value) in (int, F):
        return F(value)
    if type(value) is str and re.fullmatch(r'-?\d+(?:/[1-9]\d*)?', value):
        return F(value)
    raise ValueError('Use exact integers or rational strings, not floats or booleans')


def _integer(value, low, high, name):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(name + ' is outside its integer budget')
    return value


def normalize_function(function):
    if type(function) is not dict or function.get('kind') != 'axis_profile':
        raise ValueError('Expected an axis_profile function')
    profile = function.get('profile')
    if profile not in ('step', 'abs_power'):
        raise ValueError('Only step and absolute power -1/4 are supported')
    allowed = {'kind', 'axis', 'profile', 'amplitude', 'offset'}
    if profile == 'abs_power':
        allowed.add('exponent')
    if not set(function) <= allowed:
        raise ValueError('Unknown source fields')
    result = {'kind': 'axis_profile', 'axis': _integer(function.get('axis', 2), 0, 2, 'axis'),
              'profile': profile, 'amplitude': str(rational(function.get('amplitude', 1))),
              'offset': str(rational(function.get('offset', 0)))}
    if profile == 'abs_power':
        exponent = rational(function.get('exponent', 1))
        if exponent != F(-1, 4):
            raise ValueError('This structured model supports only exponent -1/4')
        result['exponent'] = str(exponent)
    return result


def sqrt_upper(value_squared, bits=40):
    value = rational(value_squared)
    _integer(bits, 0, 256, 'sqrt_bits')
    if value < 0:
        raise ValueError('Negative squared error')
    denominator = 1 << bits
    numerator = isqrt(value.numerator * denominator**2 // value.denominator)
    if F(numerator**2, denominator**2) < value:
        numerator += 1
    return F(numerator, denominator)


def _add(p, q):
    result = [F(0)] * max(len(p), len(q))
    for i, a in enumerate(p): result[i] += a
    for i, a in enumerate(q): result[i] += a
    return result


def _scale(p, c): return [c*a for a in p]


def _mul(p, q):
    result = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): result[i+j] += a*b
    return result


def _compose(p, affine):
    result = [F(0)]
    for a in reversed(p):
        result = _add(_mul(result, affine), [a])
    while len(result) > 1 and not result[-1]: result.pop()
    return result


def _shifted_legendre(degree):
    basis = [[F(1)]]
    if degree: basis.append([F(-1), F(2)])
    for l in range(1, degree):
        basis.append(_scale(_add(_scale(_mul([F(-1), F(2)], basis[-1]), 2*l+1),
                                 _scale(basis[-2], -l)), F(1, l+1)))
    return basis


def _integrate_polynomial(p, a, b):
    return sum((c*(b**(k+1)-a**(k+1))/F(k+1) for k, c in enumerate(p)), F(0))


def _power_moment(k, left_root, right_root):
    # t = root**4, so t**(k+3/4) = root**(4*k+3), including root=0.
    return F(4, 4*k+3) * (right_root**(4*k+3)-left_root**(4*k+3))


def _cell(left_root, right_root, degree, amplitude, offset, core=False):
    a, b = left_root**4, right_root**4
    width = b-a
    moments = [offset*(b**(k+1)-a**(k+1))/F(k+1)
               + amplitude*_power_moment(k, left_root, right_root)
               for k in range(degree+1)]
    source_squared = (offset**2*width + 2*offset*amplitude*_power_moment(0, left_root, right_root)
                      + 2*amplitude**2*(right_root**2-left_root**2))
    entries = []
    local_polynomial = [F(0)]
    approximation_squared = F(0)
    for l, local_basis in enumerate(_shifted_legendre(degree)):
        in_t = _compose(local_basis, [-a/width, 1/width])
        inner = sum((c*moments[k] for k, c in enumerate(in_t)), F(0))
        coefficient = (2*l+1)*inner/width
        mass = width/F(2*l+1)
        contribution = coefficient**2*mass
        approximation_squared += contribution
        local_polynomial = _add(local_polynomial, _scale(local_basis, coefficient))
        entries.append({'degree': l, 'shifted_legendre_in_s': list(map(str, local_basis)),
                        'source_inner_product_dt': str(inner), 'basis_mass_dt': str(mass),
                        'coefficient': str(coefficient), 'norm_squared_dt': str(contribution)})
    # Independently integrate the actual returned polynomial and cross product.
    polynomial_in_t = _compose(local_polynomial, [-a/width, 1/width])
    direct_squared = _integrate_polynomial(_mul(polynomial_in_t, polynomial_in_t), a, b)
    cross = sum((c*moments[k] for k, c in enumerate(polynomial_in_t)), F(0))
    residual = source_squared - 2*cross + direct_squared
    if cross != approximation_squared or direct_squared != approximation_squared or residual < 0:
        raise ArithmeticError('Exact local projection identity failed')
    return {'left': str(a), 'right': str(b), 'left_fourth_root': str(left_root),
            'right_fourth_root': str(right_root), 'core': core, 'degree': degree,
            'coordinate': 's=(abs(t)-left)/(right-left)',
            'polynomial_in_s': list(map(str, local_polynomial)),
            'source_monomial_moments_dt': list(map(str, moments)),
            'projection': entries, 'source_norm_squared_dt': str(source_squared),
            'approximation_norm_squared_dt': str(direct_squared), 'cross_integral_dt': str(cross),
            'error_squared_dt': str(residual)}


def _scaled_cell(template, right_root, amplitude, offset):
    """Transport one exactly integrated shape under t=right_root**4*x."""
    ratio = F(template['left_fourth_root'])
    left_root = ratio*right_root
    a, b = left_root**4, right_root**4
    width = b-a
    degree = template['degree']
    local_polynomial = _add(_scale(list(map(F, template['polynomial_in_s'])), amplitude/right_root), [offset])
    moments = [offset*(b**(k+1)-a**(k+1))/F(k+1)
               + amplitude*right_root**(4*k+3)*F(template['source_monomial_moments_dt'][k])
               for k in range(degree+1)]
    entries = []
    projection_norm = F(0)
    for row in template['projection']:
        l = row['degree']
        coefficient = amplitude/right_root*F(row['coefficient']) + (offset if l == 0 else 0)
        inner = amplitude*right_root**3*F(row['source_inner_product_dt']) + (offset*width if l == 0 else 0)
        mass = right_root**4*F(row['basis_mass_dt'])
        contribution = coefficient**2*mass
        if coefficient*mass != inner:
            raise ArithmeticError('Scaled projection coefficient identity failed')
        projection_norm += contribution
        entries.append({'degree': l, 'shifted_legendre_in_s': row['shifted_legendre_in_s'],
                        'source_inner_product_dt': str(inner), 'basis_mass_dt': str(mass),
                        'coefficient': str(coefficient), 'norm_squared_dt': str(contribution)})
    source_norm = (offset**2*width + 2*offset*amplitude*_power_moment(0, left_root, right_root)
                   + 2*amplitude**2*(right_root**2-left_root**2))
    residual = amplitude**2*right_root**2*F(template['error_squared_dt'])
    if source_norm-projection_norm != residual:
        raise ArithmeticError('Scaled full residual identity failed')
    return {'left': str(a), 'right': str(b), 'left_fourth_root': str(left_root),
            'right_fourth_root': str(right_root), 'core': False, 'degree': degree,
            'coordinate': 's=(abs(t)-left)/(right-left)',
            'polynomial_in_s': list(map(str, local_polynomial)),
            'source_monomial_moments_dt': list(map(str, moments)), 'projection': entries,
            'source_norm_squared_dt': str(source_norm),
            'approximation_norm_squared_dt': str(projection_norm), 'cross_integral_dt': str(projection_norm),
            'error_squared_dt': str(residual)}


def piecewise_model(function, levels=4, degree=3, root_ratio='1/2', sqrt_bits=40):
    """Return exact whole-sphere L2 error of a piecewise polynomial model.

    For abs_power the positive half has one constant core plus `levels`
    degree-`degree` annuli; even reflection represents the other half.
    For step, the two constant halves represent the source exactly.
    """
    function = normalize_function(function)
    _integer(levels, 1, MAX_LEVELS, 'levels')
    _integer(degree, 0, MAX_DEGREE, 'degree')
    _integer(sqrt_bits, 0, 256, 'sqrt_bits')
    ratio = rational(root_ratio)
    if not 0 < ratio < 1 or ratio.denominator > 4096:
        raise ValueError('root_ratio must lie in (0,1) with denominator<=4096')
    amplitude, offset = F(function['amplitude']), F(function['offset'])
    if function['profile'] == 'step':
        cells = [{'left': '-1', 'right': '0', 'polynomial_in_t': [str(offset)]},
                 {'left': '0', 'right': '1', 'polynomial_in_t': [str(offset+amplitude)]}]
        error = F(0)
        parameters = 2
        representation = 'two_halves_with_separate_value_at_zero'
        source_norm = (offset**2+(offset+amplitude)**2)/2
        shape_proof = None
    else:
        shape_proof = _cell(ratio, F(1), degree, F(1), F(0))
        cells = [_cell(F(0), ratio**levels, 0, amplitude, offset, core=True)]
        cells += [_scaled_cell(shape_proof, ratio**j, amplitude, offset)
                  for j in range(levels-1, -1, -1)]
        error = sum((F(cell['error_squared_dt']) for cell in cells), F(0))
        geometric_error = amplitude**2*(F(2, 9)*ratio**(2*levels)
            + F(shape_proof['error_squared_dt'])*(1-ratio**(2*levels))/(1-ratio**2))
        if error != geometric_error:
            raise ArithmeticError('Complete geometric residual sum failed')
        source_norm = offset**2 + F(8, 3)*offset*amplitude + 2*amplitude**2
        if sum((F(cell['source_norm_squared_dt']) for cell in cells), F(0)) != source_norm:
            raise ArithmeticError('Whole-domain source norm failed')
        parameters = 1 + levels*(degree+1)
        representation = 'even_reflection_of_positive_half_cells'
    upper = sqrt_upper(error, sqrt_bits)
    return {'format': FORMAT, 'function': function, 'scope': SCOPE,
            'norm': 'L2_probability', 'error_squared': str(error), 'error_upper': str(upper),
            'sqrt_upper_squared': str(upper**2), 'source_norm_squared': str(source_norm),
            'levels': levels, 'degree': degree, 'root_ratio': str(ratio), 'sqrt_bits': sqrt_bits,
            'representation': representation, 'cells': cells,
            'shape_proof': shape_proof,
            'shape_rule': 't=R^4*x scales centered profile coefficients by 1/R and squared residual integral by R^2',
            'source_value_at_axis_zero': str(offset),
            'axis_zero_convention': 'step(0)=0; singular_power(0)=0 as an a.e. representative',
            'endpoint_convention': 'at nonzero seams choose the cell starting at abs(t); final endpoint included',
            'expanded_polynomial_coefficient_slots': parameters,
            'shape_coefficient_count': 2 if function['profile'] == 'step' else degree+1,
            'bound_statement': 'probability_L2_norm(original_function-piecewise_model)<=error_upper',
            'spectral_transfer_claimed': False, 'global_polynomial': False,
            'formal_assistant_checked': False}


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False)


def verify_piecewise(certificate, expected_function=None):
    try:
        if type(certificate) is not dict or certificate.get('format') != FORMAT:
            return False
        function = normalize_function(certificate['function'])
        if expected_function is not None and function != normalize_function(expected_function):
            return False
        rebuilt = piecewise_model(function, certificate['levels'], certificate['degree'],
                                  certificate['root_ratio'], certificate['sqrt_bits'])
        return _canonical(certificate) == _canonical(rebuilt)
    except (KeyError, ValueError, TypeError, ArithmeticError, OverflowError, IndexError):
        return False


def evaluate(certificate, t):
    """Evaluate the certified representative at one exact axis coordinate."""
    if not verify_piecewise(certificate):
        raise ValueError('Invalid piecewise model')
    t = rational(t)
    if not -1 <= t <= 1:
        raise ValueError('Axis coordinate must be in [-1,1]')
    if certificate['function']['profile'] == 'step':
        return F(certificate['function']['offset']) + (F(certificate['function']['amplitude']) if t > 0 else 0)
    t = abs(t)
    for cell in certificate['cells']:
        a, b = F(cell['left']), F(cell['right'])
        if a <= t < b or t == b == 1:
            s = (t-a)/(b-a)
            value = F(0)
            for coefficient in reversed(cell['polynomial_in_s']): value = value*s + F(coefficient)
            return value
    raise ArithmeticError('A verified cover must contain this coordinate')


def budget_table(function=None, budgets=((1, 1), (2, 2), (4, 3), (5, 3), (6, 3), (8, 4))):
    """Reproducible whole-domain exact errors, with no spectral-accuracy claim."""
    if function is None:
        function = {'kind': 'axis_profile', 'profile': 'abs_power', 'exponent': '-1/4', 'amplitude': '-1'}
    rows = []
    for levels, degree in budgets:
        cert = piecewise_model(function, levels=levels, degree=degree)
        if not verify_piecewise(cert, function): raise ArithmeticError('Budget row failed replay')
        rows.append({'levels': levels, 'degree': degree,
                     'coefficient_slots': cert['expanded_polynomial_coefficient_slots'],
                     'error_squared': cert['error_squared'], 'error_upper': cert['error_upper'],
                     'certificate_bytes': len(_canonical(cert).encode('utf-8'))})
    return rows


if __name__ == '__main__':
    print(json.dumps(budget_table(), indent=2, sort_keys=True))
