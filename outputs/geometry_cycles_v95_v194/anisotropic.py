"""V175--V184: exact coupled harmonics for nonaxisymmetric potentials on S².

All integrals use normalized surface measure. Certificates cover the compressed
operator on ALL real mean-zero H1 functions, with every omitted harmonic degree
bounded. Degree and matrix caps keep this an auditable small prototype.
Analytic dependencies: complete spherical harmonics, their Laplace eigenvalues,
the min-max principle and positive-operator Schur comparison. No formal kernel.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import factorial, comb
import json
from common import exact, same, projected, base, arithmetic

FORMAT = 'anisotropic_mean_zero_sphere_v182'
ROTATION_FORMAT = 'sphere_rotation_transfer_v184'
CHEAP_FORMAT = 'anisotropic_cheap_ground_v175'
SCOPE = 'all_real_mean_zero_H1_functions_on_unit_S2'
ZERO = (0, 0, 0)
R2 = {(2, 0, 0): F(1), (0, 2, 0): F(1), (0, 0, 2): F(1)}


def polynomial(raw, maximum_degree=4):
    """V175: canonical sparse rational xyz polynomial; reject floats and bools."""
    if not isinstance(raw, dict) or len(raw) > 512:
        raise ValueError('Require a sparse exponent-to-rational dictionary')
    result = {}
    for key, value in raw.items():
        if isinstance(key, str):
            parts = key.split(',')
            if len(parts) != 3 or any(not p.isdigit() for p in parts):
                raise ValueError('Exponent keys must be a,b,c with nonnegative integers')
            exponent = tuple(int(p) for p in parts)
        else:
            exponent = key
        if not isinstance(exponent, tuple) or len(exponent) != 3 or any(
                type(e) is not int or e < 0 for e in exponent) or sum(exponent) > maximum_degree:
            raise ValueError('Unsupported exponent or polynomial degree')
        if type(value) not in (int, str, F):
            raise ValueError('Coefficients must be exact rationals')
        value = exact(value)
        if value:
            if exponent in result:
                raise ValueError('Duplicate normalized monomial')
            result[exponent] = value
    return result


def encode(p):
    return {','.join(map(str, e)): str(c) for e, c in sorted(p.items()) if c}


def add(*polynomials):
    result = {}
    for p in polynomials:
        for e, c in p.items():
            result[e] = result.get(e, F(0)) + c
    return {e: c for e, c in result.items() if c}


def scale(p, factor):
    return {e: c*factor for e, c in p.items() if c*factor}


def multiply(p, q):
    result = {}
    for e, c in p.items():
        for f, d in q.items():
            key = tuple(e[i]+f[i] for i in range(3))
            result[key] = result.get(key, F(0)) + c*d
    return {e: c for e, c in result.items() if c}


def power(p, n):
    if type(n) is not int or not 0 <= n <= 24:
        raise ValueError('Power must be an integer in 0..24')
    result = {ZERO: F(1)}
    for _ in range(n):
        result = multiply(result, p)
    return result


def derivative(p, axis):
    if type(axis) is not int or axis not in range(3):
        raise ValueError('Axis must be 0, 1 or 2')
    result = {}
    for e, c in p.items():
        if e[axis]:
            f = list(e)
            f[axis] -= 1
            result[tuple(f)] = c*e[axis]
    return result


def _double_factorial(n):
    result = 1
    while n > 0:
        result *= n
        n -= 2
    return result


@lru_cache(maxsize=4096)
def sphere_moment(exponent):
    """V176: E[x^a y^b z^c]; normalized area, not coordinate cubature."""
    if not isinstance(exponent, tuple) or len(exponent) != 3 or any(
            type(e) is not int or e < 0 for e in exponent):
        raise ValueError('Invalid monomial exponent')
    if any(e % 2 for e in exponent):
        return F(0)
    a, b, c = exponent
    return F(_double_factorial(a-1)*_double_factorial(b-1)*_double_factorial(c-1),
             _double_factorial(a+b+c+1))


def integrate(p):
    return sum((c*sphere_moment(e) for e, c in p.items()), F(0))


def inner(p, q):
    return integrate(multiply(p, q))


def solid_harmonic(l, m, part='real'):
    """V177: rational homogeneous real/imaginary solid associated Legendre basis.

    The harmless Condon--Shortley sign is omitted. On r=1, this is
    (d^m P_l/dz^m) Re/Im[(x+iy)^m]; homogenization inserts r^(2k).
    """
    if type(l) is not int or type(m) is not int or not 0 <= m <= l <= 8:
        raise ValueError('Require 0 <= m <= l <= 8')
    if part not in ('real', 'imag') or (m == 0 and part == 'imag'):
        raise ValueError('Invalid real harmonic component')
    azimuth = {}
    for j in range(m+1):
        if j % 2 == (part == 'imag'):
            azimuth[(m-j, j, 0)] = F(comb(m, j)*(-1)**(j//2))
    result = {}
    for k in range((l-m)//2+1):
        coefficient = F((-1)**k*factorial(2*l-2*k),
                        2**l*factorial(k)*factorial(l-k)*factorial(l-2*k-m))
        term = multiply(azimuth, {(0, 0, l-2*k-m): coefficient})
        result = add(result, multiply(term, power(R2, k)))
    laplacian = add(*(derivative(derivative(result, a), a) for a in range(3)))
    if laplacian:
        raise ArithmeticError('Constructed homogeneous polynomial is not harmonic')
    return result


def harmonic_basis(start, stop):
    if type(start) is not int or type(stop) is not int or not 0 <= start <= stop <= 8:
        raise ValueError('Unsupported harmonic window')
    result = []
    for l in range(start, stop+1):
        for m in range(l+1):
            for part in (('real',) if m == 0 else ('real', 'imag')):
                result.append({'degree': l, 'm': m, 'part': part,
                               'polynomial': solid_harmonic(l, m, part)})
    return result


def tangent_energy(p, q):
    """Exact integral of Euclidean gradient product minus radial components."""
    gradients = [inner(derivative(p, a), derivative(q, a)) for a in range(3)]
    radial_p = {e: sum(e)*c for e, c in p.items()}
    radial_q = {e: sum(e)*c for e, c in q.items()}
    return sum(gradients, F(0))-inner(radial_p, radial_q)


def gram_energy(L=3):
    """V178: exact mass and independently integrated tangential energy forms."""
    if type(L) is not int or not 1 <= L <= 4:
        raise ValueError('Retained harmonic degree must be in 1..4')
    basis = harmonic_basis(1, L)
    n = len(basis)
    gram = [[F(0) for _ in basis] for _ in basis]
    energy = [[F(0) for _ in basis] for _ in basis]
    for i in range(n):
        for j in range(i+1):
            a, b = basis[i], basis[j]
            mass = inner(a['polynomial'], b['polynomial'])
            value = tangent_energy(a['polynomial'], b['polynomial'])
            if (i != j and mass != 0) or (i == j and mass <= 0):
                raise ArithmeticError('Harmonics are not a positive orthogonal basis')
            if value != a['degree']*(a['degree']+1)*mass:
                raise ArithmeticError('Exact gradient calculation disagrees with Laplace eigenvalue')
            gram[i][j] = gram[j][i] = mass
            energy[i][j] = energy[j][i] = value
    return {'basis': basis, 'gram': gram, 'energy': energy,
            'mass': [gram[i][i] for i in range(n)]}


def coupled_form(q, L=3):
    """V179: retain every cross-m/cross-component entry for an xyz potential."""
    q = polynomial(q)
    data = gram_energy(L)
    basis = data['basis']
    n = len(basis)
    potential_form = [[F(0) for _ in basis] for _ in basis]
    products = [multiply(q, b['polynomial']) for b in basis]
    for i in range(n):
        for j in range(i+1):
            value = inner(basis[i]['polynomial'], products[j])
            potential_form[i][j] = potential_form[j][i] = value
    data['potential_form'] = potential_form
    data['A'] = [[data['energy'][i][j]+potential_form[i][j] for j in range(n)] for i in range(n)]
    return data


def tail_couplings(q, L=3, finite=None):
    """V180: all retained-to-discarded entries, including all 2l+1 harmonics.

    A multiplier of degree d cannot couple retained degrees <=L to l>L+d.
    The l=0 coefficient is removed by the mean-zero orthogonal projection.
    """
    q = polynomial(q)
    if type(L) is not int or not 1 <= L <= 4:
        raise ValueError('Retained harmonic degree must be in 1..4')
    finite = coupled_form(q, L) if finite is None else finite
    d = max(map(sum, q), default=0)
    tail = harmonic_basis(L+1, L+d) if d else []
    columns = []
    q_kept = [multiply(q, b['polynomial']) for b in finite['basis']]
    for b in tail:
        p = b['polynomial']
        columns.append({'degree': b['degree'], 'm': b['m'], 'part': b['part'],
                        'mass': inner(p, p),
                        'entries': [inner(p, v) for v in q_kept]})
    return columns


def sphere_reduce(q):
    """V183: exact division q-reduced=(x²+y²+z²-1)*quotient."""
    q = polynomial(q)
    current, quotient = dict(q), {}
    divisor = add(R2, {ZERO: F(-1)})
    while any(e[2] >= 2 for e in current):
        e = max((e for e in current if e[2] >= 2), key=lambda e: (e[2], e))
        factor = {(e[0], e[1], e[2]-2): current[e]}
        quotient = add(quotient, factor)
        current = add(current, scale(multiply(divisor, factor), F(-1)))
    if add(q, scale(current, -1)) != multiply(divisor, quotient):
        raise ArithmeticError('Sphere reduction identity failed')
    return {'original': encode(q), 'reduced': encode(current), 'quotient': encode(quotient),
            'identity': 'original-reduced=(x²+y²+z²-1)*quotient'}


def _coefficient_range(q):
    constant = q.get(ZERO, F(0))
    radius = sum((abs(c) for e, c in q.items() if e != ZERO), F(0))
    return constant-radius, constant+radius


def potential_range(q, reduce=True):
    if type(reduce) is not bool:
        raise ValueError('reduce must be boolean')
    q = polynomial(q)
    reduction = sphere_reduce(q)
    reduced = polynomial(reduction['reduced'])
    lo, hi = _coefficient_range(q)
    newlo, newhi = _coefficient_range(reduced)
    if reduce and newhi-newlo < hi-lo:
        used, lo, hi = reduced, newlo, newhi
    else:
        used = q
    return {'q': encode(q), 'used': encode(used), 'lower': str(lo), 'upper': str(hi),
            'reduction_enabled': reduce, 'sphere_identity': reduction,
            'rule': 'constant plus/minus sum of absolute nonconstant coefficients on unit sphere'}


class CoupledKernel:
    """V181: degree-graded Schur bounds with correct non-unit mass denominators."""
    def __init__(self, q, L=3, reduce=True):
        self.q = polynomial(q)
        self.L = L
        self.range = potential_range(self.q, reduce)
        self.used = polynomial(self.range['used'])
        self.qlo, self.qhi = F(self.range['lower']), F(self.range['upper'])
        self.finite = coupled_form(self.used, L)
        self.a, self.mass = self.finite['A'], self.finite['mass']
        self.columns = tail_couplings(self.used, L, self.finite)
        self.beta = (L+1)*(L+2)+self.qlo

    def matrix(self, x, kind):
        x = exact(x)
        if kind not in ('lower', 'upper_tail', 'upper_ritz'):
            raise ValueError('Unknown comparison kind')
        if kind != 'upper_ritz' and x >= self.beta:
            raise ValueError('All omitted degrees must have positive lower comparison')
        result = [row[:] for row in self.a]
        for i, mass in enumerate(self.mass):
            result[i][i] -= x*mass
        if kind != 'upper_ritz':
            for column in self.columns:
                l, mass = column['degree'], column['mass']
                bound = self.qlo if kind == 'lower' else self.qhi
                denominator = mass*(l*(l+1)+bound-x)
                entries = [(i, c) for i, c in enumerate(column['entries']) if c]
                for i, a in entries:
                    for j, b in entries:
                        result[i][j] -= a*b/denominator
        return result

    def count(self, x, kind, independent=False):
        matrix = self.matrix(x, kind)
        return base.inertia(matrix) if independent else arithmetic.banded_inertia(matrix)

    def lower_holds(self, x):
        return x < self.beta and self.count(x, 'lower')[0] == 0

    def upper_holds(self, x, kind):
        if kind != 'upper_ritz' and x >= self.beta:
            return False
        counts = self.count(x, kind)
        return counts[0]+counts[1] >= 1


def _matrix_json(matrix):
    return [[str(c) for c in row] for row in matrix]


def _record(kernel, lower, upper, kind, independent=False):
    lo, hi = exact(lower), exact(upper)
    if lo > hi or not kernel.lower_holds(lo) or not kernel.upper_holds(hi, kind):
        raise ValueError('Spectral endpoints fail exact Schur comparison')
    lower_counts = kernel.count(lo, 'lower', independent)
    upper_counts = kernel.count(hi, kind, independent)
    if lower_counts[0] != 0 or upper_counts[0]+upper_counts[1] < 1:
        raise ValueError('Independent inertia disagrees')
    labels = [{k: b[k] for k in ('degree', 'm', 'part')} for b in kernel.finite['basis']]
    columns = [{**{k: c[k] for k in ('degree', 'm', 'part')}, 'mass': str(c['mass']),
                'entries': [str(x) for x in c['entries']]} for c in kernel.columns]
    return {'format': FORMAT, 'geometry': 'unit_S2', 'scope': SCOPE,
            'potential': encode(kernel.q), 'potential_degree_cap': 4,
            'mean_zero': True, 'eigenvalue_index_in_constrained_space': 1,
            'retained_degree': kernel.L, 'retained_basis': labels,
            'mass': [str(x) for x in kernel.mass], 'finite_form': _matrix_json(kernel.a),
            'range_proof': kernel.range, 'tail_start': kernel.L+1,
            'tail_lower': str(kernel.beta), 'complete_tail_couplings': columns,
            'far_tail_rule': 'all l>L+degree(q_used) have zero finite coupling; all l>L bounded',
            'upper_kind': kind, 'lower_inertia': lower_counts, 'upper_inertia': upper_counts,
            'status': 'certified_enclosure', 'formal_assistant_checked': False,
            **arithmetic.enclosure_fields(lo, hi)}


def certify(q, L=3, bits=24, reduce=True):
    """V182: certify the whole mean-zero H1 space, not only a matrix eigenvalue."""
    if type(bits) is not int or not 8 <= bits <= 48:
        raise ValueError('Require bits in 8..48')
    kernel = CoupledKernel(q, L, reduce)
    if all(e == ZERO for e in kernel.used):
        value = F(2)+kernel.used.get(ZERO, F(0))
        return _record(kernel, value, value, 'upper_ritz')
    low, high = F(2)+kernel.qlo-1, F(2)+kernel.qhi+1
    kind = 'upper_tail' if high < kernel.beta else 'upper_ritz'
    for _ in range(bits):
        mid = (low+high)/2
        if kernel.upper_holds(mid, kind):
            high = mid
        else:
            low = mid
    upper = high
    low, high = F(2)+kernel.qlo-1, upper
    distance = F(1)
    while not kernel.lower_holds(low):
        distance *= 2
        low = F(2)+kernel.qlo-distance
    for _ in range(bits):
        mid = (low+high)/2
        if kernel.lower_holds(mid):
            low = mid
        else:
            high = mid
    certificate = _record(kernel, low, upper, kind)
    if not verify(certificate, expected_q=q, independent=True):
        raise ArithmeticError('Independent certificate replay failed')
    return certificate


def verify(certificate, expected_q=None, independent=True):
    try:
        if not isinstance(certificate, dict):
            return False
        if expected_q is not None and polynomial(expected_q) != polynomial(certificate['potential']):
            return False
        kernel = CoupledKernel(certificate['potential'], certificate['retained_degree'],
                               certificate['range_proof']['reduction_enabled'])
        expected = _record(kernel, exact(certificate['lower']), exact(certificate['upper']),
                           certificate['upper_kind'], independent)
        return same(certificate, expected)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, OverflowError):
        return False


def rayleigh(q, u):
    """Exact mean-zero polynomial trial upper bound, independently integrated."""
    q, u = polynomial(q), polynomial(u, maximum_degree=8)
    mass = inner(u, u)
    if not mass > 0 or integrate(u) != 0:
        raise ValueError('Trial function must be nonzero and mean zero on the sphere')
    numerator = tangent_energy(u, u)+integrate(multiply(q, multiply(u, u)))
    return {'potential': encode(q), 'trial': encode(u), 'mass': str(mass),
            'mean': '0', 'numerator': str(numerator), 'upper': str(numerator/mass)}


def cheap_ground_enclosure(q):
    """V175: immediate whole-space lower bound and nine exact linear trials.

    This intentionally uses the raw coefficient range. Later harmonic Schur and
    sphere-identity capabilities can sharpen this inexpensive initial result.
    """
    q = polynomial(q)
    qlo, qhi = _coefficient_range(q)
    coordinates = [{tuple(int(j == i) for j in range(3)): F(1)} for i in range(3)]
    directions = list(coordinates)
    for i in range(3):
        for j in range(i+1, 3):
            for sign in (-1, 1):
                directions.append(add(coordinates[i], scale(coordinates[j], sign)))
    trials = [rayleigh(q, u) for u in directions]
    winner = min(range(len(trials)), key=lambda j: F(trials[j]['upper']))
    lower, upper = F(2)+qlo, F(trials[winner]['upper'])
    return {'format': CHEAP_FORMAT, 'geometry': 'unit_S2', 'scope': SCOPE,
            'potential': encode(q), 'potential_degree_cap': 4, 'mean_zero': True,
            'eigenvalue_index_in_constrained_space': 1,
            'potential_lower': str(qlo), 'potential_upper': str(qhi),
            'lower_rule': 'mean-zero spherical Poincare 2 plus global coefficient potential lower bound',
            'trials': trials, 'selected_trial': winner,
            'formal_assistant_checked': False, **arithmetic.enclosure_fields(lower, upper)}


def verify_cheap(certificate, expected_q=None):
    try:
        if expected_q is not None and polynomial(expected_q) != polynomial(certificate['potential']):
            return False
        expected = cheap_ground_enclosure(certificate['potential'])
        return same(certificate, expected)
    except (ValueError, KeyError, TypeError, ArithmeticError, IndexError):
        return False


def _rotation(raw):
    if not isinstance(raw, list) or len(raw) != 3 or any(
            not isinstance(row, list) or len(row) != 3 for row in raw):
        raise ValueError('Require a 3 by 3 rational proper rotation')
    r = [[exact(v) for v in row] for row in raw]
    for i in range(3):
        for j in range(3):
            if sum((r[i][k]*r[j][k] for k in range(3)), F(0)) != int(i == j):
                raise ValueError('Rotation must be exactly orthogonal')
    determinant = sum((r[0][j]*(r[1][(j+1)%3]*r[2][(j+2)%3]
                      -r[1][(j+2)%3]*r[2][(j+1)%3]) for j in range(3)), F(0))
    if determinant != 1:
        raise ValueError('Rotation determinant must equal +1')
    return r


def rotated_potential(coefficients, rotation):
    if not isinstance(coefficients, list) or not 1 <= len(coefficients) <= 5:
        raise ValueError('Require a degree <=4 univariate potential')
    coefficients = [exact(v) for v in coefficients]
    r = _rotation(rotation)
    coordinate = {(1, 0, 0): r[2][0], (0, 1, 0): r[2][1], (0, 0, 1): r[2][2]}
    return add(*(scale(power(coordinate, j), c) for j, c in enumerate(coefficients)))


def transfer_rotation(q, coefficients, rotation, child=None):
    """V184: exact isometry transfers the validated V94 full-sphere result."""
    q = polynomial(q)
    rotated = rotated_potential(coefficients, rotation)
    if q != rotated:
        raise ValueError('Requested xyz potential does not equal the rotated axisymmetric one')
    r = _rotation(rotation)
    coefficients = [exact(v) for v in coefficients]
    child = projected.full_ground([str(c) for c in coefficients], mean_zero=True,
               modes=6, max_modes=6, bits=32, max_m=4, tolerance=F(1, 10**7)) if child is None else child
    if not projected.verify_full(child, expected_q=[str(c) for c in coefficients], expected_mean_zero=True):
        raise ValueError('Invalid source full-sphere certificate')
    return {'format': ROTATION_FORMAT, 'geometry': 'unit_S2', 'scope': SCOPE,
            'potential': encode(q), 'univariate_coefficients': [str(c) for c in coefficients],
            'rotation': _matrix_json(r), 'rotated_coordinate': [str(c) for c in r[2]],
            'polynomial_identity': 'q(x)=p((R*x)_3)', 'source_certificate': child,
            'analytic_rule': 'proper orthogonal isometry preserves area, Dirichlet form and mean zero',
            'formal_assistant_checked': False,
            **arithmetic.enclosure_fields(F(child['lower']), F(child['upper']))}


def verify_rotation(certificate, expected_q=None):
    try:
        if expected_q is not None and polynomial(expected_q) != polynomial(certificate['potential']):
            return False
        expected = transfer_rotation(certificate['potential'], certificate['univariate_coefficients'],
                                     certificate['rotation'], certificate['source_certificate'])
        return same(certificate, expected)
    except (ValueError, KeyError, TypeError, ArithmeticError, IndexError):
        return False
