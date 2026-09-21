"""Certified sphere spectra with mean-zero projection and azimuthal sectors.

The potential q(t) is axisymmetric, but full_ground controls ALL functions in
the indicated full-sphere function space, including nonaxisymmetric modes.
Finite-degree and infinite-angular tails are both explicitly bounded.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'geometry_v36_v90'))
import spectral_certifier as base
import error_bounds as arithmetic
import v05_range
import sphere_sectors

SECTOR_FORMAT = 'projected_sphere_sector_v91_v92'
FULL_FORMAT = 'projected_full_sphere_ground_v93'
INEQUALITY_FORMAT = 'weighted_poincare_certificate_v94'


def canonical_equal(a, b):
    # In particular, reject bools masquerading as integers in proof metadata.
    return json.dumps(a, sort_keys=True, separators=(',', ':')) == json.dumps(b, sort_keys=True, separators=(',', ':'))


def potential(q):
    if not isinstance(q, list):
        raise ValueError('Potential must be a list of rational power coefficients')
    return base.potential([str(x) for x in q])


def parameters(m, mean_zero, k, modes):
    if type(m) is not int or not 0 <= m <= 64 or type(mean_zero) is not bool:
        raise ValueError('Require integer 0 <= m <= 64 and boolean mean_zero')
    if type(k) is not int or type(modes) is not int or not 1 <= k <= modes <= 64:
        raise ValueError('Require integer 1 <= k <= modes <= 64')


def scope(mean_zero):
    return ('all_mean_zero_H1_functions_on_unit_S2' if mean_zero
            else 'all_H1_functions_on_unit_S2')


class Kernel:
    def __init__(self, q, m=0, mean_zero=True, modes=12, range_proof=None):
        parameters(m, mean_zero, 1, modes)
        self.q = potential(q)
        self.m, self.mean_zero, self.n = m, mean_zero, modes
        self.start = 1 if m == 0 and mean_zero else m
        self.tail_start = self.start+modes
        self.range_proof = v05_range.make_range(self.q) if range_proof is None else range_proof
        self.qlo, self.qhi = v05_range.verify_bounds(self.q, self.range_proof)
        self.beta = self.tail_start*(self.tail_start+1)+self.qlo
        data = sphere_sectors.assemble(self.q, m, self.start, modes)
        self.a, self.mass = data['A'], data['M']
        self.couplings = []
        # A degree-d multiplier couples at most d omitted radial degrees.
        # Projection removes l=0 only AFTER multiplication by the whole q.
        for l in range(self.tail_start, self.tail_start+len(self.q)-1):
            column = [(i-self.start, self.mass[i-self.start]*v)
                      for i, v in sphere_sectors.q_times_basis(self.q, m, l).items()
                      if self.start <= i < self.tail_start and v]
            self.couplings.append((l, sphere_sectors.basis_mass(m, l), column))

    def matrix(self, x, kind):
        x = F(x)
        if kind not in ('lower', 'upper_tail', 'upper_ritz'):
            raise ValueError('Unknown Schur comparison kind')
        if kind != 'upper_ritz' and x >= self.beta:
            raise ValueError('Infinite radial tail must be strictly positive')
        a = [row[:] for row in self.a]
        for i in range(self.n):
            a[i][i] -= x*self.mass[i]
        if kind == 'upper_ritz':
            return a
        for l, mass, column in self.couplings:
            diagonal = l*(l+1)+(self.qlo if kind == 'lower' else self.qhi)
            denominator = mass*(diagonal-x)
            for i, c in column:
                for j, d in column:
                    a[i][j] -= c*d/denominator
        return a

    def count(self, x, kind, independent=False):
        matrix = self.matrix(x, kind)
        return base.inertia(matrix) if independent else arithmetic.banded_inertia(matrix)

    def lower_holds(self, k, x):
        return x < self.beta and self.count(x, 'lower')[0] < k

    def upper_holds(self, k, x, kind):
        if kind == 'upper_tail' and x >= self.beta:
            return False
        negative, zero, _ = self.count(x, kind)
        return negative+zero >= k


def sector_certificate(kernel, k, lower, upper, upper_kind, independent=False):
    parameters(kernel.m, kernel.mean_zero, k, kernel.n)
    lo, hi = F(lower), F(upper)
    if lo > hi or upper_kind not in ('upper_tail', 'upper_ritz'):
        raise ValueError('Invalid spectral endpoints or upper comparison')
    low_count = kernel.count(lo, 'lower', independent)
    up_count = kernel.count(hi, upper_kind, independent)
    if low_count[0] >= k or up_count[0]+up_count[1] < k:
        raise ValueError('Exact spectral counts do not establish the endpoints')
    return {'format': SECTOR_FORMAT, 'geometry': 'unit_S2',
            'q_coefficients': [str(x) for x in kernel.q], 'azimuth_m': kernel.m,
            'mean_zero': kernel.mean_zero,
            'projection': 'remove_l0' if kernel.m == 0 and kernel.mean_zero else 'none',
            'scope': 'single_azimuth_sector_only',
            'degree_start': kernel.start, 'modes': kernel.n, 'eigenvalue_index': k,
            'range_proof': kernel.range_proof,
            'radial_tail_start': kernel.tail_start, 'radial_tail_lower': str(kernel.beta),
            'upper_kind': upper_kind, 'lower_inertia': low_count, 'upper_inertia': up_count,
            'formal_assistant_checked': False, **arithmetic.enclosure_fields(lo, hi)}


def verify_sector(cert, expected_q=None, expected_m=None, expected_mean_zero=None,
                  expected_k=None, independent=True):
    try:
        if not isinstance(cert, dict):
            return False
        if expected_q is not None and potential(expected_q) != potential(cert['q_coefficients']):
            return False
        if expected_m is not None and (type(expected_m) is not int or cert['azimuth_m'] != expected_m):
            return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool
                                               or cert['mean_zero'] is not expected_mean_zero):
            return False
        if expected_k is not None and (type(expected_k) is not int or cert['eigenvalue_index'] != expected_k):
            return False
        kernel = Kernel(cert['q_coefficients'], cert['azimuth_m'], cert['mean_zero'],
                        cert['modes'], cert['range_proof'])
        expected = sector_certificate(kernel, cert['eigenvalue_index'], base.rational(cert['lower']),
                    base.rational(cert['upper']), cert['upper_kind'], independent=independent)
        return canonical_equal(cert, expected)
    except (ValueError, KeyError, TypeError, ArithmeticError, IndexError):
        return False


def certify_sector(q, m=0, mean_zero=True, k=1, modes=12, bits=44, range_proof=None):
    parameters(m, mean_zero, k, modes)
    if type(bits) is not int or not 8 <= bits <= 160:
        raise ValueError('Require integer bits in 8..160')
    kernel = Kernel(q, m, mean_zero, modes, range_proof)
    if len(kernel.q) == 1:
        l = kernel.start+k-1
        exact = l*(l+1)+kernel.q[0]
        cert = sector_certificate(kernel, k, exact, exact, 'upper_ritz')
    else:
        unperturbed_low = kernel.start*(kernel.start+1)
        l = kernel.start+k-1
        low, high = unperturbed_low+kernel.qlo-1, l*(l+1)+kernel.qhi+1
        kind = 'upper_tail' if high < kernel.beta else 'upper_ritz'
        if not kernel.upper_holds(k, high, kind):
            raise ArithmeticError('Analytic Ritz upper bracket failed')
        for _ in range(bits):
            mid = (low+high)/2
            if kernel.upper_holds(k, mid, kind):
                high = mid
            else:
                low = mid
        if kind == 'upper_ritz' and high < kernel.beta:
            kind = 'upper_tail'
            low = unperturbed_low+kernel.qlo-1
            for _ in range(bits):
                mid = (low+high)/2
                if kernel.upper_holds(k, mid, kind):
                    high = mid
                else:
                    low = mid
        upper = high
        distance = F(1)
        low = unperturbed_low+kernel.qlo-distance
        # Termination follows from the mass term dominating as low -> -infty.
        while not kernel.lower_holds(k, low):
            distance *= 2
            low = unperturbed_low+kernel.qlo-distance
        high = upper
        for _ in range(bits):
            mid = (low+high)/2
            if kernel.lower_holds(k, mid):
                low = mid
            else:
                high = mid
        cert = sector_certificate(kernel, k, low, upper, kind)
    if not verify_sector(cert, independent=True):
        raise ArithmeticError('Independent dense sector verifier rejected certificate')
    return cert


def full_certificate(q, mean_zero, sectors, range_proof, tolerance, independent=True):
    q = potential(q)
    if type(mean_zero) is not bool or not isinstance(sectors, list) or not 1 <= len(sectors) <= 65:
        raise ValueError('At least one consecutive azimuth sector is required')
    tolerance = F(tolerance)
    if not F(1, 10**30) <= tolerance <= 1:
        raise ValueError('Require 1e-30 <= tolerance <= 1')
    qlo, _ = v05_range.verify_bounds(q, range_proof)
    for m, cert in enumerate(sectors):
        if not verify_sector(cert, expected_q=q, expected_m=m,
                             expected_mean_zero=mean_zero, expected_k=1, independent=independent):
            raise ValueError('Sector proof has wrong input, index, projection, or evidence')
    next_m = len(sectors)
    angular_tail = next_m*(next_m+1)+qlo
    finite_low = min(F(c['lower']) for c in sectors)
    upper = min(F(c['upper']) for c in sectors)
    lower = min(finite_low, angular_tail)
    covered = angular_tail >= upper
    return {'format': FULL_FORMAT, 'geometry': 'unit_S2', 'scope': scope(mean_zero),
            'q_coefficients': [str(x) for x in q], 'mean_zero': mean_zero,
            'potential_symmetry': 'axisymmetric_polynomial_degree_at_most_6',
            'eigenvalue_index_in_constrained_space': 1,
            'sectors': sectors, 'range_proof': range_proof,
            'omitted_azimuth_m_start': next_m, 'angular_tail_lower': str(angular_tail),
            'angular_tail_cannot_improve_best_upper': covered,
            'tolerance': str(tolerance),
            'status': 'certified_target_met' if upper-lower <= tolerance else 'certified_bound_open_gap',
            'formal_assistant_checked': False, **arithmetic.enclosure_fields(lower, upper)}


def verify_full(cert, expected_q=None, expected_mean_zero=None, independent=True):
    try:
        if not isinstance(cert, dict):
            return False
        if expected_q is not None and potential(expected_q) != potential(cert['q_coefficients']):
            return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool
                                               or cert['mean_zero'] is not expected_mean_zero):
            return False
        expected = full_certificate(cert['q_coefficients'], cert['mean_zero'], cert['sectors'],
                                     cert['range_proof'], base.rational(cert['tolerance']), independent)
        return canonical_equal(cert, expected)
    except (ValueError, KeyError, TypeError, ArithmeticError, IndexError):
        return False


def full_ground(q, mean_zero=True, modes=8, max_modes=32, bits=44, max_m=16,
                tolerance=F(1, 10**10)):
    parameters(0, mean_zero, 1, modes)
    if type(max_modes) is not int or not modes <= max_modes <= 64:
        raise ValueError('Require modes <= max_modes <= 64')
    if type(max_m) is not int or not 0 <= max_m <= 64:
        raise ValueError('Require 0 <= max_m <= 64')
    tolerance = F(tolerance)
    if not F(1, 10**30) <= tolerance <= 1:
        raise ValueError('Require 1e-30 <= tolerance <= 1')
    q = potential(q)
    proof = v05_range.make_range(q)
    qlo, _ = v05_range.verify_bounds(q, proof)
    sectors = []
    for m in range(max_m+1):
        size = modes
        while True:
            cert = certify_sector(q, m=m, mean_zero=mean_zero, modes=size, bits=bits, range_proof=proof)
            if F(cert['exact_width']) <= tolerance or size == max_modes:
                break
            size = min(max_modes, size*2)
        sectors.append(cert)
        next_m = m+1
        if next_m*(next_m+1)+qlo >= min(F(c['upper']) for c in sectors):
            break
    return full_certificate(q, mean_zero, sectors, proof, tolerance)


def inequality_certificate(spectral_certificate, threshold):
    if not verify_full(spectral_certificate, expected_mean_zero=True):
        raise ValueError('A verified FULL mean-zero spectrum certificate is required')
    threshold = base.rational(str(threshold))
    lower, upper = F(spectral_certificate['lower']), F(spectral_certificate['upper'])
    status = ('proved' if lower >= threshold else 'refuted_by_spectral_existence'
              if upper < threshold else 'undetermined')
    return {'format': INEQUALITY_FORMAT, 'geometry': 'unit_S2',
            'statement': 'integral(grad_u_squared + q*u_squared) >= threshold*integral(u_squared)',
            'quantifier': 'every_real_H1_function_u_with_integral_u_zero',
            'measure': 'standard_sphere_area', 'q_coefficients': spectral_certificate['q_coefficients'],
            'threshold': str(threshold), 'status': status,
            'evidence': spectral_certificate,
            'refutation_does_not_include_an_explicit_function': status == 'refuted_by_spectral_existence',
            'not_a_universal_Gagliardo_Nirenberg_constant_certificate': True,
            'formal_assistant_checked': False}


def verify_inequality(cert, expected_q=None, expected_threshold=None):
    try:
        if expected_q is not None and potential(expected_q) != potential(cert['q_coefficients']):
            return False
        if expected_threshold is not None and base.rational(str(expected_threshold)) != base.rational(cert['threshold']):
            return False
        return canonical_equal(cert, inequality_certificate(cert['evidence'], cert['threshold']))
    except (ValueError, KeyError, TypeError, ArithmeticError, IndexError):
        return False


def weighted_poincare(q, threshold, **kwargs):
    if 'mean_zero' in kwargs:
        raise ValueError('Weighted Poincare operation always requires mean_zero=True')
    return inequality_certificate(full_ground(q, mean_zero=True, **kwargs), threshold)
