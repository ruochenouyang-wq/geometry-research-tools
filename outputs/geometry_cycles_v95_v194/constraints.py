"""V145--V154: exact sphere spectra with harmonic and moment constraints.

Axisymmetric polynomial potentials of degree <= 6. Complete low harmonic
exclusions support full-sphere statements; arbitrary finite moments support
ONE azimuthal sector only. All constraints have support in the retained head.
This is a rational software certificate, not a formal proof-assistant theorem.
"""
from fractions import Fraction as F
import json
from common import exact, same, sectors, projected, base, arithmetic

SECTOR_FORMAT = 'harmonic_exclusion_sector_v148'
FULL_FORMAT = 'harmonic_exclusion_full_ground_v149'
MOMENT_FORMAT = 'finite_moment_sector_v153'
INEQUALITY_FORMAT = 'constrained_weighted_inequality_v154'
ANALYTIC_EXCLUSION_FORMAT = 'analytic_harmonic_exclusion_bound_v145'
ANALYTIC_MOMENT_FORMAT = 'analytic_finite_moment_bound_v150'


def integer(x, name, low, high):
    if type(x) is not int or not low <= x <= high:
        raise ValueError('%s must be an integer in %d..%d' % (name, low, high))
    return x


def canonical_exclusion(cutoff):
    """Remove ALL spherical harmonics with degrees 0..cutoff."""
    integer(cutoff, 'cutoff', -1, 8)
    return {'kind': 'all_harmonics_degrees_leq', 'cutoff': cutoff,
            'angular_coverage': 'every_azimuth_and_both_real_components'}


def analytic_exclusion_bound(q, cutoff=1):
    """V145: full-space Poincare lower bound and an actual harmonic Ritz upper."""
    q = sectors.potential(q)
    spec = canonical_exclusion(cutoff)
    proof = projected.v05_range.make_range(q)
    qlo, _ = projected.v05_range.verify_bounds(q, proof)
    degree = cutoff+1
    form = sectors.assemble(q, 0, degree, 1)
    mass, energy = form['M'][0], form['A'][0][0]
    lower, upper = degree*(degree+1)+qlo, energy/mass
    return {'format': ANALYTIC_EXCLUSION_FORMAT, 'geometry': 'unit_S2',
            'scope': 'all_H1_functions_orthogonal_to_every_harmonic_degree_leq_cutoff',
            'constraint_spec': spec, 'q_coefficients': [str(x) for x in q],
            'range_proof': proof, 'eigenvalue_index_in_constrained_space': 1,
            'laplace_floor': degree*(degree+1),
            'ritz_witness': {'azimuth_m': 0, 'harmonic_degree': degree,
                             'radial_mass': str(mass), 'energy': str(energy)},
            'entire_radial_and_angular_space_controlled': True,
            'formal_assistant_checked': False, **arithmetic.enclosure_fields(lower, upper)}


def verify_analytic_exclusion(cert, expected_q=None, expected_cutoff=None):
    """Replay uses a multiplication recurrence for the Ritz quotient, not assembly."""
    try:
        q = sectors.potential(cert['q_coefficients'])
        cutoff = cert['constraint_spec']['cutoff']
        spec = canonical_exclusion(cutoff)
        if expected_q is not None and sectors.potential(expected_q) != q:
            return False
        if expected_cutoff is not None and not same(expected_cutoff, cutoff):
            return False
        qlo, _ = projected.v05_range.verify_bounds(q, cert['range_proof'])
        degree = cutoff+1
        mass = sectors.basis_mass(0, degree)
        upper = degree*(degree+1)+sectors.q_times_basis(q, 0, degree).get(degree, F(0))
        lower = degree*(degree+1)+qlo
        expected = {'format': ANALYTIC_EXCLUSION_FORMAT, 'geometry': 'unit_S2',
            'scope': 'all_H1_functions_orthogonal_to_every_harmonic_degree_leq_cutoff',
            'constraint_spec': spec, 'q_coefficients': [str(x) for x in q],
            'range_proof': cert['range_proof'], 'eigenvalue_index_in_constrained_space': 1,
            'laplace_floor': degree*(degree+1),
            'ritz_witness': {'azimuth_m': 0, 'harmonic_degree': degree,
                             'radial_mass': str(mass), 'energy': str(mass*upper)},
            'entire_radial_and_angular_space_controlled': True,
            'formal_assistant_checked': False, **arithmetic.enclosure_fields(lower, upper)}
        return same(cert, expected)
    except (KeyError, ValueError, TypeError, IndexError, ArithmeticError):
        return False


def sector_assembly(q, m=0, cutoff=1, modes=8):
    """V146: restriction after complete multiplication, start=max(m,c+1)."""
    canonical_exclusion(cutoff)
    integer(m, 'm', 0, 16)
    integer(modes, 'modes', 1, 16)
    return sectors.assemble(sectors.potential(q), m, max(m, cutoff+1), modes)


def tail_floors(q, cutoff, m, modes, next_m):
    """V147: radial and uncomputed angular floors in the constrained space."""
    canonical_exclusion(cutoff)
    integer(m, 'm', 0, 16)
    integer(modes, 'modes', 1, 16)
    integer(next_m, 'next_m', 0, 17)
    q = sectors.potential(q)
    proof = projected.v05_range.make_range(q)
    qlo, qhi = projected.v05_range.verify_bounds(q, proof)
    radial_degree = max(m, cutoff+1)+modes
    angular_degree = max(next_m, cutoff+1)
    return {'range_proof': proof, 'q_lower': str(qlo), 'q_upper': str(qhi),
            'radial_degree': radial_degree,
            'radial_floor': str(radial_degree*(radial_degree+1)+qlo),
            'angular_minimum_degree': angular_degree,
            'angular_floor': str(angular_degree*(angular_degree+1)+qlo)}


def canonical_moments(rows, m=0, modes=8):
    """g_j = sum rows[j][i] P_(m+i)^m; require <u,g_j>=0.

    Rows are harmonic expansion coefficients, NOT already weighted functionals.
    Reject support outside the retained head rather than truncate it silently.
    """
    integer(m, 'm', 0, 16)
    integer(modes, 'modes', 1, 16)
    if not isinstance(rows, list) or not 1 <= len(rows) <= 4:
        raise ValueError('Provide 1..4 finite harmonic coefficient rows')
    result = []
    for row in rows:
        if not isinstance(row, list) or not 1 <= len(row) <= modes:
            raise ValueError('Moment support must lie entirely in the retained head')
        values = [exact(x) for x in row]
        if any(abs(x) > 10**6 or x.denominator > 10**9 for x in values):
            raise ValueError('Moment coefficient limit exceeded')
        while len(values) > 1 and not values[-1]:
            values.pop()
        result.append([str(x) for x in values])
    return {'kind': 'single_azimuth_finite_moments', 'azimuth_m': m,
            'coefficient_degree_start': m, 'rows': result,
            'real_angular_component': 'cos(m*phi); constant 1 when m=0',
            'interpretation': 'L2_inner_product_with_harmonic_expansion'}


def analytic_moment_bound(q, rows, m=0, modes=8):
    """V150: full radial variational bracket for finitely many exact moments.

    At most r constraints leave a nonzero vector in the first r+1 harmonic
    modes; min-max gives lambda_1 <= (m+r)(m+r+1)+q_upper. Redundancy may
    make this coarse, but does not invalidate it. No finite tail is discarded.
    """
    q = sectors.potential(q)
    spec = canonical_moments(rows, m, modes)
    proof = projected.v05_range.make_range(q)
    qlo, qhi = projected.v05_range.verify_bounds(q, proof)
    r = len(rows)
    degree = m+r
    return {'format': ANALYTIC_MOMENT_FORMAT, 'geometry': 'unit_S2',
            'scope': 'single_real_cosine_azimuth_component_only',
            'q_coefficients': [str(x) for x in q], 'azimuth_m': m, 'modes': modes,
            'constraint_spec': spec, 'range_proof': proof,
            'eigenvalue_index_in_constrained_sector': 1,
            'variational_trial_dimension': r+1, 'constraint_count_upper_bound': r,
            'trial_highest_harmonic_degree': degree, 'entire_radial_tail_controlled': True,
            'formal_assistant_checked': False,
            **arithmetic.enclosure_fields(m*(m+1)+qlo, degree*(degree+1)+qhi)}


def verify_analytic_moment(cert, expected_q=None, expected_rows=None, expected_m=None):
    try:
        q = sectors.potential(cert['q_coefficients'])
        m, modes = cert['azimuth_m'], cert['modes']
        spec = canonical_moments(cert['constraint_spec']['rows'], m, modes)
        if expected_q is not None and sectors.potential(expected_q) != q:
            return False
        if expected_rows is not None and not same(canonical_moments(expected_rows, m, modes), spec):
            return False
        if expected_m is not None and not same(expected_m, m):
            return False
        qlo, qhi = projected.v05_range.verify_bounds(q, cert['range_proof'])
        r, degree = len(spec['rows']), m+len(spec['rows'])
        expected = {'format': ANALYTIC_MOMENT_FORMAT, 'geometry': 'unit_S2',
            'scope': 'single_real_cosine_azimuth_component_only',
            'q_coefficients': [str(x) for x in q], 'azimuth_m': m, 'modes': modes,
            'constraint_spec': spec, 'range_proof': cert['range_proof'],
            'eigenvalue_index_in_constrained_sector': 1,
            'variational_trial_dimension': r+1, 'constraint_count_upper_bound': r,
            'trial_highest_harmonic_degree': degree, 'entire_radial_tail_controlled': True,
            'formal_assistant_checked': False,
            **arithmetic.enclosure_fields(m*(m+1)+qlo, degree*(degree+1)+qhi)}
        return same(cert, expected)
    except (KeyError, ValueError, TypeError, IndexError, ArithmeticError):
        return False


def _wire_matrix(a):
    return [[str(x) for x in row] for row in a]


def mass_nullspace(rows, masses):
    """V151: exact RREF of weighted rows and an exact full-rank nullspace."""
    n = len(masses)
    if not 1 <= n <= 16:
        raise ValueError('Mass dimension outside 1..16')
    masses = [exact(x) for x in masses]
    if any(x <= 0 for x in masses):
        raise ValueError('Strictly positive harmonic masses required')
    if not isinstance(rows, list) or len(rows) > 4:
        raise ValueError('At most four rows are supported')
    weighted = []
    for row in rows:
        if not isinstance(row, list) or not 1 <= len(row) <= n:
            raise ValueError('Row dimension incompatible with retained head')
        values = [exact(x) for x in row]+[F(0)]*(n-len(row))
        weighted.append([a*b for a, b in zip(values, masses)])
    a = [row[:] for row in weighted]
    pivots, rank = [], 0
    for column in range(n):
        pivot = next((j for j in range(rank, len(a)) if a[j][column]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][column]
        a[rank] = [x/scale for x in a[rank]]
        for j in range(len(a)):
            if j != rank and a[j][column]:
                scale = a[j][column]
                a[j] = [x-scale*y for x, y in zip(a[j], a[rank])]
        pivots.append(column)
        rank += 1
    if rank >= n:
        raise ValueError('Constraints eliminate the complete retained head; increase modes')
    free = [i for i in range(n) if i not in pivots]
    t = [[F(0) for _ in free] for _ in range(n)]
    for j, free_column in enumerate(free):
        t[free_column][j] = F(1)
        for row, pivot_column in enumerate(pivots):
            t[pivot_column][j] = -a[row][free_column]
    return {'rank': rank, 'pivots': pivots, 'free_columns': free,
            'weighted_rows': _wire_matrix(weighted),
            'rref_nonzero_rows': _wire_matrix(a[:rank]), 'T': _wire_matrix(t)}


def _congruence(a, t):
    n, dimension = len(t), len(t[0])
    at = [[sum((a[i][j]*t[j][k] for j in range(n)), F(0))
           for k in range(dimension)] for i in range(n)]
    return [[sum((t[i][j]*at[i][k] for i in range(n)), F(0))
             for k in range(dimension)] for j in range(dimension)]


def constrained_forms(q, m, modes, rows=None, cutoff=-1):
    """V152: T^T A T, full T^T M T, and transformed head-tail coupling.

    The mass is generally NON-diagonal. A diagonal-only subtraction would
    certify the wrong Rayleigh quotient after this change of coordinates.
    """
    q = sectors.potential(q)
    data = sector_assembly(q, m, cutoff, modes)
    if rows is None:
        spec = canonical_exclusion(cutoff)
        rows = []
    else:
        if cutoff != -1:
            raise ValueError('General moment mode has a fixed coefficient start m')
        spec = canonical_moments(rows, m, modes)
        rows = spec['rows']
    null = mass_nullspace(rows, data['M'])
    t = [[F(x) for x in row] for row in null['T']]
    mass_matrix = [[data['M'][i] if i == j else F(0) for j in range(modes)]
                   for i in range(modes)]
    a, mass = _congruence(data['A'], t), _congruence(mass_matrix, t)
    if base.inertia(mass) != [0, 0, len(mass)]:
        raise ArithmeticError('Reduced mass is not strictly positive')
    tail_start = data['degrees'][-1]+1
    columns = []
    for degree in range(tail_start, tail_start+len(q)-1):
        expansion = sectors.q_times_basis(q, m, degree)
        old_column = [data['M'][i]*expansion.get(l, F(0))
                      for i, l in enumerate(data['degrees'])]
        column = [sum((t[i][j]*old_column[i] for i in range(modes)), F(0))
                  for j in range(len(t[0]))]
        columns.append({'degree': degree, 'mass': str(sectors.basis_mass(m, degree)),
                        'column': [str(x) for x in column]})
    return {'constraint_spec': spec, 'degrees': data['degrees'],
            'head_mass': [str(x) for x in data['M']], 'nullspace': null,
            'reduced_A': _wire_matrix(a), 'reduced_mass': _wire_matrix(mass),
            'tail_start': tail_start, 'couplings': columns,
            'reduced_dimension': len(mass), 'mass_is_diagonal': all(
                not mass[i][j] for i in range(len(mass)) for j in range(i))}


class ConstraintKernel:
    def __init__(self, q, m=0, cutoff=1, modes=8, rows=None, range_proof=None):
        self.q = sectors.potential(q)
        self.m, self.cutoff, self.modes = m, cutoff, modes
        self.forms = constrained_forms(self.q, m, modes, rows, cutoff)
        self.spec = self.forms['constraint_spec']
        self.moment_mode = rows is not None
        self.start = self.forms['degrees'][0]
        self.tail = self.forms['tail_start']
        self.dimension = self.forms['reduced_dimension']
        self.a = [[F(x) for x in row] for row in self.forms['reduced_A']]
        self.mass = [[F(x) for x in row] for row in self.forms['reduced_mass']]
        self.range_proof = (projected.v05_range.make_range(self.q)
                            if range_proof is None else range_proof)
        self.qlo, self.qhi = projected.v05_range.verify_bounds(self.q, self.range_proof)
        self.beta = self.tail*(self.tail+1)+self.qlo

    def matrix(self, x, kind):
        x = exact(x)
        if kind not in ('lower', 'upper_tail', 'upper_ritz'):
            raise ValueError('Unknown comparison')
        if kind != 'upper_ritz' and x >= self.beta:
            raise ValueError('Schur comparison requires a strictly positive radial tail')
        matrix = [[self.a[i][j]-x*self.mass[i][j] for j in range(self.dimension)]
                  for i in range(self.dimension)]
        if kind == 'upper_ritz':
            return matrix
        for entry in self.forms['couplings']:
            l, mass = entry['degree'], F(entry['mass'])
            diagonal = l*(l+1)+(self.qlo if kind == 'lower' else self.qhi)
            denominator = mass*(diagonal-x)
            column = [F(v) for v in entry['column']]
            for i, a in enumerate(column):
                for j, b in enumerate(column):
                    matrix[i][j] -= a*b/denominator
        return matrix

    def count(self, x, kind, independent=False):
        matrix = self.matrix(x, kind)
        return base.inertia(matrix) if independent else arithmetic.banded_inertia(matrix)

    def lower_holds(self, k, x):
        return x < self.beta and self.count(x, 'lower')[0] < k

    def upper_holds(self, k, x, kind):
        if kind != 'upper_ritz' and x >= self.beta:
            return False
        negative, zero, _ = self.count(x, kind)
        return negative+zero >= k


def _sector_evidence(kernel, k, lower, upper, kind, independent=False):
    integer(k, 'k', 1, kernel.dimension)
    lower, upper = exact(lower), exact(upper)
    if lower > upper or kind not in ('upper_tail', 'upper_ritz'):
        raise ValueError('Reversed endpoints or invalid upper comparison')
    low_count = kernel.count(lower, 'lower', independent)
    up_count = kernel.count(upper, kind, independent)
    if low_count[0] >= k or up_count[0]+up_count[1] < k:
        raise ValueError('Exact inertia does not establish claimed endpoints')
    return {'format': MOMENT_FORMAT if kernel.moment_mode else SECTOR_FORMAT,
            'geometry': 'unit_S2', 'scope': ('single_real_cosine_azimuth_component_only'
                if kernel.moment_mode else 'single_azimuth_sector_only'),
            'q_coefficients': [str(x) for x in kernel.q], 'azimuth_m': kernel.m,
            'constraint_spec': kernel.spec, 'modes': kernel.modes,
            'eigenvalue_index_in_constrained_sector': k,
            'index_counts_retained_nullspace_and_entire_radial_tail': True,
            'forms': kernel.forms, 'range_proof': kernel.range_proof,
            'radial_tail_lower': str(kernel.beta), 'upper_kind': kind,
            'lower_inertia': low_count, 'upper_inertia': up_count,
            'formal_assistant_checked': False,
            **arithmetic.enclosure_fields(lower, upper)}


def verify_sector(cert, expected_q=None, expected_m=None, expected_cutoff=None,
                  expected_rows=None, expected_k=None, independent=True):
    """V148/V153 replay reconstructs rows, mass, nullspace, forms and tails."""
    try:
        if not isinstance(cert, dict):
            return False
        spec = cert['constraint_spec']
        m, modes = cert['azimuth_m'], cert['modes']
        moment_mode = cert['format'] == MOMENT_FORMAT
        if not moment_mode and cert['format'] != SECTOR_FORMAT:
            return False
        rows = spec['rows'] if moment_mode else None
        cutoff = -1 if moment_mode else spec['cutoff']
        kernel = ConstraintKernel(cert['q_coefficients'], m, cutoff, modes, rows,
                                  cert['range_proof'])
        if expected_q is not None and sectors.potential(expected_q) != kernel.q:
            return False
        if expected_m is not None and not same(expected_m, m):
            return False
        if expected_cutoff is not None and (moment_mode or not same(expected_cutoff, cutoff)):
            return False
        if expected_rows is not None and (not moment_mode or not same(
                canonical_moments(expected_rows, m, modes), spec)):
            return False
        k = cert['eigenvalue_index_in_constrained_sector']
        if expected_k is not None and not same(expected_k, k):
            return False
        return same(cert, _sector_evidence(kernel, k, exact(cert['lower']),
                    exact(cert['upper']), cert['upper_kind'], independent))
    except (KeyError, ValueError, TypeError, IndexError, ArithmeticError):
        return False


def _certify(kernel, k, bits):
    integer(k, 'k', 1, kernel.dimension)
    integer(bits, 'bits', 8, 80)
    if len(kernel.q) == 1 and not kernel.moment_mode:
        degree = kernel.start+k-1
        value = degree*(degree+1)+kernel.q[0]
        cert = _sector_evidence(kernel, k, value, value, 'upper_ritz')
    else:
        floor = kernel.start*(kernel.start+1)+kernel.qlo
        low = floor-1
        high = (kernel.tail-1)*kernel.tail+kernel.qhi+1
        kind = 'upper_tail' if high < kernel.beta else 'upper_ritz'
        if not kernel.upper_holds(k, high, kind):
            raise ArithmeticError('Finite Ritz upper bracket failed')
        for _ in range(bits):
            mid = (low+high)/2
            if kernel.upper_holds(k, mid, kind):
                high = mid
            else:
                low = mid
        if kind == 'upper_ritz' and high < kernel.beta:
            kind, low = 'upper_tail', floor-1
            for _ in range(bits):
                mid = (low+high)/2
                if kernel.upper_holds(k, mid, kind):
                    high = mid
                else:
                    low = mid
        upper = high
        distance = F(1)
        low = floor-distance
        while not kernel.lower_holds(k, low):
            distance *= 2
            low = floor-distance
        high = min(upper, kernel.beta)
        for _ in range(bits):
            mid = (low+high)/2
            if kernel.lower_holds(k, mid):
                low = mid
            else:
                high = mid
        cert = _sector_evidence(kernel, k, low, upper, kind)
    if not verify_sector(cert):
        raise ArithmeticError('Independent exact inertia replay failed')
    return cert


def certify_exclusion_sector(q, m=0, cutoff=1, k=1, modes=8, bits=40, range_proof=None):
    """V148: infinite-dimensional fixed-m spectrum after harmonic exclusion."""
    return _certify(ConstraintKernel(q, m, cutoff, modes, range_proof=range_proof), k, bits)


def certify_moment_sector(q, rows, m=0, k=1, modes=8, bits=40):
    """V153: arbitrary finite rational moment constraints in ONE m sector."""
    return _certify(ConstraintKernel(q, m, -1, modes, rows=rows), k, bits)


def _full_evidence(q, cutoff, certificates, range_proof, tolerance, independent=True):
    q = sectors.potential(q)
    spec = canonical_exclusion(cutoff)
    tolerance = exact(tolerance)
    if not F(1, 10**24) <= tolerance <= 1:
        raise ValueError('Tolerance must lie in [1e-24,1]')
    if not isinstance(certificates, list) or not 1 <= len(certificates) <= 17:
        raise ValueError('Consecutive azimuth certificates required')
    qlo, _ = projected.v05_range.verify_bounds(q, range_proof)
    for m, cert in enumerate(certificates):
        if not verify_sector(cert, q, m, cutoff, expected_k=1, independent=independent):
            raise ValueError('Sector mismatch: require the same full harmonic exclusion')
    next_m = len(certificates)
    degree = max(next_m, cutoff+1)
    angular_floor = degree*(degree+1)+qlo
    lower = min(angular_floor, min(F(c['lower']) for c in certificates))
    upper = min(F(c['upper']) for c in certificates)
    return {'format': FULL_FORMAT, 'geometry': 'unit_S2',
            'scope': 'all_H1_functions_orthogonal_to_every_harmonic_degree_leq_cutoff',
            'q_coefficients': [str(x) for x in q], 'constraint_spec': spec,
            'eigenvalue_index_in_constrained_space': 1,
            'sectors': certificates, 'range_proof': range_proof,
            'omitted_azimuth_m_start': next_m, 'angular_minimum_degree': degree,
            'angular_tail_lower': str(angular_floor),
            'angular_tail_cannot_improve_best_upper': angular_floor >= upper,
            'tolerance': str(tolerance),
            'status': 'certified_target_met' if upper-lower <= tolerance else 'certified_bound_open_gap',
            'formal_assistant_checked': False, **arithmetic.enclosure_fields(lower, upper)}


def full_exclusion_ground(q, cutoff=1, modes=8, bits=40, max_m=12, tolerance=F(1, 10**9)):
    """V149: complete sphere, including both infinite radial and angular tails."""
    canonical_exclusion(cutoff)
    integer(max_m, 'max_m', 0, 16)
    q = sectors.potential(q)
    proof = projected.v05_range.make_range(q)
    qlo, _ = projected.v05_range.verify_bounds(q, proof)
    certificates = []
    for m in range(max_m+1):
        certificates.append(certify_exclusion_sector(q, m, cutoff, 1, modes, bits, proof))
        degree = max(m+1, cutoff+1)
        if degree*(degree+1)+qlo >= min(F(c['upper']) for c in certificates):
            break
    return _full_evidence(q, cutoff, certificates, proof, tolerance)


def verify_full(cert, expected_q=None, expected_cutoff=None, independent=True):
    try:
        if not isinstance(cert, dict) or cert.get('format') != FULL_FORMAT:
            return False
        q, cutoff = cert['q_coefficients'], cert['constraint_spec']['cutoff']
        if expected_q is not None and sectors.potential(expected_q) != sectors.potential(q):
            return False
        if expected_cutoff is not None and not same(cutoff, expected_cutoff):
            return False
        return same(cert, _full_evidence(q, cutoff, cert['sectors'], cert['range_proof'],
                                         exact(cert['tolerance']), independent))
    except (KeyError, ValueError, TypeError, IndexError, ArithmeticError):
        return False


def inequality_certificate(spectral, threshold):
    """V154: bind quantifier to either full harmonic exclusion or one sector."""
    if not isinstance(spectral, dict):
        raise ValueError('A spectrum certificate is required')
    if spectral.get('format') == FULL_FORMAT and verify_full(spectral):
        quantifier = 'all_H1_functions_on_unit_S2_satisfying_complete_harmonic_exclusion'
    elif spectral.get('format') in (MOMENT_FORMAT, SECTOR_FORMAT) and verify_sector(spectral, expected_k=1):
        quantifier = ('all_H1_functions_in_this_single_azimuth_cosine_component_satisfying_constraints'
                      if spectral['format'] == MOMENT_FORMAT else
                      'all_H1_functions_in_this_single_azimuth_sector_satisfying_constraints')
    else:
        raise ValueError('A verified ground eigenvalue certificate is required')
    threshold = exact(threshold)
    lower, upper = F(spectral['lower']), F(spectral['upper'])
    status = ('proved' if threshold <= lower else 'refuted_by_spectral_existence'
              if threshold > upper else 'undetermined')
    return {'format': INEQUALITY_FORMAT, 'geometry': 'unit_S2', 'measure': 'standard_sphere_area',
            'statement': 'integral(grad_u_squared + q*u_squared) >= threshold*integral(u_squared)',
            'quantifier': quantifier, 'constraint_spec': spectral['constraint_spec'],
            'q_coefficients': spectral['q_coefficients'], 'threshold': str(threshold),
            'status': status, 'evidence': spectral, 'explicit_counterexample_included': False,
            'not_a_universal_Gagliardo_Nirenberg_constant_certificate': True,
            'formal_assistant_checked': False}


def constrained_inequality(q, threshold, cutoff=1, rows=None, m=0, **options):
    if rows is None:
        spectral = full_exclusion_ground(q, cutoff=cutoff, **options)
    else:
        if cutoff != -1:
            raise ValueError('For general moments explicitly use cutoff=-1; scope is single m')
        spectral = certify_moment_sector(q, rows, m=m, **options)
    return inequality_certificate(spectral, threshold)


def verify_inequality(cert, expected_q=None, expected_threshold=None):
    try:
        if expected_q is not None and sectors.potential(expected_q) != sectors.potential(cert['q_coefficients']):
            return False
        if expected_threshold is not None and exact(expected_threshold) != exact(cert['threshold']):
            return False
        return same(cert, inequality_certificate(cert['evidence'], exact(cert['threshold'])))
    except (KeyError, ValueError, TypeError, IndexError, ArithmeticError):
        return False


def verify(cert):
    if not isinstance(cert, dict):
        return False
    if cert.get('format') in (SECTOR_FORMAT, MOMENT_FORMAT):
        return verify_sector(cert)
    if cert.get('format') == FULL_FORMAT:
        return verify_full(cert)
    if cert.get('format') == INEQUALITY_FORMAT:
        return verify_inequality(cert)
    if cert.get('format') == ANALYTIC_EXCLUSION_FORMAT:
        return verify_analytic_exclusion(cert)
    if cert.get('format') == ANALYTIC_MOMENT_FORMAT:
        return verify_analytic_moment(cert)
    return False
