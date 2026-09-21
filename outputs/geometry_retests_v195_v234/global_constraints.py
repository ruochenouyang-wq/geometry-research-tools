"""V195--204: complete real sphere spectra with finite cross-harmonic moments.

All integrals use area/(4*pi). The head includes every real spherical harmonic
of degrees 0..L, including the constant and every sine/cosine component.
Constraints are complete xyz polynomials, checked for support in that head.
Certificates are rational software proofs, not formal-assistant theorems.
"""
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
from support import exact, same, arithmetic, base, old_anisotropic as a

SCOPE = 'all_real_H1_functions_on_unit_S2_satisfying_all_recorded_moments'
TRIAL = 'global_moment_trial_enclosure_v195'
SPECTRAL = 'global_moment_schur_v200'
TRANSFER = 'equivalent_global_moments_transfer_v202'
ADAPTIVE = 'global_moment_adaptive_v203'
INEQUALITY = 'global_moment_inequality_v204'


def _integer(value, low, high, name):
    if type(value) is not int or not low <= value <= high:
        raise ValueError('%s must be an integer in %d..%d' % (name, low, high))
    return value


def moments(raw):
    if not isinstance(raw, list) or len(raw) > 12:
        raise ValueError('Require a list of at most 12 complete polynomial moment functions')
    return [a.polynomial(row, maximum_degree=4) for row in raw]


def _wire_matrix(matrix):
    return [[str(x) for x in row] for row in matrix]


def _rref(rows, dimension):
    matrix = [[exact(x) for x in row] for row in rows]
    if any(len(row) != dimension for row in matrix):
        raise ValueError('RREF row dimension mismatch')
    transform = [[F(i == j) for j in range(len(rows))] for i in range(len(rows))]
    rank, pivots = 0, []
    for column in range(dimension):
        pivot = next((j for j in range(rank, len(matrix)) if matrix[j][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        transform[rank], transform[pivot] = transform[pivot], transform[rank]
        scale = matrix[rank][column]
        matrix[rank] = [x/scale for x in matrix[rank]]
        transform[rank] = [x/scale for x in transform[rank]]
        for j in range(len(matrix)):
            if j != rank and matrix[j][column]:
                scale = matrix[j][column]
                matrix[j] = [x-scale*y for x, y in zip(matrix[j], matrix[rank])]
                transform[j] = [x-scale*y for x, y in zip(transform[j], transform[rank])]
        pivots.append(column)
        rank += 1
    return matrix[:rank], pivots, transform[:rank]


@lru_cache(maxsize=5)
def _geometry(L):
    _integer(L, 0, 4, 'L')
    basis = a.harmonic_basis(0, L)
    n = len(basis)
    gram = [[F(0) for _ in range(n)] for _ in range(n)]
    kinetic = [[F(0) for _ in range(n)] for _ in range(n)]
    for i, left in enumerate(basis):
        for j in range(i+1):
            mass = a.inner(left['polynomial'], basis[j]['polynomial'])
            energy = a.tangent_energy(left['polynomial'], basis[j]['polynomial'])
            if (i != j and mass) or (i == j and mass <= 0):
                raise ArithmeticError('Complete real harmonic basis is not orthogonal positive')
            if energy != left['degree']*(left['degree']+1)*mass:
                raise ArithmeticError('Integrated gradient does not match the harmonic Laplacian')
            gram[i][j] = gram[j][i] = mass
            kinetic[i][j] = kinetic[j][i] = energy
    return {'basis': basis, 'gram': gram, 'kinetic': kinetic,
            'mass': [gram[i][i] for i in range(n)]}


def full_harmonics(L=2):
    """V196: all (L+1)^2 real modes, with exact Gram and gradient forms."""
    data = deepcopy(_geometry(L))
    return {'L': L, 'dimension': len(data['basis']),
            'basis': [{**{k: b[k] for k in ('degree', 'm', 'part')},
                       'polynomial': a.encode(b['polynomial'])} for b in data['basis']],
            'gram': _wire_matrix(data['gram']), 'kinetic': _wire_matrix(data['kinetic']),
            'mass': list(map(str, data['mass'])), 'includes_constant': True,
            'includes_all_real_sine_cosine_components': True}


def constraint_nullspace(constraints, L=2):
    """V197: exact global moment rows, full-support check, cross-m nullspace."""
    constraints = moments(constraints)
    data = _geometry(L)
    basis, mass = data['basis'], data['mass']
    n = len(basis)
    rows, projections = [], []
    for function in constraints:
        row = [a.inner(function, b['polynomial']) for b in basis]
        coefficients = [value/weight for value, weight in zip(row, mass)]
        residual = a.inner(function, function)-sum((v*v*w for v, w in zip(coefficients, mass)), F(0))
        if residual != 0:
            raise ValueError('Complete constraint support exceeds retained harmonic head')
        rows.append(row)
        projections.append({'coefficients': list(map(str, coefficients)),
                            'omitted_L2_norm_squared': str(residual)})
    reduced, pivots, transform = _rref(rows, n)
    free = [j for j in range(n) if j not in pivots]
    t = [[F(0) for _ in free] for _ in range(n)]
    for j, column in enumerate(free):
        t[column][j] = F(1)
        for i, pivot in enumerate(pivots):
            t[pivot][j] = -reduced[i][column]
    target = [F(1)]+[F(0)]*(n-1)
    combination = [F(0)]*len(rows)
    remainder = target[:]
    for i, pivot in enumerate(pivots):
        coefficient = target[pivot]
        remainder = [x-coefficient*y for x, y in zip(remainder, reduced[i])]
        combination = [x+coefficient*y for x, y in zip(combination, transform[i])]
    constant_present = not any(remainder)
    if constant_present:
        residual = a.add({a.ZERO: F(1)}, *(a.scale(g, -v) for g, v in zip(constraints, combination)))
        if a.inner(residual, residual) != 0:
            raise ArithmeticError('Constant row-span identity failed on the sphere')
    return {'constraints': [a.encode(g) for g in constraints], 'retained_degree': L,
            'constraint_projections': projections, 'weighted_rows': _wire_matrix(rows),
            'rref': _wire_matrix(reduced), 'rank': len(pivots), 'pivots': pivots,
            'free_columns': free, 'T': _wire_matrix(t), 'dimension': len(free),
            'constant_in_row_span': constant_present,
            'constant_span_coefficients': list(map(str, combination)) if constant_present else None}


def trial_enclosure(q, constraints, trial, support_L=4):
    """V195: validate ALL moments, actual Ritz trial and global Poincare floor."""
    q, constraints, trial = a.polynomial(q), moments(constraints), a.polynomial(trial, 8)
    system = constraint_nullspace(constraints, support_L)
    mass = a.inner(trial, trial)
    inner_products = [a.inner(trial, g) for g in constraints]
    if mass <= 0 or any(inner_products):
        raise ValueError('Trial is zero or violates a recorded full-sphere moment')
    energy = a.tangent_energy(trial, trial)+a.integrate(a.multiply(q, a.multiply(trial, trial)))
    proof = a.potential_range(q)
    floor = 2 if system['constant_in_row_span'] else 0
    lower, upper = floor+F(proof['lower']), energy/mass
    return {'format': TRIAL, 'geometry': 'unit_S2', 'scope': SCOPE,
            'potential': a.encode(q), 'constraints': [a.encode(g) for g in constraints],
            'support_L': support_L, 'trial': a.encode(trial),
            'moment_inner_products': list(map(str, inner_products)),
            'mass': str(mass), 'energy': str(energy), 'range_proof': proof,
            'constant_in_row_span': system['constant_in_row_span'],
            'constant_span_coefficients': system['constant_span_coefficients'],
            'laplace_lower': floor, 'formal_assistant_checked': False,
            **arithmetic.enclosure_fields(lower, upper)}


def _congruence(matrix, t):
    n, dim = len(t), len(t[0])
    product = [[sum((matrix[i][j]*t[j][k] for j in range(n)), F(0)) for k in range(dim)]
               for i in range(n)]
    return [[sum((t[i][j]*product[i][k] for i in range(n)), F(0)) for k in range(dim)]
            for j in range(dim)]


def reduced_forms(q, constraints, L=2):
    """V198: global T^T A T and the complete, possibly dense T^T M T."""
    q = a.polynomial(q)
    data, system = _geometry(L), constraint_nullspace(constraints, L)
    if system['dimension'] == 0:
        raise ValueError('Constraints eliminate retained head; increase L without dropping constraints')
    t = [[F(x) for x in row] for row in system['T']]
    products = [a.multiply(q, b['polynomial']) for b in data['basis']]
    matrix = [[data['kinetic'][i][j]+a.inner(b['polynomial'], products[j])
               for j in range(len(products))] for i, b in enumerate(data['basis'])]
    reduced_a, reduced_m = _congruence(matrix, t), _congruence(data['gram'], t)
    if base.inertia(reduced_m) != [0, 0, len(reduced_m)]:
        raise ArithmeticError('Constrained mass is not positive definite')
    return {'retained_degree': L, 'harmonics': full_harmonics(L), 'constraint_system': system,
            'reduced_A': _wire_matrix(reduced_a), 'reduced_mass': _wire_matrix(reduced_m),
            'mass_is_diagonal': all(not reduced_m[i][j] for i in range(len(reduced_m)) for j in range(i))}


def transformed_tail(q, constraints, L=2):
    """V199: all T^T B columns for all real harmonics L+1..L+degree(q)."""
    q = a.polynomial(q)
    data, system = _geometry(L), constraint_nullspace(constraints, L)
    if not system['dimension']:
        raise ValueError('No admissible head directions')
    t = [[F(x) for x in row] for row in system['T']]
    degree = max(map(sum, q), default=0)
    tail = a.harmonic_basis(L+1, L+degree) if degree else []
    q_head = [a.multiply(q, b['polynomial']) for b in data['basis']]
    columns = []
    for b in tail:
        old = [a.inner(b['polynomial'], p) for p in q_head]
        values = [sum((t[i][j]*old[i] for i in range(len(t))), F(0))
                  for j in range(system['dimension'])]
        columns.append({**{k: b[k] for k in ('degree', 'm', 'part')},
                        'mass': str(a.inner(b['polynomial'], b['polynomial'])),
                        'entries': list(map(str, values))})
    return {'L': L, 'potential_degree': degree, 'columns': columns,
            'far_tail_rule': 'all degrees >L+degree(q) have zero head coupling; all degrees >L remain bounded',
            'complete_real_components': True}


class Kernel:
    def __init__(self, q, constraints, L=2):
        self.q, self.constraints = a.polynomial(q), moments(constraints)
        self.L = _integer(L, 0, 4, 'L')
        self.range = a.potential_range(self.q)
        self.used = a.polynomial(self.range['used'])
        self.qlo, self.qhi = F(self.range['lower']), F(self.range['upper'])
        self.forms = reduced_forms(self.used, self.constraints, L)
        self.tail = transformed_tail(self.used, self.constraints, L)
        self.a = [[F(x) for x in row] for row in self.forms['reduced_A']]
        self.mass = [[F(x) for x in row] for row in self.forms['reduced_mass']]
        self.dimension = len(self.mass)
        self.beta = (L+1)*(L+2)+self.qlo
        self.floor = (2 if self.forms['constraint_system']['constant_in_row_span'] else 0)+self.qlo

    def matrix(self, x, kind):
        x = exact(x)
        if kind not in ('lower', 'upper_tail', 'upper_ritz'):
            raise ValueError('Invalid Schur comparison kind')
        if kind != 'upper_ritz' and x >= self.beta:
            raise ValueError('Tail must be strictly positive')
        result = [[self.a[i][j]-x*self.mass[i][j] for j in range(self.dimension)]
                  for i in range(self.dimension)]
        if kind != 'upper_ritz':
            bound = self.qlo if kind == 'lower' else self.qhi
            for column in self.tail['columns']:
                degree, mass = column['degree'], F(column['mass'])
                denominator = mass*(degree*(degree+1)+bound-x)
                values = [(i, F(v)) for i, v in enumerate(column['entries']) if F(v)]
                for i, v in values:
                    for j, w in values:
                        result[i][j] -= v*w/denominator
        return result

    def count(self, x, kind, independent=False):
        matrix = self.matrix(x, kind)
        return base.inertia(matrix) if independent else arithmetic.banded_inertia(matrix)

    def lower_holds(self, x):
        return x < self.beta and self.count(x, 'lower')[0] == 0

    def upper_holds(self, x, kind):
        if kind != 'upper_ritz' and x >= self.beta:
            return False
        count = self.count(x, kind)
        return count[0]+count[1] >= 1


def _record(kernel, lower, upper, kind, snapped=False, independent=False):
    lower, upper = exact(lower), exact(upper)
    if type(snapped) is not bool or kind not in ('upper_tail', 'upper_ritz') or lower > upper:
        raise ValueError('Invalid comparison, endpoints or exactness flag')
    low_count, up_count = kernel.count(lower, 'lower', independent), kernel.count(upper, kind, independent)
    if low_count[0] or up_count[0]+up_count[1] < 1:
        raise ValueError('Exact inertia does not establish ground enclosure')
    if snapped and (lower != upper or lower not in [
            F(l*(l+1))+kernel.used.get(a.ZERO, F(0)) for l in range(kernel.L+1)]):
        raise ValueError('Invalid exact Laplace-level claim')
    return {'format': SPECTRAL, 'geometry': 'unit_S2', 'scope': SCOPE,
            'potential': a.encode(kernel.q), 'constraints': [a.encode(g) for g in kernel.constraints],
            'eigenvalue_index_in_constrained_space': 1, 'retained_degree': kernel.L,
            'forms': kernel.forms, 'range_proof': kernel.range, 'tail': kernel.tail,
            'infinite_tail_lower': str(kernel.beta), 'upper_kind': kind,
            'lower_inertia': low_count, 'upper_inertia': up_count,
            'exact_level_verified_by_inertia': snapped, 'formal_assistant_checked': False,
            **arithmetic.enclosure_fields(lower, upper)}


def exact_level(q, constraints, L=2):
    """V201: exact snapping is permitted ONLY when both inertia tests pass."""
    kernel = Kernel(q, constraints, L)
    for l in range(L+1):
        value = F(l*(l+1))+kernel.used.get(a.ZERO, F(0))
        if kernel.lower_holds(value) and kernel.upper_holds(value, 'upper_ritz'):
            return _record(kernel, value, value, 'upper_ritz', snapped=True)
    return None


def certify(q, constraints, L=2, bits=32, snap=True):
    """V200: complete constrained full-sphere Schur enclosure, including tail."""
    _integer(bits, 8, 64, 'bits')
    if type(snap) is not bool:
        raise ValueError('snap must be boolean')
    if snap:
        snapped = exact_level(q, constraints, L)
        if snapped is not None:
            return snapped
    kernel = Kernel(q, constraints, L)
    low, high = kernel.floor-1, F(L*(L+1))+kernel.qhi+1
    kind = 'upper_tail' if high < kernel.beta else 'upper_ritz'
    if not kernel.upper_holds(high, kind):
        raise ArithmeticError('Retained maximum-degree Ritz bracket failed')
    for _ in range(bits):
        mid = (low+high)/2
        if kernel.upper_holds(mid, kind):
            high = mid
        else:
            low = mid
    upper = high
    distance = F(1)
    low = kernel.floor-distance
    while not kernel.lower_holds(low):
        distance *= 2
        low = kernel.floor-distance
    high = min(upper, kernel.beta)
    for _ in range(bits):
        mid = (low+high)/2
        if kernel.lower_holds(mid):
            low = mid
        else:
            high = mid
    result = _record(kernel, low, upper, kind)
    if not verify_spectral(result):
        raise ArithmeticError('Independent spectral replay rejected')
    return result


def verify_spectral(cert, expected_q=None, expected_constraints=None):
    try:
        kernel = Kernel(cert['potential'], cert['constraints'], cert['retained_degree'])
        if expected_q is not None and a.polynomial(expected_q) != kernel.q:
            return False
        if expected_constraints is not None and moments(expected_constraints) != kernel.constraints:
            return False
        expected = _record(kernel, exact(cert['lower']), exact(cert['upper']), cert['upper_kind'],
                           cert['exact_level_verified_by_inertia'], independent=True)
        return same(cert, expected)
    except (ValueError, KeyError, TypeError, IndexError, ArithmeticError):
        return False


def equivalent_constraints(source, target_constraints):
    """V202: reuse a spectral proof for an exactly equivalent global row space."""
    if not verify_spectral(source):
        raise ValueError('Source spectral certificate failed replay')
    target = constraint_nullspace(target_constraints, source['retained_degree'])
    original = source['forms']['constraint_system']
    if target['rref'] != original['rref']:
        raise ValueError('Constraint row spaces differ on the complete retained basis')
    return {'format': TRANSFER, 'geometry': 'unit_S2', 'scope': SCOPE,
            'potential': source['potential'], 'constraints': target['constraints'],
            'source': source, 'target_constraint_system': target,
            'same_full_constraint_rowspace': True, 'formal_assistant_checked': False,
            **arithmetic.enclosure_fields(F(source['lower']), F(source['upper']))}


def adaptive(q, constraints, tolerance=F(1, 10**8), start_L=1, max_L=4, bits=32):
    """V203: never drop constraint support; increase head or return explicit open."""
    q, constraints = a.polynomial(q), moments(constraints)
    tolerance = exact(tolerance)
    if not F(1, 10**20) <= tolerance <= 1:
        raise ValueError('Tolerance outside [1e-20,1]')
    _integer(start_L, 0, 4, 'start_L'); _integer(max_L, start_L, 4, 'max_L')
    _integer(bits, 8, 64, 'bits')
    attempts, final = [], None
    for L in range(start_L, max_L+1):
        try:
            system = constraint_nullspace(constraints, L)
        except ValueError:
            attempts.append({'L': L, 'status': 'complete_constraint_support_exceeds_head'})
            continue
        if not system['dimension']:
            attempts.append({'L': L, 'status': 'constraints_eliminate_head'})
            continue
        final = certify(q, constraints, L=L, bits=bits)
        attempts.append({'L': L, 'status': 'certified', 'certificate': final})
        if F(final['exact_width']) <= tolerance:
            break
    status = ('unsupported_within_budget' if final is None else
              'certified_target_met' if F(final['exact_width']) <= tolerance else 'certified_bound_open_gap')
    return {'format': ADAPTIVE, 'geometry': 'unit_S2', 'scope': SCOPE,
            'potential': a.encode(q), 'constraints': [a.encode(g) for g in constraints],
            'tolerance': str(tolerance), 'start_L': start_L, 'max_L': max_L, 'bits': bits,
            'attempts': attempts, 'final_certificate': final, 'status': status,
            'constraint_support_was_never_truncated': True, 'formal_assistant_checked': False}


def verify_adaptive(cert):
    try:
        q, constraints = a.polynomial(cert['potential']), moments(cert['constraints'])
        start, maximum, bits = cert['start_L'], cert['max_L'], cert['bits']
        tolerance = exact(cert['tolerance'])
        _integer(start, 0, 4, 'start_L'); _integer(maximum, start, 4, 'max_L'); _integer(bits, 8, 64, 'bits')
        if not F(1, 10**20) <= tolerance <= 1 or not isinstance(cert['attempts'], list):
            return False
        last, met = None, False
        for i, attempt in enumerate(cert['attempts']):
            L = start+i
            if L > maximum or met or not same(attempt['L'], L):
                return False
            try:
                system = constraint_nullspace(constraints, L)
            except ValueError:
                expected = {'L': L, 'status': 'complete_constraint_support_exceeds_head'}
            else:
                if not system['dimension']:
                    expected = {'L': L, 'status': 'constraints_eliminate_head'}
                else:
                    child = attempt['certificate']
                    if not verify_spectral(child, q, constraints) or child['retained_degree'] != L:
                        return False
                    last = child
                    met = F(child['exact_width']) <= tolerance
                    expected = {'L': L, 'status': 'certified', 'certificate': child}
            if not same(attempt, expected):
                return False
        if not cert['attempts'] or (not met and len(cert['attempts']) != maximum-start+1):
            return False
        status = ('unsupported_within_budget' if last is None else
                  'certified_target_met' if met else 'certified_bound_open_gap')
        expected = {'format': ADAPTIVE, 'geometry': 'unit_S2', 'scope': SCOPE,
            'potential': a.encode(q), 'constraints': [a.encode(g) for g in constraints],
            'tolerance': str(tolerance), 'start_L': start, 'max_L': maximum, 'bits': bits,
            'attempts': cert['attempts'], 'final_certificate': last, 'status': status,
            'constraint_support_was_never_truncated': True, 'formal_assistant_checked': False}
        return same(cert, expected)
    except (ValueError, KeyError, TypeError, IndexError, ArithmeticError):
        return False


def _negative_direction(matrix):
    n = len(matrix)
    if not n:
        return None
    negative = next((i for i in range(n) if matrix[i][i] < 0), None)
    if negative is not None:
        return [F(i == negative) for i in range(n)]
    pivot = next((i for i in range(n) if matrix[i][i] > 0), None)
    if pivot is None:
        pair = next(((i, j) for i in range(n) for j in range(i+1, n) if matrix[i][j]), None)
        if pair is None:
            return None
        i, j = pair
        result = [F(0)]*n
        result[i], result[j] = F(1), F(-1 if matrix[i][j] > 0 else 1)
        return result
    other = [i for i in range(n) if i != pivot]
    d = matrix[pivot][pivot]
    schur = [[matrix[i][j]-matrix[i][pivot]*matrix[pivot][j]/d for j in other] for i in other]
    child = _negative_direction(schur)
    if child is None:
        return None
    result = [F(0)]*n
    for i, v in zip(other, child):
        result[i] = v
    result[pivot] = -sum((matrix[pivot][j]*result[j] for j in other), F(0))/d
    return result


def violating_trial(q, constraints, threshold, L=2):
    """Exact negative direction, lifted through T and optionally the Schur tail."""
    threshold = exact(threshold)
    kernel = Kernel(q, constraints, L)
    kind = 'upper_ritz'
    vector = _negative_direction(kernel.matrix(threshold, kind))
    if vector is None and threshold < kernel.beta:
        kind = 'upper_tail'
        vector = _negative_direction(kernel.matrix(threshold, kind))
    if vector is None:
        return None
    t = [[F(x) for x in row] for row in kernel.forms['constraint_system']['T']]
    head = [sum((v*w for v, w in zip(row, vector)), F(0)) for row in t]
    trial = a.add(*(a.scale(b['polynomial'], v) for b, v in zip(_geometry(L)['basis'], head)))
    if kind == 'upper_tail':
        for column in kernel.tail['columns']:
            l = column['degree']
            denominator = F(column['mass'])*(l*(l+1)+kernel.qhi-threshold)
            coefficient = -sum((F(v)*w for v, w in zip(column['entries'], vector)), F(0))/denominator
            trial = a.add(trial, a.scale(a.solid_harmonic(l, column['m'], column['part']), coefficient))
    proof = trial_enclosure(q, constraints, trial, support_L=L)
    if F(proof['upper']) >= threshold:
        raise ArithmeticError('Constructed negative direction did not strictly violate threshold')
    return proof


def inequality_record(spectral, threshold, trial=None):
    if not verify_spectral(spectral):
        raise ValueError('Verified complete global spectral proof required')
    threshold = exact(threshold)
    if trial is not None:
        if not verify(trial) or trial['format'] != TRIAL or trial['potential'] != spectral['potential'] or trial['constraints'] != spectral['constraints']:
            raise ValueError('Trial does not bind the complete original problem')
    status = ('proved' if F(spectral['lower']) >= threshold else
              'refuted_with_trial' if trial is not None and F(trial['upper']) < threshold else 'undetermined')
    return {'format': INEQUALITY, 'geometry': 'unit_S2', 'scope': SCOPE,
            'potential': spectral['potential'], 'constraints': spectral['constraints'],
            'statement': 'integral(grad_u_squared+q*u_squared)>=threshold*integral(u_squared)',
            'threshold': str(threshold), 'spectral_certificate': spectral, 'trial_certificate': trial,
            'status': status, 'spectral_existence_below_threshold': F(spectral['upper']) < threshold,
            'formal_assistant_checked': False}


def inequality(q, constraints, threshold, L=2, bits=32):
    """V204: full-sphere prove / explicit refutation / honest undecided."""
    spectral = certify(q, constraints, L, bits)
    trial = None if F(spectral['lower']) >= exact(threshold) else violating_trial(q, constraints, threshold, L)
    return inequality_record(spectral, threshold, trial)


def verify(cert, expected_q=None, expected_constraints=None):
    try:
        if not isinstance(cert, dict):
            return False
        if expected_q is not None and a.polynomial(cert['potential']) != a.polynomial(expected_q):
            return False
        if expected_constraints is not None and moments(cert['constraints']) != moments(expected_constraints):
            return False
        kind = cert.get('format')
        if kind == SPECTRAL:
            return verify_spectral(cert)
        if kind == TRIAL:
            expected = trial_enclosure(cert['potential'], cert['constraints'], cert['trial'], cert['support_L'])
        elif kind == TRANSFER:
            expected = equivalent_constraints(cert['source'], cert['constraints'])
        elif kind == ADAPTIVE:
            return verify_adaptive(cert)
        elif kind == INEQUALITY:
            expected = inequality_record(cert['spectral_certificate'], cert['threshold'], cert['trial_certificate'])
        else:
            return False
        return same(cert, expected)
    except (ValueError, KeyError, TypeError, IndexError, ArithmeticError):
        return False
