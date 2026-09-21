"""V225--V234: full-sphere mean-zero spectra using exact reflection blocks.

Every accepted input is a rational xyz polynomial of degree <=4. Reflection
invariance is checked as an exact identity on S². All eight characters remain
in the final proof, either with a spectral certificate or an analytic floor.
The high-degree harmonic construction is local; previous releases stay frozen.
"""
from fractions import Fraction as F
from math import factorial, comb
from itertools import product
import copy
import time
from support import exact, same, canonical, digest, base, arithmetic, old_anisotropic as a

PARITIES = tuple(product((0, 1), repeat=3))
SCOPE = 'all_real_mean_zero_H1_functions_on_unit_S2'
BLOCK_FORMAT = 'reflection_block_ground_v230'
FULL_FORMAT = 'eight_reflection_blocks_ground_v231'
DRIVER_FORMAT = 'adaptive_reflection_ground_v234'
MAX_L = 12
MAX_HARMONIC_DEGREE = 16


def _json(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        if any(isinstance(k, tuple) for k in value):
            return a.encode(value)
        return {k: _json(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json(v) for v in value]
    return value


def parity(raw):
    if not isinstance(raw, (list, tuple)) or len(raw) != 3 or any(
            type(v) is not int or v not in (0, 1) for v in raw):
        raise ValueError('Reflection parity is exactly three integer bits')
    return tuple(raw)


def _potential(q):
    return a.polynomial(q, maximum_degree=4)


def potential_range(q):
    """V225: exact simplex Bernstein range for even quadratic/quartic potentials.

    On S², X=x²,Y=y²,Z=z² lie in the simplex X+Y+Z=1. Homogenizing to
    their total degree gives Bernstein coefficients; their convex hull bounds q.
    General inputs receive a safe coefficient bound, never an unjustified split.
    """
    q = _potential(q)
    reduction = a.sphere_reduce(q)
    p = _potential(reduction['reduced'])
    if all(all(v % 2 == 0 for v in e) for e in p):
        squared = {tuple(v//2 for v in e): c for e, c in p.items()}
        degree = max(map(sum, squared), default=0)
        total = {(1, 0, 0): F(1), (0, 1, 0): F(1), (0, 0, 1): F(1)}
        homogeneous = {}
        for e, c in squared.items():
            homogeneous = a.add(homogeneous, a.multiply({e: c}, a.power(total, degree-sum(e))))
        bernstein = []
        for i in range(degree+1):
            for j in range(degree-i+1):
                e = (i, j, degree-i-j)
                coefficient = homogeneous.get(e, F(0))*F(factorial(i)*factorial(j)*factorial(e[2]), factorial(degree))
                bernstein.append({'powers': list(e), 'coefficient': str(coefficient)})
        lower = min(F(v['coefficient']) for v in bernstein)
        upper = max(F(v['coefficient']) for v in bernstein)
        rule = 'convex_hull_of_exact_simplex_Bernstein_coefficients'
    else:
        degree, bernstein = None, []
        constant = p.get((0, 0, 0), F(0))
        radius = sum((abs(c) for e, c in p.items() if e != (0, 0, 0)), F(0))
        lower, upper = constant-radius, constant+radius
        rule = 'coefficient_absolute_value_bound_on_unit_sphere'
    return {'potential': a.encode(q), 'effective_potential': a.encode(p),
            'sphere_identity': reduction, 'rule': rule, 'simplex_degree': degree,
            'Bernstein_coefficients': bernstein, 'lower': str(lower), 'upper': str(upper)}


def verify_range(certificate, expected_q=None):
    try:
        if expected_q is not None and _potential(expected_q) != _potential(certificate['potential']):
            return False
        return same(certificate, potential_range(certificate['potential']))
    except (ValueError, TypeError, KeyError, ArithmeticError):
        return False


def reflection_symmetry(q):
    """V226: verify all eight sign flips as exact sphere identities."""
    q = _potential(q)
    effective = _potential(a.sphere_reduce(q)['reduced'])
    checks = []
    for flip in PARITIES:
        reflected = {e: c*(-1)**sum(e[j]*flip[j] for j in range(3)) for e, c in effective.items()}
        difference = a.add(effective, a.scale(reflected, -1))
        checks.append({'flip': list(flip), 'difference': a.encode(difference), 'invariant': not difference})
    return {'potential': a.encode(q), 'effective_potential': a.encode(effective),
            'all_eight_reflections': all(row['invariant'] for row in checks), 'checks': checks}


def parity_project(p, character):
    """Exact character projector: retain precisely matching monomial parities."""
    p = a.polynomial(p, maximum_degree=MAX_HARMONIC_DEGREE)
    character = parity(character)
    return {e: c for e, c in p.items() if tuple(v % 2 for v in e) == character}


def _require_even(q):
    proof = reflection_symmetry(q)
    if not proof['all_eight_reflections']:
        raise ValueError('Potential is not invariant under all coordinate reflections; eight-block splitting is invalid')
    return _potential(proof['effective_potential'])


def lowest_degree(character):
    p = parity(character)
    return sum(p) if sum(p) else 2


def next_omitted_degree(character, L):
    if type(L) is not int or not 0 <= L <= MAX_L:
        raise ValueError('Retained degree must be an integer in 0..12')
    degree = lowest_degree(character)
    while degree <= L:
        degree += 2
    return degree


def _solid(l, m, part):
    if not 0 <= m <= l <= MAX_HARMONIC_DEGREE:
        raise ValueError('Harmonic construction degree cap exceeded')
    azimuth = {(m-j, j, 0): F(comb(m, j)*(-1)**(j//2))
               for j in range(m+1) if j % 2 == (part == 'imag')}
    result = {}
    for k in range((l-m)//2+1):
        coefficient = F((-1)**k*factorial(2*l-2*k),
                        2**l*factorial(k)*factorial(l-k)*factorial(l-2*k-m))
        result = a.add(result, a.multiply(a.multiply(azimuth, {(0, 0, l-2*k-m): coefficient}),
                                           a.power(a.R2, k)))
    if a.add(*(a.derivative(a.derivative(result, j), j) for j in range(3))):
        raise ArithmeticError('Solid harmonic failed exact Laplace check')
    return result


def degree_basis(character, degree):
    """V227: every harmonic for one exact reflection character at one degree."""
    p = parity(character)
    if type(degree) is not int or not 0 <= degree <= MAX_HARMONIC_DEGREE:
        raise ValueError('Unsupported harmonic degree')
    if degree < sum(p) or (degree-sum(p)) % 2:
        return []
    result = []
    for m in range(degree+1):
        for part in (('real',) if m == 0 else ('real', 'imag')):
            expected = (m % 2, 0, (degree-m) % 2) if part == 'real' else ((m-1) % 2, 1, (degree-m) % 2)
            if expected != p:
                continue
            poly = _solid(degree, m, part)
            if parity_project(poly, p) != poly:
                raise ArithmeticError('Constructed harmonic has the wrong reflection character')
            result.append({'degree': degree, 'm': m, 'part': part, 'polynomial': poly,
                           'mass': a.inner(poly, poly)})
    expected_dimension = (degree-sum(p))//2+1
    if len(result) != expected_dimension:
        raise ArithmeticError('Parity harmonic dimension is incomplete')
    return result


def parity_basis(character, L):
    p = parity(character)
    if type(L) is not int or not 0 <= L <= MAX_L:
        raise ValueError('Retained degree must be an integer in 0..12')
    basis = [b for l in range(lowest_degree(p), L+1, 2) for b in degree_basis(p, l)]
    for i, b in enumerate(basis):
        if a.integrate(b['polynomial']) != 0 or b['mass'] <= 0:
            raise ArithmeticError('Mean-zero or positive mass check failed')
        for c in basis[:i]:
            if a.inner(b['polynomial'], c['polynomial']) != 0:
                raise ArithmeticError('Parity basis is not orthogonal')
    return basis


def block_form(q, character, L):
    """V228: exact finite block form; symmetry proves other characters decouple."""
    q = _potential(q)
    effective = _require_even(q)
    p = parity(character)
    basis = parity_basis(p, L)
    if not basis:
        raise ValueError('No retained nonconstant harmonic for the requested character')
    n = len(basis)
    products = [a.multiply(effective, b['polynomial']) for b in basis]
    matrix = [[F(0) for _ in basis] for _ in basis]
    for i in range(n):
        for j in range(i+1):
            value = a.inner(basis[i]['polynomial'], products[j])
            if i == j:
                l = basis[i]['degree']
                value += l*(l+1)*basis[i]['mass']
            matrix[i][j] = matrix[j][i] = value
    return {'potential': a.encode(q), 'parity': list(p), 'retained_degree': L,
            'basis': basis, 'A': matrix, 'mass': [b['mass'] for b in basis]}


def cross_block_form(q, first, second, L=3):
    """Independent exact cross-form check, also usable to expose invalid splitting."""
    q = _potential(q)
    first, second = parity(first), parity(second)
    if first == second:
        raise ValueError('Two different characters required')
    left, right = parity_basis(first, L), parity_basis(second, L)
    return [[a.inner(b['polynomial'], a.multiply(q, c['polynomial']))
             for c in right] for b in left]


def tail_couplings(q, character, L, finite=None):
    """V229: complete same-character tail and its actual next degree."""
    q = _potential(q)
    effective = _require_even(q)
    p = parity(character)
    finite = block_form(q, p, L) if finite is None else finite
    last = max(b['degree'] for b in finite['basis'])
    d = max(map(sum, effective), default=0)
    start = next_omitted_degree(p, L)
    columns = []
    products = [a.multiply(effective, b['polynomial']) for b in finite['basis']]
    for l in range(start, last+d+1, 2):
        for b in degree_basis(p, l):
            columns.append({'degree': l, 'm': b['m'], 'part': b['part'], 'mass': b['mass'],
                            'entries': [a.inner(b['polynomial'], v) for v in products]})
    return {'first_omitted_degree': start, 'last_retained_degree': last,
            'last_potentially_coupled_degree': last+d, 'columns': columns}


class FormCache:
    """V233: immutable-copy form/range reuse bound to original q, parity and L."""
    def __init__(self):
        self._ranges, self._forms = {}, {}
        self.stats = {'range_builds': 0, 'range_hits': 0, 'form_builds': 0, 'form_hits': 0,
                      'built_dimensions': [], 'inertia_calls': 0}

    def range(self, q):
        key = canonical(a.encode(_potential(q)))
        if key in self._ranges:
            self.stats['range_hits'] += 1
        else:
            self._ranges[key] = potential_range(q)
            self.stats['range_builds'] += 1
        return copy.deepcopy(self._ranges[key])

    def form(self, q, character, L):
        q, p = _potential(q), parity(character)
        key = (canonical(a.encode(q)), p, L)
        if key in self._forms:
            self.stats['form_hits'] += 1
        else:
            finite = block_form(q, p, L)
            tail = tail_couplings(q, p, L, finite)
            self._forms[key] = (finite, tail)
            self.stats['form_builds'] += 1
            self.stats['built_dimensions'].append(len(finite['mass']))
        finite, tail = copy.deepcopy(self._forms[key])
        if finite['potential'] != a.encode(q) or finite['parity'] != list(p) or finite['retained_degree'] != L:
            raise ValueError('Cached form has a different input, parity or degree')
        return finite, tail


class Kernel:
    def __init__(self, q, character, L=3, cache=None):
        self.q, self.parity, self.L = _potential(q), parity(character), L
        _require_even(self.q)
        self.cache = FormCache() if cache is None else cache
        if not isinstance(self.cache, FormCache):
            raise ValueError('An internal FormCache object is required')
        self.range = self.cache.range(self.q)
        self.qlo, self.qhi = F(self.range['lower']), F(self.range['upper'])
        self.finite, self.tail = self.cache.form(self.q, self.parity, L)
        self.a, self.mass = self.finite['A'], self.finite['mass']
        self.first = lowest_degree(self.parity)
        nxt = self.tail['first_omitted_degree']
        self.beta = nxt*(nxt+1)+self.qlo

    def matrix(self, x, kind):
        x = exact(x)
        if kind not in ('lower', 'upper_tail', 'upper_ritz'):
            raise ValueError('Unknown Schur comparison')
        if kind != 'upper_ritz' and x >= self.beta:
            raise ValueError('The complete same-parity tail must be positive')
        matrix = [row[:] for row in self.a]
        for i, mass in enumerate(self.mass):
            matrix[i][i] -= x*mass
        if kind != 'upper_ritz':
            for column in self.tail['columns']:
                l, mass = column['degree'], column['mass']
                bound = self.qlo if kind == 'lower' else self.qhi
                denominator = mass*(l*(l+1)+bound-x)
                entries = [(i, c) for i, c in enumerate(column['entries']) if c]
                for i, c in entries:
                    for j, d in entries:
                        matrix[i][j] -= c*d/denominator
        return matrix

    def count(self, x, kind, independent=False):
        self.cache.stats['inertia_calls'] += 1
        matrix = self.matrix(x, kind)
        return base.inertia(matrix) if independent else arithmetic.banded_inertia(matrix)

    def lower_holds(self, x):
        return x < self.beta and self.count(x, 'lower')[0] == 0

    def upper_holds(self, x, kind):
        if kind != 'upper_ritz' and x >= self.beta:
            return False
        counts = self.count(x, kind)
        return counts[0]+counts[1] >= 1


def _block_record(kernel, lo, hi, kind, independent=True):
    lo, hi = exact(lo), exact(hi)
    if lo > hi or lo >= kernel.beta:
        raise ValueError('Invalid block enclosure or tail positivity')
    lower_counts = kernel.count(lo, 'lower', independent)
    upper_counts = kernel.count(hi, kind, independent)
    if lower_counts[0] or upper_counts[0]+upper_counts[1] < 1:
        raise ValueError('Exact block inertia does not prove the endpoints')
    labels = [{k: b[k] for k in ('degree', 'm', 'part')} for b in kernel.finite['basis']]
    return {'format': BLOCK_FORMAT, 'geometry': 'unit_S2', 'scope': 'one_reflection_character_in_mean_zero_H1',
            'potential': a.encode(kernel.q), 'parity': list(kernel.parity), 'mean_zero': True,
            'projection': 'remove_l0_constant_in_000', 'lowest_degree': kernel.first,
            'retained_degree': kernel.L, 'basis': labels, 'mass': [str(v) for v in kernel.mass],
            'finite_form': _json(kernel.a), 'range_proof': kernel.range,
            'complete_tail': _json(kernel.tail), 'tail_lower': str(kernel.beta),
            'upper_kind': kind, 'lower_inertia': list(lower_counts), 'upper_inertia': list(upper_counts),
            'formal_assistant_checked': False, **arithmetic.enclosure_fields(lo, hi)}


def block_ground(q, character, L=3, bits=40, cache=None):
    """V230: a rigorously projected reflection block including its infinite tail."""
    if type(bits) is not int or not 8 <= bits <= 160:
        raise ValueError('Require bits in 8..160')
    kernel = Kernel(q, character, L, cache)
    qeff = _potential(kernel.range['effective_potential'])
    if all(e == (0, 0, 0) for e in qeff):
        val = kernel.first*(kernel.first+1)+qeff.get((0, 0, 0), F(0))
        cert = _block_record(kernel, val, val, 'upper_ritz')
    else:
        start = F(kernel.first*(kernel.first+1))
        low, high = start+kernel.qlo-1, start+kernel.qhi+1
        kind = 'upper_tail' if high < kernel.beta else 'upper_ritz'
        for _ in range(bits):
            mid = (low+high)/2
            if kernel.upper_holds(mid, kind):
                high = mid
            else:
                low = mid
        upper = high
        low, high = start+kernel.qlo-1, upper
        distance = F(1)
        while not kernel.lower_holds(low):
            distance *= 2
            low = start+kernel.qlo-distance
        for _ in range(bits):
            mid = (low+high)/2
            if kernel.lower_holds(mid):
                low = mid
            else:
                high = mid
        cert = _block_record(kernel, low, upper, kind)
    if not verify_block(cert, expected_q=q, expected_parity=character):
        raise ArithmeticError('Independent block replay failed')
    return cert


def verify_block(cert, expected_q=None, expected_parity=None):
    try:
        if expected_q is not None and _potential(expected_q) != _potential(cert['potential']):
            return False
        if expected_parity is not None and parity(expected_parity) != parity(cert['parity']):
            return False
        kernel = Kernel(cert['potential'], cert['parity'], cert['retained_degree'])
        return same(cert, _block_record(kernel, exact(cert['lower']), exact(cert['upper']), cert['upper_kind']))
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


def analytic_floor(q, character):
    q, p = _potential(q), parity(character)
    _require_even(q)
    bound = potential_range(q)
    l = lowest_degree(p)
    return {'format': 'reflection_analytic_floor_v230', 'potential': a.encode(q), 'parity': list(p),
            'mean_zero': True, 'lowest_degree': l, 'range_proof': bound,
            'lower': str(F(l*(l+1))+F(bound['lower'])),
            'rule': 'entire_character_Laplace_floor_plus_global_potential_lower'}


def _full_record(q, entries, target):
    q, target = _potential(q), exact(target)
    symmetry = reflection_symmetry(q)
    if not symmetry['all_eight_reflections'] or not isinstance(entries, list) or len(entries) != 8:
        raise ValueError('All eight invariant reflection characters must be accounted for')
    if not F(1, 10**30) <= target <= 1:
        raise ValueError('Target must lie between 1e-30 and 1')
    uppers = []
    for p, entry in zip(PARITIES, entries):
        if entry.get('format') == BLOCK_FORMAT:
            if not verify_block(entry, expected_q=q, expected_parity=p):
                raise ValueError('Wrong or invalid reflection block proof')
            uppers.append(F(entry['upper']))
        elif not same(entry, analytic_floor(q, p)):
            raise ValueError('Wrong or invalid whole-character analytic floor')
    if not uppers:
        raise ValueError('At least one actual upper-bound block certificate is required')
    lower, upper = min(F(e['lower']) for e in entries), min(uppers)
    return {'format': FULL_FORMAT, 'geometry': 'unit_S2', 'scope': SCOPE,
            'potential': a.encode(q), 'mean_zero': True, 'eigenvalue_index_in_constrained_space': 1,
            'symmetry_proof': symmetry, 'characters': entries, 'target_width': str(target),
            'status': 'target_met' if upper-lower <= target else 'certified_open_gap',
            'formal_assistant_checked': False, **arithmetic.enclosure_fields(lower, upper)}


def full_ground(q, L=3, bits=40, target=F(1, 10**8), cache=None, max_L=None):
    """V231: merge all eight blocks, safely pruning by whole-character floors."""
    q = _potential(q)
    _require_even(q)
    if type(L) is not int or not 1 <= L <= MAX_L:
        raise ValueError('Initial retained degree must be an integer in 1..12')
    if max_L is not None and (type(max_L) is not int or not L <= max_L <= MAX_L):
        raise ValueError('Initial degree cap must be an integer with L <= max_L <=12')
    cache = FormCache() if cache is None else cache
    entries = [analytic_floor(q, p) for p in PARITIES]
    upper = None
    for index in sorted(range(8), key=lambda i: F(entries[i]['lower'])):
        if upper is not None and F(entries[index]['lower']) >= upper:
            continue
        p = PARITIES[index]
        size = max(L, lowest_degree(p))
        if max_L is not None and size > max_L:
            continue
        certificate = block_ground(q, p, size, bits, cache)
        entries[index] = certificate
        value = F(certificate['upper'])
        upper = value if upper is None else min(upper, value)
    return _full_record(q, entries, target)


def verify_full(cert, expected_q=None):
    try:
        if expected_q is not None and _potential(expected_q) != _potential(cert['potential']):
            return False
        expected = _full_record(cert['potential'], cert['characters'], exact(cert['target_width']))
        return same(cert, expected)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


def error_diagnosis(q, character, L=3, bits=32, extra_bits=24):
    """V232: actual fixed-space precision and next-parity-degree experiments."""
    p = parity(character)
    cache = FormCache()
    low = block_ground(q, p, L, bits, cache)
    high = block_ground(q, p, L, bits+extra_bits, cache)
    next_L = next_omitted_degree(p, L)
    grown = block_ground(q, p, next_L, bits+extra_bits, cache)
    wlo, whi, wgrown = (F(c['exact_width']) for c in (low, high, grown))
    return {'potential': a.encode(_potential(q)), 'parity': list(p), 'retained_degree': L,
            'bits': [bits, bits+extra_bits], 'next_same_parity_degree': next_L,
            'precision_low': low, 'precision_high': high, 'larger_block': grown,
            'fixed_space_width_ratio': str(whi/wlo) if wlo else None,
            'larger_space_width_ratio': str(wgrown/whi) if whi else None,
            'diagnosis': 'tail_dominated' if whi*2 > wlo and wgrown*2 < whi else 'mixed_or_rounding_limited',
            'work': copy.deepcopy(cache.stats)}


def adaptive_ground(q, target=F(1, 10**8), start_L=3, max_L=9, bits=40,
                    wall_budget_seconds=90):
    """V234: meet the original full-space width or return an explicit open gap.

    Only characters whose current lower bound can compete with the best upper
    are refined. Other full-space lower proofs remain in the final certificate.
    Time is checked between mathematical operations; a proof is never cut short.
    """
    q, target = _potential(q), exact(target)
    if type(start_L) is not int or type(max_L) is not int or not 1 <= start_L <= max_L <= MAX_L:
        raise ValueError('Require integer 1 <= start_L <= max_L <=12')
    if type(bits) is not int or not 8 <= bits <= 136:
        raise ValueError('Adaptive initial bits must lie in 8..136')
    if type(wall_budget_seconds) not in (int, float) or not 0 < wall_budget_seconds <= 3600:
        raise ValueError('Wall budget must be positive and at most 3600 seconds')
    _require_even(q)
    began = time.perf_counter()
    cache = FormCache()
    current = full_ground(q, start_L, bits, target, cache, max_L=max_L)
    trace = [{'action': 'initial_all_character_coverage', 'lower': current['lower'],
              'upper': current['upper'], 'width': current['exact_width'],
              'elapsed_seconds': time.perf_counter()-began}]
    block_bits = {p: bits for p in PARITIES}
    limit_reason = None
    for _ in range(64):
        if F(current['exact_width']) <= target:
            break
        if time.perf_counter()-began > wall_budget_seconds:
            limit_reason = 'wall_budget_reached_between_operations'
            break
        entries = copy.deepcopy(current['characters'])
        best_upper = F(current['upper'])
        candidates = [i for i, e in enumerate(entries) if F(e['lower']) < best_upper-target]
        if not candidates:
            limit_reason = 'no_refinable_competitor'
            break
        index = min(candidates, key=lambda i: F(entries[i]['lower']))
        p, old = PARITIES[index], entries[index]
        if old['format'] != BLOCK_FORMAT:
            new_L, action = max(start_L, lowest_degree(p)), 'replace_analytic_floor_with_block'
            if new_L > max_L:
                limit_reason = 'retained_degree_cap_reached'
                break
        else:
            rounding_scale = (lowest_degree(p)*(lowest_degree(p)+1)
                              +F(old['range_proof']['upper'])-F(old['range_proof']['lower'])+2)/2**block_bits[p]
            if rounding_scale > target/16 and block_bits[p]+16 <= 160:
                block_bits[p] += 16
                new_L, action = old['retained_degree'], 'increase_bits_in_existing_block'
            else:
                new_L, action = next_omitted_degree(p, old['retained_degree']), 'increase_competing_character_degree'
                if new_L > max_L:
                    limit_reason = 'retained_degree_cap_reached'
                    break
        new = block_ground(q, p, new_L, block_bits[p], cache)
        entries[index] = new
        updated = _full_record(q, entries, target)
        trace.append({'action': action, 'parity': list(p), 'retained_degree': new_L,
                      'bits': block_bits[p], 'dimension': len(new['mass']),
                      'lower': updated['lower'], 'upper': updated['upper'], 'width': updated['exact_width'],
                      'elapsed_seconds': time.perf_counter()-began})
        current = updated
    else:
        limit_reason = 'adaptive_iteration_cap_reached'
    if not verify_full(current, expected_q=q):
        raise ArithmeticError('Final all-character certificate failed independent replay')
    elapsed = time.perf_counter()-began
    return {'format': DRIVER_FORMAT, 'potential': a.encode(q), 'scope': SCOPE,
            'requested_width': str(target), 'status': current['status'], 'certificate': current,
            'trace': trace, 'work': copy.deepcopy(cache.stats), 'elapsed_seconds': elapsed,
            'requested_limits': {'start_L': start_L, 'max_L': max_L, 'initial_bits': bits},
            'wall_budget_seconds': wall_budget_seconds, 'within_wall_budget': elapsed <= wall_budget_seconds,
            'limit_reason': limit_reason, 'formal_assistant_checked': False}


def verify_result(result, expected_q=None, expected_width=None):
    """Verify mathematical acceptance; runtime/cache counters are diagnostics."""
    try:
        if result['format'] != DRIVER_FORMAT or result['scope'] != SCOPE:
            return False
        q = _potential(result['potential'])
        if expected_q is not None and q != _potential(expected_q):
            return False
        target = exact(result['requested_width'])
        if expected_width is not None and target != exact(expected_width):
            return False
        certificate = result['certificate']
        if not verify_full(certificate, expected_q=q) or certificate['target_width'] != str(target):
            return False
        return result['status'] == certificate['status'] and result['formal_assistant_checked'] is False
    except (ValueError, TypeError, KeyError, ArithmeticError):
        return False
