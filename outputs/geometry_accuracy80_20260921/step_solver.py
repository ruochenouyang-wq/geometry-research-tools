"""Exact C1 piecewise trials for an axisymmetric step potential.

Search is rounded rational inverse iteration. Acceptance recomputes the full
strong residual by direct polynomial integration and replays the infinite
space source; no finite Ritz matrix is used as a spectral lower bound.
"""
from copy import deepcopy
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
from time import perf_counter
import hashlib
import json
import sys

sys.dont_write_bytecode = True
FROZEN = Path(__file__).resolve().parent.parent / 'geometry_precision_followup'
if str(FROZEN) not in sys.path:
    sys.path.insert(0, str(FROZEN))
import direct_moments as direct

FORMAT = 'step_c1_strong_residual_temple_v1'
RUN_FORMAT = 'step_c1_search_v1'


def _r(x):
    return direct.f.rational(x)


def _same(a, b):
    return json.dumps(a, sort_keys=True, separators=(',', ':')) == json.dumps(b, sort_keys=True, separators=(',', ':'))


def _digest(x):
    return hashlib.sha256(json.dumps(x, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def normalize(function):
    q = direct.normalize(function)
    if q['profile'] != 'step':
        raise ValueError('This solver requires an original step profile')
    return q


def basis(degree=6):
    direct.integer(degree, 2, 12, 'degree')
    return ((0, 0), (0, 1)) + tuple((side, k) for k in range(2, degree+1) for side in (-1, 1))


def _basis(raw):
    if not isinstance(raw, (list, tuple)) or not 2 <= len(raw) <= 24:
        raise ValueError('Require 2..24 C1 basis terms')
    result = []
    for term in raw:
        if not isinstance(term, (list, tuple)) or len(term) != 2:
            raise ValueError('Basis term is [side, power]')
        side, k = term
        if type(side) is not int or side not in (-1, 0, 1):
            raise ValueError('Invalid side')
        direct.integer(k, 0, 12, 'power')
        if side != 0 and k < 2:
            raise ValueError('One-sided powers below two violate C1 matching')
        result.append((side, k))
    if len(set(result)) != len(result):
        raise ValueError('Duplicate basis terms')
    # Avoid exact linear dependencies between global and one-sided powers.
    if any((0, k) in result and (-1, k) in result and (1, k) in result for _, k in result):
        raise ValueError('Linearly dependent split basis')
    return tuple(result)


def _mode(m):
    return direct.integer(m, 0, 1, 'azimuth_m')


def _poly(term, side):
    support, k = term
    return {k: F(1)} if support in (0, side) else {}


def _action(poly, m, potential):
    out = {}
    for k, a in poly.items():
        out[k] = out.get(k, F(0)) + a*((k+m)*(k+m+1)+potential)
        if k >= 2:
            out[k-2] = out.get(k-2, F(0))-a*k*(k-1)
    return {k: a for k, a in out.items() if a}


def _moment(k, side, m):
    # Exact integral of z^k (1-z^2)^m on one hemisphere (m=0 or 1).
    sign = -1 if side == -1 and k % 2 else 1
    return sign*(F(1, k+1)-(F(1, k+3) if m else 0))


def _inner(a, b, side, m):
    return sum((x*y*_moment(i+j, side, m) for i, x in a.items() for j, y in b.items()), F(0))


@lru_cache(maxsize=32)
def _matrices(amplitude, offset, m, terms):
    n = len(terms)
    M = [[F(0)]*n for _ in range(n)]
    H = [[F(0)]*n for _ in range(n)]
    R = [[F(0)]*n for _ in range(n)]
    for side in (-1, 1):
        q = offset+(amplitude if side == 1 else 0)
        ps = [_poly(t, side) for t in terms]
        hs = [_action(p, m, q) for p in ps]
        for i in range(n):
            for j in range(i, n):
                M[i][j] += _inner(ps[i], ps[j], side, m)
                H[i][j] += _inner(ps[i], hs[j], side, m)
                R[i][j] += _inner(hs[i], hs[j], side, m)
    for A in (M, H, R):
        for i in range(n):
            for j in range(i):
                A[i][j] = A[j][i]
    return tuple(tuple(tuple(row) for row in A) for A in (M, H, R))


def trial_matrices(function, terms=None, m=1, degree=6):
    q = normalize(function); m = _mode(m)
    terms = basis(degree) if terms is None else _basis(terms)
    return _matrices(F(q['amplitude']), F(q['offset']), m, terms)


def _vector(coefficients, terms):
    if not isinstance(coefficients, (list, tuple)) or len(coefficients) != len(terms):
        raise ValueError('One exact coefficient per basis term is required')
    v = tuple(_r(x) for x in coefficients)
    if not any(v):
        raise ValueError('Zero trial is not allowed')
    return v


def _statistics(mass, energy, norm):
    if mass <= 0:
        raise ValueError('Nonpositive mass')
    mu = energy/mass
    variance = norm/mass-mu*mu
    if variance < 0:
        raise ArithmeticError('Negative full residual')
    return {'mass': str(mass), 'operator_form': str(energy),
            'operator_norm_squared': str(norm), 'rayleigh': str(mu),
            'residual_squared': str(variance)}


def trial_statistics(function, terms, coefficients, m=1):
    terms = _basis(terms); v = _vector(coefficients, terms)
    matrices = trial_matrices(function, terms, m)
    values = [sum((v[i]*A[i][j]*v[j] for i in range(len(v)) for j in range(len(v))), F(0)) for A in matrices]
    return _statistics(*values)


def _direct_statistics(q, terms, coefficients, m):
    """Verifier path: combine the trial first; never form its Gram matrices."""
    v = _vector(coefficients, terms)
    mass = energy = norm = F(0)
    for side in (-1, 1):
        p = {}
        for (support, k), a in zip(terms, v):
            if support in (0, side):
                p[k] = p.get(k, F(0))+a
        h = _action(p, m, F(q['offset'])+(F(q['amplitude']) if side == 1 else 0))
        mass += _inner(p, p, side, m)
        energy += _inner(p, h, side, m)
        norm += _inner(h, h, side, m)
    return _statistics(mass, energy, norm)


def _ldl(A):
    n = len(A); L = [[F(i == j) for j in range(n)] for i in range(n)]; D = []
    for i in range(n):
        d = A[i][i]-sum((L[i][k]**2*D[k] for k in range(i)), F(0))
        if d <= 0:
            raise ValueError('Search shift does not give a positive definite matrix')
        D.append(d)
        for j in range(i+1, n):
            L[j][i] = (A[j][i]-sum((L[j][k]*L[i][k]*D[k] for k in range(i)), F(0)))/d
    def solve(rhs):
        y = []
        for i in range(n):
            y.append(rhs[i]-sum((L[i][j]*y[j] for j in range(i)), F(0)))
        x = [F(0)]*n
        for i in range(n-1, -1, -1):
            x[i] = y[i]/D[i]-sum((L[j][i]*x[j] for j in range(i+1, n)), F(0))
        return x
    return solve


def propose_trial(function, degree=6, m=1, precision_bits=96, iterations=12):
    q = normalize(function); m = _mode(m); terms = basis(degree)
    direct.integer(precision_bits, 32, 192, 'precision_bits')
    direct.integer(iterations, 0, 32, 'iterations')
    M, H, _ = trial_matrices(q, terms, m)
    qmin = min(F(q['offset']), F(q['offset'])+F(q['amplitude']))
    shift = m*(m+1)+qmin-F(1, 16)
    solve = _ldl([[H[i][j]-shift*M[i][j] for j in range(len(terms))] for i in range(len(terms))])
    v = [F(i == 0) for i in range(len(terms))]; denominator = 2**precision_bits
    for _ in range(iterations):
        rhs = [sum((M[i][j]*v[j] for j in range(len(v))), F(0)) for i in range(len(v))]
        x = solve(rhs); scale = max(map(abs, x))
        if not scale:
            raise ArithmeticError('Zero inverse iterate')
        v = [F(round(a/scale*denominator), denominator) for a in x]
    return {'basis': [list(t) for t in terms], 'coefficients': list(map(str, v)),
            'precision_bits': precision_bits, 'iterations': iterations,
            'search_shift': str(shift), 'method': 'shifted_rational_inverse_iteration'}


def _gap(q, m, mean_zero, proof):
    if proof is None:
        lower = (m+1)*(m+2)+min(F(q['offset']), F(q['offset'])+F(q['amplitude']))
        return {'format': 'step_pointwise_sector_gap_v1', 'function': q,
                'azimuth_m': m, 'mean_zero': mean_zero, 'eigenvalue_index': 2,
                'real_component': 'one_real_azimuth_component', 'lower': str(lower)}
    if isinstance(proof, dict) and proof.get('format') == 'step_pointwise_sector_gap_v1':
        if not _same(proof, _gap(q, m, mean_zero, None)):
            raise ValueError('Invalid pointwise gap')
        return deepcopy(proof)
    import spectral_gap
    if not spectral_gap.verify_gap(proof, expected_function=q, expected_m=m, expected_mean_zero=mean_zero):
        raise ValueError('Invalid complete second-sector gap proof')
    return deepcopy(proof)


def _certificate(q, terms, v, source, gap, tolerance, mean_zero, independent):
    if type(mean_zero) is not bool:
        raise ValueError('mean_zero must be boolean')
    m = 1 if mean_zero else 0
    tol = _r(tolerance)
    if not F(1, 10**30) <= tol <= 1:
        raise ValueError('Tolerance outside 1e-30..1')
    if not direct.verify_full(source, q, mean_zero):
        raise ValueError('Source must cover the original complete sphere problem')
    if len(source['sectors']) <= m:
        raise ValueError('Source must explicitly cover the trial sector')
    stats = _direct_statistics(q, terms, v, m) if independent else trial_statistics(q, terms, v, m)
    proof = _gap(q, m, mean_zero, gap)
    beta = F(proof['lower']); mu = F(stats['rayleigh']); var = F(stats['residual_squared'])
    if mu >= beta:
        raise ValueError('Temple requires strict mu < beta')
    temple = mu-var/(beta-mu)
    sector_lower = max(temple, F(source['sectors'][m]['lower']))
    others = [{'azimuth_m': c['azimuth_m'], 'lower': c['lower']} for c in source['sectors'] if c['azimuth_m'] != m]
    lower = min([sector_lower, F(source['angular_tail_lower'])]+[F(c['lower']) for c in others])
    pointwise_global = (2 if mean_zero else 0)+min(F(q['offset']), F(q['offset'])+F(q['amplitude']))
    lower = max(lower, pointwise_global)
    upper = min(mu, F(source['upper']))
    if lower > upper:
        raise ArithmeticError('Contradictory global bounds')
    return {'format': FORMAT, 'function': q, 'geometry': 'unit_S2', 'measure': source['measure'],
            'scope': source['scope'], 'mean_zero': mean_zero, 'eigenvalue_index': 1,
            'trial': {'azimuth_m': m, 'real_component': 'cos_phi' if m else 'axisymmetric',
                      'basis': [list(t) for t in terms], 'coefficients': list(map(str, v)),
                      'strong_operator_domain': 'piecewise_polynomial_C1_at_equator',
                      'gap_uses_complete_radial_component': True},
            'statistics': stats, 'gap': proof, 'temple_lower': str(temple),
            'selected_sector_lower': str(sector_lower), 'other_sector_lowers': others,
            'pointwise_global_lower': str(pointwise_global),
            'angular_tail_lower': source['angular_tail_lower'],
            'source': deepcopy(source), 'source_digest': _digest(source),
            'lower': str(lower), 'upper': str(upper), 'exact_width': str(upper-lower),
            'tolerance': str(tol), 'status': 'target_met' if upper-lower <= tol else 'certified_open',
            'full_strong_residual_not_projected': True, 'full_infinite_space_covered': True,
            'function_approximation_error': '0', 'formal_proof_assistant_checked': False}


def trial_certificate(function, terms, coefficients, source, gap=None,
                      tolerance='1/100000000', mean_zero=True):
    q = normalize(function); terms = _basis(terms); v = _vector(coefficients, terms)
    return _certificate(q, terms, v, source, gap, tolerance, mean_zero, False)


def verify(cert, expected_function=None, expected_source=None, expected_tolerance=None,
           expected_mean_zero=None):
    try:
        if not isinstance(cert, dict) or cert.get('format') != FORMAT:
            return False
        q = normalize(cert['function']); mean_zero = cert['mean_zero']
        if expected_function is not None and normalize(expected_function) != q: return False
        if expected_source is not None and not _same(expected_source, cert['source']): return False
        if expected_tolerance is not None and _r(expected_tolerance) != F(cert['tolerance']): return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool or expected_mean_zero is not mean_zero): return False
        trial = cert['trial']; terms = _basis(trial['basis']); v = _vector(trial['coefficients'], terms)
        expected = _certificate(q, terms, v, cert['source'], cert['gap'], cert['tolerance'], mean_zero, True)
        return _same(cert, expected)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, ImportError):
        return False


def solve(function, source=None, gap=None, degree=10, tolerance='1/100000000',
          mean_zero=True, precision_bits=96, iterations=12):
    """One bounded solve, including final independent verification.

    Also accepts a task {'function': ..., 'mean_zero': ..., 'tolerance': ...}.
    A valid wide interval is returned as certified_open; it is never hidden.
    """
    started = perf_counter()
    if isinstance(function, dict) and 'function' in function:
        task = function
        if task.get('kind', 'spectrum') != 'spectrum':
            raise ValueError('This entry point certifies the lowest spectral value')
        function = task['function']; mean_zero = task.get('mean_zero', mean_zero)
        tolerance = task.get('tolerance', tolerance)
    q = normalize(function)
    if type(mean_zero) is not bool: raise ValueError('mean_zero must be boolean')
    m = 1 if mean_zero else 0
    if source is None:
        source = direct.full_ground(q, mean_zero, modes=4, bits=20, max_m=2,
                                    tolerance=tolerance, near_tail=4)
        if len(source['sectors']) <= m:
            # A constant potential may close the angular comparison after m=0;
            # the Temple trial still needs the explicitly bound m=1 component.
            sectors = source['sectors']+[direct.certify_sector(q, m, mean_zero, 4, 20, 40, 4)]
            source = direct.full_certificate(q, mean_zero, sectors, tolerance)
    proposal = propose_trial(q, degree, m, precision_bits, iterations)
    cert = trial_certificate(q, proposal['basis'], proposal['coefficients'], source,
                             gap, tolerance, mean_zero)
    if not verify(cert, q, source, tolerance, mean_zero):
        raise ArithmeticError('Final independent replay failed')
    return {'format': RUN_FORMAT, 'certificate': cert, 'proposal': proposal,
            'elapsed_seconds_including_verify': perf_counter()-started,
            'status': cert['status'], 'degree': degree,
            'development_only_not_holdout': True}
