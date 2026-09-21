"""One singularity-adapted, exact-residual trial method for -|z|**(-1/4).

The only candidate search is rounded rational inverse iteration.  Acceptance
uses exact moments and Temple's inequality in a *single real* m=1 component.
The independently replayed direct-moment source covers all other components.
"""
from copy import deepcopy
from functools import lru_cache
from time import perf_counter
import json

from bridge import F, same, canonical, digest, previous_functions as f
import direct_moments as direct

FORMAT = 'singular_fractional_trial_temple_v1'
RUN_FORMAT = 'singular_fractional_trial_search_v1'
FUNCTION = {'kind': 'axis_profile', 'axis': 2, 'profile': 'abs_power',
            'exponent': '-1/4', 'amplitude': '-1', 'offset': '0'}
DEFAULT_EXPONENTS = tuple(F(x) for x in (
    '0', '7/4', '2', '7/2', '15/4', '4', '21/4', '11/2',
    '23/4', '6', '7', '29/4', '15/2', '31/4', '8', '35/4'))


def exponents(raw):
    if not isinstance(raw, (list, tuple)) or not 1 <= len(raw) <= 16:
        raise ValueError('Use 1..16 distinct fractional powers')
    ss = tuple(f.rational(s) for s in raw)
    if len(set(ss)) != len(ss) or tuple(sorted(ss)) != ss:
        raise ValueError('Powers must be distinct and increasing')
    if any(s != 0 and s <= F(3, 2) for s in ss):
        raise ValueError('Strong operator domain requires s=0 or s>3/2')
    if any(s > 64 or s.denominator > 1000000 for s in ss):
        raise ValueError('Power range or denominator exceeds budget')
    return ss


def weighted_moment(power):
    """Integral_{-1}^1 (1-z^2)|z|^power dz, not probability measure."""
    r = f.rational(power)
    if r <= -1:
        raise ValueError('Moment diverges at the equator')
    return F(4)/((r+1)*(r+3))


def action_terms(power):
    s = exponents([power])[0]
    terms = {s: (s+1)*(s+2), s-F(1, 4): F(-1)}
    if s*(s-1):
        terms[s-2] = terms.get(s-2, F(0))-s*(s-1)
    return {r: a for r, a in terms.items() if a}


def _inner(a, b):
    return sum((x*y*weighted_moment(r+s)
                for r, x in a.items() for s, y in b.items()), F(0))


@lru_cache(maxsize=32)
def _matrices(ss):
    basis = [{s: F(1)} for s in ss]
    action = [action_terms(s) for s in ss]
    M = tuple(tuple(_inner(a, b) for b in basis) for a in basis)
    H = tuple(tuple(_inner(a, b) for b in action) for a in basis)
    R = tuple(tuple(_inner(a, b) for b in action) for a in action)
    if any(H[i][j] != H[j][i] for i in range(len(ss)) for j in range(len(ss))):
        raise ArithmeticError('Operator integration-by-parts symmetry failed')
    return M, H, R


def trial_matrices(powers):
    """Exact mass, operator form, and complete strong-operator Gram matrix."""
    return deepcopy(_matrices(exponents(powers)))


def _quadratic(A, v):
    return sum((v[i]*A[i][j]*v[j] for i in range(len(v))
                for j in range(len(v))), F(0))


def trial_statistics(powers, coefficients):
    ss = exponents(powers)
    if not isinstance(coefficients, (list, tuple)) or len(coefficients) != len(ss):
        raise ValueError('A rational coefficient is required for every power')
    v = tuple(f.rational(x) for x in coefficients)
    if not any(v):
        raise ValueError('Trial function must be nonzero')
    M, H, R = _matrices(ss)
    mass = _quadratic(M, v)
    if mass <= 0:
        raise ArithmeticError('Nonpositive trial mass')
    energy = _quadratic(H, v)
    operator_norm = _quadratic(R, v)
    mu = energy/mass
    residual = operator_norm/mass-mu*mu
    if residual < 0:
        raise ArithmeticError('Negative complete residual norm')
    return {'mass': str(mass), 'operator_form': str(energy),
            'operator_norm_squared': str(operator_norm),
            'rayleigh': str(mu), 'residual_squared': str(residual)}


def second_eigenvalue_bound():
    """The second eigenvalue within one real m=1 radial component."""
    form = direct.form_bound(FUNCTION, sqrt_bits=40)
    beta = direct.tail_bound(form, 2)
    return {'form_bound': form, 'unperturbed_second_eigenvalue': '6',
            'azimuth_m': 1, 'real_component': 'cos_phi_only',
            'radial_parities': 'both_even_and_odd', 'eigenvalue_index': 2,
            'lower': str(beta),
            'sin_component_is_unitarily_equivalent': True,
            'not_the_second_eigenvalue_of_cos_plus_sin': True}


def _encoded_matrix(A):
    return [[str(x) for x in row] for row in A]


def _validated_source(source):
    if not direct.verify_full(source, FUNCTION, True):
        raise ValueError('Source must replay for the original global mean-zero problem')
    if len(source['sectors']) < 2:
        raise ValueError('Source must explicitly cover m=0 and m=1')
    return source


def trial_certificate(powers, coefficients, source, tolerance='1/100000000'):
    """Certify a supplied rational trial; this function does not search."""
    ss = exponents(powers)
    tol = f.rational(tolerance)
    if not F(1, 10**30) <= tol <= 1:
        raise ValueError('Tolerance must be in 1e-30..1')
    source = _validated_source(source)
    stats = trial_statistics(ss, coefficients)
    v = [str(f.rational(x)) for x in coefficients]
    mu, residual = F(stats['rayleigh']), F(stats['residual_squared'])
    gap = second_eigenvalue_bound()
    beta = F(gap['lower'])
    if mu >= beta:
        raise ValueError('Temple requires Rayleigh quotient below the second-sector bound')
    temple = mu-residual/(beta-mu)
    m1_lower = max(temple, F(source['sectors'][1]['lower']))
    other = [{'azimuth_m': c['azimuth_m'], 'lower': c['lower']}
             for c in source['sectors'] if c['azimuth_m'] != 1]
    angular = F(source['angular_tail_lower'])
    lower = min([m1_lower, angular]+[F(c['lower']) for c in other])
    upper = min(mu, F(source['upper']))
    if lower > upper:
        raise ArithmeticError('Global spectral contradiction')
    M, H, R = _matrices(ss)
    return {'format': FORMAT, 'function': deepcopy(FUNCTION),
            'geometry': 'unit_S2', 'measure': source['measure'],
            'scope': 'all_real_mean_zero_H1_on_unit_S2',
            'mean_zero': True, 'eigenvalue_index': 1,
            'trial': {'azimuth_m': 1, 'real_component': 'cos_phi',
                      'radial_parity': 'even',
                      'formula': 'sqrt(1-z^2)*sum(a_s*abs(z)^s)*cos(phi)',
                      'powers': list(map(str, ss)), 'coefficients': v,
                      'strong_operator_domain': 's=0_or_s>3/2',
                      'azimuth_factor_cancels': True},
            'matrices': {'M': _encoded_matrix(M), 'H': _encoded_matrix(H),
                         'R': _encoded_matrix(R)},
            'statistics': stats, 'second_eigenvalue_proof': gap,
            'temple_lower': str(temple), 'm1_lower': str(m1_lower),
            'other_sector_lowers': other, 'angular_tail_lower': str(angular),
            'source': deepcopy(source), 'source_digest': digest(source),
            'lower': str(lower), 'upper': str(upper), 'exact_width': str(upper-lower),
            'tolerance': str(tol),
            'status': 'target_met' if upper-lower <= tol else 'certified_open',
            'full_infinite_space_covered': True,
            'original_function_moments_exact': True,
            'function_approximation_error': '0',
            'full_strong_residual_not_projected_residual': True,
            'formal_proof_assistant_checked': False}


def verify(cert, expected_source=None, expected_tolerance=None):
    try:
        if not isinstance(cert, dict) or cert.get('format') != FORMAT:
            return False
        if expected_source is not None and not same(expected_source, cert['source']):
            return False
        if expected_tolerance is not None and f.rational(expected_tolerance) != F(cert['tolerance']):
            return False
        trial = cert['trial']
        return same(cert, trial_certificate(trial['powers'], trial['coefficients'],
                                           cert['source'], cert['tolerance']))
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


def _ldl_solver(A):
    n = len(A)
    L = [[F(i == j) for j in range(n)] for i in range(n)]
    D = []
    for i in range(n):
        d = A[i][i]-sum((L[i][k]**2*D[k] for k in range(i)), F(0))
        if d <= 0:
            raise ValueError('Candidate inverse-iteration matrix must be positive definite')
        D.append(d)
        for j in range(i+1, n):
            L[j][i] = (A[j][i]-sum((L[j][k]*L[i][k]*D[k]
                                     for k in range(i)), F(0)))/d

    def solve(rhs):
        y = []
        for i in range(n):
            y.append(rhs[i]-sum((L[i][j]*y[j] for j in range(i)), F(0)))
        x = [F(0)]*n
        for i in range(n-1, -1, -1):
            x[i] = y[i]/D[i]-sum((L[j][i]*x[j] for j in range(i+1, n)), F(0))
        return x
    return solve


def propose_trial(powers, precision_bits=112, iterations=14):
    """Rounded exact inverse iteration; proposals are not spectral proofs."""
    ss = exponents(powers)
    direct.integer(precision_bits, 32, 256, 'precision_bits')
    direct.integer(iterations, 0, 64, 'iterations')
    started = perf_counter()
    M, H, _ = _matrices(ss)
    solve = _ldl_solver(H)
    v = [F(i == 0) for i in range(len(ss))]
    denominator = 2**precision_bits
    for _ in range(iterations):
        rhs = [sum((M[i][j]*v[j] for j in range(len(v))), F(0))
               for i in range(len(v))]
        x = solve(rhs)
        scale = max(abs(z) for z in x)
        if not scale:
            raise ArithmeticError('Zero inverse iteration proposal')
        v = [F(round(z/scale*denominator), denominator) for z in x]
    return {'powers': list(map(str, ss)), 'coefficients': list(map(str, v)),
            'method': 'rational_inverse_iteration_with_dyadic_rounding',
            'precision_bits': precision_bits, 'iterations': iterations,
            'elapsed_seconds': perf_counter()-started}


def enrich(source, max_terms=16, tolerance='1/100000000', precision_bits=112,
           iterations=14):
    """Try one nested fractional-power basis; preserve every attempted width."""
    direct.integer(max_terms, 1, 16, 'max_terms')
    _validated_source(source)
    started = perf_counter()
    counts = sorted(set([n for n in (1, 3, 6, 10, 16) if n <= max_terms]+[max_terms]))
    attempts = []
    best = None
    for n in counts:
        begin = perf_counter()
        proposal = propose_trial(DEFAULT_EXPONENTS[:n], precision_bits, iterations)
        cert = trial_certificate(proposal['powers'], proposal['coefficients'], source, tolerance)
        attempts.append({'terms': n, 'powers': proposal['powers'],
                         'coefficient_bits': precision_bits, 'iterations': iterations,
                         'elapsed_seconds': perf_counter()-begin,
                         'proposal_seconds': proposal['elapsed_seconds'],
                         'rayleigh': cert['statistics']['rayleigh'],
                         'residual_squared': cert['statistics']['residual_squared'],
                         'lower': cert['lower'], 'upper': cert['upper'],
                         'exact_width': cert['exact_width'], 'status': cert['status']})
        if best is None or F(cert['exact_width']) < F(best['exact_width']):
            best = cert
        if best['status'] == 'target_met':
            break
    if not verify(best, expected_source=source, expected_tolerance=tolerance):
        raise ArithmeticError('Final independent certificate replay failed')
    return {'format': RUN_FORMAT, 'certificate': best, 'attempts': attempts,
            'elapsed_seconds_including_final_replay': perf_counter()-started,
            'status': best['status'], 'max_terms': max_terms,
            'precision_goal': str(f.rational(tolerance)),
            'search_is_not_a_global_optimality_claim_over_bases': True}


if __name__ == '__main__':
    from pathlib import Path
    source_path = Path(__file__).resolve().parent/'certificates'/'singular_n4.json'
    result = enrich(json.loads(source_path.read_text()))
    print(json.dumps(result, sort_keys=True, indent=2))
