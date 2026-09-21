"""Mean-zero singular ground-state certificates with a complete radial gap."""
from copy import deepcopy
from time import perf_counter

from spectral_gap import F, canonical, direct, certified_gap, verify_gap, _json_native
from backend import digest
import adaptive_singular as trial

FORMAT = 'singular_complete_gap_temple_v1'


def certificate(function, powers, coefficients, gap, m0, m1_source,
                tolerance='1/100000000'):
    q = trial.normalize(function)
    tol = direct.f.rational(tolerance)
    if not F(1, 10**30) <= tol <= 1:
        raise ValueError('Tolerance must lie in 1e-30..1')
    if not verify_gap(gap, q, 1, True):
        raise ValueError('Require a replayable complete m1 second-eigenvalue lower bound')
    if not direct.verify_sector(m0, q, 0, True):
        raise ValueError('Require a complete mean-zero m0 sector certificate')
    if not direct.verify_sector(m1_source, q, 1, True):
        raise ValueError('Require a complete matching m1 source certificate')
    ss = trial.enriched.exponents(powers)
    v = [direct.f.rational(x) for x in coefficients]
    stats = trial.statistics(q, ss, v)
    mu, variance, beta = map(F, (stats['rayleigh'], stats['residual_squared'], gap['beta']))
    if mu >= beta:
        raise ValueError('Temple is not applicable: mu >= beta')
    temple = mu - variance / (beta - mu)
    m1_lower = max(temple, F(m1_source['lower']))
    angular_lower = m1_lower + 3
    low = min(F(m0['lower']), m1_lower, angular_lower)
    high = min(mu, F(m0['upper']), F(m1_source['upper']))
    if low > high:
        raise ArithmeticError('Inconsistent global bounds')
    matrices = trial.trial_matrices(q, ss)
    return {
        'format': FORMAT, 'function': q, 'geometry': 'unit_S2', 'mean_zero': True,
        'axis_coordinate': q['axis'],
        'axis_isometry': 'coordinate_permutation_preserves_energy_mass_and_mean',
        'scope': 'all_real_mean_zero_H1_on_unit_S2', 'eigenvalue_index': 1,
        'trial_space': 'one_real_cos_m1_even_radial_trial',
        'strong_domain': 's=0_or_s>3/2', 'powers': list(map(str, ss)),
        'coefficients': list(map(str, v)), 'statistics': stats,
        'matrix_digest': digest([[[str(x) for x in row] for row in mat] for mat in matrices]),
        'gap': deepcopy(gap), 'temple_lower': str(temple), 'm1_lower': str(m1_lower),
        'm0_source': deepcopy(m0), 'm1_source': deepcopy(m1_source),
        'angular_coverage': {
            'method': 'fourier_sector_domination', 'first_omitted_m': 2,
            'reference_m': 1, 'reference_lower': str(m1_lower),
            'inequality': 'lambda1(m)>=lambda1(1)+(m*m-1)',
            'form_difference': '(m*m-1)/(1-z*z)>=m*m-1',
            'm_ge_2_domain_embeds_in_m1': True,
            'sine_cosine_unitarily_equivalent': True, 'lower': str(angular_lower),
        },
        'lower': str(low), 'upper': str(high), 'exact_width': str(high-low),
        'tolerance': str(tol), 'status': 'target_met' if high-low <= tol else 'certified_open',
        'full_infinite_space_covered': True, 'full_residual_not_projected_residual': True,
        'function_approximation_error': '0', 'formal_proof_assistant_checked': False,
    }


def verify(cert, expected_function=None, expected_mean_zero=None, expected_tolerance=None):
    try:
        if not _json_native(cert) or type(cert) is not dict or cert.get('format') != FORMAT:
            return False
        if expected_function is not None and trial.normalize(expected_function) != cert['function']:
            return False
        if expected_mean_zero is not None and expected_mean_zero is not True:
            return False
        if expected_tolerance is not None and direct.f.rational(expected_tolerance) != F(cert['tolerance']):
            return False
        rebuilt = certificate(cert['function'], cert['powers'], cert['coefficients'],
                              cert['gap'], cert['m0_source'], cert['m1_source'], cert['tolerance'])
        return canonical(cert) == canonical(rebuilt)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, OverflowError):
        return False


def solve(task):
    """Bounded generic policy; every open attempt remains in the returned trace."""
    if type(task) is not dict or task.get('kind', 'spectrum') != 'spectrum' or task.get('mean_zero') is not True:
        raise ValueError('Require a mean-zero spectral task')
    q = trial.normalize(task['function'])
    tolerance = direct.f.rational(task.get('tolerance', '1/100000000'))
    if not F(1, 10**30) <= tolerance <= 1:
        raise ValueError('Tolerance must lie in 1e-30..1')
    begin = perf_counter()
    # Neither budget nor basis depends on a case id or an amplitude lookup.
    gap = certified_gap(q, 1, True, modes=8, near_tail=12, bits=24)
    m0 = direct.certify_sector(q, 0, True, 8, 24, 40, 12)
    m1 = direct.certify_sector(q, 1, True, 8, 24, 40, 12)
    attempts, best = [], None
    shift = F(m1['lower']) - F(1, 1024)
    for n in (3, 6, 10, 16):
        started = perf_counter()
        ss = trial.generated_powers(q, n)
        proposal = trial._inverse_proposal(q, ss, 112, 14, shift=shift)
        try:
            current = certificate(q, ss, proposal['coefficients'], gap, m0, m1, tolerance)
        except ValueError as error:
            if str(error) != 'Temple is not applicable: mu >= beta':
                raise
            attempts.append({'terms': n, 'powers': list(map(str, ss)),
                             'status': 'temple_gap_unavailable', 'reason': str(error),
                             'elapsed_seconds': perf_counter()-started})
            continue
        attempts.append({'terms': n, 'powers': list(map(str, ss)),
                         'status': current['status'], 'statistics': current['statistics'],
                         'exact_width': current['exact_width'],
                         'elapsed_seconds': perf_counter()-started})
        if best is None or F(current['exact_width']) < F(best['exact_width']):
            best = current
        if best['status'] == 'target_met':
            break
    if best is None:
        return {'status': 'gap_unavailable', 'certificate': None, 'attempts': attempts,
                'gap': gap, 'm0_source': m0, 'm1_source': m1,
                'elapsed_seconds': perf_counter()-begin}
    if not verify(best, q, True, tolerance):
        raise ArithmeticError('Final complete singular certificate failed to replay')
    return {'status': best['status'], 'certificate': best, 'attempts': attempts,
            'elapsed_seconds': perf_counter()-begin,
            'policy': 'complete_radial_gap_and_nested_singular_trials',
            'parameter_scope': 'one_exact_parameter_point', 'human_designed_policy': True}
