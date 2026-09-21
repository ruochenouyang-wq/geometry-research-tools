"""Exact full-space m=0 fractional trials and replayable Temple certificates.

Candidate search is numerical research, never acceptance evidence.  Acceptance
rebuilds exact weak/strong moments, the second-eigenvalue gap, and the form bound.
The full-space reduction uses Fourier-sector domination and requires mean_zero
to be False; a mean-zero task is a different mathematical problem.
"""
from copy import deepcopy
from fractions import Fraction as F
from heapq import heappop, heappush
from pathlib import Path
from time import perf_counter
import hashlib
import json
import math
import sys

sys.dont_write_bytecode = True
PREVIOUS = Path(__file__).resolve().parents[1] / 'geometry_strategy_v238_v317'
if str(PREVIOUS) not in sys.path:
    sys.path.insert(0, str(PREVIOUS))
import backend

direct = backend.direct
FORMAT = 'fullspace_m0_fractional_temple_v1'
FORM_FORMAT = 'fullspace_negative_power_form_v1'
DOMINATION = {
    'method': 'axisymmetric_Fourier_sector_form_domination',
    'm0_contains_each_other_sector_radial_form_domain': True,
    'additional_energy': 'm^2*integral(|f(z)|^2/(1-z^2))',
    'additional_energy_nonnegative': True,
    'conclusion': 'full_space_spectral_infimum_equals_m0_spectral_infimum',
    'requires_mean_zero_false': True,
}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def normalize(function):
    q = direct.normalize(function)
    if q['profile'] != 'abs_power' or not -F(1, 2) < F(q['exponent']) < 0:
        raise ValueError('Require an abs_power potential with -1/2 < alpha < 0')
    if F(q['amplitude']) >= 0:
        raise ValueError('This route requires a strictly negative amplitude')
    return q


def form_certificate(function, method='auto', clip_base=2, sqrt_bits=40):
    """A complete form bound, with an exact clipped-potential fallback.

    For alpha=-r/s, delta=base**(-s) makes every clipping power rational.
    The choice is bounded and deterministic, not a spectral parameter search.
    """
    q = normalize(function)
    direct.integer(sqrt_bits, 8, 160, 'sqrt_bits')
    if method not in ('auto', 'centered_L2', 'clipped_negative_L2'):
        raise ValueError('Unknown form-bound method')
    a, alpha, b = F(q['amplitude']), F(q['exponent']), F(q['offset'])
    center = b+a/(1+alpha)
    variance = a*a*alpha*alpha/((1+2*alpha)*(1+alpha)**2)
    eta = direct.f.sqrt_upper(variance, sqrt_bits)
    if method == 'centered_L2' or (method == 'auto' and eta < 1):
        if eta >= 1:
            raise ValueError('Centered L2 comparison has no positive energy factor')
        proof = {'center': str(center), 'variance': str(variance),
                 'eta_upper': str(eta)}
        name, energy, mass = 'centered_L2', 1-eta, center-eta
    else:
        direct.integer(clip_base, 2, 16, 'clip_base')
        if alpha.denominator > 64:
            raise ValueError('Clipped exact-power denominator exceeds 64')
        r, s = -alpha.numerator, alpha.denominator
        bases = (2, 4, 8, 16) if method == 'auto' else (clip_base,)
        for base in bases:
            delta = F(1, base**s)
            cap = -a*base**r
            residual_squared = a*a*F(1, base**(s-2*r))*(
                1/(1+2*alpha)-2/(1+alpha)+1)
            residual_eta = direct.f.sqrt_upper(residual_squared, sqrt_bits)
            if residual_eta < 1:
                break
        else:
            raise ValueError('Bounded clipping choices did not give a positive energy factor')
        eta = residual_eta
        proof = {'clip_base': base, 'delta': str(delta), 'cap': str(cap),
                 'residual_definition': '(abs(amplitude)*abs(z)^alpha-cap)_+',
                 'residual_L2_squared': str(residual_squared), 'eta_upper': str(eta)}
        name, energy, mass = 'clipped_negative_L2', 1-eta, b-cap-eta
    return {'format': FORM_FORMAT, 'function': q, 'method': name,
            'geometry': 'unit_S2', 'measure': 'd_sigma/(4*pi)',
            'scope': 'all_real_H1_on_unit_S2', 'sqrt_bits': sqrt_bits,
            'embedding': 'norm_L4_squared <= Dirichlet_energy + L2_mass',
            'embedding_source': 'https://arxiv.org/pdf/1210.1853; equation (1), d=2,p=4',
            'proof': proof, 'form': {'energy_factor': str(energy),
                                      'mass_offset': str(mass)}}


def verify_form(cert, expected_function=None):
    try:
        import power_forms
        if isinstance(cert, dict) and cert.get('format') in (power_forms.FORMAT, power_forms.SELECTED_FORMAT):
            return power_forms.verify_form(cert, expected_function)
        if not isinstance(cert, dict) or cert.get('format') != FORM_FORMAT:
            return False
        if expected_function is not None and normalize(expected_function) != cert['function']:
            return False
        rebuilt = form_certificate(cert['function'], cert['method'],
                                   cert['proof'].get('clip_base', 2), cert['sqrt_bits'])
        return canonical(rebuilt) == canonical(cert)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, OverflowError):
        return False


def powers(raw):
    if not isinstance(raw, (list, tuple)) or not 1 <= len(raw) <= 16:
        raise ValueError('Require 1..16 distinct powers')
    ss = tuple(direct.f.rational(x) for x in raw)
    if tuple(sorted(set(ss))) != ss or ss[0] != 0:
        raise ValueError('Powers must be strictly increasing and include the constant')
    if any(s != 0 and s <= F(3, 2) for s in ss):
        raise ValueError('Strong domain requires each power to be zero or >3/2')
    if any(s > 64 or s.denominator > 1000000 for s in ss):
        raise ValueError('Power range exceeds the finite diagnostic budget')
    return ss


def generated_powers(function, count=10):
    q = normalize(function)
    direct.integer(count, 1, 16, 'count')
    increments = (F(2), F(2)+F(q['exponent']))
    heap, seen, result = [F(0)], {F(0)}, []
    while len(result) < count:
        s = heappop(heap)
        result.append(s)
        for step in increments:
            t = s+step
            if t not in seen:
                seen.add(t)
                heappush(heap, t)
    return powers(result)


def moment(power):
    r = F(power)
    if r <= -1:
        raise ValueError('Moment is not integrable at the equator')
    return F(2)/(r+1)


def operator_terms(function, power):
    q = normalize(function)
    s = F(power)
    if s != 0 and s <= F(3, 2):
        raise ValueError('Strong domain requires s=0 or s>3/2')
    terms = {}
    for r, coefficient in ((s-2, -s*(s-1)),
                            (s, s*(s+1)+F(q['offset'])),
                            (s+F(q['exponent']), F(q['amplitude']))):
        if coefficient:
            terms[r] = terms.get(r, F(0))+coefficient
    return {r: c for r, c in terms.items() if c}


def trial_matrices(function, trial_powers):
    """Rebuild exact mass, weak energy, and complete strong-operator Gram.

    Moments here use dz, not dz/2.  The common factor cancels in all quotients.
    No memoized certificate data is trusted by verification.
    """
    q, ss = normalize(function), powers(trial_powers)
    a, alpha, b = F(q['amplitude']), F(q['exponent']), F(q['offset'])
    actions = [operator_terms(q, s) for s in ss]
    M, H, R = [], [], []
    for s, action_s in zip(ss, actions):
        mass_row, energy_row, residual_row = [], [], []
        for t, action_t in zip(ss, actions):
            mass = moment(s+t)
            kinetic = F(0) if not s*t else s*t*(moment(s+t-2)-mass)
            weak = kinetic+a*moment(s+t+alpha)+b*mass
            strong = sum((c*moment(s+r) for r, c in action_t.items()), F(0))
            if weak != strong:
                raise ArithmeticError('Weak and strong operator forms disagree')
            mass_row.append(mass)
            energy_row.append(weak)
            residual_row.append(sum((c*d*moment(r+v) for r, c in action_s.items()
                                     for v, d in action_t.items()), F(0)))
        M.append(mass_row); H.append(energy_row); R.append(residual_row)
    if any(H[i][j] != H[j][i] or R[i][j] != R[j][i]
           for i in range(len(ss)) for j in range(len(ss))):
        raise ArithmeticError('Exact operator matrices must be symmetric')
    return M, H, R


def matvec(matrix, vector):
    return [sum((a*b for a, b in zip(row, vector)), F(0)) for row in matrix]


def quadratic(matrix, vector):
    return sum((a*b for a, b in zip(vector, matvec(matrix, vector))), F(0))


def congruence(matrix, transform):
    n, k = len(matrix), len(transform[0])
    ac = [[sum((matrix[i][r]*transform[r][j] for r in range(n)), F(0))
           for j in range(k)] for i in range(n)]
    return [[sum((transform[r][i]*ac[r][j] for r in range(n)), F(0))
             for j in range(k)] for i in range(k)]


def statistics(function, trial_powers, coefficients, matrices=None):
    ss = powers(trial_powers)
    if not isinstance(coefficients, (list, tuple)) or len(coefficients) != len(ss):
        raise ValueError('One rational coefficient per power is required')
    v = tuple(direct.f.rational(x) for x in coefficients)
    if not any(v):
        raise ValueError('A nonzero trial is required')
    M, H, R = trial_matrices(function, ss) if matrices is None else matrices
    mass = quadratic(M, v)
    if mass <= 0:
        raise ValueError('Trial mass is not positive')
    energy, norm = quadratic(H, v), quadratic(R, v)
    mu = energy/mass
    residual = norm/mass-mu*mu
    if residual < 0:
        raise ArithmeticError('Complete residual variance is negative')
    return {'mass': str(mass), 'operator_form': str(energy),
            'operator_norm_squared': str(norm), 'rayleigh': str(mu),
            'residual_squared': str(residual)}


def cancellation_transform(function, trial_powers):
    q, ss = normalize(function), powers(trial_powers)
    delta = 2+F(q['exponent'])
    if delta not in ss:
        return [[F(i == j) for j in range(len(ss))] for i in range(len(ss))]
    eliminated, constant = ss.index(delta), ss.index(F(0))
    free = [i for i in range(len(ss)) if i != eliminated]
    C = [[F(i == j) for j in free] for i in range(len(ss))]
    C[eliminated][free.index(constant)] = F(q['amplitude'])/(delta*(delta-1))
    return C


def proposal(function, trial_powers, form, beta, precision_bits=112,
             iterations=24, residual_steps=2):
    q, ss = normalize(function), powers(trial_powers)
    direct.integer(precision_bits, 32, 256, 'precision_bits')
    direct.integer(iterations, 1, 64, 'iterations')
    direct.integer(residual_steps, 0, 4, 'residual_steps')
    if not verify_form(form, q):
        raise ValueError('Proposal shift requires a verified matching form bound')
    M, H, R = trial_matrices(q, ss)
    C = cancellation_transform(q, ss)
    mr, hr, rr = (congruence(matrix, C) for matrix in (M, H, R))
    shift = F(form['form']['mass_offset'])-F(1, 1024)
    inverse = backend.enriched._ldl_solver(
        [[h-shift*mr[i][j] for j, h in enumerate(row)] for i, row in enumerate(hr)])
    v = [F(i == 0) for i in range(len(mr))]
    denominator = 2**precision_bits
    def rounded(x):
        scale = max(abs(t) for t in x)
        if not scale:
            raise ArithmeticError('Zero inverse proposal')
        return [F(round(t/scale*denominator), denominator) for t in x]
    for _ in range(iterations):
        v = rounded(inverse(matvec(mr, v)))
    trace = []
    beta = F(beta)
    for _ in range(residual_steps):
        before = statistics(q, ss, matvec(C, v), (M, H, R))
        mu, variance = F(before['rayleigh']), F(before['residual_squared'])
        kappa = max(F(1, 2**24), variance)
        residual_matrix = [[rr[i][j]-2*mu*hr[i][j]+(mu*mu+kappa)*mr[i][j]
                            for j in range(len(mr))] for i in range(len(mr))]
        candidate = rounded(backend.enriched._ldl_solver(residual_matrix)(matvec(mr, v)))
        after = statistics(q, ss, matvec(C, candidate), (M, H, R))
        amu, avar = F(after['rayleigh']), F(after['residual_squared'])
        accepted = amu < beta and (mu >= beta or avar/(beta-amu) < variance/(beta-mu))
        trace.append({'before_residual_squared': str(variance),
                      'after_residual_squared': str(avar), 'accepted': accepted})
        if not accepted:
            break
        v = candidate
    return {'powers': list(map(str, ss)), 'coefficients': list(map(str, matvec(C, v))),
            'method': 'shifted_rational_inverse_iteration_then_complete_residual_filter',
            'leading_singularity_cancelled': 2+F(q['exponent']) in ss,
            'precision_bits': precision_bits, 'iterations': iterations,
            'inverse_shift': str(shift), 'residual_refinement': trace}


def certificate(function, trial_powers, coefficients, gap, tolerance='1/100000000',
                mean_zero=False):
    import spectral_gap
    q, ss = normalize(function), powers(trial_powers)
    if mean_zero is not False:
        raise ValueError('Full-space certificate requires mean_zero=False')
    tol = direct.f.rational(tolerance)
    if not F(1, 10**30) <= tol <= 1:
        raise ValueError('Tolerance must lie in 1e-30..1')
    if not spectral_gap.verify_gap(gap, expected_function=q, expected_m=0,
                                   expected_mean_zero=False, form_verifier=verify_form):
        raise ValueError('Require a replayable full m=0 second-eigenvalue bound')
    matrices = trial_matrices(q, ss)
    stats = statistics(q, ss, coefficients, matrices)
    mu, residual, beta = F(stats['rayleigh']), F(stats['residual_squared']), F(gap['beta'])
    if mu >= beta:
        raise ValueError('Temple is not applicable: Rayleigh quotient >= second-eigenvalue bound')
    lower = mu-residual/(beta-mu)
    encoded = [[[str(x) for x in row] for row in matrix] for matrix in matrices]
    return {'format': FORMAT, 'function': q, 'geometry': 'unit_S2',
            'measure': 'd_sigma/(4*pi)', 'scope': 'all_real_H1_on_unit_S2',
            'mean_zero': False, 'eigenvalue_index': 1,
            'trial': {'azimuth_m': 0, 'radial_parity': 'even',
                      'formula': 'sum(a_s*abs(z)^s)', 'powers': list(map(str, ss)),
                      'coefficients': list(map(str, map(direct.f.rational, coefficients))),
                      'strong_domain': 's=0_or_s>3/2', 'matrix_measure': 'dz'},
            'matrix_digest': digest(encoded), 'statistics': stats,
            'gap': deepcopy(gap), 'gap_digest': digest(gap),
            'full_space_domination': deepcopy(DOMINATION),
            'lower': str(lower), 'upper': str(mu), 'exact_width': str(mu-lower),
            'tolerance': str(tol), 'status': 'target_met' if mu-lower <= tol else 'certified_open',
            'full_infinite_space_covered': True, 'original_function_moments_exact': True,
            'function_approximation_error': '0', 'formal_proof_assistant_checked': False}


def verify(cert, expected_function=None, expected_mean_zero=None, expected_tolerance=None):
    try:
        if not isinstance(cert, dict) or cert.get('format') != FORMAT:
            return False
        if expected_mean_zero is not None and expected_mean_zero is not False:
            return False
        if cert.get('mean_zero') is not False:
            return False
        if expected_function is not None and normalize(expected_function) != cert['function']:
            return False
        if expected_tolerance is not None and F(expected_tolerance) != F(cert['tolerance']):
            return False
        trial = cert['trial']
        rebuilt = certificate(cert['function'], trial['powers'], trial['coefficients'],
                              cert['gap'], cert['tolerance'], False)
        return canonical(rebuilt) == canonical(cert)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, OverflowError):
        return False


def solve(function, tolerance='1/100000000', mean_zero=False, max_terms=16,
          gap_modes=8, gap_near_tail=16, gap_bits=20, iterations=24,
          precision_bits=112, wall_seconds=8.0):
    """One bounded stage schedule; all attempted widths/failures are retained.

    wall_seconds is cooperative.  The evaluator must enforce the outer hard
    deadline, including the independent verification following this call.
    """
    import spectral_gap
    import power_forms
    started = perf_counter()
    q = normalize(function)
    if mean_zero is not False:
        raise ValueError('This solver must not change a mean-zero task into full space')
    direct.integer(max_terms, 1, 16, 'max_terms')
    direct.integer(gap_modes, 1, 32, 'gap_modes')
    tolerance = str(direct.f.rational(tolerance))
    if not F(1, 10**30) <= F(tolerance) <= 1:
        raise ValueError('Tolerance must lie in 1e-30..1')
    if isinstance(wall_seconds, bool) or not math.isfinite(float(wall_seconds)) or not 0 <= float(wall_seconds) <= 60:
        raise ValueError('wall_seconds must be finite and lie in 0..60')
    wall_seconds = float(wall_seconds)
    try:
        centered = form_certificate(q, method='centered_L2')
    except ValueError:
        centered = None
    gap_stages = [gap_modes] if centered is not None else sorted(set((max(12, gap_modes), max(16, gap_modes))))
    stages = sorted(set([n for n in (3, 6, 10, 14, 16) if n <= max_terms]+[max_terms]))
    attempts, gap_attempts, best, form = [], [], None, None
    reason = 'stage_budget'
    for gap_stage, modes in enumerate(gap_stages):
        if perf_counter()-started >= wall_seconds:
            reason = 'wall_budget'
            break
        begin = perf_counter()
        gap_row = {'modes': modes, 'near_tail': gap_near_tail, 'verified': False}
        try:
            form = centered if centered is not None else power_forms.select_form(q, modes)
            gap_row['form_certificate'] = deepcopy(form)
            gap = spectral_gap.certified_gap(q, m=0, mean_zero=False, modes=modes,
                                            near_tail=gap_near_tail, bits=gap_bits,
                                            form_certificate=form, form_verifier=verify_form)
            valid_gap = spectral_gap.verify_gap(gap, q, 0, False, form_verifier=verify_form)
            if not valid_gap:
                raise ArithmeticError('Independent complete gap replay failed')
            gap_row.update(status='gap_verified', verified=True, beta=gap['beta'],
                           initial_tail_lower=gap['radial_tail_lower'])
        except (ValueError, TypeError, ArithmeticError, KeyError) as error:
            gap_row.update(status='gap_failed', error=type(error).__name__, reason=str(error),
                           form_proposals=deepcopy(getattr(error, 'attempts', [])))
        gap_row['elapsed_seconds'] = perf_counter()-begin
        gap_attempts.append(gap_row)
        if not gap_row['verified']:
            continue
        # A larger gap stage revisits only the declared largest candidate,
        # rather than silently repeating every already-recorded small basis.
        counts = stages if gap_stage == 0 else [max_terms]
        for count in counts:
            if perf_counter()-started >= wall_seconds:
                reason = 'wall_budget'
                break
            begin = perf_counter()
            row = {'terms': count, 'gap_modes': modes, 'verified': False}
            try:
                candidate = proposal(q, generated_powers(q, count), form, gap['beta'],
                                     precision_bits, iterations)
                cert = certificate(q, candidate['powers'], candidate['coefficients'], gap, tolerance)
                valid = verify(cert, q, False, tolerance)
                row.update(candidate=candidate, verified=valid, lower=cert['lower'],
                           upper=cert['upper'], exact_width=cert['exact_width'], status=cert['status'])
                if not valid:
                    raise ArithmeticError('Independent exact certificate replay failed')
                if best is None or F(cert['exact_width']) < F(best['exact_width']):
                    best = cert
                if best['status'] == 'target_met':
                    reason = 'target_met'
            except (ValueError, TypeError, ArithmeticError, KeyError) as error:
                row.update(status='attempt_failed', error=type(error).__name__, reason=str(error))
            row['elapsed_seconds'] = perf_counter()-begin
            attempts.append(row)
            if reason == 'target_met':
                break
        if reason in ('target_met', 'wall_budget'):
            break
    elapsed = perf_counter()-started
    return {'status': best['status'] if best is not None else 'unsupported',
            'target_met': best is not None and best['status'] == 'target_met',
            'certificate': best, 'attempts': attempts, 'gap_attempts': gap_attempts,
            'stop_reason': reason,
            'form_certificate': best['gap']['form_certificate'] if best is not None else form,
            'elapsed_seconds': elapsed,
            'within_budget': elapsed <= wall_seconds,
            'budget': {'max_terms': max_terms, 'gap_modes_stages': gap_stages,
                       'wall_seconds': wall_seconds,
                       'wall_budget_is_cooperative': True},
            'actual_model_calls': 0, 'actual_model_tokens': None}
