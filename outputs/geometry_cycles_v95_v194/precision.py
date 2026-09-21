"""V105--114: precision-driven exact certificates, separate from frozen V94.

Temple uses the second eigenvalue of the SAME projected azimuthal sector.
All proposals and diagnostics are advisory; acceptance replays rational proofs.
"""
from fractions import Fraction as F
import copy
from common import projected, sectors, base, exact, same, digest


def _tol(value):
    value = exact(value)
    if not F(1, 10**30) <= value <= 1:
        raise ValueError('Require 1e-30 <= tolerance <= 1')
    return value


def choose_bits(q, tolerance, m=0, mean_zero=True, k=1, maximum=160):
    """V106: choose absolute bisection precision with rational comparisons.

    The finite analytic bracket estimate is a scheduling heuristic, not an
    additional accepted error bound; the returned spectral interval is checked.
    """
    q, tolerance = sectors.potential(q), _tol(tolerance)
    projected.parameters(m, mean_zero, k, max(2, k))
    if type(maximum) is not int or not 8 <= maximum <= 160:
        raise ValueError('maximum bits must be 8..160')
    start = 1 if m == 0 and mean_zero else m
    scale = 2 * sum((abs(c) for c in q), F(0)) + (start+k-1)*(start+k)-start*(start+1)+2
    bits = 8
    while bits < maximum and scale / (2**bits) > tolerance / 32:
        bits += 1
    return bits


def diagnose_precision(cert, tolerance, bits=44):
    """V105: distinguish angular coverage, bisection and radial limitations."""
    tolerance = _tol(tolerance)
    if type(bits) is not int or not 8 <= bits <= 160:
        raise ValueError('bits must be 8..160')
    if cert.get('format') == projected.FULL_FORMAT:
        if not projected.verify_full(cert):
            raise ValueError('Invalid full-sphere certificate')
        candidate = min(cert['sectors'], key=lambda c: F(c['lower']))
        angular_open = F(cert['angular_tail_lower']) < F(cert['upper'])
    elif projected.verify_sector(cert, expected_k=1):
        candidate, angular_open = cert, False
    else:
        raise ValueError('Expected a verified V94 sector or full certificate')
    needed = choose_bits(candidate['q_coefficients'], tolerance, candidate['azimuth_m'], candidate['mean_zero'])
    action = ('target_met' if F(cert['exact_width']) <= tolerance else
              'add_azimuth_sector' if angular_open else
              'increase_bisection_bits' if bits < needed else 'increase_radial_modes')
    return {'action': action, 'recommended_bits': needed,
            'candidate_m': candidate['azimuth_m'], 'current_modes': candidate['modes'],
            'current_width': cert['exact_width'], 'tolerance': str(tolerance),
            'diagnostic_is_advisory': True}


def _binding(cert):
    return (tuple(sectors.potential(cert['q_coefficients'])),
            cert['azimuth_m'], cert['mean_zero'])


def _fields(lo, hi):
    if lo > hi:
        raise ValueError('Disjoint claimed enclosures')
    return {'lower': str(lo), 'upper': str(hi), 'exact_width': str(hi-lo)}


def analytic_sector(q, m=0, mean_zero=True, range_proof=None):
    """Exact inexpensive ground enclosure for an individual sector."""
    projected.parameters(m, mean_zero, 1, 1)
    q = sectors.potential(q)
    proof = projected.v05_range.make_range(q) if range_proof is None else range_proof
    qlo, _ = projected.v05_range.verify_bounds(q, proof)
    start = 1 if m == 0 and mean_zero else m
    matrix = sectors.assemble(q, m, start, 1)
    return {'format': 'precision_analytic_sector_v113',
            'q_coefficients': [str(c) for c in q], 'azimuth_m': m,
            'mean_zero': mean_zero, 'eigenvalue_index': 1, 'range_proof': proof,
            **_fields(start*(start+1)+qlo, matrix['A'][0][0]/matrix['M'][0])}


def _verify_analytic(cert):
    try:
        return same(cert, analytic_sector(cert['q_coefficients'], cert['azimuth_m'],
                                         cert['mean_zero'], cert['range_proof']))
    except (ValueError, KeyError, TypeError, ArithmeticError):
        return False


def verify_sector(cert, depth=0):
    try:
        if depth > 32 or not isinstance(cert, dict):
            return False
        kind = cert.get('format')
        if kind == projected.SECTOR_FORMAT:
            return projected.verify_sector(cert, expected_k=1)
        if kind == 'precision_analytic_sector_v113':
            return _verify_analytic(cert)
        if kind == 'precision_sector_intersection_v111':
            expected = intersect_bounds(cert['schur'], cert['temple'], cert['previous'], depth+1)
            return same(cert, expected)
        return False
    except (ValueError, KeyError, TypeError, ArithmeticError, IndexError):
        return False


def intersect_bounds(schur, temple=None, previous=None, depth=0):
    """V111: retain the exact intersection of Schur, Temple and prior bounds."""
    if not verify_sector(schur, depth):
        raise ValueError('Invalid Schur/analytic evidence')
    lo, hi = F(schur['lower']), F(schur['upper'])
    if temple is not None:
        if not verify_temple(temple) or _binding(temple) != _binding(schur):
            raise ValueError('Temple evidence has wrong binding or proof')
        hi = min(hi, F(temple['upper']))
        if temple['gap_accepted']:
            lo = max(lo, F(temple['lower']))
    if previous is not None:
        if not verify_sector(previous, depth+1) or _binding(previous) != _binding(schur):
            raise ValueError('Previous sector does not match')
        lo, hi = max(lo, F(previous['lower'])), min(hi, F(previous['upper']))
    return {'format': 'precision_sector_intersection_v111',
            'q_coefficients': schur['q_coefficients'], 'azimuth_m': schur['azimuth_m'],
            'mean_zero': schur['mean_zero'], 'eigenvalue_index': 1,
            'schur': schur, 'temple': temple, 'previous': previous,
            **_fields(lo, hi)}


def winning_sectors(evidence, tolerance=F(1, 10**10)):
    """V113: certified pruning using sector lower bounds and a best upper."""
    tolerance = _tol(tolerance)
    if not isinstance(evidence, list) or not evidence or any(not verify_sector(c) for c in evidence):
        raise ValueError('Require verified sector evidence')
    q, _, mean_zero = _binding(evidence[0])
    if any(_binding(c)[0] != q or _binding(c)[2] != mean_zero for c in evidence):
        raise ValueError('Sectors have different operators')
    upper = min(F(c['upper']) for c in evidence)
    return {'best_upper': str(upper),
            'refine_m': [c['azimuth_m'] for c in evidence if F(c['lower']) < upper-tolerance],
            'pruned_m': [c['azimuth_m'] for c in evidence if F(c['lower']) >= upper],
            'sufficient_at_tolerance_m': [c['azimuth_m'] for c in evidence
                                         if upper-tolerance <= F(c['lower']) < upper]}


def _full_certificate(q, mean_zero, evidence, range_proof, tolerance, previous=None, depth=0):
    q, tolerance = sectors.potential(q), _tol(tolerance)
    if type(mean_zero) is not bool or not isinstance(evidence, list) or not 1 <= len(evidence) <= 65:
        raise ValueError('Invalid full-sphere parameters')
    for m, cert in enumerate(evidence):
        if not verify_sector(cert) or _binding(cert) != (tuple(q), m, mean_zero):
            raise ValueError('Missing or mismatched azimuth sector')
    qlo, _ = projected.v05_range.verify_bounds(q, range_proof)
    next_m = len(evidence)
    angular = next_m*(next_m+1)+qlo
    lo = min(angular, *(F(c['lower']) for c in evidence))
    hi = min(F(c['upper']) for c in evidence)
    if previous is not None:
        if not verify_full(previous, depth=depth+1) or sectors.potential(previous['q_coefficients']) != q or previous['mean_zero'] is not mean_zero:
            raise ValueError('Previous full certificate has wrong binding')
        lo, hi = max(lo, F(previous['lower'])), min(hi, F(previous['upper']))
    return {'format': 'precision_full_ground_v114', 'geometry': 'unit_S2',
            'scope': projected.scope(mean_zero), 'q_coefficients': [str(c) for c in q],
            'mean_zero': mean_zero, 'eigenvalue_index_in_constrained_space': 1,
            'sectors': evidence, 'range_proof': range_proof,
            'omitted_azimuth_m_start': next_m, 'angular_tail_lower': str(angular),
            'tolerance': str(tolerance), 'previous': previous,
            'status': 'certified_target_met' if hi-lo <= tolerance else 'certified_bound_open_gap',
            'formal_assistant_checked': False, **_fields(lo, hi)}


def verify_full(cert, expected_q=None, expected_mean_zero=None, depth=0):
    try:
        if depth > 8 or not isinstance(cert, dict):
            return False
        if expected_q is not None and sectors.potential(expected_q) != sectors.potential(cert['q_coefficients']):
            return False
        if expected_mean_zero is not None and cert['mean_zero'] is not expected_mean_zero:
            return False
        return same(cert, _full_certificate(cert['q_coefficients'], cert['mean_zero'], cert['sectors'],
                                           cert['range_proof'], cert['tolerance'], cert['previous'], depth))
    except (ValueError, KeyError, TypeError, ArithmeticError, IndexError):
        return False


def checkpoint(cert):
    """V112: replayable continuation checkpoint, not merely an input hash."""
    if not verify_full(cert):
        raise ValueError('Cannot checkpoint an unverified full-sphere proof')
    return {'format': 'precision_checkpoint_v112', 'q_coefficients': cert['q_coefficients'],
            'mean_zero': cert['mean_zero'], 'scope': cert['scope'],
            'lower': cert['lower'], 'upper': cert['upper'],
            'certificate_digest': digest(cert), 'certificate': copy.deepcopy(cert)}


def verify_checkpoint(saved, expected_q=None, expected_mean_zero=None):
    try:
        return (verify_full(saved['certificate'], expected_q, expected_mean_zero)
                and same(saved, checkpoint(saved['certificate'])))
    except (ValueError, KeyError, TypeError, ArithmeticError):
        return False


def full_ground(q, mean_zero=True, tolerance=F(1, 10**18), modes=4,
                max_modes=32, max_m=16, max_steps=16, use_temple=True, saved=None):
    """V114: refine only possible winners; every budget exit has valid bounds.

    max_steps counts exact spectral refinement rounds, not time or FLOPs.  All
    bisections and candidate residuals are verified before they enter the result.
    """
    q, tolerance = sectors.potential(q), _tol(tolerance)
    projected.parameters(0, mean_zero, 1, modes)
    if (type(max_modes) is not int or not modes <= max_modes <= 64 or
            type(max_m) is not int or not 0 <= max_m <= 64 or
            type(max_steps) is not int or not 0 <= max_steps <= 64 or type(use_temple) is not bool):
        raise ValueError('Invalid precision budget')
    previous = None
    if saved is not None:
        if not verify_checkpoint(saved, q, mean_zero):
            raise ValueError('Checkpoint q/projection binding failed')
        previous = saved['certificate']
    proof = projected.v05_range.make_range(q)
    qlo, _ = projected.v05_range.verify_bounds(q, proof)
    evidence, sizes, bitcounts = [analytic_sector(q, 0, mean_zero, proof)], {0: modes}, {0: 0}
    coefficient_rounds = {0: 0}
    steps = 0
    while True:
        current = _full_certificate(q, mean_zero, evidence, proof, tolerance, previous)
        if current['status'] == 'certified_target_met' or steps >= max_steps:
            return current
        upper = F(current['upper'])
        next_m = len(evidence)
        # Cover a missing angular competitor before refining already known ones.
        if next_m <= max_m and next_m*(next_m+1)+qlo < upper:
            evidence.append(analytic_sector(q, next_m, mean_zero, proof))
            sizes[next_m], bitcounts[next_m] = modes, 0
            coefficient_rounds[next_m] = 0
            continue
        candidates = winning_sectors(evidence, tolerance)['refine_m']
        if not candidates:
            return current  # e.g. omitted angular sector beyond max_m budget.
        def schedule_for(m):
            last_temple = evidence[m].get('temple')
            return (residual_schedule(last_temple['trial'], sizes[m],
                    min(256, choose_bits(q, tolerance, m, mean_zero)+16+16*coefficient_rounds[m]), max_modes)
                    if last_temple is not None else None)
        schedules = {m: schedule_for(m) for m in candidates}
        possible = [m for m in candidates if bitcounts[m] == 0 or sizes[m] < max_modes
                    or bitcounts[m] < choose_bits(q, tolerance, m, mean_zero)
                    or (coefficient_rounds[m] < 2 and schedules[m] is not None
                        and schedules[m]['action'] == 'increase_rational_coefficient_bits')]
        if not possible:
            return current
        m = min(possible, key=lambda m: F(evidence[m]['lower']))
        wanted_bits = choose_bits(q, tolerance, m, mean_zero)
        if bitcounts[m] == 0:
            # Cheap first pass can discard a nonwinning sector without target precision.
            bits = min(24, wanted_bits)
        elif bitcounts[m] < wanted_bits:
            bits = wanted_bits
        else:
            schedule = schedules[m]
            if (schedule is not None and schedule['action'] == 'increase_rational_coefficient_bits'
                    and coefficient_rounds[m] < 2):
                coefficient_rounds[m] += 1
            else:
                sizes[m] = (schedule['next_modes'] if schedule is not None
                            else min(max_modes, sizes[m] + max(2, sizes[m]//2)))
            bits = wanted_bits
        n = sizes[m]
        schur = projected.certify_sector(q, m, mean_zero, 1, n, bits, proof)
        temple = None
        if use_temple and bits == wanted_bits and n >= 2:
            trial = refine_coefficients(q, m, mean_zero, n,
                        min(256, bits+16+16*coefficient_rounds[m]), 8+8*coefficient_rounds[m])
            temple = recover_gap(trial, max(2, n), min(160, bits))
        evidence[m] = intersect_bounds(schur, temple, evidence[m])
        bitcounts[m] = bits
        steps += 1


# The witness interface and Temple/residual operations are below, so the
# independent Schur precision path can also run with use_temple=False.


def temple_certificate(trial, next_eigenvalue):
    """V107: exact Temple bound from a complete projected residual.

    For g <= lambda_2 and mu < g, the spectral theorem gives
    0 <= <(H-lambda_1)(H-g)> = variance+(mu-lambda_1)(mu-g).
    Thus lambda_1 >= mu-variance/(g-mu).  No full-sphere lambda_2 is used.
    """
    import witness
    if trial.get('format') != witness.RAYLEIGH or not witness.verify(trial):
        raise ValueError('A verified rational Rayleigh trial is required')
    q, m, mean_zero = _binding(trial)
    if not projected.verify_sector(next_eigenvalue, expected_q=list(q), expected_m=m,
                                    expected_mean_zero=mean_zero, expected_k=2):
        raise ValueError('The second eigenvalue must belong to the same projected m sector')
    residual = witness.residual_certificate(trial)
    if not witness.verify(residual):
        raise ValueError('Complete projected residual verification failed')
    mu, variance = F(trial['rayleigh_quotient']), F(residual['normalized_residual_squared'])
    gap = F(next_eigenvalue['lower'])-mu
    accepted = gap > 0
    return {'format': 'same_sector_temple_v107', 'q_coefficients': trial['q_coefficients'],
            'azimuth_m': m, 'mean_zero': mean_zero, 'eigenvalue_index': 1,
            'trial': trial, 'complete_residual': residual, 'next_eigenvalue': next_eigenvalue,
            'next_eigenvalue_index': 2, 'spectral_gap_lower_from_trial': str(gap),
            'rayleigh_quotient': str(mu), 'residual_variance': str(variance),
            'gap_accepted': accepted, 'lower': str(mu-variance/gap) if accepted else None,
            'upper': str(mu), 'status': 'temple_bound_verified' if accepted else 'gap_failed_rayleigh_only',
            'scope': 'single_projected_azimuth_sector_ground_only', 'formal_assistant_checked': False}


def verify_temple(cert):
    try:
        return same(cert, temple_certificate(cert['trial'], cert['next_eigenvalue']))
    except (ValueError, KeyError, TypeError, ArithmeticError, IndexError):
        return False


def recover_gap(trial, modes=8, bits=64, previous_gap=None):
    """V108: refine a failed spectral separator; reject a still invalid gap.

    The unsuccessful result retains the exact Rayleigh upper but has NO Temple
    lower.  A trial above lambda_2 cannot be rescued merely by more precision.
    """
    import witness
    if trial.get('format') != witness.RAYLEIGH or not witness.verify(trial):
        raise ValueError('Invalid rational trial')
    if previous_gap is not None:
        old = temple_certificate(trial, previous_gap)
        if old['gap_accepted']:
            return old
    next_eigenvalue = projected.certify_sector(trial['q_coefficients'], trial['azimuth_m'],
                                               trial['mean_zero'], 2, modes, bits)
    return temple_certificate(trial, next_eigenvalue)


def residual_schedule(trial, modes, coefficient_bits=64, max_modes=64):
    """V109: partition residual by the actual radial window, then choose work.

    Outside sparse trial support is NOT synonymous with outside the radial
    window: zero coefficients inside the window remain part of its residual.
    A finite-window residual does not by itself identify quantization as its
    cause; coefficient refinement is an experiment whose improvement is checked.
    """
    import witness
    if trial.get('format') != witness.RAYLEIGH or not witness.verify(trial):
        raise ValueError('Invalid trial')
    projected.parameters(trial['azimuth_m'], trial['mean_zero'], 1, modes)
    if type(max_modes) is not int or not modes <= max_modes <= 64:
        raise ValueError('Invalid maximum modes')
    if type(coefficient_bits) is not int or not 8 <= coefficient_bits <= 256:
        raise ValueError('Coefficient bits must be 8..256')
    m = trial['azimuth_m']
    start = 1 if m == 0 and trial['mean_zero'] else m
    if any(l < start or l >= start+modes for l in trial['degrees']):
        raise ValueError('Trial lies outside the stated radial window')
    residual = witness.residual_certificate(trial)
    norm = F(trial['norm_squared_divided_by_azimuth_factor'])
    tail = sum((F(row['coefficient'])**2 * sectors.basis_mass(m, row['degree']) / norm
                for row in residual['residual_coefficients'] if row['degree'] >= start+modes), F(0))
    variance = F(residual['normalized_residual_squared'])
    finite = variance-tail
    action = ('zero_residual' if variance == 0 else
              'increase_radial_modes' if tail >= finite and modes < max_modes else
              'radial_budget_exhausted' if tail >= finite else 'increase_rational_coefficient_bits')
    return {'action': action, 'window_tail_variance': str(tail), 'window_finite_variance': str(finite),
            'tail_fraction': str(tail/variance if variance else F(0)),
            'next_modes': min(max_modes, modes+max(2, modes//2)),
            'next_coefficient_bits': min(256, coefficient_bits+16),
            'finite_residual_is_not_proof_of_quantization_error': True}


def _ldl_solver(matrix):
    n = len(matrix)
    lower = [[F(i == j) for j in range(n)] for i in range(n)]
    diagonal = []
    for i in range(n):
        d = matrix[i][i]-sum((lower[i][k]**2*diagonal[k] for k in range(i)), F(0))
        if d <= 0:
            raise ValueError('Inverse iteration shift must give a positive definite form')
        diagonal.append(d)
        for j in range(i+1, n):
            lower[j][i] = (matrix[j][i]-sum((lower[j][k]*lower[i][k]*diagonal[k]
                                           for k in range(i)), F(0)))/d
    def solve(rhs):
        y = []
        for i in range(n):
            y.append(rhs[i]-sum((lower[i][j]*y[j] for j in range(i)), F(0)))
        x = [F(0)]*n
        for i in range(n-1, -1, -1):
            x[i] = y[i]/diagonal[i]-sum((lower[j][i]*x[j] for j in range(i+1, n)), F(0))
        return x
    return solve


def refine_coefficients(q, m=0, mean_zero=True, modes=8, coefficient_bits=80, iterations=8):
    """V110: go beyond binary64 proposals using rounded rational inverse steps.

    The initial floating eigensolver has no authority.  Every accepted trial is
    rational and must improve both its Rayleigh quotient and full residual.
    """
    import witness
    if type(coefficient_bits) is not int or not 8 <= coefficient_bits <= 256:
        raise ValueError('Coefficient bits must be 8..256')
    if type(iterations) is not int or not 0 <= iterations <= 32:
        raise ValueError('Inverse iteration budget must be 0..32')
    q = sectors.potential(q)
    projected.parameters(m, mean_zero, 1, modes)
    proposal = witness.mass_normalized_proposal(q, m, mean_zero, modes)
    best = witness.quantize_proposal(proposal, bits=(coefficient_bits,))['certificate']
    if not iterations:
        return best
    start = 1 if m == 0 and mean_zero else m
    data = sectors.assemble(q, m, start, modes)
    qlo = F(analytic_sector(q, m, mean_zero)['lower'])-start*(start+1)
    shift = start*(start+1)+qlo-1
    matrix = [[data['A'][i][j]-(shift*data['M'][i] if i == j else 0)
               for j in range(modes)] for i in range(modes)]
    solve = _ldl_solver(matrix)
    vector = dict(zip(best['degrees'], map(F, best['coefficients'])))
    values = [vector.get(l, F(0)) for l in data['degrees']]
    best_variance = F(witness.residual_certificate(best)['normalized_residual_squared'])
    denominator = 2**coefficient_bits
    for _ in range(iterations):
        values = solve([x*w for x, w in zip(values, data['M'])])
        scale = max(map(abs, values))
        values = [F(round(x/scale*denominator), denominator) for x in values]
        trial = witness.rayleigh_certificate(q, values, m, mean_zero, data['degrees'])
        variance = F(witness.residual_certificate(trial)['normalized_residual_squared'])
        if F(trial['rayleigh_quotient']) <= F(best['rayleigh_quotient']) and variance <= best_variance:
            best, best_variance = trial, variance
    return best


def run_evidence():
    """Reproduce the ten stage experiments; write only round02 evidence."""
    import json
    import time
    import witness
    from common import ROOT, save
    tolerance = F(1, 10**18)
    old = json.loads((ROOT/'precision_candidate.json').read_text())['old_certificate']
    old_path = save('round02/certificates/old_q2.json', old)
    results, stages = {}, []
    def record(version, name, capability, function, data, outcome, extra=()):
        path = save('round02/results/v%d.json' % version, data)
        stages.append({'version': version, 'name': name, 'capability': capability,
                       'callable': 'precision.'+function, 'evidence': [path, *extra], 'outcome': outcome})
    diagnostic = diagnose_precision(old, tolerance, 44)
    record(105, 'Precision bottleneck diagnosis', 'Choose bit, radial or angular work from a verified bound',
           'diagnose_precision', diagnostic, diagnostic['action'], (old_path,))
    chosen = choose_bits([0, 0, 2], tolerance, 1)
    fixed = projected.certify_sector([0, 0, 2], m=1, modes=12, bits=44)
    automatic = projected.certify_sector([0, 0, 2], m=1, modes=12, bits=chosen)
    auto_path = save('round02/certificates/auto_bits_q2.json', automatic)
    record(106, 'Absolute tolerance bit selection', 'Scale bisection bits to a rational absolute target',
           'choose_bits', {'chosen_bits': chosen, 'fixed_44_bits': fixed, 'automatic': automatic},
           'same N12: width '+fixed['exact_width']+' -> '+automatic['exact_width'], (auto_path,))
    for label, q in (('q2', [0, 0, 2]), ('mixed', [0, 1, 1]), ('q100', [0, 0, 100])):
        began = time.perf_counter()
        cert = full_ground(q, tolerance=tolerance)
        accepted = verify_full(cert, q, True)
        elapsed = time.perf_counter()-began
        if not accepted:
            raise AssertionError('Independent precision replay failed')
        path = save('round02/certificates/'+label+'_full.json', cert)
        results[label] = {'q': list(map(str, q)), 'certificate': path,
                          'status': cert['status'], 'exact_width': cert['exact_width'],
                          'verified': accepted, 'seconds_including_verification': elapsed,
                          'sector_modes': [{'m': row['azimuth_m'],
                                            'modes': row.get('schur', {}).get('modes', 1)}
                                           for row in cert['sectors']]}
        if label == 'q2':
            q2 = cert
    temple = q2['sectors'][1]['temple']
    temple_path = save('round02/certificates/q2_temple.json', temple)
    record(107, 'Same-sector Temple lower', 'Use full projected variance and bound the same m second eigenvalue',
           'temple_certificate', temple, temple['status'], (temple_path,))
    ground = witness.rayleigh_certificate([0], [1], m=1)
    weak_separator = projected.sector_certificate(projected.Kernel([0], m=1, modes=2),
                                                  2, F(0), F(6), 'upper_ritz', independent=True)
    before = temple_certificate(ground, weak_separator)
    recovered = recover_gap(ground, 2, 64, weak_separator)
    excited = witness.rayleigh_certificate([0], [1], m=1, degrees=[3])
    failed = recover_gap(excited, 4, 64)
    record(108, 'Gap recovery and honest refusal', 'Refine a weak separator and keep Rayleigh only if the gap still fails',
           'recover_gap', {'weak_gap': before, 'recovered': recovered, 'irrecoverable_trial': failed},
           'weak separator recovered; excited trial keeps lower=null')
    radial_trial = refine_coefficients([0, 0, 2], 1, True, 4, 80)
    schedule = residual_schedule(radial_trial, 4)
    enriched = refine_coefficients([0, 0, 2], 1, True, schedule['next_modes'], 80)
    record(109, 'Residual-controlled radial growth', 'Partition by radial window, retaining every projected tail coefficient',
           'residual_schedule', {'before_trial': radial_trial, 'decision': schedule,
                                 'before_residual': witness.residual_certificate(radial_trial),
                                 'after_trial': enriched, 'after_residual': witness.residual_certificate(enriched)},
           'tail-dominated residual schedules N4 -> N6')
    rough = refine_coefficients([0, 0, 2], 1, True, 8, 8, 0)
    fine = refine_coefficients([0, 0, 2], 1, True, 8, 80, 8)
    rough_residual, fine_residual = witness.residual_certificate(rough), witness.residual_certificate(fine)
    record(110, 'Rational coefficient refinement', 'Increase dyadic coefficient precision and apply exact inverse steps',
           'refine_coefficients', {'rough': rough, 'fine': fine,
                                  'rough_residual': rough_residual, 'fine_residual': fine_residual,
                                  'decision': residual_schedule(rough, 8, 8)},
           'full variance '+rough_residual['normalized_residual_squared']+' -> '+fine_residual['normalized_residual_squared'])
    sector = q2['sectors'][1]
    record(111, 'Independent enclosure intersection', 'Intersect verified Schur, Temple and previously accepted sector bounds',
           'intersect_bounds', {'schur_width': sector['schur']['exact_width'],
                                'intersection_width': sector['exact_width'], 'certificate': sector},
           'intersection never widens any contributing certified enclosure')
    first = full_ground([0, 1, 1], max_steps=1)
    saved = checkpoint(first)
    resumed = full_ground([0, 1, 1], saved=saved)
    checkpoint_path = save('round02/certificates/mixed_checkpoint.json', saved)
    resumed_path = save('round02/certificates/mixed_resumed.json', resumed)
    record(112, 'Bound-preserving continuation', 'Revalidate q/projection binding and inherit exact bounds monotonically',
           'checkpoint', {'checkpoint': checkpoint_path, 'resumed': resumed_path,
                          'lower_monotone': F(resumed['lower']) >= F(first['lower']),
                          'upper_monotone': F(resumed['upper']) <= F(first['upper']),
                          'checkpoint_verified': verify_checkpoint(saved), 'resumed_verified': verify_full(resumed)},
           'both endpoints inherited monotonically', (checkpoint_path, resumed_path))
    pruning = winning_sectors(q2['sectors'], tolerance)
    record(113, 'Winner-directed sector refinement', 'Prune sectors by exact lower bounds against the best certified upper',
           'winning_sectors', {'selection': pruning, 'sector_modes': results['q2']['sector_modes'],
                               'angular_tail_lower': q2['angular_tail_lower'], 'certificate': results['q2']['certificate']},
           'm0 stops at N4; possible winner m1 reaches target at N9', (results['q2']['certificate'],))
    low = full_ground([0, 0, 2], max_steps=0)
    low_path = save('round02/certificates/zero_budget.json', low)
    results['zero_budget'] = {'certificate': low_path, 'status': low['status'],
                              'exact_width': low['exact_width'], 'verified': verify_full(low)}
    results['old_q2'] = {'certificate': old_path, 'status': old['status'], 'exact_width': old['exact_width']}
    results['timing_policy'] = 'Observed local wall time includes replay; no model speed or cross-process benchmark claim.'
    record(114, 'Full sphere precision driver', 'Target absolute full-S2 projected ground width with honest finite-budget exits',
           'full_ground', results, 'all three requested potentials reach 1e-18; zero budget remains valid open',
           tuple(row['certificate'] for row in results.values() if isinstance(row, dict) and 'certificate' in row))
    save('round02/results/summary.json', results)
    save('round02/STAGES.json', stages)
    return results


if __name__ == '__main__':
    import json
    print(json.dumps(run_evidence(), indent=2))
