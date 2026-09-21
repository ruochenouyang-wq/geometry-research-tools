"""V95--V104: explicit rational spherical-harmonic violation witnesses.

Real basis: P_l^m(t) cos(m phi), with cos(0 phi)=1 and the standard
Condon--Shortley associated Legendre convention. Only rational certificates
have mathematical authority; floating eigensolvers merely propose trials.
"""
from math import sqrt, isfinite, isqrt
import re
from common import ROOT, F, projected, sectors, base, arithmetic, save, exact, same

SCOPE = 'explicit_real_spherical_harmonic_trial_on_unit_S2'
RAYLEIGH = 'explicit_sector_rayleigh_v95'
RESIDUAL = 'complete_projected_residual_v100'
RITZ = 'two_trial_ritz_v102'
FINAL = 'explicit_sphere_counterexample_v104'


def _strings(v):
    return list(map(str, v))


def _coefficient(value):
    """Exact rational coefficients have no inherited small-denominator cap."""
    if type(value) in (int, F):
        return F(value)
    if isinstance(value, str) and re.fullmatch(r'-?\d+(?:/[1-9]\d*)?', value):
        return F(value)
    raise ValueError('Exact integer or rational trial coefficient required')


def _input(q, coefficients, m, mean_zero, degrees=None):
    q = sectors.potential(q)
    if type(m) is not int or not 0 <= m <= 64 or type(mean_zero) is not bool:
        raise ValueError('Integer m in 0..64 and boolean mean_zero required')
    if not isinstance(coefficients, list) or not 1 <= len(coefficients) <= 256:
        raise ValueError('One to 256 exact trial coefficients required')
    coefficients = [_coefficient(x) for x in coefficients]
    start = 1 if m == 0 and mean_zero else m
    degrees = list(range(start, start+len(coefficients))) if degrees is None else degrees
    if not isinstance(degrees, list) or len(degrees) != len(coefficients):
        raise ValueError('Degree and coefficient dimensions differ')
    if any(type(l) is not int or not m <= l <= 250 for l in degrees) or len(set(degrees)) != len(degrees):
        raise ValueError('Distinct integer harmonic degrees m<=l<=250 required')
    series = {l: x for l, x in zip(degrees, coefficients) if x}
    if not series:
        raise ValueError('Trial must be nonzero')
    if mean_zero and m == 0 and series.get(0):
        raise ValueError('Nonzero l=0 violates mean-zero constraint')
    return q, dict(sorted(series.items()))


def _dot(u, v, m):
    return sum((x*v.get(l, F(0))*sectors.basis_mass(m, l) for l, x in u.items()), F(0))


def _apply(q, u, m, mean_zero):
    acted = {}
    for l, value in u.items():
        acted[l] = acted.get(l, F(0))+l*(l+1)*value
        for k, factor in sectors.q_times_basis(q, m, l).items():
            acted[k] = acted.get(k, F(0))+value*factor
    removed = acted.pop(0, F(0)) if m == 0 and mean_zero else F(0)
    return {l: value for l, value in sorted(acted.items()) if value}, removed


def _series(cert):
    return dict(zip(cert['degrees'], map(F, cert['coefficients'])))


def _encode_series(series):
    return [{'degree': l, 'coefficient': str(x)} for l, x in sorted(series.items()) if x]


def rayleigh_certificate(q, coefficients, m=0, mean_zero=True, degrees=None, threshold=None):
    """V95: exact Rayleigh quotient of any nonzero rational real trial."""
    q, u = _input(q, coefficients, m, mean_zero, degrees)
    h, _ = _apply(q, u, m, mean_zero)
    mass = _dot(u, u, m)
    kinetic = sum((l*(l+1)*x*x*sectors.basis_mass(m, l) for l, x in u.items()), F(0))
    energy = _dot(u, h, m)
    quotient = energy/mass
    threshold = None if threshold is None else exact(threshold)
    return {'format': RAYLEIGH, 'geometry': 'unit_S2', 'scope': SCOPE,
            'q_coefficients': _strings(q), 'azimuth_m': m, 'mean_zero': mean_zero,
            'degrees': list(u), 'coefficients': _strings(u.values()),
            'harmonic_terms': [{'l': l, 'm': m, 'coefficient': str(x), 'real_component': 'cos'} for l, x in u.items()],
            'basis_convention': 'standard_Condon_Shortley_P_l^m(t)*cos(m*phi)',
            'azimuth_integral_factor': '2*pi' if m == 0 else 'pi',
            'norm_squared_divided_by_azimuth_factor': str(mass),
            'kinetic_divided_by_azimuth_factor': str(kinetic),
            'potential_divided_by_azimuth_factor': str(energy-kinetic),
            'energy_divided_by_azimuth_factor': str(energy), 'rayleigh_quotient': str(quotient),
            'full_space_spectral_upper': str(quotient),
            'full_space': projected.scope(mean_zero),
            'threshold': None if threshold is None else str(threshold),
            'strict_margin': None if threshold is None else str(threshold-quotient),
            'status': 'strict_counterexample' if threshold is not None and quotient < threshold else
                      'no_strict_violation' if threshold is not None else 'rayleigh_certified',
            'formal_assistant_checked': False}


def _ray_verify(cert):
    return same(cert, rayleigh_certificate(cert['q_coefficients'], cert['coefficients'],
                cert['azimuth_m'], cert['mean_zero'], cert['degrees'], cert['threshold']))


def _compatible(left, right):
    return all(left[key] == right[key] for key in ['q_coefficients', 'azimuth_m', 'mean_zero'])


def select_sector(full_certificate, threshold=None):
    """V96: use a verified full-spectrum certificate to choose a trial sector."""
    if not projected.verify_full(full_certificate):
        raise ValueError('A complete valid full-sphere certificate is required')
    threshold = None if threshold is None else exact(threshold)
    ranked = sorted(full_certificate['sectors'], key=lambda row: (F(row['upper']), row['azimuth_m']))
    return {'format': 'candidate_sector_selection_v96', 'source': full_certificate,
            'q_coefficients': full_certificate['q_coefficients'], 'mean_zero': full_certificate['mean_zero'],
            'ranked_sectors': [{'m': row['azimuth_m'], 'lower': row['lower'], 'upper': row['upper']} for row in ranked],
            'selected_m': ranked[0]['azimuth_m'],
            'threshold': None if threshold is None else str(threshold),
            'spectral_existence_below_threshold': threshold is not None and F(ranked[0]['upper']) < threshold,
            'limitation': 'Sector selection does not yet exhibit a function.', 'formal_assistant_checked': False}


def _proposal(q, m, mean_zero, degrees):
    q, _ = _input(q, [1]*len(degrees), m, mean_zero, degrees)
    if len(degrees) > 64:
        raise ValueError('Floating proposal dimension at most 64')
    mass = [sectors.basis_mass(m, l) for l in degrees]
    columns = [_apply(q, {l: F(1)}, m, mean_zero)[0] for l in degrees]
    matrix = [[mass[i]*columns[j].get(l, F(0)) for j in range(len(degrees))]
              for i, l in enumerate(degrees)]
    # Convert only dimensionless rational ratios, avoiding huge associated-
    # Legendre normalization constants before the floating eigensolver.
    normalized = [[(1 if matrix[i][j] >= 0 else -1)*sqrt(float(matrix[i][j]**2/(mass[i]*mass[j])))
                   for j in range(len(degrees))] for i in range(len(degrees))]
    vector = arithmetic.jacobi_lowest(normalized)
    length = sqrt(sum(x*x for x in vector))
    vector = [x/length for x in vector]
    reference = min(mass)
    coefficients = [x/sqrt(float(w/reference)) for x, w in zip(vector, mass)]
    scale = max(map(abs, coefficients))
    coefficients = [x/scale for x in coefficients]
    return {'format': 'mass_normalized_proposal_v97', 'q_coefficients': _strings(q),
            'azimuth_m': m, 'mean_zero': mean_zero, 'degrees': degrees,
            'mass_normalized_eigenvector': vector, 'coefficients_float': coefficients,
            'numerical_rayleigh': sum(vector[i]*normalized[i][j]*vector[j]
                                      for i in range(len(vector)) for j in range(len(vector))),
            'authority': 'untrusted proposal; no spectral certificate'}


def mass_normalized_proposal(q, m=0, mean_zero=True, modes=6):
    """V97: generalized finite eigenproblem in unit-mass coordinates."""
    if type(modes) is not int or not 1 <= modes <= 64:
        raise ValueError('Proposal mode budget in 1..64')
    start = 1 if m == 0 and mean_zero else m
    return _proposal(q, m, mean_zero, list(range(start, start+modes)))


def quantize_proposal(proposal, bits=(2, 4, 8, 16, 24), threshold=None):
    """V98: exact acceptance of a denominator ladder, retaining the best trial."""
    if not isinstance(bits, (list, tuple)) or not 1 <= len(bits) <= 32:
        raise ValueError('Quantization ladder size')
    if any(type(b) is not int or not 0 <= b <= 512 for b in bits):
        raise ValueError('Quantization bits in 0..512')
    floating = proposal['coefficients_float']
    if not all(type(x) in (int, float) and isfinite(x) for x in floating):
        raise ValueError('Finite floating proposal required')
    amplitude = max(map(abs, floating), default=0)
    if not amplitude:
        raise ValueError('Zero floating proposal')
    floating = [x/amplitude for x in floating]
    trace, best = [], None
    for bit in bits:
        values = [F(round(x*2**bit), 2**bit) for x in floating]
        if not any(values):
            trace.append({'bits': bit, 'status': 'rejected_zero_trial'})
            continue
        cert = rayleigh_certificate(proposal['q_coefficients'], values, proposal['azimuth_m'],
                        proposal['mean_zero'], proposal['degrees'], threshold)
        accepted = best is None or F(cert['rayleigh_quotient']) < F(best['rayleigh_quotient'])
        if accepted:
            best = cert
        trace.append({'bits': bit, 'status': 'accepted_improvement' if accepted else 'retained_old_trial', 'certificate': cert})
    if best is None:
        raise ValueError('Every quantization vanished')
    return {'format': 'quantization_ladder_v98', 'certificate': best, 'trace': trace}


def sparsify(candidate, keep_counts=(1, 2, 3)):
    """V99: truncate by exact mass contribution; reject any Rayleigh increase."""
    if not verify(candidate):
        raise ValueError('Invalid trial certificate')
    u = _series(candidate)
    ranking = sorted(u, key=lambda l: (-u[l]**2*sectors.basis_mass(candidate['azimuth_m'], l), l))
    best, trace = candidate, []
    for count in keep_counts:
        if type(count) is not int or not 1 <= count <= 256:
            raise ValueError('Positive sparse support size')
        support = sorted(ranking[:count])
        cert = rayleigh_certificate(candidate['q_coefficients'], [u[l] for l in support],
                candidate['azimuth_m'], candidate['mean_zero'], support, candidate['threshold'])
        quotient, old = F(cert['rayleigh_quotient']), F(best['rayleigh_quotient'])
        accepted = quotient < old or (quotient == old and len(support) < len(best['degrees']))
        if accepted:
            best = cert
        trace.append({'keep_terms': count, 'status': 'accepted' if accepted else 'rejected_no_gain',
                      'certificate': cert})
    return {'format': 'sparse_trial_acceptance_v99', 'baseline': candidate, 'certificate': best, 'trace': trace}


def residual_certificate(candidate):
    """V100: full projected residual, including every generated high degree."""
    if not _ray_verify(candidate):
        raise ValueError('Invalid rational trial')
    q, u = sectors.potential(candidate['q_coefficients']), _series(candidate)
    m, mean_zero = candidate['azimuth_m'], candidate['mean_zero']
    h, removed = _apply(q, u, m, mean_zero)
    mu = F(candidate['rayleigh_quotient'])
    residual = dict(h)
    for l, x in u.items():
        residual[l] = residual.get(l, F(0))-mu*x
    residual = {l: x for l, x in sorted(residual.items()) if x}
    norm = _dot(u, u, m)
    residual_norm = _dot(residual, residual, m)
    finite = _dot({l: x for l, x in residual.items() if l in u}, residual, m)
    return {'format': RESIDUAL, 'scope': 'complete_residual_of_the_stated_projected_operator',
            'trial': candidate, 'q_coefficients': candidate['q_coefficients'],
            'azimuth_m': m, 'mean_zero': mean_zero, 'projection': 'remove_l0' if m == 0 and mean_zero else 'none',
            'acted_coefficients': _encode_series(h), 'residual_coefficients': _encode_series(residual),
            'removed_l0_coefficient': str(removed), 'orthogonality_to_trial': str(_dot(u, residual, m)),
            'residual_norm_squared_divided_by_azimuth_factor': str(residual_norm),
            'normalized_residual_squared': str(residual_norm/norm),
            'in_support_residual_squared': str(finite/norm),
            'outside_support_residual_squared': str((residual_norm-finite)/norm),
            'limitation': 'Residual norm alone gives no ground-energy lower bound without spectral separation.',
            'formal_assistant_checked': False}


def residual_enrich(candidate, max_new_degrees=6, bits=(8, 16, 24)):
    """V101: grow support in degrees detected by the COMPLETE residual."""
    residual = residual_certificate(candidate)
    if type(max_new_degrees) is not int or not 1 <= max_new_degrees <= 32:
        raise ValueError('Support growth budget in 1..32')
    terms = {r['degree']: F(r['coefficient']) for r in residual['residual_coefficients']}
    new = [l for l in terms if l not in candidate['degrees'] and l <= 250]
    new.sort(key=lambda l: (-terms[l]**2*sectors.basis_mass(candidate['azimuth_m'], l), l))
    available = max(0, 64-len(candidate['degrees']))
    if new and not available:
        return {'format': 'residual_support_enrichment_v101', 'certificate': candidate,
                'residual': residual, 'status': 'support_dimension_budget', 'added_degrees': []}
    support = sorted(candidate['degrees']+new[:min(max_new_degrees, available)])
    if support == candidate['degrees']:
        return {'format': 'residual_support_enrichment_v101', 'certificate': candidate,
                'residual': residual, 'status': 'support_did_not_grow', 'added_degrees': []}
    proposed = _proposal(candidate['q_coefficients'], candidate['azimuth_m'], candidate['mean_zero'], support)
    quantized = quantize_proposal(proposed, bits, candidate['threshold'])
    accepted = F(quantized['certificate']['rayleigh_quotient']) < F(candidate['rayleigh_quotient'])
    return {'format': 'residual_support_enrichment_v101', 'certificate': quantized['certificate'] if accepted else candidate,
            'residual': residual, 'quantization': quantized, 'added_degrees': sorted(set(support)-set(candidate['degrees'])),
            'status': 'accepted_improvement' if accepted else 'rejected_no_gain'}


def _combine(left, right, weights):
    result = {}
    for series, factor in zip((left, right), weights):
        for l, x in series.items():
            result[l] = result.get(l, F(0))+factor*x
    return {l: x for l, x in sorted(result.items()) if x}


def _two_forms(first, second):
    if not _ray_verify(first) or not _ray_verify(second) or not _compatible(first, second):
        raise ValueError('Two valid compatible trials required')
    u, v = _series(first), _series(second)
    m, q = first['azimuth_m'], first['q_coefficients']
    h = [_apply(q, w, m, first['mean_zero'])[0] for w in (u, v)]
    mass = [[_dot(x, y, m) for y in (u, v)] for x in (u, v)]
    form = [[_dot(x, y, m) for y in h] for x in (u, v)]
    if base.inertia(mass) != [0, 0, 2]:
        raise ValueError('Trial span must have dimension two')
    return (u, v), mass, form


def _negative_direction(matrix):
    x, b, y = matrix[0][0], matrix[0][1], matrix[1][1]
    if x <= 0:
        return [F(1), F(0)]
    if y <= 0:
        return [F(0), F(1)]
    if x*y-b*b <= 0:
        return [-b, x]
    raise ValueError('No nonpositive rational direction')


def _ritz_record(first, second, lower, weights):
    (u, v), mass, form = _two_forms(first, second)
    lower = _coefficient(lower)
    weights = [_coefficient(x) for x in weights]
    if len(weights) != 2 or not any(weights):
        raise ValueError('Nonzero Ritz mixing coefficients required')
    count = base.inertia([[form[i][j]-lower*mass[i][j] for j in range(2)] for i in range(2)])
    if count[0]:
        raise ValueError('Claimed two-dimensional Ritz lower is invalid')
    trial = _combine(u, v, weights)
    candidate = rayleigh_certificate(first['q_coefficients'], list(trial.values()), first['azimuth_m'],
                first['mean_zero'], list(trial), first['threshold'])
    upper = F(candidate['rayleigh_quotient'])
    return {'format': RITZ, 'scope': 'minimum_Rayleigh_quotient_in_the_span_of_the_two_stated_trials_only',
            'first': first, 'second': second, 'mass_matrix': [_strings(row) for row in mass],
            'form_matrix': [_strings(row) for row in form], 'mixing_coefficients': _strings(weights),
            'ritz_lower': str(lower), 'actual_trial_upper': str(upper), 'ritz_gap': str(upper-lower),
            'lower_inertia': count, 'certificate': candidate,
            'limitation': 'ritz_lower is not a lower bound for the full sphere or the full projected sector.',
            'formal_assistant_checked': False}


def two_trial_ritz(first, second, bits=32):
    """V102: exact 2x2 Ritz lower plus an explicit rational trial upper."""
    (u, v), mass, form = _two_forms(first, second)
    if type(bits) is not int or not 1 <= bits <= 256:
        raise ValueError('Ritz bisection bits in 1..256')
    q = sectors.potential(first['q_coefficients'])
    lower = min(l*(l+1) for l in set(u)|set(v))-sum(map(abs, q))
    upper = min(F(first['rayleigh_quotient']), F(second['rayleigh_quotient']))
    for _ in range(bits):
        middle = (lower+upper)/2
        inertia = base.inertia([[form[i][j]-middle*mass[i][j] for j in range(2)] for i in range(2)])
        if inertia[0]:
            upper = middle
        else:
            lower = middle
            if inertia[1]:
                upper = middle
                break
    shifted = [[form[i][j]-upper*mass[i][j] for j in range(2)] for i in range(2)]
    return _ritz_record(first, second, lower, _negative_direction(shifted))


def _rational_sqrt(x):
    if x < 0:
        return None
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    return F(a, b) if a*a == x.numerator and b*b == x.denominator else None


def residual_line_search(candidate, bits=24):
    """V103: rational stationary-root / dyadic search along negative residual."""
    if type(bits) is not int or not 1 <= bits <= 128:
        raise ValueError('Line-search quantization bits in 1..128')
    residual = residual_certificate(candidate)
    direction = {r['degree']: -F(r['coefficient']) for r in residual['residual_coefficients']}
    if not direction or any(l > 250 for l in direction):
        return {'format': 'residual_line_search_v103', 'certificate': candidate,
                'residual': residual, 'status': 'zero_residual_or_degree_budget', 'trace': []}
    dcert = rayleigh_certificate(candidate['q_coefficients'], list(direction.values()),
                candidate['azimuth_m'], candidate['mean_zero'], list(direction), candidate['threshold'])
    (u, d), mass, form = _two_forms(candidate, dcert)
    n, s, t = mass[0][0], mass[0][1], mass[1][1]
    a, b, c = form[0][0], form[0][1], form[1][1]
    polynomial = [b*n-a*s, c*n-a*t, c*s-b*t]
    scale = max(map(abs, polynomial), default=F(0))
    normalized = [x/scale for x in polynomial] if scale else polynomial
    roots, roots_exact = [], True
    root_failure = None
    if normalized[2]:
        discriminant = normalized[1]**2-4*normalized[2]*normalized[0]
        square = _rational_sqrt(discriminant)
        if square is not None:
            roots = [(-normalized[1]+sign*square)/(2*normalized[2]) for sign in (-1, 1)]
        elif discriminant >= 0:
            roots_exact = False
            try:
                square_float = sqrt(float(discriminant))
                roots = [F(round((float(-normalized[1])+sign*square_float)/
                                  float(2*normalized[2])*2**bits), 2**bits) for sign in (-1, 1)]
            except (OverflowError, ValueError, ZeroDivisionError) as exc:
                roots = []
                root_failure = str(exc)
    elif normalized[1]:
        roots = [-normalized[0]/normalized[1]]
    steps = sorted(set([F(0)]+[F(1, 2**k) for k in range(21)]+roots))
    best, trace = candidate, []
    for step in steps:
        trial = _combine(u, d, [F(1), step])
        cert = rayleigh_certificate(candidate['q_coefficients'], list(trial.values()), candidate['azimuth_m'],
                    candidate['mean_zero'], list(trial), candidate['threshold'])
        accepted = F(cert['rayleigh_quotient']) < F(best['rayleigh_quotient'])
        if accepted:
            best = cert
        trace.append({'step': str(step), 'accepted': accepted, 'certificate': cert})
    return {'format': 'residual_line_search_v103', 'certificate': best, 'baseline': candidate,
            'residual': residual, 'derivative_numerator_divided_by_two': _strings(polynomial),
            'stationary_roots_exact': roots_exact, 'stationary_root_proposals': _strings(roots),
            'root_proposal_failure': root_failure,
            'status': 'accepted_improvement' if F(best['rayleigh_quotient']) < F(candidate['rayleigh_quotient']) else 'no_improvement',
            'trace': trace, 'limitation': 'Only evaluated rational line trials are accepted; rounded roots do not prove exact line optimality.'}


def find_counterexample(q, threshold, mean_zero=True, full_certificate=None, modes=4, max_steps=2):
    """V104: explicitly exhibit R<threshold, or report honestly not found."""
    q, threshold = sectors.potential(q), exact(threshold)
    if type(max_steps) is not int or not 0 <= max_steps <= 8:
        raise ValueError('At most eight enrichment steps')
    if full_certificate is None:
        full_certificate = projected.full_ground(q, mean_zero=mean_zero, modes=max(2, modes),
                            max_modes=max(8, modes), bits=36, tolerance=F(1, 10**8))
    if not projected.verify_full(full_certificate, expected_q=q, expected_mean_zero=mean_zero):
        raise ValueError('Full spectral source has wrong q, mean constraint, or proof')
    selection = select_sector(full_certificate, threshold)
    proposal = mass_normalized_proposal(q, selection['selected_m'], mean_zero, modes)
    quantized = quantize_proposal(proposal, threshold=threshold)
    sparse = sparsify(quantized['certificate'])
    best, trace = sparse['certificate'], []
    for _ in range(max_steps):
        line = residual_line_search(best)
        if F(line['certificate']['rayleigh_quotient']) < F(best['rayleigh_quotient']):
            best = line['certificate']
        enriched = residual_enrich(best)
        if F(enriched['certificate']['rayleigh_quotient']) < F(best['rayleigh_quotient']):
            best = enriched['certificate']
        trace.append({'line_search': line, 'support_enrichment': enriched})
        if F(best['rayleigh_quotient']) < threshold:
            break
    return {'format': FINAL, 'geometry': 'unit_S2', 'scope': projected.scope(mean_zero),
            'q_coefficients': _strings(q), 'mean_zero': mean_zero, 'threshold': str(threshold),
            'status': 'explicit_counterexample_verified' if F(best['rayleigh_quotient']) < threshold else 'not_found',
            'certificate': best, 'selection': selection, 'quantization': quantized, 'sparsification': sparse,
            'trace': trace, 'explicit_real_harmonic_expansion': best['harmonic_terms'],
            'not_found_is_not_a_proof_of_the_inequality': True, 'formal_assistant_checked': False}


def verify(cert, expected_q=None, expected_m=None, expected_mean_zero=None,
           expected_coefficients=None, expected_threshold=None):
    """Independent exact replay with explicit input/scope binding."""
    try:
        kind = cert['format']
        if kind == RAYLEIGH:
            valid, trial = _ray_verify(cert), cert
        elif kind == RESIDUAL:
            valid, trial = same(cert, residual_certificate(cert['trial'])), cert['trial']
        elif kind == RITZ:
            valid = same(cert, _ritz_record(cert['first'], cert['second'], cert['ritz_lower'], cert['mixing_coefficients']))
            trial = cert['certificate']
        elif kind == 'candidate_sector_selection_v96':
            valid = same(cert, select_sector(cert['source'], cert['threshold']))
            if expected_q is not None and sectors.potential(expected_q) != sectors.potential(cert['q_coefficients']):
                return False
            if expected_m is not None and (type(expected_m) is not int or expected_m != cert['selected_m']):
                return False
            if expected_mean_zero is not None and (type(expected_mean_zero) is not bool or expected_mean_zero is not cert['mean_zero']):
                return False
            if expected_threshold is not None and exact(expected_threshold) != F(cert['threshold']):
                return False
            return valid and expected_coefficients is None
        elif kind == FINAL:
            trial = cert['certificate']
            valid = (_ray_verify(trial) and verify(cert['selection'], expected_q=cert['q_coefficients'],
                         expected_m=trial['azimuth_m'], expected_mean_zero=cert['mean_zero'], expected_threshold=cert['threshold'])
                     and trial['q_coefficients'] == cert['q_coefficients'] and trial['mean_zero'] is cert['mean_zero']
                     and trial['threshold'] == cert['threshold'] and cert['scope'] == projected.scope(cert['mean_zero'])
                     and cert['geometry'] == 'unit_S2' and cert['explicit_real_harmonic_expansion'] == trial['harmonic_terms']
                     and cert['status'] == ('explicit_counterexample_verified' if F(trial['rayleigh_quotient']) < F(cert['threshold']) else 'not_found')
                     and cert['not_found_is_not_a_proof_of_the_inequality'] is True and cert['formal_assistant_checked'] is False)
        else:
            return False
        if not valid:
            return False
        if expected_q is not None and sectors.potential(expected_q) != sectors.potential(trial['q_coefficients']):
            return False
        if expected_m is not None and (type(expected_m) is not int or expected_m != trial['azimuth_m']):
            return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool or expected_mean_zero is not trial['mean_zero']):
            return False
        if expected_coefficients is not None and [_coefficient(x) for x in expected_coefficients] != list(map(F, trial['coefficients'])):
            return False
        if expected_threshold is not None and (trial['threshold'] is None or exact(expected_threshold) != F(trial['threshold'])):
            return False
        return True
    except (ValueError, TypeError, KeyError, IndexError, ArithmeticError, OverflowError):
        return False


def run_round01():
    """Reproduce every V95--V104 example and save its original certificates."""
    import json
    from time import perf_counter
    source = json.loads((ROOT/'round01/baseline.json').read_text())['old_certificate']['evidence']
    q, threshold = [0, 0, 2], F(12, 5)
    one = rayleigh_certificate(q, [1], m=1, threshold=threshold)
    stages = []

    def record(version, name, function, test, example, result, certificates, gain):
        paths = []
        for label, certificate in certificates:
            if not verify(certificate):
                raise ArithmeticError('Stage '+str(version)+' certificate failed: '+label)
            paths.append(save('round01/certs/V'+str(version)+'_'+label+'.json', certificate))
        result_path = save('round01/results/V'+str(version)+'.json', result)
        stages.append({'version': version, 'name': name, 'function': 'witness.'+function,
                       'test': 'test_witness.WitnessStages.'+test, 'example': example,
                       'evidence_path': result_path, 'evidence_paths': [result_path]+paths,
                       'certificate_paths': paths, 'actual_outcome': gain,
                       'status': 'implemented_and_example_verified', 'new_mathematical_theorem_claimed': False})

    start = perf_counter()
    simple = rayleigh_certificate(q, [1, F(-1, 36)], m=1, degrees=[1, 3], threshold=threshold)
    record(95, '精确任意有理球谐Rayleigh证书', 'rayleigh_certificate', 'test_v95_exact_trial',
           'q=2t²，u=(P1^1−P3^1/36)cosφ，等价于x1(5−x3²)的非零倍数',
           simple, [('rayleigh', simple)], '商722/303，比阈值12/5严格小26/1515；系数分母不受旧小预算限制。')
    selection = select_sector(source, threshold)
    record(96, '从完整谱证书选择候选方位扇区', 'select_sector', 'test_v96_sector_selection',
           '对baseline完整零均值谱证书按扇区上界排序', selection, [('selection', selection)],
           '选m=1；m=0不能给本阈值的违反者，显式保留选择尚未展示函数的限制。')
    proposal = mass_normalized_proposal(q, selection['selected_m'], True, modes=4)
    proposed_trial = rayleigh_certificate(q, [F(round(x*2**16), 2**16) for x in proposal['coefficients_float']],
                                          m=1, degrees=proposal['degrees'], threshold=threshold)
    record(97, '质量归一化最低本征向量提议', 'mass_normalized_proposal', 'test_v97_mass_normalized_proposal',
           'm=1、连续四个degree的有限形式矩阵；只将量化后的函数用于证明',
           {'proposal': proposal, 'certificate': proposed_trial}, [('rationalized_trial', proposed_trial)],
           '浮点单位质量向量给R约2.38266866的四模候选，低于单模2.4；浮点数本身不作证。')
    quantized = quantize_proposal(proposal, threshold=threshold)
    record(98, '逐级有理量化且保留最佳商', 'quantize_proposal', 'test_v98_quantization_keeps_best',
           '同一提议依次使用2、4、8、16、24bits', quantized,
           [('bits'+str(row['bits']), row['certificate']) for row in quantized['trace'] if 'certificate' in row],
           '8/16/24bits改善；4bits没有改善而保留旧候选，不假设分母更大必更好。')
    noisy = rayleigh_certificate(q, [1, F(-1, 36), F(1, 1000)], m=1, degrees=[1, 3, 7], threshold=threshold)
    sparse = sparsify(noisy, (1, 2))
    record(99, '质量贡献排序的稀疏截断与精确接受', 'sparsify', 'test_v99_sparse_acceptance_and_rejection',
           '三项含P7^1噪声试探，分别保留一项、两项', sparse,
           [('noisy_baseline', noisy)]+[('keep'+str(row['keep_terms']), row['certificate']) for row in sparse['trace']],
           '保留一项使商回到2.4，明确拒绝；删除P7^1噪声保留两项严格改善。')
    residual = residual_certificate(quantized['certificate'])
    calibration = residual_certificate(rayleigh_certificate([0, 1], [1], m=0))
    record(100, '完整投影残差与所有势乘法高阶', 'residual_certificate', 'test_v100_complete_projected_residual',
           '量化后q=2t²试探及q=t,u=P1的投影校准',
           {'candidate_residual': residual, 'projection_calibration': calibration},
           [('candidate_residual', residual), ('projection_calibration', calibration)],
           '校准移除l0系数1/3，完整归一化残差平方4/15；不误用未投影的3/5。')
    enriched = residual_enrich(one)
    no_growth = residual_enrich(rayleigh_certificate([0], [1], m=1))
    record(101, '残差驱动新增试探支持', 'residual_enrich', 'test_v101_residual_driven_support_growth',
           'q=2t²从单个P1^1出发，由非零残差检测并新增degree3',
           {'improvement': enriched, 'no_growth': no_growth},
           [('enriched_trial', enriched['certificate']), ('zero_residual', no_growth['residual'])],
           '自动新增degree3并严格降低R；q=0精确本征函数无新增degree，保留support_did_not_grow。')
    second = rayleigh_certificate(q, [1], m=1, degrees=[3], threshold=threshold)
    ritz = two_trial_ritz(one, second)
    record(102, '二维Ritz精确下界及实际试探上界', 'two_trial_ritz', 'test_v102_two_trial_ritz_scope_and_actual_upper',
           'span(P1^1 cosφ,P3^1 cosφ)的2×2广义本征问题', ritz,
           [('ritz', ritz), ('actual_trial', ritz['certificate'])],
           '二维最小商误差小于1e-8；Ritz下界仅限此二维空间，禁止充当完整球面谱下界。')
    line = residual_line_search(one)
    exact_line = residual_line_search(rayleigh_certificate([0], [1, 1], m=0, degrees=[1, 2]))
    record(103, '负残差方向的驻点及有理线搜索', 'residual_line_search', 'test_v103_negative_residual_line_search',
           'q=2t²从P1^1沿负残差下降；另用q=0,P1+P2演示有理驻点',
           {'irrational_root_proposals': line, 'exact_stationary_case': exact_line},
           [('line_trial', line['certificate']), ('exact_root_trial', exact_line['certificate'])],
           '前者从2.4降低到约2.3826686613；后者有理驻点精确获得R=2。根转浮点前按有理系数归一化，避免缩放溢出。')
    final = find_counterexample(q, threshold, full_certificate=source)
    failed = find_counterexample([0], F(1), modes=2, max_steps=0)
    record(104, '自动构造显式球谐违反者', 'find_counterexample', 'test_v104_explicit_refutation_and_honest_not_found',
           '输入q=2t²、threshold=12/5和baseline源证书；另测试q=0、threshold=1',
           {'successful': final, 'not_found': failed}, [('counterexample', final), ('honest_not_found', failed)],
           '自动选m=1并给严格R<threshold及明确球谐系数；q=0阈值1未找到，未伪称证明或反例。')
    if [row['version'] for row in stages] != list(range(95, 105)):
        raise ArithmeticError('Stage coverage mismatch')
    save('round01/STAGES.json', stages)
    return {'stage_count': len(stages), 'all_example_certificates_verified': True,
            'seconds_including_verification': perf_counter()-start,
            'automatic_rayleigh': final['certificate']['rayleigh_quotient'],
            'automatic_strict_margin': final['certificate']['strict_margin'],
            'simple_rayleigh': simple['rayleigh_quotient'], 'simple_margin': simple['strict_margin'],
            'not_found_status': failed['status']}


if __name__ == '__main__':
    import json
    result = run_round01()
    save('round01/results/summary.json', result)
    print(json.dumps(result, ensure_ascii=False, indent=2))
