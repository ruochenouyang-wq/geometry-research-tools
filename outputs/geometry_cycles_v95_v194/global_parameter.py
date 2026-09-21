"""V135--V144: globally bounded affine families of projected sphere operators.

C(q) is the infimum of the Rayleigh quotient on the FULL mean-zero H1(S2)
space. Its concavity in q and its L-infinity Lipschitz estimate are analytic
dependencies. Rational certificates replay finite/infinite spectral tails,
polynomial extrema, the complete parameter cover, and all bound directions.
"""
from fractions import Fraction as F
from pathlib import Path
import json
from common import (ROOT, arithmetic as A, projected, exact_extrema as E,
                    exact, rational_vector, same, digest, save)

FORMAT = 'global_affine_projected_sphere_v144'


def canonical_objective(raw):
    """Bind the geometry, constraint, affine family, penalty and domain."""
    if not isinstance(raw, dict) or set(raw) != {'q0', 'direction', 'penalty', 'domain'}:
        raise ValueError('Expected q0, direction, penalty, domain only')
    q0 = projected.potential(raw['q0'])
    direction = projected.potential(raw['direction'])
    penalty = A.trim(rational_vector(raw['penalty'], maximum=5))
    domain = rational_vector(raw['domain'], maximum=2)
    if len(domain) != 2 or not domain[0] < domain[1]:
        raise ValueError('A nonempty closed rational interval is required')
    # Validate supported potential inputs at the endpoints. Coefficient norms
    # are convex, hence endpoint bounds cover the whole affine segment.
    for s in domain:
        projected.potential([str(v) for v in A.add(q0, A.scale(direction, s))])
    return {k: [str(v) for v in p] for k, p in
            [('q0', q0), ('direction', direction), ('penalty', penalty), ('domain', domain)]}


def _interval_product(a, b):
    products = [x*y for x in a for y in b]
    return min(products), max(products)


def _power_interval(a, b, degree):
    if degree == 0:
        return F(1), F(1)
    values = [a**degree, b**degree]
    if degree % 2 == 0 and a <= 0 <= b:
        values.append(F(0))
    return min(values), max(values)


def _cheap_assemble(obj, coefficients, potential_terms, penalty_terms, candidates):
    low_q = sum((v[0] for v in potential_terms), F(0))
    low_p = sum((v[0] for v in penalty_terms), F(0))
    lower = 2+low_q+low_p
    chosen = min(candidates, key=lambda name: (candidates[name], name))
    upper = candidates[chosen]
    if lower > upper:
        raise ArithmeticError('Contradictory cheap family bounds')
    strings = lambda rows: [[str(x) for x in row] for row in rows]
    return {'format': 'cheap_global_affine_sphere_v135', 'objective': obj,
        'objective_sha256': digest(obj), 'geometry': 'unit_S2',
        'function_space': 'all_mean_zero_H1_functions',
        'potential_coefficient_intervals': strings(coefficients),
        'potential_term_intervals': strings(potential_terms),
        'penalty_term_intervals': strings(penalty_terms),
        'potential_lower': str(low_q), 'penalty_lower': str(low_p),
        'lower': str(lower), 'upper': str(upper), 'gap': str(upper-lower),
        'feasible_parameter': str(sum(map(F, obj['domain']), F(0))/2),
        'trial_function': chosen, 'trial_squared_L2_probability_norm': '1/3',
        'trial_Dirichlet_probability_energy': '2/3',
        'trial_objective_upper_candidates': {k: str(v) for k, v in candidates.items()},
        'point_spectrum_solver_calls': 0, 'formal_assistant_checked': False,
        'analytic_dependencies': ['zero_mean_unit_sphere_Poincare_constant_2',
             'Rayleigh_principle_on_full_mean_zero_space', 'exact_sphere_coordinate_moments']}


def cheap_objective_enclosure(objective):
    """V135: a whole-family bound using interval arithmetic and l=1 trials.

    No matrix assembly, eigenvalue calculation or polynomial extremum solver
    is called. Correlations are deliberately discarded only in the lower
    bound; the upper bound uses an actual midpoint and actual zero-mean trial.
    """
    obj = canonical_objective(objective)
    a, b = map(F, obj['domain'])
    q0, direction = list(map(F, obj['q0'])), list(map(F, obj['direction']))
    size = max(len(q0), len(direction))
    q0 += [F(0)]*(size-len(q0))
    direction += [F(0)]*(size-len(direction))
    coefficients = [tuple(sorted([c+a*d, c+b*d])) for c, d in zip(q0, direction)]
    potential_terms = [_interval_product(c, _power_interval(F(-1), F(1), i))
                       for i, c in enumerate(coefficients)]
    penalty_terms = [_interval_product((c, c), _power_interval(a, b, i))
                     for i, c in enumerate(map(F, obj['penalty']))]
    middle = (a+b)/2
    q = _q_at(obj, middle)
    penalty = A.evaluate(list(map(F, obj['penalty'])), middle)
    candidates = {
        'x1': 2+penalty+sum((F(3)*q[i]/((i+1)*(i+3)) for i in range(0, len(q), 2)), F(0)),
        'x3': 2+penalty+sum((F(3)*q[i]/(i+3) for i in range(0, len(q), 2)), F(0))}
    return _cheap_assemble(obj, coefficients, potential_terms, penalty_terms, candidates)


def verify_cheap_objective(cert, expected_objective=None):
    """Recompute from the bound original input, independently of generation.

    In particular the upper trial values are recalculated as integrals of
    q*x3^2 and q*(1-x3^2)/2, rather than trusting the generator's closed forms.
    """
    try:
        if not isinstance(cert, dict) or cert.get('format') != 'cheap_global_affine_sphere_v135':
            return False
        obj = canonical_objective(cert['objective'])
        if expected_objective is not None and not same(obj, canonical_objective(expected_objective)):
            return False
        a, b = map(F, obj['domain'])
        q0, direction = list(map(F, obj['q0'])), list(map(F, obj['direction']))
        coefficients, potential_terms, penalty_terms = [], [], []
        for i in range(max(len(q0), len(direction))):
            c = q0[i] if i < len(q0) else F(0)
            d = direction[i] if i < len(direction) else F(0)
            lo, hi = min(c+a*d, c+b*d), max(c+a*d, c+b*d)
            coefficients.append((lo, hi))
            if i == 0:
                potential_terms.append((lo, hi))
            elif i % 2:
                radius = max(abs(lo), abs(hi))
                potential_terms.append((-radius, radius))
            else:
                potential_terms.append((min(F(0), lo), max(F(0), hi)))
        for i, c in enumerate(map(F, obj['penalty'])):
            powers = [a**i, b**i]
            if i > 0 and i % 2 == 0 and a <= 0 <= b:
                powers.append(F(0))
            values = [c*v for v in powers]
            penalty_terms.append((min(values), max(values)))
        middle = (a+b)/2
        q = []
        for i in range(max(len(q0), len(direction))):
            q.append((q0[i] if i < len(q0) else 0)+middle*(direction[i] if i < len(direction) else 0))
        def moment(degree):
            return F(0) if degree % 2 else F(1, degree+1)
        potential_x3 = sum((c*moment(i+2) for i, c in enumerate(q)), F(0))
        potential_x1 = sum((c*(moment(i)-moment(i+2))/2 for i, c in enumerate(q)), F(0))
        penalty = sum((c*middle**i for i, c in enumerate(map(F, obj['penalty']))), F(0))
        candidates = {'x1': (F(2, 3)+potential_x1)/F(1, 3)+penalty,
                      'x3': (F(2, 3)+potential_x3)/F(1, 3)+penalty}
        return same(cert, _cheap_assemble(obj, coefficients, potential_terms, penalty_terms, candidates))
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


def _poly(proof, polynomial, interval):
    p = [str(v) for v in A.trim(polynomial)]
    return (isinstance(proof, dict) and proof.get('polynomial') == p
            and proof.get('interval') == [str(v) for v in interval] and E.verify(proof))


def penalty_enclosure(penalty, interval, evidence=None):
    """V136: exact global polynomial extrema, including nonconvex penalties."""
    p, domain = A.trim(rational_vector(penalty, maximum=5)), list(map(exact, interval))
    if len(domain) != 2 or not domain[0] < domain[1]:
        raise ValueError('Invalid penalty interval')
    if evidence is None:
        evidence = {'negative_maximum': E.maximize(A.scale(p, -1), domain),
                    'positive_maximum': E.maximize(p, domain)}
    if (not _poly(evidence['negative_maximum'], A.scale(p, -1), domain)
            or not _poly(evidence['positive_maximum'], p, domain)):
        raise ValueError('Penalty proof not bound to input and interval')
    return {'lower': str(-F(evidence['negative_maximum']['maximum_upper'])),
            'upper': evidence['positive_maximum']['maximum_upper'],
            'negative_maximum': evidence['negative_maximum'],
            'positive_maximum': evidence['positive_maximum']}


def _q_at(obj, s):
    return A.add(list(map(F, obj['q0'])), A.scale(list(map(F, obj['direction'])), s))


def _point_from_spectrum(obj, s, spectrum):
    q = _q_at(obj, s)
    if not projected.verify_full(spectrum, expected_q=q, expected_mean_zero=True):
        raise ValueError('Point must carry the full mean-zero spectral proof')
    penalty = A.evaluate(list(map(F, obj['penalty'])), s)
    return {'parameter': str(s), 'penalty': str(penalty),
            'lower': str(F(spectrum['lower'])+penalty),
            'upper': str(F(spectrum['upper'])+penalty), 'spectrum': spectrum}


def feasible_point(objective, parameter, bits=32, modes=4):
    """V138: a rational feasible parameter supplies a certified incumbent."""
    obj, s = canonical_objective(objective), exact(parameter)
    if not F(obj['domain'][0]) <= s <= F(obj['domain'][1]):
        raise ValueError('Incumbent is outside the parameter domain')
    spectrum = projected.full_ground(_q_at(obj, s), mean_zero=True, modes=modes,
          max_modes=max(8, modes), bits=bits, max_m=12, tolerance=F(1, 10**8))
    return _point_from_spectrum(obj, s, spectrum)


def cached_point(objective, parameter, cache, bits=32, modes=4):
    """V141: only input-bound, replayed point evidence may enter the cache."""
    obj, s = canonical_objective(objective), exact(parameter)
    if not F(obj['domain'][0]) <= s <= F(obj['domain'][1]):
        raise ValueError('Cache point outside domain')
    key = str(s)
    if key in cache:
        expected = _point_from_spectrum(obj, s, cache[key]['spectrum'])
        if not same(expected, cache[key]):
            raise ValueError('Stale or corrupt cached point')
    else:
        cache[key] = feasible_point(obj, s, bits=bits, modes=modes)
    return cache[key]


def lipschitz_cell(objective, interval, points, penalty=None):
    """V137: every point yields a lower bound over the entire cell."""
    obj = canonical_objective(objective)
    left, right = map(exact, interval)
    if not F(obj['domain'][0]) <= left < right <= F(obj['domain'][1]):
        raise ValueError('Cell outside domain')
    penalty = penalty_enclosure(obj['penalty'], [left, right], penalty)
    constant = sum(map(abs, map(F, obj['direction'])), F(0))
    lower = []
    for point in points:
        s = F(point['parameter'])
        expected = _point_from_spectrum(obj, s, point['spectrum'])
        if not left <= s <= right or not same(expected, point):
            raise ValueError('Unbound point in Lipschitz cell')
        lower.append(F(point['spectrum']['lower'])-constant*max(s-left, right-s)
                     + F(penalty['lower']))
    if not lower:
        raise ValueError('A point proof is required')
    return {'lower': str(max(lower)), 'lipschitz_constant': str(constant),
            'penalty_bounds': penalty}


def concavity_cell(objective, interval, left_point, right_point, proof=None):
    """V140: C(q_s) lies above the chord of its endpoint LOWER bounds."""
    obj = canonical_objective(objective)
    a, b = map(exact, interval)
    if not F(obj['domain'][0]) <= a < b <= F(obj['domain'][1]):
        raise ValueError('Invalid chord interval')
    for s, point in [(a, left_point), (b, right_point)]:
        if not same(_point_from_spectrum(obj, s, point['spectrum']), point):
            raise ValueError('Chord endpoint is not input-bound')
    ca, cb = F(left_point['spectrum']['lower']), F(right_point['spectrum']['lower'])
    slope = (cb-ca)/(b-a)
    polynomial = A.add(list(map(F, obj['penalty'])), [ca-a*slope, slope])
    if proof is None:
        proof = E.maximize(A.scale(polynomial, -1), [a, b])
    if not _poly(proof, A.scale(polynomial, -1), [a, b]):
        raise ValueError('Chord-polynomial proof invalid')
    return {'lower': str(-F(proof['maximum_upper'])),
            'chord_plus_penalty': [str(v) for v in polynomial], 'negative_maximum': proof}


def _cell(obj, a, b, points, strategy, evidence=None):
    keys = [str(a), str((a+b)/2), str(b)]
    lip = lipschitz_cell(obj, [a, b], [points[k] for k in keys],
                         None if evidence is None else evidence['lipschitz']['penalty_bounds'])
    chord = (concavity_cell(obj, [a, b], points[keys[0]], points[keys[2]],
                           None if evidence is None else evidence['concavity']['negative_maximum'])
             if strategy == 'concavity' else None)
    return {'interval': [str(a), str(b)], 'point_keys': keys,
            'lower': str(max(F(lip['lower']), F(chord['lower'])) if chord else F(lip['lower'])),
            'lipschitz': lip, 'concavity': chord}


def _assemble(obj, points, cells, strategy, tolerance, budget, seed=None):
    cells = sorted(cells, key=lambda c: F(c['interval'][0]))
    a, b = map(F, obj['domain'])
    if (not cells or F(cells[0]['interval'][0]) != a or F(cells[-1]['interval'][1]) != b
            or any(cells[i]['interval'][1] != cells[i+1]['interval'][0]
                   for i in range(len(cells)-1))):
        raise ValueError('Incomplete or overlapping domain cover')
    incumbent = min(points, key=lambda k: (F(points[k]['upper']), F(k)))
    upper = F(points[incumbent]['upper'])
    lower = min(F(c['lower']) for c in cells)
    if seed is not None:
        lower = max(lower, F(seed['lower']))
    if lower > upper:
        raise ArithmeticError('Contradictory bounds')
    decorated = [dict(c, status='excluded' if F(c['lower']) > upper else 'live') for c in cells]
    return {'format': FORMAT, 'objective': obj, 'objective_sha256': digest(obj),
            'geometry': 'unit_S2', 'function_space': 'all_mean_zero_H1_functions',
            'strategy': strategy, 'points': points, 'cells': decorated,
            'lower': str(lower), 'upper': str(upper), 'gap': str(upper-lower),
            'incumbent_parameter': incumbent, 'tolerance': str(tolerance),
            'max_leaves': budget, 'partition_leaves': len(cells),
            'status': 'certified_target_met' if upper-lower <= tolerance else 'certified_bound_open_gap',
            'minimizer_intervals': [c['interval'] for c in decorated if c['status'] == 'live'],
            'uniqueness_claimed': False, 'seed': seed, 'formal_assistant_checked': False,
            'analytic_dependencies': ['minmax_full_mean_zero_sphere',
              'concavity_of_infimum_of_affine_Rayleigh_forms', 'bounded_multiplier_L_infinity_perturbation']}


def branch_and_bound(objective, max_leaves=16, tolerance=F(1, 10000),
                     strategy='concavity', bits=32, modes=4, seed=None):
    """V139: deterministic best-lower-first search retaining a full cover."""
    obj, tolerance = canonical_objective(objective), exact(tolerance)
    if type(max_leaves) is not int or not 1 <= max_leaves <= 128 or tolerance <= 0:
        raise ValueError('Positive tolerance and 1..128 leaves required')
    if strategy not in ('lipschitz', 'concavity'):
        raise ValueError('Unknown lower-bound strategy')
    if seed is not None:
        if not verify(seed, expected_objective=obj) or seed['partition_leaves'] > max_leaves:
            raise ValueError('Checkpoint has wrong input, invalid proof or excessive leaves')
        ancestor, depth = seed, 0
        while ancestor['seed'] is not None:
            ancestor, depth = ancestor['seed'], depth+1
        if depth >= 8:
            raise ValueError('At most eight nested checkpoint continuations supported')
        points = json.loads(json.dumps(seed['points']))
        intervals = [list(map(F, c['interval'])) for c in seed['cells']]
    else:
        points = {}
        intervals = [list(map(F, obj['domain']))]
    def new_cell(a, b):
        for s in [a, (a+b)/2, b]:
            cached_point(obj, s, points, bits=bits, modes=modes)
        return _cell(obj, a, b, points, strategy)
    cells = [new_cell(a, b) for a, b in intervals]
    while True:
        cert = _assemble(obj, points, cells, strategy, tolerance, max_leaves, seed)
        if F(cert['gap']) <= tolerance or len(cells) >= max_leaves:
            return cert
        selected = min(range(len(cells)), key=lambda i:
             (F(cells[i]['lower']), F(cells[i]['interval'][0])))
        old = cells.pop(selected)
        a, b = map(F, old['interval'])
        middle = (a+b)/2
        cells.extend([new_cell(a, middle), new_cell(middle, b)])


def minimizer_enclosure(certificate):
    """V142: enclose every minimizer; exclusion is strict, never uniqueness."""
    if not verify(certificate):
        raise ValueError('Invalid global certificate')
    return {'all_minimizers_in': certificate['minimizer_intervals'],
            'excluded_intervals': [c['interval'] for c in certificate['cells'] if c['status'] == 'excluded'],
            'uniqueness_claimed': False}


def resume(objective, checkpoint, max_leaves=32, tolerance=F(1, 100000), **kwargs):
    """V143: input-bound checkpoints preserve previously certified bounds."""
    return branch_and_bound(objective, seed=checkpoint, max_leaves=max_leaves,
                            tolerance=tolerance, **kwargs)


def verify(cert, expected_objective=None, _depth=0):
    """Replay evidence and coverage without rerunning the optimization search."""
    try:
        if not isinstance(cert, dict) or _depth > 8 or cert.get('format') != FORMAT:
            return False
        obj = canonical_objective(cert['objective'])
        if expected_objective is not None and not same(obj, canonical_objective(expected_objective)):
            return False
        budget, tolerance, strategy = cert['max_leaves'], exact(cert['tolerance']), cert['strategy']
        if type(budget) is not int or not 1 <= len(cert['cells']) <= budget <= 128 or tolerance <= 0:
            return False
        if strategy not in ('lipschitz', 'concavity'):
            return False
        points = cert['points']
        if not isinstance(points, dict) or not 1 <= len(points) <= 2048:
            return False
        for key, point in points.items():
            s = exact(key)
            if str(s) != key or not F(obj['domain'][0]) <= s <= F(obj['domain'][1]):
                return False
            if not same(point, _point_from_spectrum(obj, s, point['spectrum'])):
                return False
        seed = cert['seed']
        if seed is not None:
            if not verify(seed, obj, _depth+1):
                return False
            if any(k not in points or not same(v, points[k]) for k, v in seed['points'].items()):
                return False
        cells = []
        for cell in cert['cells']:
            a, b = map(exact, cell['interval'])
            expected = _cell(obj, a, b, points, strategy, evidence=cell)
            cells.append(expected)
        expected = _assemble(obj, points, cells, strategy, tolerance, budget, seed)
        return same(expected, cert)
    except (ValueError, KeyError, TypeError, ArithmeticError, IndexError, OverflowError):
        return False


def threshold_certificate(global_certificate, threshold):
    """V144: uniform proof, spectral-existence refutation, or an honest gap."""
    if not verify(global_certificate):
        raise ValueError('Invalid global certificate')
    threshold = exact(threshold)
    lo, hi = F(global_certificate['lower']), F(global_certificate['upper'])
    return {'format': 'global_affine_threshold_v144', 'threshold': str(threshold),
            'status': 'proved' if lo >= threshold else
              'refuted_by_feasible_parameter_spectral_existence' if hi < threshold else 'undetermined',
            'quantifier': 'every_parameter_in_domain_and_every_mean_zero_H1_sphere_function',
            'statement': 'integral(grad_u_squared + (q_s + penalty(s))*u_squared) >= threshold*integral(u_squared)',
            'evidence': global_certificate, 'explicit_function_witness_included': False}


def verify_threshold(cert, expected_objective=None, expected_threshold=None):
    try:
        if expected_objective is not None and not verify(cert['evidence'], expected_objective):
            return False
        if expected_threshold is not None and exact(expected_threshold) != exact(cert['threshold']):
            return False
        return same(cert, threshold_certificate(cert['evidence'], cert['threshold']))
    except (ValueError, KeyError, TypeError, ArithmeticError):
        return False


def run_trials():
    import time
    start = time.perf_counter()
    quartic = {'q0': ['0'], 'direction': ['0'], 'penalty': ['0', '0', '-1', '0', '1'],
               'domain': ['-1', '1']}
    sphere = {'q0': ['0', '0', '2'], 'direction': ['0', '1'],
              'penalty': ['0', '0', '1/10'], 'domain': ['-2', '2']}
    boundary = {'q0': ['0'], 'direction': ['1'], 'penalty': ['0'], 'domain': ['-1', '1']}
    stages = []
    def stage(version, capability, call, result, outcome):
        path = save(f'round05/stage_{version}.json', result)
        stages.append({'version': version, 'capability': capability, 'callable': call,
                       'evidence': [path], 'outcome': outcome})
    obj = canonical_objective(sphere)
    cheap = cheap_objective_enclosure(obj)
    stage(135, 'Cheap whole-family interval lower bound and feasible spherical-harmonic upper bound',
          'global_parameter.cheap_objective_enclosure', cheap,
          'Certified complete-domain enclosure [0,12/5] with no spectral solver calls')
    penalty = penalty_enclosure(quartic['penalty'], quartic['domain'])
    stage(136, 'Nonconvex polynomial penalty interval extrema', 'global_parameter.penalty_enclosure', penalty,
          'Penalty minimum -1/4 enclosed independently of endpoint samples')
    cache = {}
    for s in [F(-2), F(0), F(2)]:
        cached_point(obj, s, cache)
    lip = lipschitz_cell(obj, [-2, 2], list(cache.values()))
    stage(137, 'Certified spectral Lipschitz lower bound on a complete interval',
          'global_parameter.lipschitz_cell', lip, 'Uniform lower bound includes uncovered inter-sample regions')
    stage(138, 'Feasible rational parameter with full-sphere spectral incumbent',
          'global_parameter.feasible_point', cache['0'], 'Upper bound is tied to an actual parameter')
    rough = branch_and_bound(quartic, max_leaves=4, strategy='lipschitz', tolerance=F(1, 10**8))
    stage(139, 'Best-lower-first complete-domain branch and bound', 'global_parameter.branch_and_bound', rough,
          'Three-point sampling misses quartic minima; certified global bound covers them')
    chord = concavity_cell(obj, [-2, 2], cache['-2'], cache['2'])
    stage(140, 'Concavity-chord spectral bound plus globally minimized penalty',
          'global_parameter.concavity_cell', chord, 'Chord bound tightens the whole-cell spectral lower bound')
    before = len(cache)
    hit = cached_point(obj, F(0), cache)
    stage(141, 'Validated cache and endpoint/midpoint reuse', 'global_parameter.cached_point',
          {'cache_before': before, 'cache_after': len(cache), 'reused_point': hit},
          'Reused full-spectrum proof without a repeated solve')
    well = branch_and_bound(quartic, max_leaves=12, tolerance=F(1, 10**5))
    stage(142, 'Strict cell exclusion and enclosure of every minimizer',
          'global_parameter.minimizer_enclosure', {'result': minimizer_enclosure(well), 'certificate': well},
          'Both symmetric minimizers retained; no uniqueness claim')
    resumed = resume(quartic, rough, max_leaves=16, tolerance=F(1, 10**6))
    stage(143, 'Input-bound continuation with monotone global bounds', 'global_parameter.resume', resumed,
          'Old lower and incumbent upper preserved; remaining gap reported')
    actual = branch_and_bound(sphere, max_leaves=16, tolerance=F(1, 10000))
    open_budget = branch_and_bound(sphere, max_leaves=1, tolerance=F(1, 10**8))
    bound = branch_and_bound(boundary, max_leaves=4, tolerance=F(1, 10**8))
    transitions = []
    for alpha in ['1/100', '1/10']:
        family = {'q0': ['0'], 'direction': ['0', '1'], 'penalty': ['0', '0', alpha],
                  'domain': ['-2', '2']}
        proof = branch_and_bound(family, max_leaves=16, tolerance=F(1, 10000))
        path = save('round05/transition_alpha_'+alpha.replace('/', '_')+'.json', proof)
        transitions.append({'alpha': alpha, 'certificate_path': path,
            **{k: proof[k] for k in ['lower', 'upper', 'gap', 'incumbent_parameter', 'minimizer_intervals']}})
    thresholds = [threshold_certificate(actual, F(actual['lower'])),
                  threshold_certificate(actual, F(actual['upper'])+1),
                  threshold_certificate(actual, (F(actual['lower'])+F(actual['upper']))/2)]
    stage(144, 'Uniform threshold proof/refutation/undetermined and full replay',
          'global_parameter.threshold_certificate', {'certificates': thresholds},
          'All three outcomes supported without equating samples with a proof')
    baseline = {'method': 'three_point_sampling_only', 'certifies_global_optimum': False,
         'quartic_parameters': ['-1', '0', '1'], 'sampled_values': ['2', '2', '2'],
         'actual_global_minimum': '7/4', 'uncovered_improvement': '1/4',
         'reason': 'The penalty s^4-s^2 has minima at +-1/sqrt(2) between sampled points.'}
    cert_paths = [save('round05/global_quartic.json', well),
                  save('round05/global_sphere.json', actual), save('round05/global_boundary.json', bound),
                  save('round05/resumed.json', resumed), save('round05/open_budget.json', open_budget)]
    results = {'baseline': baseline, 'quartic': {k: well[k] for k in ['lower', 'upper', 'gap', 'status', 'minimizer_intervals']},
       'cheap_objective': {k: cheap[k] for k in ['lower', 'upper', 'gap', 'trial_function', 'point_spectrum_solver_calls']},
       'cheap_objective_replay_passed': verify_cheap_objective(cheap, sphere),
       'sphere': {k: actual[k] for k in ['lower', 'upper', 'gap', 'status', 'minimizer_intervals']},
       'boundary': {k: bound[k] for k in ['lower', 'upper', 'gap', 'status', 'incumbent_parameter']},
       'lipschitz_whole_interval_lower': lip['lower'], 'chord_whole_interval_lower': chord['lower'],
       'finite_budget_open': {k: open_budget[k] for k in ['lower', 'upper', 'gap', 'status', 'minimizer_intervals']},
       'penalty_transition_cases': transitions,
       'checkpoint_monotone': F(resumed['lower']) >= F(rough['lower']) and F(resumed['upper']) <= F(rough['upper']),
       'threshold_statuses': [c['status'] for c in thresholds], 'certificate_paths': cert_paths,
       'all_certificate_replays_passed': all(verify(c) for c in [rough, well, resumed, actual, bound, open_budget])
             and all(verify(json.loads((ROOT/t['certificate_path']).read_text())) for t in transitions),
       'elapsed_seconds': time.perf_counter()-start}
    save('round05/baseline.json', baseline)
    save('round05/STAGES.json', stages)
    save('round05/results.json', results)
    return results


if __name__ == '__main__':
    print(json.dumps(run_trials(), ensure_ascii=False, indent=2))
