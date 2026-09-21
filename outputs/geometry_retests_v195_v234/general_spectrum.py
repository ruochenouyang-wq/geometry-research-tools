"""V205--V214: full and mean-zero sphere spectra, with trace-aware search.

The unprojected operator is -Delta+q. The mean-zero operator is the form
compression P(-Delta+q)P. For nonconstant q these are different operators;
no constant eigenvalue is inserted to translate between their spectra.
"""
from copy import deepcopy
from support import F, exact, same, digest, sectors, projected, arithmetic

INITIAL = 'two_scope_initial_spectrum_v205'
SECTOR = 'two_scope_sector_enumeration_v206'
MODEL = 'two_scope_verified_model_v206'
ORDERED = 'two_scope_complete_order_statistics_v207'
COUNT = 'two_scope_strict_spectral_count_v208'
TRACE = 'two_scope_complete_negative_trace_v209'
MULTI = 'trace_aware_multiquantity_certificate_v210'
SHIFT = 'projection_preserving_shifted_count_and_trace_v211'
LT = 'two_scope_fixed_potential_trace_quotient_v212'
INTERLACE = 'codimension_one_projection_interlacing_v213'
DRIVER = 'input_bound_multispace_research_v214'


def _r(x):
    if type(x) not in (int, str, F):
        raise ValueError('Use exact rational input, not float or bool')
    return exact(x)


def _i(x, low, high, name):
    if type(x) is not int or not low <= x <= high:
        raise ValueError('%s must be an integer in %d..%d' % (name, low, high))
    return x


def _scope(mean_zero):
    if type(mean_zero) is not bool:
        raise ValueError('mean_zero must be a boolean')
    return {'mean_zero': mean_zero,
            'space': 'full_mean_zero' if mean_zero else 'full_unprojected',
            'operator': 'P_mean_zero*(-Delta+q)*P_mean_zero' if mean_zero else '-Delta+q',
            'geometry': 'unit_S2', 'measure': 'standard_sphere_area',
            'real_angular_multiplicity': {'m=0': 1, 'm>0': 2}}


def _start(m, mean_zero):
    return max(1, m) if mean_zero else m


def _q(q):
    return list(map(str, sectors.potential(q)))


def _range(q):
    proof = projected.v05_range.make_range(sectors.potential(q))
    low, high = projected.v05_range.verify_bounds(sectors.potential(q), proof)
    return {'proof': proof, 'lower': str(low), 'upper': str(high)}


def _tolerance(tolerance):
    tol = _r(tolerance)
    if not F(1, 10**25) <= tol <= 1:
        raise ValueError('Tolerance must lie in [1e-25,1]')
    return tol


def laplace_values(k, mean_zero):
    _i(k, 1, 128, 'k'); _scope(mean_zero)
    values, l = [], 1 if mean_zero else 0
    while len(values) < k:
        values.extend([F(l*(l+1))]*(2*l+1))
        l += 1
    return values[:k]


def analytic_bounds(q, k=9, mean_zero=False):
    """V205: global ordered min--max bounds with the correct l=0 convention."""
    q = _q(q); scope = _scope(mean_zero)
    bounds = _range(q)
    lo, hi = F(bounds['lower']), F(bounds['upper'])
    intervals = [{'index': i+1, 'laplace_level': str(x), 'lower': str(x+lo),
                  'upper': str(x+hi), 'width': str(hi-lo)}
                 for i, x in enumerate(laplace_values(k, mean_zero))]
    return {'format': INITIAL, 'q': q, **scope, 'k': k, 'range': bounds,
            'ordered': intervals, 'l0_multiplicity': 0 if mean_zero else 1,
            'rule': 'complete Laplace ordered spectrum plus inf(q)/sup(q)',
            'schur_computations': 0, 'formal_assistant_checked': False}


def enumerate_sector(q, m=0, radial_count=1, mean_zero=False, modes=8, bits=36):
    """V206: solve the correct projected or unprojected radial operator."""
    q = _q(q); scope = _scope(mean_zero)
    _i(m, 0, 16, 'm'); _i(modes, 2, 32, 'modes')
    _i(radial_count, 1, modes, 'radial_count'); _i(bits, 8, 100, 'bits')
    certs = [projected.certify_sector(q, m=m, mean_zero=mean_zero,
                                    k=j, modes=modes, bits=bits)
             for j in range(1, radial_count+1)]
    return sector_bundle(q, m, mean_zero, certs)


def sector_bundle(q, m, mean_zero, eigenvalues):
    q = _q(q); scope = _scope(mean_zero); _i(m, 0, 16, 'm')
    if not isinstance(eigenvalues, list) or not 1 <= len(eigenvalues) <= 32:
        raise ValueError('Consecutive radial eigenvalues are required')
    for j, cert in enumerate(eigenvalues, 1):
        if not projected.verify_sector(cert, expected_q=q, expected_m=m,
                expected_mean_zero=mean_zero, expected_k=j, independent=True):
            raise ValueError('Sector proof has changed q, projection, m, k or evidence')
    return {'format': SECTOR, 'q': q, **scope, 'm': m,
            'degree_start': _start(m, mean_zero), 'multiplicity': 1 if m == 0 else 2,
            'eigenvalues': deepcopy(eigenvalues)}


def assemble(q, groups, mean_zero=False):
    q = _q(q); scope = _scope(mean_zero)
    if not isinstance(groups, list) or len(groups) > 17:
        raise ValueError('At most 17 consecutive angular sectors are supported')
    entries = []
    for m, group in enumerate(groups):
        rebuilt = sector_bundle(q, m, mean_zero, group['eigenvalues'])
        if not same(group, rebuilt):
            raise ValueError('Missing, reordered or scope-swapped angular sector')
        for j, cert in enumerate(group['eigenvalues'], 1):
            entries.append({'m': m, 'j': j, 'l': _start(m, mean_zero)+j-1,
                            'multiplicity': 1 if m == 0 else 2,
                            'lower': cert['lower'], 'upper': cert['upper']})
    return {'format': MODEL, 'q': q, **scope, 'groups': deepcopy(groups), 'entries': entries}


def _model(model):
    rebuilt = assemble(model['q'], model['groups'], model['mean_zero'])
    if not same(model, rebuilt):
        raise ValueError('Model metadata or eigenvalue multiplicity changed')
    return rebuilt


def _slots(model, cutoff):
    model = _model(model); cutoff = _r(cutoff)
    bounds = _range(model['q']); lo, hi = F(bounds['lower']), F(bounds['upper'])
    minimum = 1 if model['mean_zero'] else 0
    last = minimum-1
    while (last+1)*(last+2)+lo <= cutoff:
        last += 1
        if last > 128:
            raise ValueError('Threshold requires more than the supported 128 angular degrees')
    known = {(e['m'], e['j']): e for e in model['entries']}
    keys = set(known)
    for l in range(minimum, last+1):
        for m in range(l+1):
            keys.add((m, l-_start(m, model['mean_zero'])+1))
    slots = []
    for m, j in sorted(keys):
        l = _start(m, model['mean_zero'])+j-1
        lower, upper = l*(l+1)+lo, l*(l+1)+hi
        missing = (m, j) not in known
        if not missing:
            lower = max(lower, F(known[(m, j)]['lower']))
            upper = min(upper, F(known[(m, j)]['upper']))
        if lower > upper:
            raise ArithmeticError('Certified sector and min-max intervals disagree')
        slots.append({'m': m, 'j': j, 'l': l, 'multiplicity': 1 if m == 0 else 2,
                      'lower': str(lower), 'upper': str(upper), 'uncomputed': missing})
    radial = []
    for m, group in enumerate(model['groups']):
        l = _start(m, model['mean_zero'])+len(group['eigenvalues'])
        radial.append({'m': m, 'next_l': l, 'lower': str(l*(l+1)+lo)})
    next_m = len(model['groups'])
    first_l = _start(next_m, model['mean_zero'])
    tails = {'potential_range': bounds, 'radial_floors': radial,
             'next_uncomputed_m': next_m, 'angular_floor': str(first_l*(first_l+1)+lo),
             'degree_cutoff': last, 'cutoff': str(cutoff),
             'every_unlisted_eigenvalue_strictly_above_cutoff': True,
             'unlisted_floor': str((last+1)*(last+2)+lo),
             'rule': 'sector minmax: l(l+1)+inf(q) <= lambda_m,j <= l(l+1)+sup(q)'}
    return slots, tails


def ordered_merge(model, k=9, tolerance='1/100000000'):
    """V207: certified full-space merge including all omitted-mode floors."""
    model = _model(model); tol = _tolerance(tolerance)
    initial = analytic_bounds(model['q'], k, model['mean_zero'])
    universal = [F(e['upper']) for e in initial['ordered']]
    slots, tails = _slots(model, universal[-1])
    low, high = [], []
    for s in slots:
        low.extend([F(s['lower'])]*s['multiplicity'])
        high.extend([F(s['upper'])]*s['multiplicity'])
    low.sort(); high.sort()
    intervals = []
    for i in range(k):
        upper = min(high[i], universal[i]); lower = low[i]
        if lower > upper:
            raise ArithmeticError('Ordered interval is empty')
        intervals.append({'index': i+1, 'lower': str(lower), 'upper': str(upper),
                          'width': str(upper-lower)})
    return {'format': ORDERED, 'model': model, 'k': k, 'tolerance': str(tol),
            'ordered': intervals, 'slots': slots, 'tails': tails,
            'target_met': all(F(x['width']) <= tol for x in intervals),
            'formal_assistant_checked': False}


def strict_count(model, threshold=0):
    """V208: strict counting includes l0 when appropriate, excludes all ties."""
    model = _model(model); threshold = _r(threshold)
    slots, tails = _slots(model, threshold)
    lower = sum(s['multiplicity'] for s in slots if F(s['upper']) < threshold)
    upper = sum(s['multiplicity'] for s in slots if F(s['lower']) < threshold)
    ties = sum(s['multiplicity'] for s in slots if F(s['lower']) == threshold == F(s['upper']))
    unresolved = [s for s in slots if F(s['lower']) < threshold <= F(s['upper'])]
    return {'format': COUNT, 'model': model, 'threshold': str(threshold),
            'comparison': 'lambda < threshold', 'count_lower': lower, 'count_upper': upper,
            'exact_tie_multiplicity': ties, 'unresolved': unresolved,
            'slots': slots, 'tails': tails, 'target_met': lower == upper,
            'formal_assistant_checked': False}


def negative_trace(model, tolerance='1/100000000'):
    """V209: complete negative-part trace, with explicit omitted contribution."""
    model = _model(model); tol = _tolerance(tolerance)
    slots, tails = _slots(model, 0)
    def contribution(items, endpoint):
        return sum((s['multiplicity']*max(F(0), -F(s[endpoint])) for s in items), F(0))
    lower, upper = contribution(slots, 'upper'), contribution(slots, 'lower')
    omitted = [s for s in slots if s['uncomputed'] and F(s['lower']) < 0]
    return {'format': TRACE, 'model': model, 'tolerance': str(tol),
            'lower': str(lower), 'upper': str(upper), 'width': str(upper-lower),
            'omitted_negative_trace_lower': str(contribution(omitted, 'upper')),
            'omitted_negative_trace_upper': str(contribution(omitted, 'lower')),
            'slots': slots, 'tails': tails, 'target_met': upper-lower <= tol,
            'definition': 'sum_with_multiplicity max(0,-lambda)',
            'formal_assistant_checked': False}


def request(q, k=9, mean_zero=False, quantities=('ordered', 'count', 'trace'),
            threshold=0, tolerance='1/100000000', modes=6, max_modes=16,
            max_m=6, max_radial=8, bits=36, max_steps=128):
    q = _q(q); scope = _scope(mean_zero); _i(k, 1, 128, 'k')
    if not isinstance(quantities, (list, tuple)) or not quantities or len(set(quantities)) != len(quantities):
        raise ValueError('Distinct requested quantities are required')
    if any(x not in ('ordered', 'count', 'trace') for x in quantities):
        raise ValueError('Unknown requested spectral quantity')
    _i(modes, 2, 32, 'modes'); _i(max_modes, modes, 32, 'max_modes')
    _i(max_m, 0, 16, 'max_m'); _i(max_radial, 1, max_modes, 'max_radial')
    _i(bits, 8, 100, 'bits'); _i(max_steps, 0, 1024, 'max_steps')
    return {'q': q, **scope, 'k': k, 'quantities': sorted(quantities),
            'threshold': str(_r(threshold)), 'tolerance': str(_tolerance(tolerance)),
            'modes': modes, 'max_modes': max_modes, 'max_m': max_m,
            'max_radial': max_radial, 'bits': bits, 'max_steps': max_steps}


def _request(raw):
    result = request(raw['q'], raw['k'], raw['mean_zero'], raw['quantities'], raw['threshold'],
        raw['tolerance'], raw['modes'], raw['max_modes'], raw['max_m'], raw['max_radial'],
        raw['bits'], raw['max_steps'])
    if not same(raw, result):
        raise ValueError('Noncanonical request or changed operator scope')
    return result


def multiquantity(req, model):
    req = _request(req); model = _model(model)
    if req['q'] != model['q'] or req['mean_zero'] is not model['mean_zero']:
        raise ValueError('Request and model have different q or projection')
    answers = {}
    for quantity in req['quantities']:
        if quantity == 'ordered': answers[quantity] = ordered_merge(model, req['k'], req['tolerance'])
        if quantity == 'count': answers[quantity] = strict_count(model, req['threshold'])
        if quantity == 'trace': answers[quantity] = negative_trace(model, req['tolerance'])
    decisions = {name: 'certified_target_met' if proof['target_met'] else 'certified_bound_open'
                 for name, proof in answers.items()}
    return {'format': MULTI, 'request': req, 'model': model, 'answers': answers,
            'decisions': decisions, 'all_requested_targets_met': all(p['target_met'] for p in answers.values()),
            'formal_assistant_checked': False}


def _priority_candidates(cert):
    """Weight unresolved modes by their contribution to requested quantities."""
    scores = {}
    tol = F(cert['request']['tolerance'])
    for name, proof in cert['answers'].items():
        if proof['target_met']:
            continue
        for slot in proof['slots']:
            lo, hi = F(slot['lower']), F(slot['upper'])
            if lo == hi:
                continue
            score = F(0)
            if name == 'trace':
                score = slot['multiplicity']*(max(F(0), -lo)-max(F(0), -hi))/tol
            elif name == 'count':
                if lo < F(proof['threshold']) <= hi:
                    score = F(slot['multiplicity'], 1)/tol
            elif lo <= F(proof['ordered'][-1]['upper']):
                score = (hi-lo)/tol
            if score:
                key = slot['m'], slot['j']
                scores[key] = scores.get(key, F(0))+score
    return sorted(scores, key=lambda x: (-scores[x], x))


def adaptive_spectrum(q, **options):
    """V210: scheduling stops only when ALL requested quantities close.

    In particular a sharp first eigenvalue does not terminate a requested
    negative trace whose other radial/angular contributions remain open.
    """
    req = request(q, **options); groups = []
    cert = multiquantity(req, assemble(req['q'], groups, req['mean_zero']))
    history = [{'step': 0, 'decisions': cert['decisions'], 'action': 'analytic_initial_bounds'}]
    for step in range(1, req['max_steps']+1):
        if cert['all_requested_targets_met']:
            break
        action = None
        for target_m, target_j in _priority_candidates(cert):
            if target_m > req['max_m'] or target_j > req['max_radial']:
                continue
            # Fill the necessary consecutive prefix, never forge an index gap.
            m = min(target_m, len(groups))
            if m == len(groups):
                j = 1; size = req['modes']; action = (m, j, size, 'new_angular_sector'); break
            group = groups[m]['eigenvalues']
            j = min(target_j, len(group)+1)
            if j > len(group):
                size = max(req['modes'], j)
                action = (m, j, size, 'new_radial_eigenvalue'); break
            current_size = group[j-1]['modes']
            if current_size < req['max_modes']:
                size = min(req['max_modes'], current_size*2)
                action = (m, j, size, 'refine_radial_certificate'); break
        if action is None:
            break
        m, j, size, reason = action
        proof = projected.certify_sector(req['q'], m=m, mean_zero=req['mean_zero'],
                                        k=j, modes=size, bits=req['bits'])
        values = [] if m == len(groups) else list(groups[m]['eigenvalues'])
        if j > len(values): values.append(proof)
        else: values[j-1] = proof
        group = sector_bundle(req['q'], m, req['mean_zero'], values)
        if m == len(groups): groups.append(group)
        else: groups[m] = group
        cert = multiquantity(req, assemble(req['q'], groups, req['mean_zero']))
        history.append({'step': step, 'action': reason, 'm': m, 'j': j, 'modes': size,
                        'decisions': cert['decisions']})
    return {'certificate': cert, 'search_history': history,
            'budget_or_precision_exhausted': not cert['all_requested_targets_met'],
            'sector_solves': len(history)-1}


def spectral_shift(model, threshold):
    """V211: preserve projection under q->q-T and compute shifted count/trace."""
    model = _model(model); threshold = _r(threshold)
    shifted = list(map(F, model['q'])); shifted[0] -= threshold
    shifted = _q(shifted)
    count = strict_count(model, threshold)
    slots, tails = _slots(model, threshold)
    lower = sum((s['multiplicity']*max(F(0), threshold-F(s['upper'])) for s in slots), F(0))
    upper = sum((s['multiplicity']*max(F(0), threshold-F(s['lower'])) for s in slots), F(0))
    return {'format': SHIFT, 'model': model, 'threshold': str(threshold),
            'shifted_q': shifted, 'shifted_mean_zero': model['mean_zero'],
            'strict_negative_count_for_shifted_operator': count,
            'shifted_negative_trace_lower': str(lower), 'shifted_negative_trace_upper': str(upper),
            'slots': slots, 'tails': tails,
            'rule': 'P(H-T*I)P=PHP-T*I_on_range(P); projection remains unchanged',
            'formal_assistant_checked': False}


def fixed_potential_quotient(trace, potential):
    """V212: fixed-V quotient for either scope, with the constant-mode obstruction."""
    if trace.get('format') != TRACE or not verify(trace):
        raise ValueError('A verified complete negative trace is required')
    v = sectors.potential(potential)
    if _q([-x for x in v]) != trace['model']['q']:
        raise ValueError('The trace must use exactly q=-V')
    bounds = _range(v)
    if F(bounds['lower']) < 0:
        raise ValueError('V is not certified nonnegative')
    square = arithmetic.mul(v, v)
    mean_square = sum((x/F(i+1) for i, x in enumerate(square) if i % 2 == 0), F(0))
    if mean_square <= 0:
        raise ValueError('V must be nonzero')
    divisor = 4*mean_square
    lower, upper = F(trace['lower'])/divisor, F(trace['upper'])/divisor
    mean_zero = trace['model']['mean_zero']
    obstruction = None if mean_zero else {
        'potential': 'V=epsilon constant with 0<epsilon<2',
        'negative_count': 1, 'negative_trace': 'epsilon',
        'pi_scaled_trace_quotient': '1/(4*epsilon)',
        'limit_as_epsilon_to_zero_positive': 'positive_infinity',
        'reason': 'the unprojected constant spherical harmonic is an eigenfunction'}
    return {'format': LT, 'trace': trace, 'V': list(map(str, v)), 'nonnegative_range': bounds,
            'mean_V_squared': str(mean_square), 'rational_divisor': str(divisor),
            'quantity': 'pi*negative_trace/integral_S2(V^2)',
            'area_normalization': 'integral_S2(V^2)=4*pi*mean(V^2)',
            'lower': str(lower), 'upper': str(upper),
            'claim': 'fixed_potential_ratio_only',
            'finite_universal_L2_trace_inequality_compatible_with_scope': mean_zero,
            'unprojected_constant_mode_obstruction': obstruction,
            'not_a_universal_upper_bound': True, 'formal_assistant_checked': False}


def projection_interlacing(full, constrained, k=None):
    """V213: codimension-one min--max interlacing for the SAME potential."""
    if full.get('format') != ORDERED or constrained.get('format') != ORDERED:
        raise ValueError('Two complete ordered spectral proofs are required')
    if not verify(full, expected_mean_zero=False) or not verify(constrained, expected_mean_zero=True):
        raise ValueError('Interlacing requires verified unprojected and mean-zero proofs')
    if full['model']['q'] != constrained['model']['q']:
        raise ValueError('Interlacing must use exactly the same potential')
    maximum = min(full['k']-1, constrained['k'])
    k = maximum if k is None else _i(k, 1, maximum, 'k')
    if k < 1:
        raise ValueError('At least two unprojected eigenvalue bounds are required')
    rows = []
    for i in range(k):
        a, b, c = full['ordered'][i], constrained['ordered'][i], full['ordered'][i+1]
        lower = max(F(a['lower']), F(b['lower']))
        upper = min(F(c['upper']), F(b['upper']))
        if lower > upper or F(a['lower']) > F(b['upper']) or F(b['lower']) > F(c['upper']):
            raise ArithmeticError('Numerical bounds contradict codimension-one interlacing')
        rows.append({'index': i+1, 'full_k': a, 'projected_k': b, 'full_k_plus_1': c,
                     'tightened_projected_lower': str(lower), 'tightened_projected_upper': str(upper)})
    return {'format': INTERLACE, 'full': full, 'projected': constrained, 'k': k, 'rows': rows,
            'relation': 'lambda_k(full) <= lambda_k(mean_zero) <= lambda_(k+1)(full)',
            'proof_rule': 'minmax restriction to the codimension-one orthogonal complement of constants',
            'nonconstant_q_does_not_have_an_added_constant_eigenvalue': True,
            'formal_assistant_checked': False}


def driver_certificate(requests, results):
    if not isinstance(requests, list) or not 1 <= len(requests) <= 2 or len(requests) != len(results):
        raise ValueError('One or two scope requests and results required')
    if len({r['mean_zero'] for r in requests}) != len(requests):
        raise ValueError('Repeated scope request')
    if len({tuple(_q(r['q'])) for r in requests}) != 1:
        raise ValueError('Both spaces must use the same potential')
    for req, result in zip(requests, results):
        _request(req)
        if not isinstance(result, dict) or result.get('format') != MULTI or not verify(result) or not same(result['request'], req):
            raise ValueError('Result or cached certificate is bound to different inputs')
    decisions = {r['request']['space']: r['decisions'] for r in results}
    return {'format': DRIVER, 'requests': requests, 'results': results, 'decisions': decisions,
            'all_requested_targets_met': all(x['all_requested_targets_met'] for x in results),
            'formal_assistant_checked': False}


def research_driver(q, spaces=('full_mean_zero', 'full_unprojected'), cache=None, **options):
    """V214: requested-quantity decisions and verified full-input cache reuse."""
    if not isinstance(spaces, (list, tuple)) or not 1 <= len(spaces) <= 2 or len(set(spaces)) != len(spaces):
        raise ValueError('One or two distinct spaces are required')
    if any(s not in ('full_mean_zero', 'full_unprojected') for s in spaces):
        raise ValueError('Unknown function space')
    if 'mean_zero' in options:
        raise ValueError('Use spaces to select the projections')
    if cache is None: cache = {}
    if not isinstance(cache, dict):
        raise ValueError('Cache must be a dictionary')
    requests, results, histories = [], [], []
    hits, rejects, solves = 0, 0, 0
    for space in spaces:
        mean_zero = space == 'full_mean_zero'
        req = request(q, mean_zero=mean_zero, **options); key = digest(req)
        previous = cache.get(key)
        if isinstance(previous, dict) and previous.get('format') == MULTI and verify(previous) and same(previous['request'], req):
            proof = deepcopy(previous); hits += 1; history = []
        else:
            if previous is not None: rejects += 1
            solved = adaptive_spectrum(q, mean_zero=mean_zero, **options)
            proof = solved['certificate']; history = solved['search_history']; solves += solved['sector_solves']
            cache[key] = deepcopy(proof)
        requests.append(req); results.append(proof); histories.append(history)
    return {'certificate': driver_certificate(requests, results),
            'execution': {'cache_hits': hits, 'cache_rejections': rejects, 'sector_solves': solves},
            'search_histories': histories}


def _proof_inputs(cert):
    fmt = cert['format']
    if fmt in (INITIAL, SECTOR, MODEL): return cert['q'], cert['mean_zero'], cert.get('k')
    if fmt in (ORDERED, COUNT, TRACE, SHIFT): return cert['model']['q'], cert['model']['mean_zero'], cert.get('k')
    if fmt == MULTI: return cert['request']['q'], cert['request']['mean_zero'], cert['request']['k']
    if fmt == LT: return _proof_inputs(cert['trace'])
    if fmt == INTERLACE: return cert['full']['model']['q'], None, cert['k']
    if fmt == DRIVER: return cert['requests'][0]['q'], None, None
    raise ValueError('Unknown certificate format')


def verify(cert, expected_q=None, expected_mean_zero=None, expected_k=None):
    """Rebuild every mathematical field; bind optional external inputs."""
    try:
        q, mean_zero, k = _proof_inputs(cert)
        if expected_q is not None and _q(expected_q) != q: return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool or expected_mean_zero is not mean_zero): return False
        if expected_k is not None and (type(expected_k) is not int or expected_k != k): return False
        fmt = cert['format']
        if fmt == INITIAL: rebuilt = analytic_bounds(q, cert['k'], mean_zero)
        elif fmt == SECTOR: rebuilt = sector_bundle(q, cert['m'], mean_zero, cert['eigenvalues'])
        elif fmt == MODEL: rebuilt = assemble(q, cert['groups'], mean_zero)
        elif fmt == ORDERED: rebuilt = ordered_merge(cert['model'], cert['k'], cert['tolerance'])
        elif fmt == COUNT: rebuilt = strict_count(cert['model'], cert['threshold'])
        elif fmt == TRACE: rebuilt = negative_trace(cert['model'], cert['tolerance'])
        elif fmt == MULTI: rebuilt = multiquantity(cert['request'], cert['model'])
        elif fmt == SHIFT: rebuilt = spectral_shift(cert['model'], cert['threshold'])
        elif fmt == LT: rebuilt = fixed_potential_quotient(cert['trace'], cert['V'])
        elif fmt == INTERLACE: rebuilt = projection_interlacing(cert['full'], cert['projected'], cert['k'])
        elif fmt == DRIVER: rebuilt = driver_certificate(cert['requests'], cert['results'])
        else: return False
        return same(cert, rebuilt)
    except (ValueError, TypeError, KeyError, IndexError, ArithmeticError, OverflowError):
        return False
