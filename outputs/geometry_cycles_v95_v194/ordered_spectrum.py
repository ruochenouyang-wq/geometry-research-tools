"""V155--V164: ordered full sphere spectra and negative spectral moments.

The operator is the quadratic-form compression P(-Delta+q)P to mean-zero
real H1(S2). q depends only on t=x3. Angular modes m>0 have multiplicity 2.
Uncomputed radial and angular modes retain independent min--max bounds;
they are never deleted from an order statistic or negative trace.

Certificates are rational software proofs with stated analytical rules,
not formal-assistant proofs. The LT operation gives a fixed-potential
witness, never a universal upper bound for a best constant.
"""
from copy import deepcopy
from common import F, sectors, projected, arithmetic, exact, same, save

REQUEST = 'full_mean_zero_ordered_request_v155'
INITIAL = 'full_ordered_initial_minmax_v155'
ASSEMBLY = 'angular_multiplicity_assembly_v157'
FLOORS = 'omitted_sector_floors_v158'
ORDERED = 'full_ordered_spectrum_v159'
COUNT = 'strict_full_spectral_count_v161'
TRACE = 'full_negative_trace_v162'
LT = 'fixed_potential_LT_witness_v163'
BUNDLE = 'full_spectral_research_bundle_v164'


def _rational(value):
    if type(value) not in (int, str, F):
        raise ValueError('An exact integer or rational string is required')
    return exact(value)


def _int(value, lo, hi, label):
    if type(value) is not int or not lo <= value <= hi:
        raise ValueError('%s must be an integer in %d..%d' % (label, lo, hi))
    return value


def normalize_request(q, k=8, mean_zero=True, tolerance='1/100000000',
                      modes=8, max_modes=16, max_m=4, max_radial=6, bits=36):
    """V155: canonical scope, input validation and resource limits."""
    q = sectors.potential(q)
    if mean_zero is not True:
        raise ValueError('This full spectrum interface requires mean_zero=True')
    _int(k, 1, 128, 'k')
    _int(modes, 2, 32, 'modes')
    _int(max_modes, modes, 32, 'max_modes')
    _int(max_m, 0, 16, 'max_m')
    _int(max_radial, 1, max_modes, 'max_radial')
    _int(bits, 8, 100, 'bits')
    tol = _rational(tolerance)
    if not F(1, 10**25) <= tol <= 1:
        raise ValueError('Tolerance must be in [1e-25,1]')
    return {'format': REQUEST, 'geometry': 'unit_S2',
            'operator': 'P_mean_zero*(-Delta+q)*P_mean_zero',
            'function_space': 'all_real_mean_zero_H1_S2',
            'q_coefficients': list(map(str, q)), 'k': k, 'mean_zero': True,
            'multiplicity': {'m=0': 1, 'm>0': 2},
            'tolerance': str(tol), 'modes': modes, 'max_modes': max_modes,
            'max_m': max_m, 'max_radial': max_radial, 'bits': bits}


def _request(raw):
    expected = normalize_request(raw['q_coefficients'], raw['k'], raw['mean_zero'],
                                 raw['tolerance'], raw['modes'], raw['max_modes'],
                                 raw['max_m'], raw['max_radial'], raw['bits'])
    if not same(raw, expected):
        raise ValueError('Noncanonical request or changed scope')
    return expected


def laplace_spectrum(k=8, shift=0):
    """V156: exact first k levels of -Delta+constant, with multiplicities."""
    _int(k, 1, 128, 'k')
    shift = _rational(shift)
    levels = []
    l = 1
    while len(levels) < k:
        levels.extend([str(l*(l+1)+shift)]*(2*l+1))
        l += 1
    return {'format': 'exact_mean_zero_laplace_calibrator_v156',
            'mean_zero': True, 'removed_l0_multiplicity': 1,
            'shift': str(shift), 'k': k, 'eigenvalues': levels[:k]}


def initial_order_bounds(request):
    """V155: full ordered min--max bounds without any Schur computation.

    The complete mean-zero Laplace spectrum carries the real multiplicity
    2*l+1. A bounded multiplication form shifts every ordered eigenvalue by
    an amount between inf(q) and sup(q), including repeated eigenvalues.
    """
    request = _request(request)
    q = sectors.potential(request['q_coefficients'])
    proof = projected.v05_range.make_range(q)
    qlo, qhi = projected.v05_range.verify_bounds(q, proof)
    laplace = laplace_spectrum(request['k'])
    intervals = []
    for i, raw in enumerate(laplace['eigenvalues'], 1):
        lower, upper = F(raw)+qlo, F(raw)+qhi
        intervals.append({'index': i, 'laplace_eigenvalue': raw,
                          'lower': str(lower), 'upper': str(upper),
                          'width': str(upper-lower),
                          'target_met': upper-lower <= F(request['tolerance'])})
    return {'format': INITIAL, 'request': request,
            'potential_range_proof': proof, 'q_lower': str(qlo), 'q_upper': str(qhi),
            'ordered_intervals': intervals,
            'rule': 'lambda_k(-Delta)+inf(q) <= lambda_k(P*(-Delta+q)*P) <= lambda_k(-Delta)+sup(q)',
            'projection': 'mean_zero', 'laplace_multiplicity': '2*l+1, l>=1',
            'schur_computations': 0,
            'status': 'certified_target_met' if all(x['target_met'] for x in intervals)
                      else 'certified_bound_open_gap',
            'formal_assistant_checked': False}


def verify_initial(cert, expected_q=None, expected_k=None, expected_mean_zero=True):
    """Replay initial bounds and optionally bind all external problem inputs."""
    try:
        if cert.get('format') != INITIAL:
            return False
        request = cert['request']
        if expected_q is not None and list(map(str, sectors.potential(expected_q))) != request['q_coefficients']:
            return False
        if expected_k is not None and (type(expected_k) is not int or expected_k != request['k']):
            return False
        if expected_mean_zero is not True or request['mean_zero'] is not True:
            return False
        return same(cert, initial_order_bounds(request))
    except (ValueError, TypeError, KeyError, IndexError, ArithmeticError, OverflowError):
        return False


def assemble_sectors(q, sector_lists):
    """V157: verify consecutive sectors/radial indices and angular weights."""
    q = sectors.potential(q)
    if not isinstance(sector_lists, list) or len(sector_lists) > 17:
        raise ValueError('A list of consecutive m=0..M sectors is required')
    entries = []
    for m, group in enumerate(sector_lists):
        if not isinstance(group, list) or not 1 <= len(group) <= 32:
            raise ValueError('Each sector needs consecutive radial indices')
        for j, cert in enumerate(group, 1):
            if not projected.verify_sector(cert, expected_q=q, expected_m=m,
                    expected_mean_zero=True, expected_k=j, independent=True):
                raise ValueError('Sector certificate has wrong q, m, k or projection')
            entries.append({'m': m, 'radial_index': j,
                            'laplace_degree': max(1, m)+j-1,
                            'multiplicity': 1 if m == 0 else 2,
                            'lower': cert['lower'], 'upper': cert['upper']})
    return {'format': ASSEMBLY, 'q_coefficients': list(map(str, q)),
            'mean_zero': True, 'sector_lists': deepcopy(sector_lists),
            'entries': entries,
            'represented_eigenvalues_with_multiplicity': sum(e['multiplicity'] for e in entries)}


def _assembly(raw):
    expected = assemble_sectors(raw['q_coefficients'], raw['sector_lists'])
    if not same(raw, expected):
        raise ValueError('Assembly or angular multiplicity altered')
    return expected


def omitted_floors(assembly):
    """V158: independent min--max floors for both unenumerated directions."""
    assembly = _assembly(assembly)
    q = sectors.potential(assembly['q_coefficients'])
    proof = projected.v05_range.make_range(q)
    lo, hi = projected.v05_range.verify_bounds(q, proof)
    radial = []
    for m, group in enumerate(assembly['sector_lists']):
        l = max(1, m)+len(group)
        radial.append({'m': m, 'next_radial_index': len(group)+1,
                       'next_laplace_degree': l, 'multiplicity': 1 if m == 0 else 2,
                       'all_unenumerated_radial_eigenvalues_at_least': str(l*(l+1)+lo)})
    next_m = len(assembly['sector_lists'])
    first_l = max(1, next_m)
    return {'format': FLOORS, 'q_coefficients': assembly['q_coefficients'],
            'potential_range_proof': proof, 'q_lower': str(lo), 'q_upper': str(hi),
            'radial': radial, 'omitted_azimuth_m_start': next_m,
            'all_unenumerated_angular_eigenvalues_at_least': str(first_l*(first_l+1)+lo),
            'rule': 'minmax: l(l+1)+inf(q) <= lambda_m,j <= l(l+1)+sup(q)'}


def _slots(assembly, cutoff, floors):
    """All modes whose universal floor is <= cutoff, plus enumerated modes.

    Every absent mode beyond these slots is strictly above cutoff. The
    finite cutoff is derived from the requested order statistic/threshold,
    never from a matrix-size guess.
    """
    cutoff = F(cutoff)
    lo, hi = F(floors['q_lower']), F(floors['q_upper'])
    known = {(e['m'], e['radial_index']): e for e in assembly['entries']}
    last = 0
    while (last+1)*(last+2)+lo <= cutoff:
        last += 1
        if last > 512:
            raise ValueError('Analytical cutoff exceeds supported 512 degrees')
    keys = set(known)
    for l in range(1, last+1):
        for m in range(l+1):
            keys.add((m, l-max(1, m)+1))
    result = []
    for m, j in sorted(keys):
        l = max(1, m)+j-1
        lower, upper = l*(l+1)+lo, l*(l+1)+hi
        source = 'uncomputed_minmax_placeholder'
        if (m, j) in known:
            e = known[(m, j)]
            lower, upper = max(lower, F(e['lower'])), min(upper, F(e['upper']))
            source = 'certified_sector_and_minmax'
        if lower > upper:
            raise ArithmeticError('Sector and min--max enclosures disagree')
        result.append({'m': m, 'radial_index': j, 'laplace_degree': l,
                       'multiplicity': 1 if m == 0 else 2,
                       'lower': str(lower), 'upper': str(upper), 'source': source})
    return result, {'degree_cutoff': last,
                    'every_unlisted_mode_above_cutoff': True,
                    'unlisted_floor': str((last+1)*(last+2)+lo),
                    'cutoff': str(cutoff)}


def ordered_bounds(request, assembly):
    """V159: rigorous first k order statistics with missing-mode placeholders."""
    request, assembly = _request(request), _assembly(assembly)
    if request['q_coefficients'] != assembly['q_coefficients']:
        raise ValueError('Request and spectral evidence use different potentials')
    floors = omitted_floors(assembly)
    universal = list(map(F, laplace_spectrum(request['k'], floors['q_upper'])['eigenvalues']))
    slots, cutoff = _slots(assembly, universal[-1], floors)
    lows, highs = [], []
    for slot in slots:
        lows.extend([F(slot['lower'])]*slot['multiplicity'])
        highs.extend([F(slot['upper'])]*slot['multiplicity'])
    lows.sort()
    highs.sort()
    if len(lows) < request['k']:
        raise ArithmeticError('Universal cutoff failed to cover the kth eigenvalue')
    intervals = []
    for i in range(request['k']):
        lower, upper = lows[i], min(highs[i], universal[i])
        if lower > upper:
            raise ArithmeticError('Invalid ordered spectral bounds')
        intervals.append({'index': i+1, 'lower': str(lower), 'upper': str(upper),
                          'width': str(upper-lower),
                          'target_met': upper-lower <= F(request['tolerance'])})
    kth_upper = F(intervals[-1]['upper'])
    uncomputed = [s for s in slots if s['source'] == 'uncomputed_minmax_placeholder'
                  and F(s['lower']) <= kth_upper]
    all_met = all(x['target_met'] for x in intervals)
    return {'format': ORDERED, 'request': request, 'assembly': assembly,
            'floors': floors, 'slots': slots, 'cutoff_proof': cutoff,
            'ordered_intervals': intervals,
            'uncomputed_modes_may_touch_first_k': uncomputed,
            'first_k_resolved_without_enumerating_all_modes': all_met and bool(uncomputed),
            'status': 'certified_target_met' if all_met else 'certified_bound_open_gap',
            'analytical_dependencies': ['complete real spherical harmonic decomposition',
                'quadratic form minmax under bounded multiplication',
                'verified sector radial Schur comparison'],
            'formal_assistant_checked': False}


def adaptive_spectrum(q, k=8, **options):
    """V160: grow the needed angular/radial list; exhaustion returns safe bounds."""
    request = normalize_request(q, k=k, **options)
    assembly = assemble_sectors(request['q_coefficients'], [])
    current = ordered_bounds(request, assembly)
    history = [{'sectors': 0, 'radial_eigenvalues': 0, 'status': current['status']}]
    groups = []
    # Constants can already be exact by min--max; deliberately retain a
    # first sector certificate so the research bundle also tests the kernel.
    for m in range(request['max_m']+1):
        kth_upper = F(current['ordered_intervals'][-1]['upper'])
        qlo = F(current['floors']['q_lower'])
        first_l = max(1, m)
        if m and first_l*(first_l+1)+qlo > kth_upper:
            break
        group = []
        for j in range(1, request['max_radial']+1):
            l = first_l+j-1
            if j > 1 and l*(l+1)+qlo > kth_upper:
                break
            size = max(request['modes'], j)
            while True:
                cert = projected.certify_sector(request['q_coefficients'], m=m,
                    mean_zero=True, k=j, modes=size, bits=request['bits'])
                if F(cert['exact_width']) <= F(request['tolerance']) or size == request['max_modes']:
                    break
                size = min(request['max_modes'], 2*size)
            group.append(cert)
        groups.append(group)
        assembly = assemble_sectors(request['q_coefficients'], groups)
        current = ordered_bounds(request, assembly)
        history.append({'sectors': len(groups), 'radial_eigenvalues': sum(map(len, groups)),
                        'status': current['status']})
        if current['status'] == 'certified_target_met':
            break
    # History is exploratory metadata, stored separately from proof authority.
    return {'certificate': current, 'search_history': history,
            'budget_exhausted_before_target': current['status'] != 'certified_target_met'}


def count_below(assembly, threshold=0):
    """V161: count lambda<threshold. Equal endpoints are NOT negative."""
    assembly = _assembly(assembly)
    threshold = _rational(threshold)
    floors = omitted_floors(assembly)
    slots, cutoff = _slots(assembly, threshold, floors)
    lower = sum(s['multiplicity'] for s in slots if F(s['upper']) < threshold)
    upper = sum(s['multiplicity'] for s in slots if F(s['lower']) < threshold)
    unresolved = [s for s in slots if F(s['lower']) < threshold <= F(s['upper'])]
    ties = [s for s in slots if F(s['lower']) == threshold == F(s['upper'])]
    return {'format': COUNT, 'assembly': assembly, 'threshold': str(threshold),
            'comparison': 'strictly_less_than', 'floors': floors, 'slots': slots,
            'cutoff_proof': cutoff, 'count_lower': lower, 'count_upper': upper,
            'exact_tie_multiplicity': sum(s['multiplicity'] for s in ties),
            'unresolved_threshold_clusters': unresolved,
            'status': 'exact_count' if lower == upper else 'count_open',
            'formal_assistant_checked': False}


def negative_trace(assembly, tolerance='1/100000000'):
    """V162: sum max(0,-lambda), including bounded omitted negative modes."""
    assembly = _assembly(assembly)
    tol = _rational(tolerance)
    if not F(1, 10**25) <= tol <= 1:
        raise ValueError('Tolerance must be in [1e-25,1]')
    floors = omitted_floors(assembly)
    slots, cutoff = _slots(assembly, 0, floors)
    lower = sum((s['multiplicity']*max(F(0), -F(s['upper'])) for s in slots), F(0))
    upper = sum((s['multiplicity']*max(F(0), -F(s['lower'])) for s in slots), F(0))
    omitted = [s for s in slots if s['source'] == 'uncomputed_minmax_placeholder'
               and F(s['lower']) < 0]
    tail_lower = sum((s['multiplicity']*max(F(0), -F(s['upper'])) for s in omitted), F(0))
    tail_upper = sum((s['multiplicity']*max(F(0), -F(s['lower'])) for s in omitted), F(0))
    return {'format': TRACE, 'assembly': assembly, 'floors': floors, 'slots': slots,
            'cutoff_proof': cutoff, 'sum_convention': 'sum_multiplicity*max(0,-lambda)',
            'lower': str(lower), 'upper': str(upper), 'width': str(upper-lower),
            'tolerance': str(tol), 'omitted_negative_modes': omitted,
            'omitted_negative_trace_lower': str(tail_lower),
            'omitted_negative_trace_upper': str(tail_upper),
            'omitted_negative_tail_resolved': tail_lower == tail_upper,
            'status': 'certified_trace_target_met' if upper-lower <= tol else 'certified_trace_open_gap',
            'formal_assistant_checked': False}


def _mean(poly):
    return sum((x/F(i+1) for i, x in enumerate(poly) if i % 2 == 0), F(0))


def lt_quotient(trace, potential):
    """V163: rational pi-scaled trace quotient for verified V>=0, q=-V."""
    if not verify(trace):
        raise ValueError('A valid full negative trace proof is required')
    if trace.get('format') != TRACE:
        raise ValueError('A negative trace certificate is required')
    v = sectors.potential(potential)
    q = [-x for x in v]
    if list(map(str, q)) != trace['assembly']['q_coefficients']:
        raise ValueError('The trace operator must use q=-V with the same V')
    proof = projected.v05_range.make_range(v)
    lo, _ = projected.v05_range.verify_bounds(v, proof)
    if lo < 0:
        raise ValueError('Cannot certify V>=0 on the sphere')
    mean_square = _mean(arithmetic.mul(v, v))
    if mean_square <= 0:
        raise ValueError('V must be nonzero')
    divisor = 4*mean_square
    lower, upper = F(trace['lower'])/divisor, F(trace['upper'])/divisor
    return {'format': LT, 'potential_coefficients': list(map(str, v)),
            'potential_nonnegative_proof': proof, 'trace_certificate': trace,
            'mean_V_squared': str(mean_square), 'rational_divisor': str(divisor),
            'quantity': 'pi * Tr(P_mean_zero*(-Delta-V)*P_mean_zero)_- / integral_S2(V^2)',
            'measure': 'standard_area; integral(V^2)=4*pi*mean(V^2)',
            'lower': str(lower), 'upper': str(upper), 'width': str(upper-lower),
            'universal_trace_constant_pi_lower_bound': str(lower),
            'claim': 'fixed_potential_witness_lower_bound_only',
            'orthonormal_constant_relation': 'k_LT=4*L_trace in dimension two',
            'not_an_upper_bound_for_any_universal_constant': True,
            'formal_assistant_checked': False}


def bundle(ordered, threshold=0, assertion=None, potential=None):
    """V164: bind all subproofs and a tri-state assertion to the same input."""
    if not verify(ordered) or ordered.get('format') != ORDERED:
        raise ValueError('A valid full ordered spectrum certificate is required')
    threshold = _rational(threshold)
    assembly = ordered['assembly']
    count = count_below(assembly, threshold)
    trace = negative_trace(assembly, ordered['request']['tolerance'])
    result = {'format': BUNDLE, 'ordered': ordered, 'count': count, 'negative_trace': trace,
              'threshold': str(threshold), 'assertion': assertion,
              'assertion_result': None, 'lt_quotient': None,
              'unresolved_threshold_clusters': count['unresolved_threshold_clusters'],
              'formal_assistant_checked': False}
    if assertion is not None:
        if not isinstance(assertion, dict) or set(assertion) != {'index', 'lower_bound'}:
            raise ValueError('Assertion requires exactly index and lower_bound')
        idx = _int(assertion['index'], 1, ordered['request']['k'], 'index')
        bound = _rational(assertion['lower_bound'])
        interval = ordered['ordered_intervals'][idx-1]
        result['assertion'] = {'index': idx, 'lower_bound': str(bound)}
        result['assertion_result'] = ('proved' if F(interval['lower']) >= bound else
                                      'refuted' if F(interval['upper']) < bound else 'undetermined')
    if potential is not None:
        result['lt_quotient'] = lt_quotient(trace, potential)
    return result


def verify(cert, expected_q=None):
    """Strict replay; rejects stale cross-potential, count and multiplicity swaps."""
    try:
        fmt = cert['format']
        if fmt == INITIAL:
            return verify_initial(cert, expected_q=expected_q)
        if fmt == ORDERED:
            rebuilt = ordered_bounds(cert['request'], cert['assembly'])
            q = cert['request']['q_coefficients']
        elif fmt == COUNT:
            rebuilt = count_below(cert['assembly'], cert['threshold'])
            q = cert['assembly']['q_coefficients']
        elif fmt == TRACE:
            rebuilt = negative_trace(cert['assembly'], cert['tolerance'])
            q = cert['assembly']['q_coefficients']
        elif fmt == LT:
            rebuilt = lt_quotient(cert['trace_certificate'], cert['potential_coefficients'])
            q = cert['trace_certificate']['assembly']['q_coefficients']
        elif fmt == BUNDLE:
            potential = None if cert['lt_quotient'] is None else cert['lt_quotient']['potential_coefficients']
            rebuilt = bundle(cert['ordered'], cert['threshold'], cert['assertion'], potential)
            q = cert['ordered']['request']['q_coefficients']
        else:
            return False
        if expected_q is not None and list(map(str, sectors.potential(expected_q))) != q:
            return False
        return same(cert, rebuilt)
    except (ValueError, TypeError, KeyError, IndexError, ArithmeticError, OverflowError):
        return False
