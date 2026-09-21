"""V125--V134: rational continuum certificates for affine potential boxes.

The constrained lowest eigenvalue is an infimum of affine Rayleigh forms:
therefore concave and L-infinity Lipschitz in the potential. Analytic principles
are explicit dependencies, not formal-assistant proofs. Point proofs include
all azimuthal and radial tails through the frozen V94 spectral verifier.
"""
from itertools import product, combinations
from common import F, exact, arithmetic, sectors, projected, exact_extrema, same, save
import witness

PREFIX = 'affine_sphere_family_'


def _strings(v):
    return list(map(str, v))


def _box(raw, dimension):
    if not isinstance(raw, list) or len(raw) != dimension:
        raise ValueError('Box dimension mismatch')
    result = []
    for interval in raw:
        if not isinstance(interval, list) or len(interval) != 2:
            raise ValueError('Each interval requires two rational endpoints')
        lo, hi = map(exact, interval)
        if not lo < hi:
            raise ValueError('Nondegenerate increasing intervals required')
        result.append([lo, hi])
    return result


def validate_family(q0, directions, box, mean_zero=True):
    """V125: canonical affine q(s,t) with one to three parameters."""
    q0 = sectors.potential(q0)
    if not isinstance(directions, list) or not 1 <= len(directions) <= 3:
        raise ValueError('One to three potential directions required')
    directions = [sectors.potential(q) for q in directions]
    if type(mean_zero) is not bool:
        raise ValueError('Boolean mean-zero flag required')
    box = _box(box, len(directions))
    # All real parameter values are covered; the rational solver coefficient
    # representation limit applies only to subsequently sampled anchors.
    return {'format': PREFIX+'representation_v125', 'geometry': 'unit_S2',
            'q0': _strings(q0), 'directions': [_strings(q) for q in directions],
            'box': [_strings(interval) for interval in box], 'mean_zero': mean_zero,
            'function_space': projected.scope(mean_zero), 'parameter_dimension': len(directions),
            'potential_degree_limit': 6, 'parameters_quantifier': 'every_real_point_in_closed_box',
            'formal_assistant_checked': False}


def _family(family):
    rebuilt = validate_family(family['q0'], family['directions'], family['box'], family['mean_zero'])
    if not same(family, rebuilt):
        raise ValueError('Invalid canonical family or scope')
    return family


def _cell(family, cell=None):
    _family(family)
    box = _box(family['box'] if cell is None else cell, family['parameter_dimension'])
    domain = _box(family['box'], len(box))
    if any(lo < a or hi > b for (lo, hi), (a, b) in zip(box, domain)):
        raise ValueError('Cell must be contained in family domain')
    return box


def _point(family, point, cell=None):
    box = _cell(family, cell)
    if not isinstance(point, list) or len(point) != len(box):
        raise ValueError('Point dimension mismatch')
    point = list(map(exact, point))
    if any(not lo <= x <= hi for x, (lo, hi) in zip(point, box)):
        raise ValueError('Point outside closed parameter box')
    return point


def potential_at(family, point):
    point = _point(family, point)
    q = list(map(F, family['q0']))
    for value, direction in zip(point, family['directions']):
        q = arithmetic.add(q, arithmetic.scale(list(map(F, direction)), value))
    return q



def cheap_family_enclosure(family):
    """V125: a search-free enclosure valid for every real parameter in a box.

    Exact coefficient intervals and t^j ranges bound the entire potential.
    The lowest admissible unperturbed harmonic (t or 1) supplies an affine
    Rayleigh upper, maximized over the parameter box by exact sign choices.
    """
    _family(family)
    box = _cell(family)
    q0 = list(map(F, family['q0']))
    directions = [list(map(F, q)) for q in family['directions']]
    degree = max(map(len, [q0]+directions))-1
    coefficient_intervals, monomial_intervals = [], []
    for j in range(degree+1):
        lo = hi = q0[j] if j < len(q0) else F(0)
        for direction, (a, b) in zip(directions, box):
            value = direction[j] if j < len(direction) else F(0)
            lo += min(a*value, b*value)
            hi += max(a*value, b*value)
        coefficient_intervals.append([lo, hi])
        monomial_intervals.append([lo, hi] if j == 0 else
            [-max(abs(lo), abs(hi)), max(abs(lo), abs(hi))] if j % 2 else
            [min(lo, F(0)), max(hi, F(0))])
    qlo = sum((row[0] for row in monomial_intervals), F(0))
    qhi = sum((row[1] for row in monomial_intervals), F(0))
    harmonic_degree = 1 if family['mean_zero'] else 0
    laplace_energy = F(harmonic_degree*(harmonic_degree+1))
    def expectation(q):
        return sum((c*F(2*harmonic_degree+1, j+2*harmonic_degree+1)
                    for j, c in enumerate(q) if j % 2 == 0), F(0))
    intercept = laplace_energy+expectation(q0)
    slopes = list(map(expectation, directions))
    maximizer = [hi if slope >= 0 else lo for slope, (lo, hi) in zip(slopes, box)]
    upper = intercept+sum((a*x for a, x in zip(slopes, maximizer)), F(0))
    lower = laplace_energy+qlo
    return {'format': PREFIX+'cheap_enclosure_v125', 'family': family,
            'coefficient_intervals': [_strings(row) for row in coefficient_intervals],
            'monomial_contribution_intervals': [_strings(row) for row in monomial_intervals],
            'potential_lower': str(qlo), 'potential_upper': str(qhi),
            'unperturbed_spectral_floor': str(laplace_energy),
            'fixed_trial': 't' if harmonic_degree else '1',
            'trial_mass_dt_over_2': str(F(1, 2*harmonic_degree+1)),
            'rayleigh_intercept': str(intercept), 'rayleigh_slopes': _strings(slopes),
            'rayleigh_maximizing_parameter': _strings(maximizer),
            'lower': str(lower), 'upper': str(upper), 'gap': str(upper-lower),
            'claim': 'every parameter has its lowest constrained eigenvalue in [lower,upper]',
            'rules': ['coefficient_box_and_monomial_ranges_on_minus1_plus1',
                      'unperturbed_spectral_floor_plus_bounded_multiplier',
                      'fixed_admissible_harmonic_Rayleigh_affine_box_maximum'],
            'formal_assistant_checked': False}


def direction_norms(family):
    """V126: certified upper bounds for each multiplier L-infinity norm."""
    _family(family)
    proofs = []
    for raw in family['directions']:
        q = list(map(F, raw))
        positive = exact_extrema.maximize(q, tolerance=F(1, 10**12))
        negative = exact_extrema.maximize(arithmetic.scale(q, -1), tolerance=F(1, 10**12))
        proofs.append({'positive': positive, 'negative': negative,
                       'linf_upper': str(max(F(positive['maximum_upper']), F(negative['maximum_upper'])))})
    return {'format': PREFIX+'direction_norms_v126', 'family': family, 'proofs': proofs,
            'norm_upper': [p['linf_upper'] for p in proofs]}


def _verify_norms(cert, family):
    if cert.get('format') != PREFIX+'direction_norms_v126' or not same(cert['family'], family):
        return False
    if len(cert['proofs']) != len(family['directions']):
        return False
    expected = []
    for raw, row in zip(family['directions'], cert['proofs']):
        q = list(map(F, raw))
        p, n = row['positive'], row['negative']
        if not exact_extrema.verify(p) or not exact_extrema.verify(n):
            return False
        if list(map(F, p['polynomial'])) != q or list(map(F, n['polynomial'])) != arithmetic.scale(q, -1):
            return False
        if p['interval'] != ['-1', '1'] or n['interval'] != ['-1', '1']:
            return False
        value = str(max(F(p['maximum_upper']), F(n['maximum_upper'])))
        if not same(row, {'positive': p, 'negative': n, 'linf_upper': value}):
            return False
        expected.append(value)
    return same(cert, {'format': PREFIX+'direction_norms_v126', 'family': family,
                       'proofs': cert['proofs'], 'norm_upper': expected})


def _anchor(family, point, cache=None, modes=6, bits=32):
    point = _point(family, point)
    key = tuple(_strings(point))
    q = potential_at(family, point)
    cached = None if cache is None else cache.get(key)
    if cached is not None:
        if not projected.verify_full(cached, expected_q=q, expected_mean_zero=family['mean_zero']):
            raise ValueError('Cached anchor has wrong potential, scope, or proof')
        return cached
    proof = projected.full_ground(q, mean_zero=family['mean_zero'], modes=modes,
                     max_modes=modes, bits=bits, tolerance=F(1, 10**7))
    if cache is not None:
        cache[key] = proof
    return proof


def _anchor_cell_record(family, cell, anchor, proof, norms):
    box, point = _cell(family, cell), _point(family, anchor, cell)
    if not _verify_norms(norms, family) or not projected.verify_full(proof,
               expected_q=potential_at(family, point), expected_mean_zero=family['mean_zero']):
        raise ValueError('Wrong anchor or direction proof')
    radius = sum((F(norm)*max(x-lo, hi-x) for norm, x, (lo, hi) in
                  zip(norms['norm_upper'], point, box)), F(0))
    return {'format': PREFIX+'anchor_cell_v127', 'family': family,
            'cell': [_strings(v) for v in box], 'anchor': _strings(point), 'anchor_proof': proof,
            'direction_norms': norms, 'lipschitz_radius': str(radius),
            'lower': str(F(proof['lower'])-radius), 'upper': str(F(proof['upper'])+radius),
            'scope': 'every_parameter_in_cell_on_stated_full_function_space'}


def anchor_cell(family, cell=None, anchor=None, cache=None, modes=6, bits=32, norms=None):
    """V127: a point spectrum plus certified multiplier perturbation radius."""
    box = _cell(family, cell)
    point = [(lo+hi)/2 for lo, hi in box] if anchor is None else _point(family, anchor, cell)
    return _anchor_cell_record(family, box, point, _anchor(family, point, cache, modes, bits),
                               direction_norms(family) if norms is None else norms)


def partition_certificate(family, cells):
    """V128: exact finite closed-box coverage, no missing/overlapping interior."""
    domain = _cell(family)
    if not isinstance(cells, list) or not 1 <= len(cells) <= 512:
        raise ValueError('One to 512 cells required')
    boxes = [_cell(family, cell) for cell in cells]
    for left, right in combinations(boxes, 2):
        if all(max(a, c) < min(b, d) for (a, b), (c, d) in zip(left, right)):
            raise ValueError('Cell interiors overlap')
    def volume(box):
        result = F(1)
        for lo, hi in box:
            result *= hi-lo
        return result
    if sum(map(volume, boxes), F(0)) != volume(domain):
        raise ValueError('Cell volumes do not cover the entire parent box')
    return {'format': PREFIX+'box_partition_v128', 'family': family,
            'cells': [[_strings(v) for v in box] for box in boxes],
            'cell_volumes': _strings(map(volume, boxes)), 'parent_volume': str(volume(domain)),
            'coverage': 'exact_closed_box_union_with_disjoint_interiors'}


def vertices(cell):
    return [list(v) for v in product(*_box(cell, len(cell)))]


def _vertex_record(family, cell, records):
    box = _cell(family, cell)
    points = vertices(box)
    if len(points) != len(records):
        raise ValueError('Every vertex required')
    for point, record in zip(points, records):
        if record['point'] != _strings(point) or not projected.verify_full(record['spectrum'],
                   expected_q=potential_at(family, point), expected_mean_zero=family['mean_zero']):
            raise ValueError('Missing, reordered, or wrongly bound vertex')
        if set(record) != {'point', 'spectrum'}:
            raise ValueError('Unknown vertex fields')
    lo = min(F(row['spectrum']['lower']) for row in records)
    hi = min(F(row['spectrum']['upper']) for row in records)
    return {'format': PREFIX+'vertex_concavity_v129', 'family': family,
            'cell': [_strings(v) for v in box], 'vertices': records,
            'lower': str(lo), 'minimum_upper': str(hi), 'minimum_gap': str(hi-lo),
            'maximum_lower': str(max(F(row['spectrum']['lower']) for row in records)),
            'claim': 'uniform_lower_for_all_parameters_and_enclosure_of_parameter_minimum',
            'rule': 'lowest_constrained_eigenvalue_is_infimum_of_affine_Rayleigh_forms_hence_concave',
            'warning': 'The maximum of vertex values is NOT a bound from above on the concave family maximum.'}


def cell_lower_bound(family, cell=None, cache=None, modes=6, bits=32):
    """V129: concavity certifies the entire affine box from its vertices."""
    box = _cell(family, cell)
    rows = [{'point': _strings(point), 'spectrum': _anchor(family, point, cache, modes, bits)}
            for point in vertices(box)]
    return _vertex_record(family, box, rows)


def affine_majorant(family, trial):
    """V130: one fixed admissible trial gives an exact affine spectral upper."""
    _family(family)
    if not witness.verify(trial, expected_mean_zero=family['mean_zero']) or trial.get('format') != witness.RAYLEIGH:
        raise ValueError('A valid explicit Rayleigh trial in the family function space is required')
    args = (trial['coefficients'], trial['azimuth_m'], trial['mean_zero'], trial['degrees'])
    baseline = witness.rayleigh_certificate(family['q0'], *args)
    rows = [witness.rayleigh_certificate(q, *args) for q in family['directions']]
    kinetic = F(baseline['kinetic_divided_by_azimuth_factor'])/F(baseline['norm_squared_divided_by_azimuth_factor'])
    return {'format': PREFIX+'affine_rayleigh_v130', 'family': family,
            'trial_at_q0': baseline, 'direction_trials': rows,
            'intercept': baseline['rayleigh_quotient'],
            'slopes': [str(F(row['rayleigh_quotient'])-kinetic) for row in rows],
            'claim': 'lambda(parameter) <= intercept + sum(slopes_i * parameter_i) for every parameter'}


def _majorants(family, majorants):
    _family(family)
    if not isinstance(majorants, list) or not 1 <= len(majorants) <= 128:
        raise ValueError('One to 128 affine majorants required')
    for row in majorants:
        if not same(row['family'], family) or not same(row, affine_majorant(family, row['trial_at_q0'])):
            raise ValueError('Invalid or wrongly bound affine majorant')
    return majorants


def _affine_at(row, point):
    return F(row['intercept'])+sum((F(a)*x for a, x in zip(row['slopes'], point)), F(0))


def envelope_1d(family, majorants, cell=None):
    """V131: exact supremum of a piecewise-affine 1D upper envelope."""
    box = _cell(family, cell)
    if len(box) != 1:
        raise ValueError('Line intersections require a one-dimensional family')
    _majorants(family, majorants)
    lo, hi = box[0]
    candidates = {lo, hi}
    for a, b in combinations(majorants, 2):
        slope = F(a['slopes'][0])-F(b['slopes'][0])
        if slope:
            x = (F(b['intercept'])-F(a['intercept']))/slope
            if lo <= x <= hi:
                candidates.add(x)
    rows = [{'parameter': str(x), 'envelope_value': str(min(_affine_at(row, [x]) for row in majorants))}
            for x in sorted(candidates)]
    best = max(rows, key=lambda row: F(row['envelope_value']))
    return {'format': PREFIX+'line_envelope_v131', 'family': family,
            'cell': [_strings(v) for v in box], 'majorants': majorants, 'candidates': rows,
            'maximum_upper': best['envelope_value'], 'envelope_maximizer': best['parameter'],
            'claim': 'supremum_of_lambda_over_cell_is_at_most_the_exact_envelope_maximum'}


def envelope_box(family, majorants, cell=None):
    """V132: safe box upper, min of affine maxima (not vertex max of min)."""
    box = _cell(family, cell)
    _majorants(family, majorants)
    maxima = []
    for row in majorants:
        point = [hi if F(a) >= 0 else lo for a, (lo, hi) in zip(row['slopes'], box)]
        maxima.append({'maximizing_corner': _strings(point), 'affine_maximum': str(_affine_at(row, point))})
    return {'format': PREFIX+'box_envelope_v132', 'family': family, 'cell': [_strings(v) for v in box],
            'majorants': majorants, 'individual_maxima': maxima,
            'maximum_upper': str(min(F(row['affine_maximum']) for row in maxima)),
            'claim': 'safe_supremum_upper_from_minimum_of_individual_affine_box_maxima',
            'limitation': 'Not generally the exact maximum of the minimum affine envelope in dimension above one.'}


def _adaptive_record(family, partition, anchors, tolerance, budget, majorants, sample_count):
    _family(family)
    tolerance = exact(tolerance)
    if tolerance <= 0 or type(budget) is not int or not 1 <= budget <= 128:
        raise ValueError('Positive exact tolerance and cell budget in 1..128 required')
    if family['parameter_dimension'] != 1:
        raise ValueError('Adaptive interval subdivision currently one-dimensional')
    if not same(partition, partition_certificate(family, partition['cells'])) or len(anchors) != len(partition['cells']):
        raise ValueError('Incomplete adaptive domain coverage')
    if len(anchors) > budget:
        raise ValueError('Cell budget exceeded')
    if majorants:
        _majorants(family, majorants)
    lowers, uppers, cells = [], [], []
    for cell, anchor in zip(partition['cells'], anchors):
        if not verify(anchor, expected_family=family) or not same(cell, anchor['cell']):
            raise ValueError('Cell anchor mismatch')
        lo, hi = F(anchor['anchor_proof']['lower']), F(anchor['upper'])
        envelope = None
        if majorants:
            envelope = envelope_1d(family, majorants, cell)
            hi = min(hi, F(envelope['maximum_upper']))
        lowers.append(lo)
        uppers.append(hi)
        cells.append({'anchor': anchor, 'envelope': envelope, 'maximum_upper': str(hi)})
    lo, hi = max(lowers), max(uppers)
    if hi < lo:
        raise ValueError('Inconsistent continuum maximum bounds')
    return {'format': PREFIX+'adaptive_supremum_v133', 'family': family, 'partition': partition,
            'cells': cells, 'majorants': majorants, 'tolerance': str(tolerance), 'cell_budget': budget,
            'unique_anchor_count': sample_count, 'maximum_lower': str(lo), 'maximum_upper': str(hi),
            'maximum_gap': str(hi-lo), 'status': 'target_met' if hi-lo <= tolerance else 'budget_open',
            'claim': 'certified_enclosure_of_parameter_supremum_over_the_entire_interval',
            'point_proof_reuse': 'anchor proofs retained and validated through the parameter-keyed cache'}


def adaptive_family(family, tolerance=F(1, 100), max_cells=8, majorants=None, modes=6, bits=32, cache=None):
    """V133: refine cells threatening the global supremum, with honest budget."""
    _family(family)
    if family['parameter_dimension'] != 1 or type(max_cells) is not int or not 1 <= max_cells <= 128:
        raise ValueError('One-dimensional family and cell budget in 1..128 required')
    tolerance = exact(tolerance)
    if tolerance <= 0:
        raise ValueError('Positive rational tolerance required')
    cache = {} if cache is None else cache
    norms = direction_norms(family)
    majorants = [] if majorants is None else majorants
    anchors = [anchor_cell(family, cache=cache, modes=modes, bits=bits, norms=norms)]
    while True:
        partition = partition_certificate(family, [a['cell'] for a in anchors])
        # Trace-free acceptance relies only on the final covering cells. The
        # exact number of unique final proof anchors is independently visible.
        result = _adaptive_record(family, partition, anchors, tolerance, max_cells, majorants,
                                  len({tuple(a['anchor']) for a in anchors}))
        if result['status'] == 'target_met' or len(anchors) >= max_cells:
            return result
        index = max(range(len(anchors)), key=lambda i: F(result['cells'][i]['maximum_upper']))
        old = anchors.pop(index)
        lo, hi = map(F, old['cell'][0])
        middle = (lo+hi)/2
        children = [anchor_cell(family, [[a, b]], cache=cache, modes=modes, bits=bits, norms=norms)
                    for a, b in ((lo, middle), (middle, hi))]
        anchors[index:index] = children


def _universal_record(family, threshold, lower_proof, parameter, trial):
    _family(family)
    if not family['mean_zero']:
        raise ValueError('Weighted Poincare family requires mean-zero constraint')
    threshold = exact(threshold)
    if not verify(lower_proof, expected_family=family) or lower_proof['format'] != PREFIX+'vertex_concavity_v129' or lower_proof['cell'] != family['box']:
        raise ValueError('A complete-box concavity lower proof is required')
    counterexample = None
    if parameter is not None or trial is not None:
        point = _point(family, parameter)
        if not witness.verify(trial, expected_q=potential_at(family, point), expected_mean_zero=True) or trial.get('format') != witness.RAYLEIGH:
            raise ValueError('Invalid explicit parameter and trial binding')
        counterexample = {'parameter': _strings(point), 'trial': trial}
    status = 'proved' if F(lower_proof['lower']) >= threshold else 'undetermined'
    if counterexample is not None and F(trial['rayleigh_quotient']) < threshold:
        if status == 'proved':
            raise ValueError('Contradictory certificates')
        status = 'refuted_with_explicit_parameter_and_function'
    return {'format': PREFIX+'universal_inequality_v134', 'family': family, 'threshold': str(threshold),
            'statement': 'For every parameter in the box and every mean-zero real H1(S2) function u, integral(|grad u|^2+q_parameter*u^2)>=threshold*integral(u^2).',
            'lower_proof': lower_proof, 'counterexample': counterexample, 'status': status,
            'limitation': 'This is a fixed affine potential family statement, not the universal GN sharp constant.',
            'formal_assistant_checked': False}


def universal_inequality(family, threshold, modes=6, bits=32, cache=None, lower_proof=None):
    """V134: prove universally, refute by explicit parameter+u, or unresolved."""
    lower_proof = cell_lower_bound(family, cache=cache, modes=modes, bits=bits) if lower_proof is None else lower_proof
    if (not verify(lower_proof, expected_family=family)
            or lower_proof.get('format') != PREFIX+'vertex_concavity_v129'
            or lower_proof.get('cell') != family['box']):
        raise ValueError('A whole-box vertex-concavity lower certificate is required')
    parameter = trial = None
    if F(lower_proof['lower']) < exact(threshold):
        row = min(lower_proof['vertices'], key=lambda row: F(row['spectrum']['upper']))
        parameter = row['point']
        m = min(row['spectrum']['sectors'], key=lambda c: F(c['upper']))['azimuth_m']
        proposal = witness.mass_normalized_proposal(potential_at(family, parameter), m=m, mean_zero=True, modes=modes)
        trial = witness.quantize_proposal(proposal)['certificate']
    return _universal_record(family, threshold, lower_proof, parameter, trial)


def verify(cert, expected_family=None, expected_threshold=None):
    """Exact replay: primitive searches are never rerun by this verifier."""
    try:
        kind = cert['format']
        family = cert if kind == PREFIX+'representation_v125' else cert['family']
        _family(family)
        if expected_family is not None and not same(family, expected_family):
            return False
        if kind == PREFIX+'representation_v125':
            valid = True
        elif kind == PREFIX+'cheap_enclosure_v125':
            valid = same(cert, cheap_family_enclosure(family))
        elif kind == PREFIX+'direction_norms_v126':
            valid = _verify_norms(cert, family)
        elif kind == PREFIX+'anchor_cell_v127':
            valid = same(cert, _anchor_cell_record(family, cert['cell'], cert['anchor'], cert['anchor_proof'], cert['direction_norms']))
        elif kind == PREFIX+'box_partition_v128':
            valid = same(cert, partition_certificate(family, cert['cells']))
        elif kind == PREFIX+'vertex_concavity_v129':
            valid = same(cert, _vertex_record(family, cert['cell'], cert['vertices']))
        elif kind == PREFIX+'affine_rayleigh_v130':
            valid = same(cert, affine_majorant(family, cert['trial_at_q0']))
        elif kind == PREFIX+'line_envelope_v131':
            valid = same(cert, envelope_1d(family, cert['majorants'], cert['cell']))
        elif kind == PREFIX+'box_envelope_v132':
            valid = same(cert, envelope_box(family, cert['majorants'], cert['cell']))
        elif kind == PREFIX+'adaptive_supremum_v133':
            anchors = [row['anchor'] for row in cert['cells']]
            valid = same(cert, _adaptive_record(family, cert['partition'], anchors, cert['tolerance'], cert['cell_budget'],
                    cert['majorants'], len({tuple(a['anchor']) for a in anchors})))
        elif kind == PREFIX+'universal_inequality_v134':
            counter = cert['counterexample']
            valid = same(cert, _universal_record(family, cert['threshold'], cert['lower_proof'],
                    None if counter is None else counter['parameter'], None if counter is None else counter['trial']))
        else:
            return False
        if expected_threshold is not None:
            return valid and kind == PREFIX+'universal_inequality_v134' and F(cert['threshold']) == exact(expected_threshold)
        return valid
    except (ValueError, TypeError, KeyError, IndexError, ArithmeticError, OverflowError):
        return False
