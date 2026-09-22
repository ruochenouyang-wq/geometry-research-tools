"""Whole-amplitude-interval certificates from concavity and fixed Rayleigh trials.

Ordinary solve is explicitly the unchanged baseline. The extension prepares a
self-contained family proof; its fast context is created only after replay.
"""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from math import isfinite
from pathlib import Path
import sys
from time import perf_counter, process_time

sys.dont_write_bytecode = True
PARENT = Path(__file__).resolve().parents[1]
if str(PARENT) not in sys.path:
    sys.path.insert(0, str(PARENT))
import shared_baseline as shared

VERSION = 'b_parameter_concavity_v1'
FAMILY = 'complete_amplitude_family_concavity_v1'
QUERY = 'amplitude_family_point_query_v1'
MAX_SAMPLES = 65


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def native(value, depth=0):
    if depth > 64:
        return False
    if type(value) in (str, int, bool, type(None)):
        return True
    if type(value) is list:
        return all(native(v, depth+1) for v in value)
    if type(value) is dict:
        return all(type(k) is str and native(v, depth+1) for k, v in value.items())
    return False


def normalize_family(spec):
    allowed = {'profile', 'exponent', 'axis', 'amplitude_interval', 'offset', 'mean_zero', 'tolerance'}
    if type(spec) is not dict or set(spec)-allowed:
        raise ValueError('Malformed family fields')
    if spec.get('profile', 'abs_power') != 'abs_power' or spec.get('mean_zero') is not False:
        raise ValueError('This extension requires full-space abs_power families')
    raw = spec['amplitude_interval']
    if type(raw) is not list or len(raw) != 2:
        raise ValueError('Two exact amplitude endpoints are required')
    lo, hi = map(shared.baseline.rational, raw)
    if not -1000 <= lo < hi < 0 or max(lo.denominator, hi.denominator) > 10**6:
        raise ValueError('Require two distinct strictly negative rational amplitudes')
    task = shared.baseline.normalize_task({
        'kind': 'spectrum', 'mean_zero': False, 'tolerance': spec.get('tolerance', '1/1000000'),
        'function': {'kind': 'axis_profile', 'profile': 'abs_power', 'axis': spec.get('axis', 2),
                     'exponent': spec['exponent'], 'amplitude': str(lo), 'offset': spec.get('offset', '0')}})
    q = task['function']
    if not -F(1, 2) < F(q['exponent']) < 0:
        raise ValueError('Require -1/2 < exponent < 0')
    return {'profile': 'abs_power', 'exponent': q['exponent'], 'axis': q['axis'],
            'amplitude_interval': [str(lo), str(hi)], 'offset': q['offset'],
            'mean_zero': False, 'tolerance': task['tolerance']}


def endpoint_task(spec, amplitude, tolerance=None, wall_seconds=10):
    return {'kind': 'spectrum', 'mean_zero': False,
            'function': {'kind': 'axis_profile', 'profile': spec['profile'], 'axis': spec['axis'],
                         'exponent': spec['exponent'], 'amplitude': str(amplitude), 'offset': spec['offset']},
            'tolerance': tolerance or spec['tolerance'], 'budget': {'wall_seconds': wall_seconds}}


def _rayleigh_affine(cert, spec):
    """Independently integrate M, kinetic energy and W moment (no Hu needed)."""
    import fullspace_singular as fs
    if cert.get('format') != fs.FORMAT:
        raise ValueError('An endpoint must carry the supported full-space trial certificate')
    trial = cert['trial']
    powers = fs.powers(trial['powers'])
    values = tuple(shared.baseline.rational(x) for x in trial['coefficients'])
    if len(values) != len(powers):
        raise ValueError('Trial coefficient dimension mismatch')
    mass, kinetic, potential = F(0), F(0), F(0)
    alpha = F(spec['exponent'])
    for s, v in zip(powers, values):
        for t, w in zip(powers, values):
            m = F(2)/(s+t+1)
            mass += v*w*m
            if s*t:
                kinetic += v*w*s*t*(F(2)/(s+t-1)-m)
            potential += v*w*F(2)/(s+t+alpha+1)
    if mass <= 0:
        raise ValueError('Nonpositive trial mass')
    result = {'slope': str(potential/mass),
              'intercept': str(kinetic/mass+F(spec['offset'])),
              'mass': str(mass), 'kinetic_form': str(kinetic),
              'potential_form': str(potential),
              'scope': 'fixed_admissible_m0_trial_in_the_full_space'}
    a = F(cert['function']['amplitude'])
    if mass != F(cert['statistics']['mass']) or _line(result, a) != F(cert['statistics']['rayleigh']):
        raise ArithmeticError('Independent Rayleigh affine disagrees with the endpoint proof')
    return result


def _line(line, a):
    return F(line['slope'])*a+F(line['intercept'])


def _sample(spec, amplitude, cert, replay=True):
    a = F(amplitude)
    if replay:
        task = endpoint_task(spec, a, cert['tolerance'])
        judgment = shared.verify(cert, task)
        if judgment.get('certificate_valid') is not True:
            raise ValueError('Endpoint complete certificate did not replay')
    return {'amplitude': str(a), 'certificate': deepcopy(cert),
            'rayleigh_affine': _rayleigh_affine(cert, spec)}


def _cell(left, right):
    a, b = F(left['amplitude']), F(right['amplitude'])
    if a >= b:
        raise ValueError('Empty or reversed amplitude cell')
    low_a, low_b = F(left['certificate']['lower']), F(right['certificate']['lower'])
    slope = (low_b-low_a)/(b-a)
    lower = {'slope': str(slope), 'intercept': str(low_a-slope*a)}
    lines = [deepcopy(left['rayleigh_affine']), deepcopy(right['rayleigh_affine'])]
    points = {a, b}
    slope_difference = F(lines[0]['slope'])-F(lines[1]['slope'])
    if slope_difference:
        crossing = (F(lines[1]['intercept'])-F(lines[0]['intercept']))/slope_difference
        if a < crossing < b:
            points.add(crossing)
    checks = []
    for x in sorted(points):
        low, high = _line(lower, x), min(_line(line, x) for line in lines)
        if low > high:
            raise ArithmeticError('Inconsistent family lower and upper bounds')
        checks.append({'amplitude': str(x), 'lower': str(low), 'upper': str(high),
                       'exact_width': str(high-low)})
    return {'amplitude_interval': [str(a), str(b)], 'lower_affine': lower,
            'upper_affines': lines, 'critical_points': checks,
            'uniform_width_upper': str(max(F(c['exact_width']) for c in checks)),
            'maximum_rule': 'endpoints_and_all_intersections_of_two_affine_upper_bounds'}


def _bank(spec, samples):
    samples = sorted(samples, key=lambda x: F(x['amplitude']))
    if not 2 <= len(samples) <= MAX_SAMPLES:
        raise ValueError('Family sample budget is 2..65')
    points = [F(s['amplitude']) for s in samples]
    if points != sorted(set(points)) or [str(points[0]), str(points[-1])] != spec['amplitude_interval']:
        raise ValueError('Samples do not cover exactly the declared interval')
    cells = [_cell(a, b) for a, b in zip(samples, samples[1:])]
    width = max(F(c['uniform_width_upper']) for c in cells)
    if width > F(spec['tolerance']):
        raise ValueError('The whole interval has not reached its declared uniform tolerance')
    return {'format': FAMILY, 'family': deepcopy(spec), 'samples': deepcopy(samples), 'cells': cells,
            'uniform_width_upper': str(width), 'tolerance': spec['tolerance'], 'status': 'target_met',
            'proof': {'operator': 'minus_Laplacian_plus_a_times_abs_axis_to_alpha_plus_offset',
                      'common_form_domain': 'H1_unit_S2_for_fixed_alpha_greater_than_minus_half',
                      'lower_rule': 'concavity_of_infimum_of_affine_Rayleigh_forms',
                      'upper_rule': 'Rayleigh_principle_for_each_fixed_endpoint_trial',
                      'covers_every_real_amplitude_in_interval': True,
                      'parameter_grid_is_not_the_proof': True,
                      'full_infinite_space_covered': True},
            'formal_proof_assistant_checked': False}


def verify_family(bank, expected_spec=None):
    started, cpu = perf_counter(), process_time()
    result = {'certificate_valid': False, 'target_met': False}
    try:
        if not native(bank) or type(bank) is not dict or bank.get('format') != FAMILY:
            raise ValueError('Malformed family certificate')
        spec = normalize_family(bank['family'])
        if expected_spec is not None and normalize_family(expected_spec) != spec:
            raise ValueError('Family does not match the requested parameter box')
        if type(bank['samples']) is not list or not 2 <= len(bank['samples']) <= MAX_SAMPLES:
            raise ValueError('Invalid family sample count')
        samples = [_sample(spec, s['amplitude'], s['certificate']) for s in bank['samples']]
        rebuilt = _bank(spec, samples)
        valid = canonical(rebuilt) == canonical(bank)
        result.update(certificate_valid=valid, target_met=valid,
                      uniform_width_upper=rebuilt['uniform_width_upper'],
                      status='target_met' if valid else 'verification_failed')
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, RecursionError):
        result['status'] = 'verification_failed'
    result.update(verification_wall_seconds=perf_counter()-started,
                  verification_cpu_seconds=process_time()-cpu)
    return result


def prepare_family(spec, preparation_seconds=30):
    started, cpu = perf_counter(), process_time()
    spec = normalize_family(spec)
    if isinstance(preparation_seconds, bool) or not isfinite(float(preparation_seconds)) or not 0 <= float(preparation_seconds) <= 120:
        raise ValueError('Preparation budget must be finite and between 0 and 120 seconds')
    deadline = started+float(preparation_seconds)
    endpoint_tol = str(max(F(1, 10**30), F(spec['tolerance'])/16))
    samples, attempts, bank, status = {}, [], None, 'preparation_failed'

    def add_sample(amplitude):
        nonlocal status
        if perf_counter() >= deadline:
            status = 'preparation_budget_exceeded'
            return False
        task = endpoint_task(spec, amplitude, endpoint_tol, max(0.0, deadline-perf_counter()))
        begin, begin_cpu = perf_counter(), process_time()
        response = shared.solve(task)
        cert = response.get('certificate')
        entry = {'phase': 'endpoint', 'amplitude': str(amplitude),
                 'solver_status': response.get('status'), 'solver_attempts': response.get('attempts', [])}
        try:
            if cert is None:
                raise ValueError('No complete endpoint certificate')
            sample = _sample(spec, amplitude, cert, replay=True)
            entry.update(certificate_valid=True, exact_width=cert['exact_width'])
            samples[str(amplitude)] = sample
            success = True
        except (ValueError, TypeError, KeyError, ArithmeticError, IndexError) as error:
            entry.update(certificate_valid=False, error=str(error))
            success = False
        entry.update(wall_seconds=perf_counter()-begin, cpu_seconds=process_time()-begin_cpu)
        attempts.append(entry)
        if perf_counter() > deadline:
            status = 'preparation_budget_exceeded'
            return False
        return success

    endpoints = list(map(F, spec['amplitude_interval']))
    if add_sample(endpoints[0]) and add_sample(endpoints[1]):
        while True:
            ordered = sorted(samples.values(), key=lambda s: F(s['amplitude']))
            cells = [_cell(a, b) for a, b in zip(ordered, ordered[1:])]
            worst = max(cells, key=lambda c: F(c['uniform_width_upper']))
            attempts.append({'phase': 'uniform_interval_check', 'sample_count': len(ordered),
                             'uniform_width_upper': worst['uniform_width_upper'],
                             'target_met': F(worst['uniform_width_upper']) <= F(spec['tolerance'])})
            if F(worst['uniform_width_upper']) <= F(spec['tolerance']):
                candidate = _bank(spec, ordered)
                judgment = verify_family(candidate, spec)
                attempts.append({'phase': 'final_complete_family_replay', **judgment})
                if judgment['certificate_valid'] and perf_counter() <= deadline:
                    bank, status = candidate, 'target_met'
                else:
                    status = 'preparation_budget_exceeded' if perf_counter() > deadline else 'verification_failed'
                break
            if len(samples) >= MAX_SAMPLES:
                status = 'preparation_sample_budget_exhausted'
                break
            a, b = map(F, worst['amplitude_interval'])
            if not add_sample((a+b)/2):
                break
    return {'version': VERSION, 'status': status, 'bank': bank, 'attempts': attempts,
            'preparation_wall_seconds': perf_counter()-started,
            'preparation_cpu_seconds': process_time()-cpu,
            'preparation_budget_seconds': float(preparation_seconds),
            'budget_scope': 'whole_preparation_including_endpoint_and_final_family_replays; cooperative',
            'query_results_precomputed': False}


def _bound_task(spec, task):
    task = shared.baseline.normalize_task(task)
    if task['kind'] != 'spectrum' or task['mean_zero'] is not False:
        raise ValueError('A family query must keep the full spectral space')
    q = task['function']
    a = F(q['amplitude'])
    expected = endpoint_task(spec, a)['function']
    if q != expected or not F(spec['amplitude_interval'][0]) <= a <= F(spec['amplitude_interval'][1]):
        raise ValueError('Query is outside the certified family')
    return task


def _query_certificate(bank, task):
    task = _bound_task(bank['family'], task)
    a = F(task['function']['amplitude'])
    index = next(i for i, c in enumerate(bank['cells'])
                 if F(c['amplitude_interval'][0]) <= a <= F(c['amplitude_interval'][1]))
    cell = bank['cells'][index]
    low = _line(cell['lower_affine'], a)
    high = min(_line(line, a) for line in cell['upper_affines'])
    return {'format': QUERY, 'function': task['function'], 'mean_zero': False,
            'geometry': 'unit_S2', 'scope': 'all_real_H1_on_unit_S2', 'eigenvalue_index': 1,
            'family_digest': digest(bank), 'family_certificate': deepcopy(bank), 'cell_index': index,
            'lower': str(low), 'upper': str(high), 'exact_width': str(high-low),
            'tolerance': task['tolerance'], 'status': 'target_met' if high-low <= F(task['tolerance']) else 'certified_open',
            'full_infinite_space_covered': True, 'query_is_proof_evaluation_not_endpoint_cache_hit': True}


class PreparedContext:
    """Owns a serialized snapshot verified before publication to any query."""
    def __init__(self, bank):
        begin, cpu = perf_counter(), process_time()
        if not native(bank):
            raise ValueError('A family must use exact JSON-native data')
        # Snapshot BEFORE replay, so concurrent caller mutation cannot replace
        # evidence between its verification and our retained owned copy.
        serialized = canonical(bank)
        owned = json.loads(serialized)
        judgment = verify_family(owned)
        if not judgment['certificate_valid']:
            raise ValueError('Cannot create a context from an unverified family')
        self._serialized = serialized
        self._digest = hashlib.sha256(serialized.encode()).hexdigest()
        self.verification = judgment
        self.setup_wall_seconds = perf_counter()-begin
        self.setup_cpu_seconds = process_time()-cpu

    def query(self, task):
        begin, cpu = perf_counter(), process_time()
        bank = json.loads(self._serialized)
        normalized = _bound_task(bank['family'], task)
        if normalized['budget']['wall_seconds'] == 0:
            return {'version': VERSION, 'certificate': None, 'certificate_valid': False,
                    'target_met': False, 'status': 'budget_exceeded', 'fallback': False}
        cert = _query_certificate(bank, normalized)
        judgment = self.verify(cert, normalized)
        elapsed = perf_counter()-begin
        within = elapsed <= normalized['budget']['wall_seconds']
        return {'version': VERSION, 'certificate': cert, **judgment, 'fallback': False,
                'status': cert['status'] if within else 'budget_exceeded',
                'target_met': judgment['target_met'] and within, 'within_budget': within,
                'query_wall_seconds': elapsed, 'query_cpu_seconds': process_time()-cpu,
                'route': 'uniform_parameter_certificate_evaluation',
                'verification_scope': 'previously_verified_owned_family_plus_exact_query_binding'}

    def verify(self, cert, task):
        try:
            if not native(cert) or type(cert) is not dict or cert.get('format') != QUERY:
                raise ValueError('Malformed query certificate')
            bank = cert['family_certificate']
            if cert['family_digest'] != self._digest or digest(bank) != self._digest:
                raise ValueError('Query family differs from the verified owned family')
            rebuilt = _query_certificate(json.loads(self._serialized), task)
            valid = canonical(rebuilt) == canonical(cert)
            return {'certificate_valid': valid,
                    'target_met': valid and F(cert['exact_width']) <= F(rebuilt['tolerance']),
                    'metric': rebuilt['exact_width'], 'tolerance': rebuilt['tolerance']}
        except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, RecursionError):
            return {'certificate_valid': False, 'target_met': False}


def prepare_context(bank):
    return PreparedContext(bank)


def query(bank, task):
    """Cold query includes a complete family replay; use context for a stream."""
    begin, cpu = perf_counter(), process_time()
    normalized = shared.baseline.normalize_task(task)
    context = prepare_context(bank)
    result = context.query(normalized)
    total = perf_counter()-begin
    within = total <= normalized['budget']['wall_seconds']
    result.update(total_elapsed_seconds=total, total_cpu_seconds=process_time()-cpu,
                  family_verification=context.verification, within_budget=within,
                  target_met=result['target_met'] and within)
    if not within:
        result['status'] = 'budget_exceeded'
    return result


def verify(cert, task):
    """Self-contained cold replay; no trust in a previous process's context."""
    begin, cpu = perf_counter(), process_time()
    try:
        if type(cert) is dict and cert.get('format') == QUERY:
            context = prepare_context(cert['family_certificate'])
            result = context.verify(cert, task)
        else:
            result = shared.verify(cert, task)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, RecursionError):
        result = {'certificate_valid': False, 'target_met': False}
    result.update(verification_wall_seconds=perf_counter()-begin,
                  verification_cpu_seconds=process_time()-cpu)
    return result


def solve(task):
    result = shared.solve(task)
    return {**result, 'variant': VERSION, 'fallback': True, 'fallback_reason': 'ordinary_single_point_uses_unchanged_baseline'}
