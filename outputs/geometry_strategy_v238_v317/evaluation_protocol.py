"""Independent, predeclared local evaluation; no model/network calls.

Only complete, replayed and task-bound certificates count. This evaluator is
outside v238-v317: it does not represent an algorithm iteration.
"""
from copy import deepcopy
from datetime import datetime, timezone
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import platform
import secrets
import signal
import sys
import threading
import time

import backend

ROOT = Path(__file__).resolve().parent
DATA = ROOT / 'evaluation'
VERSION = 'independent_geometry_evaluation_v1'
TASK_SECONDS = 20.0
REPLAY_SECONDS = 20.0
CONFIGS = [
    {'name': 'small', 'approximation': {'levels': 1, 'degree': 1, 'sqrt_bits': 24},
     'spectrum': {'modes': 2, 'bits': 12, 'max_m': 1, 'sqrt_bits': 24, 'near_tail': 0}},
    {'name': 'medium', 'approximation': {'levels': 3, 'degree': 3, 'sqrt_bits': 28},
     'spectrum': {'modes': 3, 'bits': 18, 'max_m': 2, 'sqrt_bits': 28, 'near_tail': 0}},
    {'name': 'large', 'approximation': {'levels': 6, 'degree': 4, 'sqrt_bits': 32},
     'spectrum': {'modes': 5, 'bits': 26, 'max_m': 3, 'sqrt_bits': 32, 'near_tail': 0}},
]
DIMENSIONS = ['certified_accuracy', 'reliability', 'automatic_selection',
              'new_parameter_transfer', 'runtime_resources', 'visible_tokens', 'usability']


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode('utf-8')).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def q(profile='step', amplitude='1', offset='0', exponent=None):
    function = {'kind': 'axis_profile', 'profile': profile, 'amplitude': amplitude, 'offset': offset, 'axis': 2}
    if exponent is not None:
        function['exponent'] = exponent
    return backend.direct.normalize(function)


def case(identifier, kind, function, tolerance, mean_zero=None, label='development'):
    task = {'kind': kind, 'function': function, 'tolerance': tolerance}
    if kind == 'spectrum':
        task['mean_zero'] = True if mean_zero is None else mean_zero
    item = {'id': identifier, 'task': task, 'classification': label}
    if kind == 'spectrum' and F(function['amplitude']) == 0:
        item['independent_exact_value'] = str(F(function['offset']) + (2 if task['mean_zero'] else 0))
    return item


def development_cases():
    return [
        case('D01', 'spectrum', q(amplitude='0', offset='1/3'), '1/1000', False, 'analytic_constant_regression'),
        case('D02', 'spectrum', q(amplitude='0', offset='1/3'), '1/1000', True, 'analytic_constant_regression'),
        case('D03', 'approximation', q(), '1/1000000', label='step_regression'),
        case('D04', 'approximation', q('abs_power', '-1', exponent='-1/4'), '1/20', label='fixed_singular_regression_not_holdout'),
        case('D05', 'spectrum', q(amplitude='1/4'), '1/1000', label='step_development'),
        case('D06', 'spectrum', q('abs_power', '-1/4', exponent='-1/5'), '1/1000', label='abs_power_development'),
    ]


def validation_cases():
    return [
        case('V01', 'spectrum', q(amplitude='0', offset='-2/5'), '1/1000', False, 'validation_not_holdout'),
        case('V02', 'approximation', q(amplitude='3/4', offset='1/5'), '1/1000000', label='validation_not_holdout'),
        case('V03', 'approximation', q('abs_power', '-1/2', exponent='-1/4'), '1/50', label='fixed_singular_validation_not_holdout'),
        case('V04', 'spectrum', q('abs_power', '-1/3', exponent='-1/6'), '1/1000', label='validation_not_holdout'),
        case('V05', 'spectrum', q(amplitude='1/3'), '1/1000', label='validation_not_holdout'),
    ]


def public_rules():
    return {
        'version': VERSION, 'dimensions': DIMENSIONS,
        'development': development_cases(), 'validation': validation_cases(), 'configurations': CONFIGS,
        'budget': {'online_solver_wall_seconds_per_task': TASK_SECONDS, 'replay_wall_seconds_per_certificate': REPLAY_SECONDS},
        'holdout': {
            'count': 6, 'fixed_axis': 2, 'scope': 'within_family_new_parameters_only',
            'layout': ['constant_full', 'constant_mean_zero', 'step_approximation', 'step_spectrum', 'abs_power_spectrum', 'abs_power_spectrum'],
            'constant_offsets': ['-5/7', '-3/7', '2/7', '4/7'],
            'step_amplitudes': ['1/7', '2/7', '3/7', '4/7'],
            'abs_power_exponents': ['-1/7', '-1/8', '-2/7', '1/5', '2/5'],
            'abs_power_amplitudes': ['-1/10', '-1/8', '-1/6', '1/10', '1/8', '1/6'],
            'spectral_tolerance': '1/1000', 'approximation_tolerance': '1/1000000',
            'sampling': 'SHA256(seed || counter_as_8_byte_big_endian), modulo list length; two distinct abs_power exponents; stable form_bound required',
            'no_fixed_minus_quarter_singular_holdout': True,
        },
        'acceptance': {
            'accuracy': 'Exact replay, original function/kind/space binding, nonnegative bound, requested tolerance, and known constant eigenvalue inclusion where applicable.',
            'reliability': 'All mutation probes must be rejected. Any falsely accepted mutation blocks a reliability improvement claim.',
            'automatic_selection': 'Compare all three baselines. No manual per-task settings for service; validation selection preparation costs are included.',
            'transfer': 'At least one extra certified holdout success over the strongest baseline without a reliability regression is an accuracy benefit. Equal success with >=20% lower total CPU including charged validation setup is only a preliminary resource benefit.',
            'new_mathematical_family': 'Not tested; constant offsets and rotations never count as new mathematical families.',
            'claims': 'No aggregate score; no model-vs-model or 30-minute/8-hour claim. Small sample and concurrent-run timing are explicit limitations.',
        },
    }


def initialize():
    """Create public protocol and a private seed without exposing seed contents."""
    DATA.mkdir(parents=True, exist_ok=True)
    rules = public_rules()
    public_path = DATA / 'protocol.json'
    if public_path.exists() and json.loads(public_path.read_text()) != rules:
        raise RuntimeError('Protocol changed: use a new evaluation version and new untouched holdout')
    write_json(public_path, rules)
    private = DATA / 'private_seed.json'
    commitment = DATA / 'seed_commitment.json'
    if not private.exists():
        if commitment.exists() or (DATA/'holdout_instances.json').exists():
            raise RuntimeError('Refusing to replace a missing committed private seed')
        seed = secrets.token_bytes(32).hex()
        descriptor = os.open(private, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        with os.fdopen(descriptor, 'w') as stream:
            json.dump({'seed_hex': seed}, stream)
            stream.write('\n')
    seed = json.loads(private.read_text())['seed_hex']
    value = {'sha256': hashlib.sha256(bytes.fromhex(seed)).hexdigest()}
    if commitment.exists() and json.loads(commitment.read_text()) != value:
        raise RuntimeError('Seed commitment mismatch')
    write_json(commitment, value)
    return {'protocol_digest': digest(rules), 'commitment_file': str(commitment), 'seed_disclosed': False}


def write_freeze_manifest(note='Code and policy frozen before holdout reveal'):
    """Call only after development completes; refuses to re-freeze used holdout."""
    initialize()
    if (DATA/'holdout_instances.json').exists():
        raise RuntimeError('Holdout already revealed; a new protocol and fresh holdout are required after changes')
    if not (DATA/'validation_selection.json').exists():
        raise RuntimeError('Run a development evaluation to freeze baseline selection first')
    files = list(sorted(ROOT.glob('*.py'))) + [DATA/'protocol.json', DATA/'validation_selection.json']
    manifest = {'status': 'frozen', 'note': note, 'created_utc': datetime.now(timezone.utc).isoformat(),
                'protocol_digest': digest(public_rules()),
                'files': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}}
    write_json(DATA/'freeze_manifest.json', manifest)
    return {'manifest': str(DATA/'freeze_manifest.json'), 'manifest_digest': digest(manifest), 'files': len(files)}


def checked_manifest():
    path = DATA/'freeze_manifest.json'
    if not path.exists():
        raise RuntimeError('Holdout requires evaluation/freeze_manifest.json before any instances are generated')
    manifest = json.loads(path.read_text())
    if manifest.get('status') != 'frozen' or manifest.get('protocol_digest') != digest(public_rules()):
        raise RuntimeError('Invalid freeze manifest or protocol mismatch')
    files = manifest.get('files', {})
    required = {p.name for p in ROOT.glob('*.py')} | {'evaluation/protocol.json', 'evaluation/validation_selection.json'}
    if not required <= set(files):
        raise RuntimeError('Freeze manifest must cover all current root Python modules and baseline selection')
    for relative, expected in files.items():
        source = (ROOT/relative).resolve()
        if not source.is_relative_to(ROOT.resolve()) or not source.is_file():
            raise RuntimeError('Invalid freeze file path')
        if hashlib.sha256(source.read_bytes()).hexdigest() != expected:
            raise RuntimeError('Source changed after freeze: ' + relative)
    return manifest


def holdout_cases():
    """No hidden instances are constructed until a complete freeze is verified."""
    manifest = checked_manifest()
    initialize()
    path = DATA/'holdout_instances.json'
    release_path = DATA/'holdout_release.json'
    if path.exists():
        release = json.loads(release_path.read_text())
        cases = json.loads(path.read_text())
        if release['freeze_digest'] != digest(manifest) or release['suite_digest'] != digest(cases):
            raise RuntimeError('Holdout release does not match frozen sources')
        return cases
    seed = bytes.fromhex(json.loads((DATA/'private_seed.json').read_text())['seed_hex'])
    rules = public_rules()['holdout']
    counter = 0
    def choose(values):
        nonlocal counter
        value = int.from_bytes(hashlib.sha256(seed + counter.to_bytes(8, 'big')).digest(), 'big')
        counter += 1
        return values[value % len(values)]
    c = choose(rules['constant_offsets'])
    a = choose(rules['step_amplitudes'])
    b = choose(rules['step_amplitudes'])
    exponents = list(rules['abs_power_exponents'])
    instances = [
        case('H01', 'spectrum', q(amplitude='0', offset=c), '1/1000', False, 'within_family_parameter_transfer_constant'),
        case('H02', 'spectrum', q(amplitude='0', offset=c), '1/1000', True, 'within_family_parameter_transfer_constant'),
        case('H03', 'approximation', q(amplitude=a), '1/1000000', label='within_family_parameter_transfer_step'),
        case('H04', 'spectrum', q(amplitude=b), '1/1000', label='within_family_parameter_transfer_step'),
    ]
    for identifier in ('H05', 'H06'):
        exponent = choose(exponents)
        exponents.remove(exponent)
        function = q('abs_power', choose(rules['abs_power_amplitudes']), exponent=exponent)
        # Domain/stability check is prescribed, not a success-dependent filter.
        backend.direct.form_bound(function)
        instances.append(case(identifier, 'spectrum', function, '1/1000', label='within_family_parameter_transfer_abs_power'))
    known = {canonical(c['task']['function']) for c in development_cases()+validation_cases()}
    if any(canonical(c['task']['function']) in known for c in instances):
        raise RuntimeError('Holdout parameter overlaps public suite')
    write_json(path, instances)
    write_json(release_path, {'freeze_digest': digest(manifest), 'suite_digest': digest(instances),
                              'released_utc': datetime.now(timezone.utc).isoformat(), 'one_shot_interpretation': True})
    return instances


class EvaluationTimeout(BaseException):
    pass


def timed_call(function, seconds):
    """A hard wall limit on the main Unix thread; report fallback explicitly."""
    start, cpu = time.perf_counter(), time.process_time()
    supported = hasattr(signal, 'SIGALRM') and threading.current_thread() is threading.main_thread()
    previous_handler = None
    previous_timer = None
    try:
        if supported:
            previous_handler = signal.getsignal(signal.SIGALRM)
            previous_timer = signal.setitimer(signal.ITIMER_REAL, 0)
            def alarm_handler(*_):
                raise EvaluationTimeout('Evaluator wall budget exhausted')
            signal.signal(signal.SIGALRM, alarm_handler)
            signal.setitimer(signal.ITIMER_REAL, max(0.000001, seconds))
        result = function()
        # Confirm the recorded envelope is fully serializable, not a summary object.
        canonical(result)
        return {'execution_ok': True, 'response': result,
                'wall_seconds': time.perf_counter()-start, 'cpu_seconds': time.process_time()-cpu,
                'hard_timeout_enforced': supported}
    except (Exception, EvaluationTimeout) as error:
        return {'execution_ok': False, 'error_type': type(error).__name__, 'reason': str(error),
                'wall_seconds': time.perf_counter()-start, 'cpu_seconds': time.process_time()-cpu,
                'hard_timeout_enforced': supported, 'timed_out': isinstance(error, EvaluationTimeout)}
    finally:
        if supported:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, previous_handler)
            old_delay, old_interval = previous_timer
            if old_delay:
                signal.setitimer(signal.ITIMER_REAL, max(0.000001, old_delay-(time.perf_counter()-start)), old_interval)


def default_verifier(certificate, expected_function=None, expected_mean_zero=None, expected_tolerance=None):
    """Replay frozen legacy formats; never accept a solver's own status flags."""
    try:
        fmt = certificate.get('format')
        if expected_function is not None and certificate.get('function') != backend.direct.normalize(expected_function):
            return False
        if expected_mean_zero is not None and certificate.get('mean_zero') is not expected_mean_zero:
            return False
        if fmt == backend.pieces.FORMAT:
            return backend.pieces.verify_piecewise(certificate, expected_function) is True
        if fmt == backend.direct.FULL:
            return backend.direct.verify_full(certificate, expected_function, expected_mean_zero, expected_tolerance) is True
        if fmt == backend.enriched.FORMAT:
            return backend.enriched.verify(certificate, expected_tolerance=expected_tolerance) is True
        return False
    except (TypeError, ValueError, KeyError, ArithmeticError):
        return False


def cert_kind(certificate):
    fmt = certificate.get('format')
    if fmt == backend.pieces.FORMAT:
        return 'approximation'
    if fmt in (backend.direct.FULL, backend.enriched.FORMAT, 'adaptive_singular_temple_v1'):
        return 'spectrum'
    if fmt == 'geometry_proof_transport_v1':
        return certificate.get('kind')
    return None


def assess(certificate, item, verifier):
    task = item['task']
    initial = {'certificate_valid': False, 'target_met': False, 'quantity':
               'function_L2_error_upper' if task['kind']=='approximation' else 'spectral_interval_width'}
    if type(certificate) is not dict or cert_kind(certificate) != task['kind']:
        return {**initial, 'reason': 'missing_or_wrong_conclusion_kind', 'replay_wall_seconds': 0.0, 'replay_cpu_seconds': 0.0}
    scope = task.get('mean_zero', True) if task['kind']=='spectrum' else None
    if certificate.get('function') != backend.direct.normalize(task['function']):
        return {**initial, 'reason': 'function_binding_mismatch', 'replay_wall_seconds': 0.0, 'replay_cpu_seconds': 0.0}
    if scope is not None and certificate.get('mean_zero') is not scope:
        return {**initial, 'reason': 'space_binding_mismatch', 'replay_wall_seconds': 0.0, 'replay_cpu_seconds': 0.0}
    replay = timed_call(lambda: verifier(deepcopy(certificate), expected_function=deepcopy(task['function']),
        expected_mean_zero=scope, expected_tolerance=task['tolerance']), REPLAY_SECONDS)
    initial.update(replay_wall_seconds=replay['wall_seconds'], replay_cpu_seconds=replay['cpu_seconds'],
                   replay_hard_timeout_enforced=replay['hard_timeout_enforced'])
    if not replay['execution_ok'] or replay.get('response') is not True:
        return {**initial, 'reason': 'independent_replay_rejected', 'replay_error': replay.get('reason')}
    try:
        value = (F(certificate['upper'])-F(certificate['lower']) if task['kind']=='spectrum' else F(certificate['error_upper']))
        if value < 0:
            raise ValueError('Negative certified bound')
        if 'independent_exact_value' in item:
            known = F(item['independent_exact_value'])
            if not F(certificate['lower']) <= known <= F(certificate['upper']):
                raise ValueError('Certificate excludes independent constant-potential exact value')
        return {**initial, 'certificate_valid': True, 'target_met': value<=F(task['tolerance']),
                'value': str(value), 'tolerance': task['tolerance'], 'reason': 'replayed'}
    except (ValueError, TypeError, KeyError, ZeroDivisionError) as error:
        return {**initial, 'reason': str(error)}


def extract_certificate(response):
    if type(response) is not dict:
        return None
    return response.get('certificate') if 'certificate' in response else response if 'format' in response else None


def text_cost(value):
    content = canonical(value)
    return {'characters': len(content), 'utf8_bytes': len(content.encode('utf-8')),
            'visible_token_count': None, 'tokenizer': None,
            'token_count_status': 'not_measured_by_this_evaluator; characters and bytes are not tokens',
            'actual_model_usage': None}


def execute_attempt(service, item, verifier, seconds=TASK_SECONDS, configuration=None):
    task = deepcopy(item['task'])
    call = timed_call(lambda: service(task), seconds)
    response = call.get('response')
    judgment = assess(extract_certificate(response), item, verifier)
    return {**call, 'configuration': configuration, 'assessment': judgment,
            'request_text': text_cost(task), 'response_text': text_cost(response) if response is not None else None,
            'accounted_wall_seconds': call['wall_seconds']+judgment['replay_wall_seconds'],
            'accounted_cpu_seconds': call['cpu_seconds']+judgment['replay_cpu_seconds']}


def legacy_service(configuration=None):
    def run(task):
        request = {'op': 'approximate' if task['kind']=='approximation' else 'spectrum', 'function': deepcopy(task['function'])}
        if task['kind']=='spectrum':
            request.update(mean_zero=task.get('mean_zero', True), tolerance=task['tolerance'])
        # The legacy approximation entry point DOES NOT accept tolerance.
        if configuration is not None:
            request.update(CONFIGS[configuration][task['kind']])
        return backend.legacy.handle(request)
    return run


def completed_row(system, item, attempts):
    valid = [i for i,a in enumerate(attempts) if a['assessment']['certificate_valid']]
    chosen = min(valid, key=lambda i: F(attempts[i]['assessment']['value'])) if valid else None
    assessment = attempts[chosen]['assessment'] if chosen is not None else {'certificate_valid': False, 'target_met': False}
    solver_seconds = sum(a['wall_seconds'] for a in attempts)
    return {'system': system, 'case_id': item['id'], 'classification': item['classification'], 'task': item['task'],
            'attempts': attempts, 'selected_attempt': chosen, 'assessment': assessment,
            'solver_wall_seconds': solver_seconds,
            'accounted_wall_seconds': sum(a['accounted_wall_seconds'] for a in attempts),
            'accounted_cpu_seconds': sum(a['accounted_cpu_seconds'] for a in attempts),
            'budget_met': solver_seconds <= TASK_SECONDS,
            'success': assessment['target_met'] and solver_seconds <= TASK_SECONDS,
            'independent_exact_value': item.get('independent_exact_value')}


def fixed_row(item, config, system):
    attempt = execute_attempt(legacy_service(config), item, default_verifier, configuration='default' if config is None else CONFIGS[config]['name'])
    return completed_row(system, item, [attempt])


def increasing_row(item):
    attempts = []
    for index in range(len(CONFIGS)):
        remaining = TASK_SECONDS-sum(a['wall_seconds'] for a in attempts)
        if remaining<=0:
            break
        attempt = execute_attempt(legacy_service(index), item, default_verifier, remaining, CONFIGS[index]['name'])
        attempts.append(attempt)
        if attempt['assessment']['target_met']:
            break
    return completed_row('fixed_increasing', item, attempts)


def validation_selection():
    path = DATA/'validation_selection.json'
    if path.exists():
        result = json.loads(path.read_text())
        if result['protocol_digest'] != digest(public_rules()):
            raise RuntimeError('Stored baseline selection belongs to another protocol')
        return result
    candidates = []
    for index, configuration in enumerate(CONFIGS):
        rows = [fixed_row(item, index, 'validation_candidate_'+configuration['name']) for item in validation_cases()]
        candidates.append({'config_index': index, 'name': configuration['name'], 'rows': rows,
                           'successes': sum(r['success'] for r in rows),
                           'valid': sum(r['assessment']['certificate_valid'] for r in rows)})
    # Success, then valid evidence, then a smaller declared configuration. No noisy timing tie-break.
    selected = max(candidates, key=lambda c: (c['successes'], c['valid'], -c['config_index']))
    result = {'protocol_digest': digest(public_rules()), 'selected_config_index': selected['config_index'],
              'selected_name': selected['name'], 'selection_rule': 'maximize validation target successes, then valid evidence, then choose smallest declared configuration',
              'all_candidates': candidates, 'preparation_wall_seconds': sum(r['accounted_wall_seconds'] for c in candidates for r in c['rows']),
              'preparation_cpu_seconds': sum(r['accounted_cpu_seconds'] for c in candidates for r in c['rows']),
              'all_failed_and_successful_branches_retained': True}
    write_json(path, result)
    return result


def mutation_probes(rows, cases, verifier):
    probes = []
    by_id = {item['id']: item for item in cases}
    # One representative valid certificate per system/conclusion kind limits replay expense.
    seen = set()
    for row in rows:
        key = (row['system'], row['task']['kind'])
        if key in seen or row['selected_attempt'] is None:
            continue
        seen.add(key)
        cert = extract_certificate(row['attempts'][row['selected_attempt']]['response'])
        item = by_id[row['case_id']]
        use_verifier = verifier if row['system']=='new_service' else default_verifier
        variants = []
        wrong_function = deepcopy(cert)
        wrong_function['function']['offset'] = str(F(cert['function']['offset'])+1)
        variants.append(('wrong_original_function', wrong_function))
        wrong_bound = deepcopy(cert)
        if row['task']['kind']=='spectrum':
            wrong_bound['lower'] = str(F(cert['upper'])+1)
            wrong_space = deepcopy(cert)
            wrong_space['mean_zero'] = not row['task']['mean_zero']
            variants.append(('wrong_function_space', wrong_space))
        else:
            wrong_bound['error_upper'] = '-1'
        variants.append(('impossible_bound', wrong_bound))
        variants.append(('status_only_forgery', {'format': cert['format'], 'function': cert['function'],
                         'mean_zero': row['task'].get('mean_zero'), 'verified': True, 'status': 'target_met'}))
        for name, mutation in variants:
            judgment = assess(mutation, item, use_verifier)
            probes.append({'system': row['system'], 'kind': row['task']['kind'], 'case_id': row['case_id'],
                           'mutation': name, 'accepted': judgment['certificate_valid'], 'assessment': judgment})
    return probes


def summaries(rows, probes, selection, phase):
    output = {}
    for system in dict.fromkeys(r['system'] for r in rows):
        group = [r for r in rows if r['system']==system]
        checks = [p for p in probes if p['system']==system]
        setup_wall = selection['preparation_wall_seconds'] if system=='validation_selected_fixed' else 0.0
        setup_cpu = selection['preparation_cpu_seconds'] if system=='validation_selected_fixed' else 0.0
        replay_wall = sum(p['assessment']['replay_wall_seconds'] for p in checks)
        replay_cpu = sum(p['assessment']['replay_cpu_seconds'] for p in checks)
        quantities = {}
        for kind in ('approximation','spectrum'):
            subset = [r for r in group if r['task']['kind']==kind]
            quantities[kind] = {'count': len(subset), 'successes': sum(r['success'] for r in subset),
                                'valid_certificates': sum(r['assessment']['certificate_valid'] for r in subset),
                                'values': [{'case_id': r['case_id'], 'value': r['assessment'].get('value'), 'tolerance': r['task']['tolerance']} for r in subset]}
        characters = sum(a['request_text']['characters']+(a['response_text']['characters'] if a['response_text'] else 0) for r in group for a in r['attempts'])
        output[system] = {
            'certified_accuracy': quantities,
            'reliability': {'probes': len(checks), 'falsely_accepted': sum(p['accepted'] for p in checks), 'replay_is_independent_of_solver_flags': True,
                            'mathematical_verifier_code_shared_with_solver_backend': True, 'formal_proof_assistant': False},
            'automatic_selection': {'successes': sum(r['success'] for r in group), 'tasks': len(group),
                                    'solver_attempts': sum(len(r['attempts']) for r in group),
                                    'validation_chosen_config': selection['selected_name'] if system=='validation_selected_fixed' else None,
                                    'retrospective_per_task_oracle_used': False, 'new_service_internal_attempts_retained_in_response': True},
            'new_parameter_transfer': {'phase': phase, 'evaluated': phase=='holdout', 'new_mathematical_families_tested': 0,
                                       'scope': 'within_family_only', 'fixed_minus_quarter_case_is_regression': True,
                                       'holdout_successes': sum(r['success'] for r in group) if phase=='holdout' else None},
            'runtime_resources': {'online_wall_seconds': sum(r['accounted_wall_seconds'] for r in group),
                                  'online_cpu_seconds': sum(r['accounted_cpu_seconds'] for r in group),
                                  'charged_selection_wall_seconds': setup_wall, 'charged_selection_cpu_seconds': setup_cpu,
                                  'reliability_probe_wall_seconds': replay_wall, 'reliability_probe_cpu_seconds': replay_cpu,
                                  'total_accounted_wall_seconds': setup_wall+replay_wall+sum(r['accounted_wall_seconds'] for r in group),
                                  'total_accounted_cpu_seconds': setup_cpu+replay_cpu+sum(r['accounted_cpu_seconds'] for r in group),
                                  'peak_memory_bytes': None, 'development_model_usage': None,
                                  'timing_warning': 'Other development agents may be active; timing is descriptive and not an isolated speed benchmark.',
                                  'cold_start_or_cache': 'one process, deterministic interleaving; imports excluded; no warm-cache speed claim'},
            'visible_tokens': {'token_count': None, 'tokenizer': None, 'actual_model_input_output_reasoning_usage': None,
                               'serialized_request_response_characters': characters,
                               'scope': 'full task-adapter evidence, not summary-view UI or actual model consumption',
                               'status': 'unknown; byte/character counts are reported separately, never treated as tokens'},
            'usability': {'execution_failures': sum(not a['execution_ok'] for r in group for a in r['attempts']),
                          'missing_or_invalid_final_certificates': sum(not r['assessment']['certificate_valid'] for r in group),
                          'timeouts': sum(a.get('timed_out', False) for r in group for a in r['attempts']),
                          'manual_per_task_configuration': False, 'certificate_data_retained': True},
        }
    return output


def run_evaluation(service, phase='development', verifier=None):
    """Evaluate task -> full certificate/envelope. None runs baseline-only.

Use ResearchService.evaluate_task directly with verification.verify_any. The
optional verifier must support expected_function/mean_zero/tolerance keywords.
The default intentionally only verifies frozen legacy certificate formats.
"""
    if phase not in ('development','validation','holdout'):
        raise ValueError('phase must be development, validation or holdout')
    if service is not None and not callable(service):
        raise TypeError('service must be a callable task adapter or None for baseline-only')
    initialize()
    if phase=='holdout':
        checked_manifest()
    selection = validation_selection()
    cases = holdout_cases() if phase=='holdout' else development_cases() if phase=='development' else validation_cases()
    verifier = verifier or default_verifier
    rows = []
    # Case-level interleaving is declared in advance; no hidden retries or dropping failures.
    for item in cases:
        rows.append(fixed_row(item, None, 'old_default'))
        rows.append(increasing_row(item))
        rows.append(fixed_row(item, selection['selected_config_index'], 'validation_selected_fixed'))
        if service is not None:
            rows.append(completed_row('new_service', item, [execute_attempt(service, item, verifier)]))
    probes = mutation_probes(rows, cases, verifier)
    report = {'version': VERSION, 'phase': phase, 'created_utc': datetime.now(timezone.utc).isoformat(),
              'protocol_digest': digest(public_rules()), 'environment': {'python': sys.version, 'platform': platform.platform(),
                'processor': platform.processor(), 'logical_cpus': os.cpu_count(), 'active_other_agents': 'unknown',
                'network_calls': 0, 'external_model_calls': 0, 'actual_model_usage': None},
              'cases': cases, 'rows': rows, 'reliability_probes': probes,
              'validation_selection_reference': str(DATA/'validation_selection.json'),
              'dimensions': summaries(rows, probes, selection, phase),
              'limitations': ['Shared directory offers discipline isolation, not permission isolation.',
                'Evaluator is independent of candidate selection but mathematical replay shares backend formulas/code.',
                'Full development and maintenance costs, model consumption and peak memory are unknown.',
                'A small deterministic suite cannot establish broad generalization or cross-model superiority.',
                'Changing code after using holdout turns these instances into development data.']}
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    path = DATA/(phase+'_'+('baselines' if service is None else 'comparison')+'_'+stamp+'.json')
    write_json(path, report)
    report['report_path'] = str(path)
    write_json(DATA/('latest_'+phase+'.json'), {'report_path': str(path), 'phase': phase,
        'protocol_digest': report['protocol_digest'], 'dimensions': report['dimensions']})
    return report


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--phase', choices=['development','validation','holdout'], default='development')
    parser.add_argument('--new-service', action='store_true')
    parser.add_argument('--initialize-only', action='store_true')
    args = parser.parse_args()
    if args.initialize_only:
        print(canonical(initialize()))
        return
    service, verifier = None, None
    if args.new_service:
        from research_service import ResearchService
        from verification import verify_any
        service = ResearchService(evidence_root=DATA/'service_evidence', cache=False).evaluate_task
        verifier = verify_any
    result = run_evaluation(service, args.phase, verifier)
    print(json.dumps({'report_path': result['report_path'], 'phase': result['phase'], 'dimensions': result['dimensions']}, ensure_ascii=False, indent=2))


if __name__=='__main__':
    main()
