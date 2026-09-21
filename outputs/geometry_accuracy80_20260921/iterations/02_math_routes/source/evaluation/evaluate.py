"""Independent accuracy evaluation; draft only, with no seed/holdout access.

Run with the documented Python and -B. Imports here are standard library only.
One pristine worker handles one method/case/repetition, and all runs are serial.
"""
from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
from datetime import datetime, timezone
from fractions import Fraction
import argparse
import hashlib
import importlib
import json
import os
from pathlib import Path
import resource
import selectors
import signal
import statistics
import subprocess
import sys
import tempfile
import threading
import time

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / 'evaluation'
RUNTIME = ROOT.parent / 'geometry_strategy_v317_runtime'
OLD = ROOT.parent / 'geometry_strategy_v238_v317'
FIRST = ROOT.parent / 'geometry_strategy_v317_release'
VERSION = 'accuracy80_evaluation_v1_parameters_approved'
METHODS = ('current', 'frozen_runtime')
TASK_SECONDS = 10.0
SETUP_SECONDS = 30.0


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False)


def forbidden(path):
    path = Path(path).resolve()
    if 'seed' in path.name.lower():
        raise RuntimeError('This draft evaluator never reads or writes seed files')
    return path


def read_json(path):
    return json.loads(forbidden(path).read_text())


def write_json(path, value):
    path = forbidden(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(canonical(value) + '\n')
    temporary.replace(path)


def file_digest(path):
    return hashlib.sha256(forbidden(path).read_bytes()).hexdigest()


def rational(value):
    if type(value) not in (int, str):
        raise ValueError('Bounds and task parameters must be exact integer/rational strings')
    return Fraction(value)


def normalize_function(function):
    """Evaluator-owned normalization, independent of solver-provided helpers."""
    if type(function) is not dict or function.get('kind') != 'axis_profile':
        raise ValueError('This preregistered family requires axis_profile')
    profile = function.get('profile')
    if profile not in ('abs_power', 'step') or type(function.get('axis', 2)) is not int:
        raise ValueError('Unsupported profile or axis')
    if function.get('axis', 2) not in (0, 1, 2):
        raise ValueError('Invalid axis')
    expected = {'kind', 'axis', 'profile', 'amplitude', 'offset'}
    if profile == 'abs_power':
        expected.add('exponent')
    if set(function) - expected:
        raise ValueError('Unknown function fields')
    result = {'kind': 'axis_profile', 'axis': function.get('axis', 2), 'profile': profile,
              'amplitude': str(rational(function.get('amplitude', '1'))),
              'offset': str(rational(function.get('offset', '0')))}
    if profile == 'abs_power':
        result['exponent'] = str(rational(function['exponent']))
    return result


def validate_task(task):
    if type(task) is not dict or task.get('kind') not in ('spectrum', 'approximation'):
        raise ValueError('Task must specify spectrum or approximation')
    normalize_function(task['function'])
    if not 0 < rational(task['tolerance']) <= 1:
        raise ValueError('Invalid target tolerance')
    if task['kind'] == 'spectrum' and type(task.get('mean_zero')) is not bool:
        raise ValueError('Spectrum must explicitly bind full or mean-zero space')
    if task['kind'] == 'approximation' and 'mean_zero' in task:
        raise ValueError('Approximation tasks do not have a spectral space')


def instance_key(task):
    """A changed tolerance never makes an old mathematical instance new."""
    validate_task(task)
    return canonical({'kind': task['kind'], 'function': normalize_function(task['function']),
                      'mean_zero': task.get('mean_zero') if task['kind'] == 'spectrum' else None})


def eta_squared(function):
    q = normalize_function(function)
    if q['profile'] != 'abs_power':
        return None
    a, alpha = rational(q['amplitude']), rational(q['exponent'])
    if not Fraction(-1, 2) < alpha < 0:
        raise ValueError('Power family requires -1/2 < alpha < 0')
    return a*a*alpha*alpha / ((1+2*alpha)*(1+alpha)**2)


def public_catalog():
    """Explicit public paths only; never enumerate or inspect private files."""
    paths = [OLD/'evaluation/holdout_instances.json',
             FIRST/'evaluation/holdout_instances.json',
             RUNTIME/'evaluation/holdout_instances.json',
             RUNTIME/'evaluation/validation_cases.json']
    result = []
    for path in paths:
        if path.exists():
            for item in read_json(path):
                result.append({'source': str(path), 'case_id': item['id'], 'task': item['task']})
    seen_path = DATA/'seen_development_cases.json'
    if seen_path.exists():
        result.extend(read_json(seen_path))
    return result


def classify_cases(cases, known=None):
    known = public_catalog() if known is None else known
    by_key, within, result = {}, {}, []
    for item in known:
        by_key.setdefault(instance_key(item['task']), []).append(
            {'source': item['source'], 'case_id': item['case_id']})
    identifiers = set()
    for source in cases:
        item = deepcopy(source)
        if not isinstance(item.get('id'), str) or item['id'] in identifiers:
            raise ValueError('Case IDs must be unique strings')
        identifiers.add(item['id'])
        key = instance_key(item['task'])
        item['budget_seconds'] = TASK_SECONDS
        item['prior_public_overlap'] = by_key.get(key, [])
        item['duplicate_of_current_draw'] = within.get(key)
        item['novel_relative_to_catalog'] = key not in by_key and key not in within
        # This command accepts disclosed development/regression tasks only.
        item['independent_new_instance'] = False
        item['classification'] = 'disclosed_development_or_regression'
        if item['task']['function'].get('profile') == 'abs_power':
            item['centered_eta_squared'] = str(eta_squared(item['task']['function']))
            item['strong_eta_ge_one'] = eta_squared(item['task']['function']) >= 1
        within.setdefault(key, item['id'])
        result.append(item)
    return result


def regression_cases():
    source = RUNTIME/'evaluation/holdout_instances.json'
    originals = read_json(source)
    wanted = {'H'+str(i).zfill(2) for i in range(5, 15)}
    cases = [deepcopy(c) for c in originals if c['id'] in wanted]
    if len(cases) != 10 or {c['id'] for c in cases} != wanted:
        raise RuntimeError('All original H05-H14 are required, without substitutions')
    for case in cases:
        case['historical_source'] = str(source)
        case['historical_source_sha256'] = file_digest(source)
        case['budget_seconds'] = TASK_SECONDS
    return classify_cases(cases)


def draft_protocol():
    return {
        'status': 'LOCKED_RULES_NO_SEED_OR_FORMAL_HOLDOUT', 'version': VERSION,
        'goal_pending_user_choice': True, 'default_metric': 'at least 16/20 independently new difficult instances succeed in all five cold repeats',
        'root_parameter_approval': '2026-09-21: fixed pools, layout and pre-sampling identity-only exclusion approved; preserve strong distribution regardless of failures',
        'parameter_generalization_scope': 'same known mathematical families, new parameters; no new mathematical structures or open-problem resolution claim',
        'provisional_layout': {
            'total': 20,
            'power_spectrum': {'count': 12, 'layout': '3 tolerances x 2 spaces x 2 centered-L2 strength strata',
                               'spaces': ['full', 'mean_zero'], 'strengths': ['eta_squared<1', 'eta_squared>=1']},
            'general_power_approximation': {'count': 6, 'layout': '2 draws at each tolerance', 'exponent': '-1/2<alpha<0, alpha != -1/4'},
            'step_spectrum': {'count': 2, 'spaces': ['full', 'mean_zero'], 'tolerance': '1/100000000'},
            'tolerances': ['1/1000000', '1/100000000', '1/10000000000'],
            'parameter_lists_and_sampler': parameter_rules_proposal(),
        },
        'methods': list(METHODS), 'task_seconds': TASK_SECONDS, 'initialization_limit_seconds': SETUP_SECONDS,
        'repetitions_provisional': 5, 'development_repetitions_allowed': [1, 3, 5],
        'single_repeat_rule': 'One run is a development smoke test only, never stable success or formal acceptance',
        'success': 'Every independent repeat must finish within the fixed deadline, serialize complete evidence, bind the original function/kind/space, independently replay, and meet the exact <= tolerance bound',
        'budget_scope': 'solve + first-use lazy imports + complete JSON encode/decode + independent verify + final raw-row serialization/write; initialization separately bounded and charged to complete cold wall/CPU',
        'hard_deadline': 'Unix main-thread SIGALRM raises BaseException; parent watchdog uses worker-start monotonic timestamp and kills the entire process group at the same deadline',
        'cold_design': 'fresh interpreter for every method/case/repeat; serial runs; alternating method order, rotating case order; no caches shared across calls',
        'failure_denominator': 'unsupported, exception, timeout, invalid evidence, loose bound, initialization or worker failure all remain; no filtering on success',
        'duplicate_identity': 'normalized exact function + conclusion kind + spectral space; tolerance ignored; preserve duplicates and show catalog overlaps',
        'primary_metric': '20 new unique nonoverlapping instances required; at least 16 stable successes; separately disclose old10 rate, per-family/scope/strength/precision rate',
        'repeat_statistics': 'all raw runs and 5/5 stable case success; repeats are not new mathematical instances; min/median/max wall/CPU and paired wins/losses descriptive only',
        'cold_accounting': 'parent-observed startup/import/setup/solve/replay/serialization/shutdown; OS child CPU; separate initialization/online stage CPU and wall; no warm/cold mixing',
        'seen_data': 'Any revealed or tuned-on instance becomes development/regression permanently, including different tolerances; source identity before/after each run must match',
        'stagnation': '5 consecutive genuine code iterations without improved stable development success trigger documented route review; no reset by version renaming; partial widths are diagnostics, not success gains',
        'tokens': {'actual_model_identity': None, 'actual_input_output_reasoning_tokens': None,
                   'actual_monetary_cost': None, 'character_or_byte_token_proxy': False},
        'formal_release_state': 'No seed access, generation, freeze or formal holdout execution is implemented by this draft',
    }


def parameter_rules_proposal():
    weak = ['-1/7', '-3/13', '-4/15', '-5/17', '-5/14', '-7/18']
    strong = ['-5/12', '-7/16', '-9/20', '-11/24', '-12/25', '-13/27']
    return {'status': 'ROOT_APPROVED_FIXED_RULES_NO_SEED_OR_DRAW',
        'axis': 2, 'offset': '0',
        'power_weak': {'exponents': weak, 'amplitudes': [str(Fraction(-k, 23)) for k in range(6, 14)]},
        'power_strong': {'exponents': strong, 'amplitudes': [str(Fraction(-k, 23)) for k in range(14, 23)]},
        'power_approximation': {'exponents': weak+strong, 'amplitudes': [str(Fraction(-k, 23)) for k in range(6, 23)]},
        'step_spectrum': {'amplitudes': [str(Fraction(-k, 23)) for k in range(6, 23)]},
        'layout_order': 'For tolerances 1e-6,1e-8,1e-10: power full weak,strong, mean_zero weak,strong; then 2 approximation draws per tolerance; then step full,mean_zero at 1e-8',
        'sampler': 'Freeze public seen-task catalog before commitment. Enumerate sorted Cartesian pools; remove only identities matching prior public function+kind+space (ignore tolerance). SHA256(seed || counter_as_8_byte_big_endian) rejection sampling gives an unbiased pool index; increment counter for every word; remove selected identity from every same-kind/same-space pool regardless of tolerance. No solver calls, support checks, outcome checks or post-reveal replacement.',
        'duplicate_policy': 'Pure novelty exclusion occurs before any solve. Audit every excluded identity and the fixed catalog hash; preserve all selected draws including every failure. The result must contain 20 distinct new identities or no formal acceptance is possible.',
        'no_new_mathematical_structure_claim': True}


class EvaluationTimeout(BaseException):
    pass


def timed_call(function, seconds):
    if not hasattr(signal, 'SIGALRM') or threading.current_thread() is not threading.main_thread():
        raise RuntimeError('Unix main-thread SIGALRM required; no cooperative fallback')
    previous_handler = signal.getsignal(signal.SIGALRM)
    old_timer = signal.setitimer(signal.ITIMER_REAL, 0)
    if old_timer[0]:
        signal.setitimer(signal.ITIMER_REAL, *old_timer)
        raise RuntimeError('Nested deadlines forbidden')
    def alarm(*_):
        raise EvaluationTimeout('hard wall deadline exhausted')
    signal.signal(signal.SIGALRM, alarm)
    try:
        if seconds <= 0:
            raise EvaluationTimeout('no remaining budget')
        signal.setitimer(signal.ITIMER_REAL, seconds)
        value = function()
        return {'ok': True, 'value': value, 'timed_out': False}
    except (Exception, EvaluationTimeout) as error:
        return {'ok': False, 'timed_out': isinstance(error, EvaluationTimeout),
                'error_type': type(error).__name__, 'reason': str(error)}
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous_handler)


def extract_certificate(response):
    if type(response) is not dict:
        return None
    return response.get('certificate', response if 'format' in response else None)


def assess(certificate, item, verifier):
    task = deepcopy(item['task'])
    invalid = {'certificate_valid': False, 'target_met': False}
    if type(certificate) is not dict:
        return {**invalid, 'reason': 'missing_certificate'}
    if certificate.get('function') != normalize_function(task['function']):
        return {**invalid, 'reason': 'original_function_mismatch'}
    scope = task.get('mean_zero') if task['kind'] == 'spectrum' else None
    if scope is not None and certificate.get('mean_zero') is not scope:
        return {**invalid, 'reason': 'original_space_mismatch'}
    if certificate.get('kind', task['kind']) != task['kind']:
        return {**invalid, 'reason': 'original_kind_mismatch'}
    if 'tolerance' in certificate and rational(certificate['tolerance']) != rational(task['tolerance']):
        return {**invalid, 'reason': 'embedded_target_mismatch'}
    if verifier(deepcopy(certificate), expected_function=deepcopy(task['function']),
                expected_mean_zero=scope, expected_tolerance=task['tolerance'],
                expected_kind=task['kind']) is not True:
        return {**invalid, 'reason': 'independent_replay_rejected'}
    try:
        if task['kind'] == 'spectrum':
            lower, upper = rational(certificate['lower']), rational(certificate['upper'])
            value = upper-lower
            if 'independent_exact_value' in item and not lower <= rational(item['independent_exact_value']) <= upper:
                raise ValueError('excludes_independent_exact_value')
        else:
            value = rational(certificate['error_upper'])
            if 'error_squared' in certificate:
                squared = rational(certificate['error_squared'])
                if squared < 0 or value*value < squared:
                    raise ValueError('L2_error_upper_does_not_cover_squared_error')
        if value < 0:
            raise ValueError('negative_bound')
        return {'certificate_valid': True, 'target_met': value <= rational(task['tolerance']),
                'value': str(value), 'tolerance': str(rational(task['tolerance'])), 'reason': 'independently_replayed',
                'quantity': 'spectral_interval_width' if task['kind'] == 'spectrum' else 'function_L2_error_upper'}
    except (ValueError, TypeError, KeyError, ArithmeticError) as error:
        return {**invalid, 'reason': str(error)}


def load_adapter(method, evidence, selftest=False):
    sys.dont_write_bytecode = True
    if method.startswith('_selftest_'):
        if not selftest:
            raise RuntimeError('Synthetic adapters allowed only in boundary selftests')
        from test_boundaries import synthetic_adapter
        return synthetic_adapter(method)
    # Prevent same-name release modules from crossing versions.
    roots = {ROOT.resolve(), RUNTIME.resolve(), FIRST.resolve(), OLD.resolve()}
    sys.path[:] = [p for p in sys.path if Path(p or os.getcwd()).resolve() not in roots]
    if method == 'current':
        sys.path.insert(0, str(ROOT))
        module = importlib.import_module('solver')
        if Path(module.__file__).resolve() != (ROOT/'solver.py').resolve():
            raise RuntimeError('Wrong current solver module imported')
        return module.solve, module.verify_certificate, {'solver': str(Path(module.__file__).resolve())}
    if method == 'frozen_runtime':
        sys.path[:0] = [str(RUNTIME), str(OLD)]
        module = importlib.import_module('release_service')
        if Path(module.__file__).resolve() != (RUNTIME/'release_service.py').resolve():
            raise RuntimeError('Wrong frozen runtime module imported')
        service = module.ResearchService(evidence_root=evidence, cache=False)
        from verification import verify_any
        return service.evaluate_task, verify_any, {'solver': str(Path(module.__file__).resolve())}
    raise ValueError('Unknown method: ' + method)


def evaluate_worker(spec):
    """Single cold case; stdout carries only parent-watchdog protocol events."""
    event_stream = sys.stdout
    def event(value):
        event_stream.write(canonical(value)+'\n')
        event_stream.flush()
    item = spec['case']
    validate_task(item['task'])
    if not spec.get('selftest') and item['budget_seconds'] != TASK_SECONDS:
        raise ValueError('Production cases always receive exactly ten seconds')
    row = {'method': spec['method'], 'case_id': item['id'], 'task': item['task'],
           'family': item.get('family', 'unspecified'), 'repetition': spec['repetition'],
           'success': False, 'assessment': {'certificate_valid': False, 'target_met': False},
           'stage': 'initialization', 'timings': {}, 'response': None}
    wall, cpu = time.monotonic(), time.process_time()
    with tempfile.TemporaryDirectory(prefix='geometry-accuracy80-') as evidence, open(spec['log'], 'w') as logs:
        with redirect_stdout(logs), redirect_stderr(logs):
            setup = timed_call(lambda: load_adapter(spec['method'], evidence, spec.get('selftest', False)), SETUP_SECONDS)
            row['timings']['initialization'] = {'wall_seconds': time.monotonic()-wall, 'cpu_seconds': time.process_time()-cpu}
            if not setup['ok']:
                row['execution'] = setup
                write_json(spec['output'], row)
                event({'event': 'setup_failed'})
                return
            adapter, verifier, implementation = setup['value']
            row['implementation'] = implementation
            started = time.monotonic()
            event({'event': 'task_started', 'monotonic': started})
            online_cpu = time.process_time()
            def stage(name, function):
                row['stage'] = name
                write_json(spec['checkpoint'], {k:v for k,v in row.items() if k != 'response'})
                before, cpu_before = time.monotonic(), time.process_time()
                try:
                    return function()
                finally:
                    row['timings'][name] = {'wall_seconds': time.monotonic()-before,
                                            'cpu_seconds': time.process_time()-cpu_before}
            def work():
                response = stage('solve', lambda: adapter(deepcopy(item['task'])))
                # Transport full response through strict JSON before verification.
                row['response'] = stage('response_json_roundtrip', lambda: json.loads(canonical(response)))
                row['assessment'] = stage('independent_replay',
                    lambda: assess(extract_certificate(row['response']), item, verifier))
                row['success'] = row['assessment']['target_met']
                row['execution'] = {'ok': True, 'timed_out': False}
                row['stage'] = 'final_raw_serialization'
                write_json(spec['output'], row)
            execution = timed_call(work, item['budget_seconds']-(time.monotonic()-started))
            finished = time.monotonic()
            event({'event': 'task_finished', 'monotonic': finished,
                   'online_wall_seconds': finished-started, 'online_cpu_seconds': time.process_time()-online_cpu,
                   'execution_ok': execution['ok'], 'timed_out': execution['timed_out']})
            if not execution['ok']:
                row['success'] = False
                row['execution'] = execution
                # Error-path output is post-deadline administrative cost; never a success.
                try:
                    write_json(spec['output'], row)
                except (Exception, EvaluationTimeout):
                    row['response'] = None
                    write_json(spec['output'], row)


def run_worker(method, item, repetition, destination, selftest=False):
    """The parent enforces the task deadline even if a child catches BaseException."""
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    stem = method+'_r'+str(repetition)+'_'+item['id']
    paths = {k: destination/(stem+'_'+k+suffix) for k,suffix in
             (('spec','.json'), ('output','.json'), ('checkpoint','.json'), ('log','.txt'))}
    if paths['spec'].exists():
        raise RuntimeError('Never overwrite or silently reuse a measured run')
    spec = {'method': method, 'case': item, 'repetition': repetition,
            'output': str(paths['output']), 'checkpoint': str(paths['checkpoint']),
            'log': str(paths['log']), 'selftest': selftest}
    write_json(paths['spec'], spec)
    start = time.monotonic()
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0')
    events, output_buffer, protocol_errors = [], b'', []
    watchdog = None
    with open(destination/(stem+'_worker_stderr.txt'), 'wb') as errors:
        process = subprocess.Popen([sys.executable, '-B', str(Path(__file__).resolve()),
                                    'worker', '--spec', str(paths['spec'])],
            cwd=str(ROOT), env=environment, stdout=subprocess.PIPE, stderr=errors, start_new_session=True)
        deadline, stage = start+SETUP_SECONDS, 'initialization'
        selector = selectors.DefaultSelector()
        selector.register(process.stdout, selectors.EVENT_READ)
        while True:
            remaining = deadline-time.monotonic()
            if remaining <= 0:
                watchdog = stage
                os.killpg(process.pid, signal.SIGKILL)
                break
            ready = selector.select(min(remaining, 0.1))
            for key, _ in ready:
                chunk = os.read(key.fileobj.fileno(), 65536)
                if not chunk:
                    selector.unregister(key.fileobj)
                    continue
                output_buffer += chunk
                while b'\n' in output_buffer:
                    line, output_buffer = output_buffer.split(b'\n', 1)
                    try:
                        message = json.loads(line)
                        events.append(message)
                        if message.get('event') == 'task_started':
                            stamp = message['monotonic']
                            if stage != 'initialization' or not start <= stamp <= time.monotonic():
                                raise ValueError('Invalid worker task-start event')
                            stage, deadline = 'online_task', stamp+item['budget_seconds']
                        elif message.get('event') == 'task_finished':
                            if stage != 'online_task':
                                raise ValueError('Unexpected task-finished event')
                            if message['monotonic'] > deadline and message.get('execution_ok'):
                                protocol_errors.append('successful_completion_after_deadline')
                            stage, deadline = 'shutdown', time.monotonic()+5.0
                        elif message.get('event') == 'setup_failed':
                            stage, deadline = 'shutdown', time.monotonic()+5.0
                        else:
                            raise ValueError('Unknown worker event')
                    except (ValueError, KeyError, TypeError) as error:
                        protocol_errors.append(str(error))
            if process.poll() is not None and not selector.get_map():
                break
        selector.close()
        process.wait()
        process.stdout.close()
    end = time.monotonic()
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = read_json(paths['output']) if paths['output'].exists() else (
        read_json(paths['checkpoint']) if paths['checkpoint'].exists() else
        {'method': method, 'case_id': item['id'], 'task': item['task'], 'repetition': repetition,
         'assessment': {'certificate_valid': False, 'target_met': False}, 'timings': {}})
    finished_events = [e for e in events if e.get('event') == 'task_finished']
    done = finished_events[-1] if finished_events else None
    within = bool(done and done['execution_ok'] and not done['timed_out'] and
                  done['online_wall_seconds'] <= item['budget_seconds'] and not watchdog and
                  process.returncode == 0 and not protocol_errors)
    result['within_budget'] = within
    result['success'] = bool(result.get('success') and within)
    result['complete'] = within
    result['budget_seconds'] = item['budget_seconds']
    result['process'] = {'wall_seconds': end-start,
        'cpu_seconds': usage_after.ru_utime+usage_after.ru_stime-usage_before.ru_utime-usage_before.ru_stime,
        'returncode': process.returncode, 'watchdog_stage': watchdog, 'protocol_errors': protocol_errors,
        'events': events, 'unparsed_stdout_bytes': len(output_buffer), 'cold_interpreter': True}
    result['online_wall_seconds'] = done['online_wall_seconds'] if done else None
    result['online_cpu_seconds'] = done['online_cpu_seconds'] if done else None
    if not within:
        result['failure_reason'] = ('parent_watchdog_'+watchdog if watchdog else
                                   result.get('execution', {}).get('reason', 'worker_failed_or_incomplete'))
    elif not result['assessment']['target_met']:
        result['failure_reason'] = result['assessment'].get('reason', 'target_not_met')
    result['raw_paths'] = {k:str(v.relative_to(ROOT)) if v.is_relative_to(ROOT) else str(v) for k,v in paths.items()}
    write_json(destination/(stem+'_result.json'), result)
    return result


def distribution(values):
    return {'raw': values, 'median': statistics.median(values), 'min': min(values), 'max': max(values)} if values else None


def summarize(rows, cases, repetitions):
    summary = {}
    for method in dict.fromkeys(r['method'] for r in rows):
        method_rows = [r for r in rows if r['method'] == method]
        per_case = []
        for item in cases:
            group = [r for r in method_rows if r['case_id'] == item['id']]
            complete_success = len(group) == repetitions and len({r['repetition'] for r in group}) == repetitions and all(r['success'] for r in group)
            robust = repetitions >= 3 and complete_success
            per_case.append({'case_id': item['id'], 'family': item.get('family', 'unspecified'),
                'mean_zero': item['task'].get('mean_zero'), 'tolerance': item['task']['tolerance'],
                'strong_eta_ge_one': item.get('strong_eta_ge_one'), 'runs': len(group),
                'successful_runs': sum(r['success'] for r in group), 'robust_success': robust,
                'all_observed_runs_successful': complete_success,
                'cold_wall_seconds': distribution([r['process']['wall_seconds'] for r in group]),
                'cold_cpu_seconds': distribution([r['process']['cpu_seconds'] for r in group]),
                'online_wall_seconds': distribution([r['online_wall_seconds'] for r in group if r['online_wall_seconds'] is not None]),
                'online_cpu_seconds': distribution([r['online_cpu_seconds'] for r in group if r['online_cpu_seconds'] is not None]),
                'missing_online_timing_runs': sum(r['online_wall_seconds'] is None for r in group),
                'failure_reasons': [r.get('failure_reason') for r in group if not r['success']],
                'values': [r.get('assessment', {}).get('value') for r in group]})
        successes = sum(c['robust_success'] for c in per_case)
        breakdown = {}
        for field in ('family', 'mean_zero', 'tolerance', 'strong_eta_ge_one'):
            breakdown[field] = {}
            for value in dict.fromkeys(str(c[field]) for c in per_case):
                group = [c for c in per_case if str(c[field]) == value]
                breakdown[field][value] = {'cases': len(group), 'robust_successes': sum(c['robust_success'] for c in group)}
        totals = []
        for repetition in range(repetitions):
            group = [r for r in method_rows if r['repetition'] == repetition]
            totals.append({'repetition': repetition, 'cold_wall_seconds': sum(r['process']['wall_seconds'] for r in group),
                          'cold_cpu_seconds': sum(r['process']['cpu_seconds'] for r in group),
                          'online_wall_seconds_known': sum(r['online_wall_seconds'] or 0 for r in group),
                          'online_cpu_seconds_known': sum(r['online_cpu_seconds'] or 0 for r in group),
                          'unknown_online_timing_runs': sum(r['online_wall_seconds'] is None for r in group)})
        summary[method] = {'cases': len(cases), 'robust_successes': successes,
            'stability_evaluated': repetitions >= 3,
            'all_observed_runs_successful_cases': sum(c['all_observed_runs_successful'] for c in per_case),
            'robust_success_fraction': str(Fraction(successes, len(cases))) if cases else None,
            'observed_rate_at_least_80_percent': bool(cases) and 5*successes >= 4*len(cases),
            'formal_acceptance': False, 'reason': 'disclosed development/regression only; independent holdout not yet authorized',
            'per_case': per_case, 'breakdown': breakdown, 'repeat_totals': totals}
    if all(m in summary for m in METHODS):
        current = {c['case_id']: c['robust_success'] for c in summary['current']['per_case']}
        baseline = {c['case_id']: c['robust_success'] for c in summary['frozen_runtime']['per_case']}
        summary['paired_comparison'] = {'gained_cases': [k for k in current if current[k] and not baseline[k]],
                                       'lost_cases': [k for k in current if baseline[k] and not current[k]]}
    return summary


def source_identity():
    # Only code files, never seed material. Capture before any measured execution.
    roots = (ROOT, RUNTIME, OLD)
    files = {}
    for root in roots:
        for path in sorted(root.rglob('*.py')):
            files[str(path.relative_to(ROOT.parent))] = file_digest(path)
    return files


def run_evaluation(cases, phase, repetitions=5, methods=METHODS, label='run'):
    if phase not in ('regression', 'development'):
        raise RuntimeError('Formal holdout is unavailable until protocol/authorization are added')
    if repetitions not in (1, 3, 5):
        raise ValueError('Use one development smoke run, or three/five independent repeats')
    if not methods or any(m not in METHODS for m in methods):
        raise ValueError('Invalid method list')
    cases = classify_cases(cases)
    if phase == 'regression':
        originals = regression_cases()
        if {c['id']:c['task'] for c in cases} != {c['id']:c['task'] for c in originals}:
            raise ValueError('Regression must retain exactly the original ten tasks')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    destination = DATA/(phase+'_'+stamp+'_'+''.join(c for c in label if c.isalnum() or c in '-_'))
    destination.mkdir(parents=True, exist_ok=False)
    identity = source_identity()
    write_json(destination/'source_identity_before.json', identity)
    write_json(destination/'cases.json', cases)
    seen = read_json(DATA/'seen_development_cases.json') if (DATA/'seen_development_cases.json').exists() else []
    for item in cases:
        seen.append({'source': str(destination/'cases.json'), 'case_id': item['id'], 'task': item['task']})
    write_json(DATA/'seen_development_cases.json', seen)
    rows = []
    for repetition in range(repetitions):
        ordered = cases[repetition % len(cases):]+cases[:repetition % len(cases)] if cases else []
        method_order = list(methods[repetition % len(methods):])+list(methods[:repetition % len(methods)])
        for item in ordered:
            for method in method_order:
                rows.append(run_worker(method, item, repetition, destination))
                write_json(destination/'checkpoint.json', {'phase': phase, 'rows': rows, 'complete': False})
    after = source_identity()
    write_json(destination/'source_identity_after.json', after)
    report = {'version': VERSION, 'phase': phase, 'repetitions': repetitions, 'cases': cases,
        'source_unchanged_during_run': identity == after, 'summary': summarize(rows, cases, repetitions),
        'rows': rows, 'complete': True, 'formal_acceptance': False,
        'statistical_scope': 'descriptive repeated measurements on a shared host; no significance claim',
        'actual_model_identity': None, 'actual_model_tokens': None, 'token_proxy_used': False}
    if identity != after:
        report['invalid_for_comparison_reason'] = 'source_changed_during_run'
    write_json(destination/'report.json', report)
    return destination/'report.json'


def prepare_draft():
    if (DATA/'protocol.json').exists() and read_json(DATA/'protocol.json') != draft_protocol():
        raise RuntimeError('Published locked protocol differs; an explicit reviewed revision is required')
    write_json(DATA/'protocol_draft.json', draft_protocol())
    write_json(DATA/'protocol.json', draft_protocol())
    write_json(DATA/'regression_cases.json', regression_cases())
    write_json(DATA/'known_public_cases.json', public_catalog())
    return {'status': 'draft_prepared', 'formal_holdout_generated': False, 'seed_accessed': False,
            'regression_cases': 10}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    commands.add_parser('prepare-draft')
    commands.add_parser('selftest')
    worker_parser = commands.add_parser('worker')
    worker_parser.add_argument('--spec', required=True)
    run_parser = commands.add_parser('run')
    run_parser.add_argument('--phase', choices=('development', 'regression'), required=True)
    run_parser.add_argument('--cases')
    run_parser.add_argument('--repetitions', type=int, choices=(1, 3, 5), default=5)
    run_parser.add_argument('--methods', nargs='+', choices=METHODS, default=list(METHODS))
    run_parser.add_argument('--label', default='run')
    args = parser.parse_args()
    if args.command == 'worker':
        evaluate_worker(read_json(args.spec))
    elif args.command == 'prepare-draft':
        print(canonical(prepare_draft()))
    elif args.command == 'selftest':
        from test_boundaries import run_tests
        result = run_tests()
        write_json(DATA/'boundary_selftest.json', result)
        print(canonical(result))
        if not result['passed']:
            raise SystemExit(1)
    elif args.command == 'run':
        cases = read_json(args.cases) if args.cases else regression_cases()
        print(run_evaluation(cases, args.phase, args.repetitions, tuple(args.methods), args.label))


if __name__ == '__main__':
    main()
