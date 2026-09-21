"""Preregistered release evaluation. No holdout generation before explicit release.

Only the Python standard library is imported before a worker's setup timer.
Run with the documented Python and -B; all old files are read-only inputs.
"""
from copy import deepcopy
from datetime import datetime, timezone
from fractions import Fraction as F
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import secrets
import signal
import statistics
import subprocess
import sys
import tempfile
import threading
import time

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'geometry_strategy_v238_v317'
DATA = ROOT / 'evaluation'
VERSION = 'geometry_release_preregistered_v2'
REPETITIONS = 5
METHODS = ['old_default', 'fixed_increasing', 'validation_selected_fixed', 'frozen_v317', 'release']
CONFIGS = [
    {'name': 'small', 'approximation': {'levels': 1, 'degree': 1, 'sqrt_bits': 24},
     'spectrum': {'modes': 2, 'bits': 12, 'max_m': 1, 'sqrt_bits': 24, 'near_tail': 0}},
    {'name': 'medium', 'approximation': {'levels': 3, 'degree': 3, 'sqrt_bits': 28},
     'spectrum': {'modes': 3, 'bits': 18, 'max_m': 2, 'sqrt_bits': 28, 'near_tail': 0}},
    {'name': 'large', 'approximation': {'levels': 6, 'degree': 4, 'sqrt_bits': 32},
     'spectrum': {'modes': 5, 'bits': 26, 'max_m': 3, 'sqrt_bits': 32, 'near_tail': 0}},
]


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def file_digest(path):
    if Path(path).name == 'private_seed.json' and Path(path).resolve().is_relative_to(OLD.resolve()):
        raise RuntimeError('The old private seed must never be read')
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2, allow_nan=False) + '\n')
    temporary.replace(path)


def read_json(path):
    if Path(path).name == 'private_seed.json' and Path(path).resolve().is_relative_to(OLD.resolve()):
        raise RuntimeError('The old private seed must never be read')
    return json.loads(Path(path).read_text())


def now():
    return datetime.now(timezone.utc).isoformat()


def q(profile='step', amplitude='1', offset='0', exponent=None):
    result = {'kind': 'axis_profile', 'axis': 2, 'profile': profile,
              'amplitude': str(F(amplitude)), 'offset': str(F(offset))}
    if exponent is not None:
        result['exponent'] = str(F(exponent))
    return result


def case(identifier, kind, function, tolerance, mean_zero=True, family='easy', classification='public_validation'):
    task = {'kind': kind, 'function': function, 'tolerance': tolerance}
    if kind == 'spectrum':
        task['mean_zero'] = mean_zero
    result = {'id': identifier, 'task': task, 'family': family, 'classification': classification,
              'budget_seconds': 3.0 if family == 'easy' else 10.0}
    if kind == 'spectrum' and F(function['amplitude']) == 0:
        result['independent_exact_value'] = str(F(function['offset']) + (2 if mean_zero else 0))
    return result


def validation_cases():
    cases = [
        case('V01', 'spectrum', q(amplitude='0', offset='-2/5'), '1/10000000000', False),
        case('V02', 'spectrum', q(amplitude='0', offset='3/5'), '1/10000000000'),
        case('V03', 'spectrum', q(amplitude='1/3'), '1/1000'),
        case('V04', 'approximation', q(amplitude='3/4', offset='1/5'), '1/10000000000'),
        case('V05', 'spectrum', q('abs_power', '-1', exponent='-1/8'), '1/1000000', family='negative_power_spectrum'),
        case('V06', 'spectrum', q('abs_power', '-1/2', exponent='-1/3'), '1/100000000', family='negative_power_spectrum'),
        case('V07', 'spectrum', q('abs_power', '-1/2', exponent='-3/8'), '1/10000000000', family='negative_power_spectrum', classification='known_failed_structure_regression'),
        case('V08', 'spectrum', q('abs_power', '-1/3', exponent='-1/5'), '1/100000000', False, 'negative_power_spectrum', 'full_space_support_challenge'),
        case('V09', 'spectrum', q(amplitude='-1/2'), '1/100000000', family='step_spectrum_hard', classification='known_failed_structure_regression'),
        case('V10', 'approximation', q('abs_power', '-1', exponent='-1/4'), '1/1000000', family='negative_power_approximation', classification='fixed_minus_quarter_regression'),
        case('V11', 'approximation', q('abs_power', '-1', exponent='-1/4'), '1/100000000', family='negative_power_approximation', classification='fixed_minus_quarter_regression'),
        case('V12', 'approximation', q('abs_power', '-1', exponent='-1/4'), '1/10000000000', family='negative_power_approximation', classification='fixed_minus_quarter_regression'),
        case('V13', 'approximation', q('abs_power', '-1/3', exponent='-1/5'), '1/100000000', family='negative_power_approximation', classification='outside_current_approximation_support'),
    ]
    return cases


def regression_cases():
    """Only the already disclosed old instances, never the old private seed."""
    path = OLD / 'evaluation' / 'holdout_instances.json'
    cases = read_json(path) if path.exists() else []
    result = []
    for item in cases:
        item = deepcopy(item)
        item.update(id='OLD_' + item['id'], family='old_disclosed_regression',
                    classification='old_disclosed_holdout_now_regression', budget_seconds=3.0)
        result.append(item)
    return result


def public_rules():
    return {
        'version': VERSION, 'methods': METHODS, 'configurations': CONFIGS,
        'validation': validation_cases(), 'repetitions': REPETITIONS,
        'budgets': {'easy_seconds': 3, 'hard_seconds': 10, 'scope': 'solver + exact independent replay + attempt selection; same whole-task deadline for every method',
                    'setup_seconds_per_session': 30, 'mutation_seconds_each': 3,
                    'timeout': 'required Unix SIGALRM on main thread; BaseException + parent subprocess watchdog; no cooperative fallback'},
        'holdout': {
            'count': 14, 'parameter_scope': 'within known mathematical families only',
            'layout': ['constant_full', 'constant_mean_zero', 'easy_step_spectrum', 'step_approximation',
                       'negative_power_spectrum_1e-6_mean_zero', 'negative_power_spectrum_1e-6_full',
                       'negative_power_spectrum_1e-8_mean_zero', 'negative_power_spectrum_1e-8_full',
                       'negative_power_spectrum_1e-10_mean_zero', 'negative_power_spectrum_1e-10_full',
                       'negative_power_approximation_1e-6', 'negative_power_approximation_1e-8',
                       'negative_power_approximation_1e-10', 'negative_step_spectrum_1e-8'],
            'negative_exponents': ['-2/9', '-2/11', '-3/10', '-4/11', '-2/5', '-3/7', '-4/9'],
            'negative_amplitudes': ['-1/3', '-2/5', '-3/5'],
            'positive_step_amplitudes': ['1/7', '2/7', '3/7', '4/7'],
            'constant_offsets': ['-5/11', '-3/11', '2/11', '7/11'],
            'tolerances': ['1/1000000', '1/100000000', '1/10000000000'],
            'sampler': 'SHA256(seed || 8-byte counter), modulo each predeclared list; no success/stability/support rejection or resampling',
            'fixed_minus_quarter': 'regression only, never a new holdout',
            'failures': 'Unsupported exponents/spaces, invalid responses, timeout and loose certificates remain in denominators and full raw results',
        },
        'selection': '5 independent cold sessions/configuration on all 13 validation cases; maximize total budget-valid successes, then replay-valid certificates, then smallest predeclared configuration; never timing or per-task oracle',
        'session_design': 'new subprocess per method/repetition; methods rotated by repetition; same task order rotated by repetition; session-local caches allowed, no cache across sessions; all imports/constructor setup and first-use lazy imports charged',
        'replay': 'frozen v317 verification.verify_any with original function, conclusion kind, space and tolerance bound; independent constant eigenvalue inclusion; never trust returned status flags',
        'reliability': 'wrong original function, wrong space, impossible bound and status-only certificate rejected; reliability failure blocks improvement claim',
        'statistics': 'report 5 raw session totals, median/min/max; paired wins and losses against every baseline, per-case successes and 5/5 robust successes, family and scope breakdown; repeats are not independent mathematical tasks; small-sample timings descriptive',
        'reward_rules': {
            'accuracy': 'An additional robust heldout success over the strongest robust-success baseline, without lost successes or reliability failure, supports a limited within-family accuracy benefit.',
            'resources': 'At equal robust successes and no per-case success loss, >=20% median full-session wall and CPU reduction versus each equal-success baseline is a preliminary runtime benefit; also show online-only results and validation preparation charged to selected baseline once.',
            'no_aggregate': 'Report accuracy, reliability, runtime, automation, support, interaction and usability separately; no aggregate reward hides regressions.',
            'model_claims': 'actual model input/output/reasoning tokens, model identity, monetary cost and 30-min/8-hour comparison unknown; characters/bytes are not token proxies',
        },
    }


def initialize():
    """Generate only a private seed and public commitment; never read it back."""
    DATA.mkdir(parents=True, exist_ok=True)
    rules = public_rules()
    protocol_path = DATA / 'protocol.json'
    if protocol_path.exists() and read_json(protocol_path) != rules:
        raise RuntimeError('Published protocol differs; do not silently alter this evaluation')
    write_json(protocol_path, rules)
    write_json(DATA / 'validation_cases.json', validation_cases())
    private, commitment = DATA / 'private_seed.json', DATA / 'seed_commitment.json'
    if not private.exists():
        if commitment.exists() or (DATA / 'holdout_instances.json').exists():
            raise RuntimeError('Committed seed is missing; refusing replacement')
        seed = secrets.token_bytes(32)
        fd = os.open(private, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        with os.fdopen(fd, 'w') as stream:
            json.dump({'seed_hex': seed.hex()}, stream)
            stream.write('\n')
        write_json(commitment, {'sha256': hashlib.sha256(seed).hexdigest(), 'created_utc': now(),
                                'seed_read_before_explicit_release': False})
        del seed
    elif not commitment.exists():
        raise RuntimeError('Private seed exists without commitment; refusing to read or replace')
    return {'protocol_digest': digest(rules), 'commitment': read_json(commitment),
            'holdout_generated': False, 'private_seed_disclosed': False}


class EvaluationTimeout(BaseException):
    pass


def timed_call(function, seconds):
    """No silent timing fallback, and Exception handlers cannot swallow alarm."""
    if not hasattr(signal, 'SIGALRM') or threading.current_thread() is not threading.main_thread():
        raise RuntimeError('Evaluation requires Unix signals on the main thread')
    if seconds <= 0:
        return {'execution_ok': False, 'timed_out': True, 'reason': 'budget_exhausted',
                'wall_seconds': 0.0, 'cpu_seconds': 0.0, 'hard_timeout_enforced': True}
    start, cpu = time.perf_counter(), time.process_time()
    handler = signal.getsignal(signal.SIGALRM)
    old_timer = signal.setitimer(signal.ITIMER_REAL, 0)
    if old_timer[0]:
        signal.setitimer(signal.ITIMER_REAL, *old_timer)
        raise RuntimeError('Nested deadlines are forbidden')
    def alarm(*_):
        raise EvaluationTimeout('whole-stage wall budget exhausted')
    signal.signal(signal.SIGALRM, alarm)
    try:
        signal.setitimer(signal.ITIMER_REAL, seconds)
        response = function()
        result = {'execution_ok': True, 'response': response, 'timed_out': False}
    except (Exception, EvaluationTimeout) as error:
        result = {'execution_ok': False, 'timed_out': isinstance(error, EvaluationTimeout),
                  'error_type': type(error).__name__, 'reason': str(error)}
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, handler)
    return {**result, 'wall_seconds': time.perf_counter()-start,
            'cpu_seconds': time.process_time()-cpu, 'hard_timeout_enforced': True}


def extract_certificate(response):
    if type(response) is not dict:
        return None
    return response.get('certificate', response if 'format' in response else None)


def assess(certificate, item, verifier, normalize):
    """Original-task replay is separately orchestrated, never solver self-report."""
    task = item['task']
    invalid = {'certificate_valid': False, 'target_met': False}
    if type(certificate) is not dict:
        return {**invalid, 'reason': 'missing_certificate'}
    if certificate.get('function') != normalize(task['function']):
        return {**invalid, 'reason': 'original_function_mismatch'}
    scope = task.get('mean_zero') if task['kind'] == 'spectrum' else None
    if scope is not None and certificate.get('mean_zero') is not scope:
        return {**invalid, 'reason': 'original_space_mismatch'}
    if verifier(deepcopy(certificate), expected_function=deepcopy(task['function']),
                expected_mean_zero=scope, expected_tolerance=task['tolerance'], expected_kind=task['kind']) is not True:
        return {**invalid, 'reason': 'independent_replay_rejected'}
    try:
        value = F(certificate['upper'])-F(certificate['lower']) if task['kind']=='spectrum' else F(certificate['error_upper'])
        if value < 0:
            raise ValueError('negative_bound')
        if 'independent_exact_value' in item and not F(certificate['lower']) <= F(item['independent_exact_value']) <= F(certificate['upper']):
            raise ValueError('excludes_independent_exact_value')
        return {'certificate_valid': True, 'target_met': value <= F(task['tolerance']),
                'value': str(value), 'tolerance': task['tolerance'], 'reason': 'replayed',
                'quantity': 'spectral_interval_width' if task['kind']=='spectrum' else 'function_L2_error_upper'}
    except (ValueError, KeyError, TypeError, ArithmeticError) as error:
        return {**invalid, 'reason': str(error)}


def load_adapter(method, evidence_root, configuration=None):
    """Invoked inside the common setup timer, in a pristine interpreter."""
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(OLD))
    import backend
    from verification import verify_any
    if method in ('old_default', 'fixed_increasing', 'validation_selected_fixed', 'validation_candidate'):
        def adapter(task, index=None):
            request = {'op': 'approximate' if task['kind']=='approximation' else 'spectrum',
                       'function': deepcopy(task['function'])}
            if task['kind'] == 'spectrum':
                request.update(mean_zero=task['mean_zero'], tolerance=task['tolerance'])
            if index is not None:
                request.update(CONFIGS[index][task['kind']])
            return backend.legacy.handle(request)
    else:
        if method == 'frozen_v317':
            from research_service import ResearchService
        elif method == 'release':
            sys.path.insert(0, str(ROOT))
            from release_service import ResearchService
        else:
            raise ValueError('Unknown method: ' + method)
        service = ResearchService(evidence_root=evidence_root, cache=False)
        def adapter(task, index=None):
            return service.evaluate_task(deepcopy(task))
    return adapter, verify_any, backend.direct.normalize


def evaluate_row(method, item, adapter, verifier, normalize, configuration=None):
    attempts = []
    state = {'active_stage': 'before_solver'}
    def work():
        configurations = range(len(CONFIGS)) if method == 'fixed_increasing' else [configuration]
        for index in configurations:
            attempt = {'configuration': CONFIGS[index]['name'] if index is not None else 'method_default',
                       'assessment': {'certificate_valid': False, 'target_met': False, 'reason': 'incomplete'}}
            attempts.append(attempt)
            state['active_stage'] = 'solver'
            started, cpu = time.perf_counter(), time.process_time()
            try:
                response = adapter(deepcopy(item['task']), index)
                attempt.update(response=response, execution_ok=True)
            except Exception as error:
                attempt.update(execution_ok=False, error_type=type(error).__name__, reason=str(error))
            finally:
                attempt.update(solver_wall_seconds=time.perf_counter()-started,
                               solver_cpu_seconds=time.process_time()-cpu)
            state['active_stage'] = 'independent_replay'
            started, cpu = time.perf_counter(), time.process_time()
            try:
                attempt['assessment'] = assess(extract_certificate(attempt.get('response')), item, verifier, normalize)
            finally:
                attempt.update(replay_wall_seconds=time.perf_counter()-started,
                               replay_cpu_seconds=time.process_time()-cpu)
            state['active_stage'] = 'serialization_and_selection'
            canonical(attempt)  # full evidence must be serializable within deadline
            if attempt['assessment']['target_met']:
                break
        valid = [i for i,a in enumerate(attempts) if a['assessment']['certificate_valid']]
        best = min(valid, key=lambda i: F(attempts[i]['assessment']['value'])) if valid else None
        return {'selected_attempt':best}
    execution = timed_call(work, item['budget_seconds'])
    # Interrupted work retains its already replayed evidence, but cannot succeed.
    if execution['execution_ok']:
        best = execution['response']['selected_attempt']
    else:
        valid = [i for i,a in enumerate(attempts) if a['assessment']['certificate_valid']]
        best = min(valid, key=lambda i: F(attempts[i]['assessment']['value'])) if valid else None
    judgment = attempts[best]['assessment'] if best is not None else {'certificate_valid': False, 'target_met': False}
    within_budget = execution['execution_ok'] and execution['wall_seconds'] <= item['budget_seconds']
    return {'method': method, 'case_id': item['id'], 'task': item['task'], 'family': item['family'],
            'classification': item['classification'], 'budget_seconds': item['budget_seconds'],
            'attempts': attempts, 'selected_attempt': best, 'assessment': judgment,
            'execution': {k:v for k,v in execution.items() if k != 'response'},
            'last_stage': state['active_stage'], 'within_budget': within_budget,
            'success': judgment['target_met'] and within_budget,
            'accounted_wall_seconds': execution['wall_seconds'], 'accounted_cpu_seconds': execution['cpu_seconds']}


def mutation_probes(rows, cases, verifier, normalize):
    by_id, seen, probes = {c['id']: c for c in cases}, set(), []
    for row in rows:
        kind = row['task']['kind']
        if kind in seen or row['selected_attempt'] is None:
            continue
        seen.add(kind)
        cert = extract_certificate(row['attempts'][row['selected_attempt']]['response'])
        item = by_id[row['case_id']]
        variants = []
        changed = deepcopy(cert)
        changed['function']['offset'] = str(F(changed['function']['offset'])+1)
        variants.append(('wrong_original_function', changed))
        impossible = deepcopy(cert)
        if kind == 'spectrum':
            impossible['lower'] = str(F(impossible['upper'])+1)
            wrong_space = deepcopy(cert)
            wrong_space['mean_zero'] = not cert['mean_zero']
            variants.append(('wrong_original_space', wrong_space))
        else:
            impossible['error_upper'] = '-1'
        variants.append(('impossible_bound', impossible))
        forged = {k:deepcopy(cert[k]) for k in ('format','function','mean_zero','kind') if k in cert}
        forged.update(certificate_valid=True, target_met=True, verified=True, status='target_met')
        variants.append(('status_only_forgery', forged))
        for name, changed in variants:
            call = timed_call(lambda: assess(changed, item, verifier, normalize), 3.0)
            judgment = call.get('response', {})
            probes.append({'case_id': row['case_id'], 'kind': kind, 'mutation': name,
                           'accepted': judgment.get('certificate_valid', False), 'execution': call})
    return probes


def worker(spec):
    start, cpu = time.perf_counter(), time.process_time()
    with tempfile.TemporaryDirectory(prefix='geometry-evaluation-') as evidence:
        setup = timed_call(lambda: load_adapter(spec['method'], evidence, spec.get('configuration')), 30.0)
        rows, probes = [], []
        if setup['execution_ok']:
            adapter, verifier, normalize = setup.pop('response')
            for item in spec['cases']:
                row = evaluate_row(spec['method'], item, adapter, verifier, normalize, spec.get('configuration'))
                rows.append(row)
                # Checkpoint includes every attempted task, including failures/timeouts.
                write_json(spec['checkpoint'], {'rows': rows, 'setup': setup, 'complete': False})
            probes = mutation_probes(rows, spec['cases'], verifier, normalize)
        report = {'method': spec['method'], 'repetition': spec['repetition'],
                  'configuration': spec.get('configuration'), 'setup': setup, 'rows': rows,
                  'reliability_probes': probes, 'session_wall_seconds': time.perf_counter()-start,
                  'session_cpu_seconds': time.process_time()-cpu,
                  'max_rss_native_units': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                  'max_rss_unit': 'bytes' if sys.platform=='darwin' else 'KiB',
                  'complete': setup['execution_ok'] and len(rows)==len(spec['cases']),
                  'first_task_latency_seconds': setup['wall_seconds'] + (rows[0]['accounted_wall_seconds'] if rows else 0),
                  'first_task_case_id': rows[0]['case_id'] if rows else None,
                  'cache_condition': 'cold interpreter per session, session-local caches permitted',
                  'actual_model_usage': None, 'actual_model_identity': None, 'token_proxy_used': False}
        write_json(spec['output'], report)


def run_session(method, cases, repetition, destination, configuration=None):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    stem = method + ('_'+str(configuration) if configuration is not None else '') + '_r'+str(repetition)
    paths = {name: destination/(stem+'_'+name+'.json') for name in ('input','output','checkpoint')}
    shift = repetition % len(cases) if cases else 0
    ordered = cases[shift:] + cases[:shift]
    spec = {'method': method, 'cases': ordered, 'repetition': repetition, 'configuration': configuration,
            'output': str(paths['output']), 'checkpoint': str(paths['checkpoint'])}
    if paths['output'].exists():
        if not paths['input'].exists() or read_json(paths['input']) != spec:
            raise RuntimeError('Refusing to reuse a session result with a different specification')
        prior = read_json(paths['output'])
        if 'process' not in prior:
            raise RuntimeError('Incomplete parent accounting; preserve this run and start a new evaluation')
        return prior
    write_json(paths['input'], spec)
    start = time.perf_counter()
    child_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    limit = 60 + sum(c['budget_seconds'] for c in cases) + 30
    environment = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0')
    try:
        completed = subprocess.run([sys.executable, '-B', str(Path(__file__).resolve()), 'worker', '--spec', str(paths['input'])],
                                   cwd=str(ROOT), env=environment, capture_output=True, text=True, timeout=limit)
        process = {'returncode': completed.returncode, 'stdout': completed.stdout, 'stderr': completed.stderr,
                   'parent_watchdog_timeout': False}
    except subprocess.TimeoutExpired as error:
        process = {'returncode': None, 'stdout': str(error.stdout or ''), 'stderr': str(error.stderr or ''),
                   'parent_watchdog_timeout': True}
    process['wall_seconds'] = time.perf_counter()-start
    child_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    process['cpu_seconds'] = child_after.ru_utime+child_after.ru_stime-child_before.ru_utime-child_before.ru_stime
    result = read_json(paths['output']) if paths['output'].exists() else read_json(paths['checkpoint']) if paths['checkpoint'].exists() else {'rows': [], 'setup': {'execution_ok': False}, 'complete': False}
    result.update(method=method, repetition=repetition, configuration=configuration, process=process, expected_cases=[c['id'] for c in cases],
                  result_file=str(paths['output']))
    returned = {r['case_id'] for r in result['rows']}
    for item in cases:
        if item['id'] not in returned:
            result['rows'].append({'method': method, 'case_id': item['id'], 'task': item['task'],
                'family': item['family'], 'classification': item['classification'], 'attempts': [],
                'selected_attempt': None, 'assessment': {'certificate_valid': False, 'target_met': False},
                'success': False, 'within_budget': False, 'accounted_wall_seconds': 0.0,
                'accounted_cpu_seconds': 0.0, 'execution': {'execution_ok': False, 'reason': 'worker_failed_or_session_not_reached'}})
    write_json(paths['output'], result)
    return result


def distribution(values):
    return {'values': values, 'median': statistics.median(values), 'min': min(values), 'max': max(values)} if values else None


def summarize(sessions, cases, preparation=None):
    summaries = {}
    for method in dict.fromkeys(s['method'] for s in sessions):
        group = [s for s in sessions if s['method']==method]
        rows = [r for s in group for r in s['rows']]
        per_case = []
        for item in cases:
            subset = [r for r in rows if r['case_id']==item['id']]
            per_case.append({'case_id': item['id'], 'family': item['family'], 'classification': item['classification'],
                'mean_zero': item['task'].get('mean_zero'), 'kind': item['task']['kind'],
                'successes': sum(r['success'] for r in subset), 'runs': len(subset),
                'robust_success': len(subset)==REPETITIONS and all(r['success'] for r in subset),
                'valid_certificates': sum(r['assessment']['certificate_valid'] for r in subset),
                'wall_seconds': distribution([r['accounted_wall_seconds'] for r in subset]),
                'values': [r['assessment'].get('value') for r in subset]})
        charge_wall = preparation['preparation_wall_seconds'] if method=='validation_selected_fixed' and preparation else 0
        charge_cpu = preparation['preparation_cpu_seconds'] if method=='validation_selected_fixed' and preparation else 0
        family_breakdown = {}
        for family in dict.fromkeys(c['family'] for c in cases):
            subset = [c for c in per_case if c['family']==family]
            family_breakdown[family] = {'cases': len(subset), 'robust_successes': sum(c['robust_success'] for c in subset),
                                      'successful_runs': sum(c['successes'] for c in subset)}
        scope_breakdown = {}
        for scope, label in ((True,'mean_zero_spectrum'),(False,'full_space_spectrum'),(None,'function_approximation')):
            subset = [c for c in per_case if c['mean_zero'] is scope]
            scope_breakdown[label] = {'cases':len(subset),'robust_successes':sum(c['robust_success'] for c in subset),
                                      'successful_runs':sum(c['successes'] for c in subset)}
        probes = [p for s in group for p in s.get('reliability_probes', [])]
        summaries[method] = {
            'per_case': per_case, 'family_breakdown': family_breakdown, 'scope_breakdown':scope_breakdown,
            'robust_successes': sum(c['robust_success'] for c in per_case), 'cases': len(cases),
            'successful_runs': sum(r['success'] for r in rows), 'attempted_task_runs': len(rows),
            'full_session_wall_seconds': distribution([s['process']['wall_seconds'] for s in group]),
            'full_session_cpu_seconds': distribution([s['process']['cpu_seconds'] for s in group]),
            'online_wall_seconds': distribution([sum(r['accounted_wall_seconds'] for r in s['rows']) for s in group]),
            'setup_wall_seconds': distribution([s['setup'].get('wall_seconds',0) for s in group]),
            'first_task_latency_seconds': distribution([s.get('first_task_latency_seconds', s['process']['wall_seconds']) for s in group]),
            'preparation_wall_seconds_charged_once': charge_wall, 'preparation_cpu_seconds_charged_once': charge_cpu,
            'total_experiment_wall_seconds': charge_wall+sum(s['process']['wall_seconds'] for s in group),
            'total_experiment_cpu_seconds': charge_cpu+sum(s['process']['cpu_seconds'] for s in group),
            'probes': len(probes), 'false_accepts': sum(p['accepted'] for p in probes),
            'probe_execution_failures': sum(not p['execution']['execution_ok'] for p in probes),
            'timeouts': sum(r['execution'].get('timed_out',False) for r in rows),
            'worker_failures': sum(s['process']['returncode'] != 0 for s in group),
            'actual_model_usage': None, 'model_tokenizer': None,
        }
    comparisons = {}
    if 'release' in summaries:
        new = summaries['release']
        for method in METHODS[:-1]:
            if method not in summaries:
                continue
            old = summaries[method]
            pairs = list(zip(new['per_case'], old['per_case']))
            comparisons[method] = {'robust_success_gains': [a['case_id'] for a,b in pairs if a['robust_success'] and not b['robust_success']],
                'robust_success_losses': [a['case_id'] for a,b in pairs if b['robust_success'] and not a['robust_success']],
                'median_full_session_wall_ratio_release_over_baseline': new['full_session_wall_seconds']['median']/old['full_session_wall_seconds']['median'],
                'median_full_session_cpu_ratio_release_over_baseline': new['full_session_cpu_seconds']['median']/old['full_session_cpu_seconds']['median'],
                'total_experiment_wall_ratio_including_preparation':new['total_experiment_wall_seconds']/old['total_experiment_wall_seconds'],
                'total_experiment_cpu_ratio_including_preparation':new['total_experiment_cpu_seconds']/old['total_experiment_cpu_seconds'],
                'median_online_wall_ratio':new['online_wall_seconds']['median']/old['online_wall_seconds']['median'] if old['online_wall_seconds']['median'] else None,
                'warning': 'Ratio includes all successes and failures. Faster failure is not a runtime benefit. Repeated timings are descriptive on a shared host.'}
    return {'methods': summaries, 'release_against_every_baseline': comparisons}


def select_validation():
    initialize()
    output = DATA / 'validation_selection.json'
    if output.exists():
        value = read_json(output)
        if value['protocol_digest'] != digest(public_rules()):
            raise RuntimeError('Selection protocol mismatch')
        return value
    if (DATA/'freeze_manifest.json').exists():
        raise RuntimeError('Cannot select a configuration after freeze')
    sessions = []
    for repetition in range(REPETITIONS):
        for index in [(i+repetition)%len(CONFIGS) for i in range(len(CONFIGS))]:
            sessions.append(run_session('validation_candidate', validation_cases(), repetition,
                                        DATA/'selection_sessions', index))
    candidates = []
    for index, config in enumerate(CONFIGS):
        group = [s for s in sessions if s.get('configuration')==index]
        rows = [r for s in group for r in s['rows']]
        candidates.append({'index': index, 'name': config['name'], 'successes': sum(r['success'] for r in rows),
                           'valid': sum(r['assessment']['certificate_valid'] for r in rows),
                           'session_files': [s['result_file'] for s in group]})
    winner = max(candidates, key=lambda c:(c['successes'],c['valid'],-c['index']))
    result = {'protocol_digest': digest(public_rules()), 'selected_config_index': winner['index'],
              'selected_name': winner['name'], 'candidates': candidates, 'created_utc': now(),
              'preparation_wall_seconds': sum(s['process']['wall_seconds'] for s in sessions),
              'preparation_cpu_seconds': sum(s['process']['cpu_seconds'] for s in sessions),
              'selection_rule': public_rules()['selection'], 'all_failed_branches_retained': True}
    write_json(output, result)
    return result


def verify_old_preservation():
    manifest = read_json(ROOT/'PRESERVED_V317.json')
    checked, private_references = 0, []
    for relative, expected in manifest['files'].items():
        path = (OLD/relative).resolve()
        if not path.is_relative_to(OLD.resolve()):
            raise RuntimeError('Invalid old preservation path')
        if path.name == 'private_seed.json':
            private_references.append({'path': relative, 'preexisting_hash': expected, 'not_read': True})
            continue
        if file_digest(path) != expected:
            raise RuntimeError('Frozen old file changed: '+relative)
        checked += 1
    return {'checked': checked, 'old_private_seed_references_without_reading': private_references}


def write_freeze_manifest():
    initialize()
    if (DATA/'holdout_instances.json').exists() or (DATA/'holdout_release.json').exists():
        raise RuntimeError('Already released; cannot refreeze used holdout')
    selection = read_json(DATA/'validation_selection.json')
    if selection['protocol_digest'] != digest(public_rules()):
        raise RuntimeError('Validation selection must precede freeze')
    preserved = verify_old_preservation()
    files = sorted(p for p in ROOT.rglob('*.py') if '__pycache__' not in p.parts)
    files += [ROOT/'PRESERVED_V317.json', ROOT/'EVALUATION.md', DATA/'protocol.json', DATA/'validation_cases.json',
              DATA/'validation_selection.json', DATA/'seed_commitment.json']
    manifest = {'status': 'frozen', 'created_utc': now(), 'protocol_digest': digest(public_rules()),
                'files': {str(p.relative_to(ROOT)): file_digest(p) for p in files}, 'old_preservation': preserved}
    destination = DATA/'freeze_manifest.json'
    if destination.exists():
        raise RuntimeError('Freeze manifest already exists; refusing silent replacement')
    write_json(destination, manifest)
    return {'manifest_digest': digest(manifest), 'files': len(files), 'holdout_released': False}


def checked_manifest():
    manifest = read_json(DATA/'freeze_manifest.json')
    if manifest.get('status') != 'frozen' or manifest.get('protocol_digest') != digest(public_rules()):
        raise RuntimeError('Invalid freeze manifest')
    required = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*.py') if '__pycache__' not in p.parts}
    required |= {'PRESERVED_V317.json','EVALUATION.md','evaluation/protocol.json','evaluation/validation_cases.json',
                 'evaluation/validation_selection.json','evaluation/seed_commitment.json'}
    if not required <= set(manifest['files']):
        raise RuntimeError('Every new Python module, old preservation manifest, public protocol and baseline selection must be frozen')
    for relative, expected in manifest['files'].items():
        path = (ROOT/relative).resolve()
        if not path.is_relative_to(ROOT.resolve()) or file_digest(path) != expected:
            raise RuntimeError('Frozen release file changed: '+relative)
    verify_old_preservation()
    return manifest


def holdout_cases(release_authorization=None):
    """The only code that reads the NEW seed; explicitly gated before any read."""
    manifest = checked_manifest()
    path, release_path = DATA/'holdout_instances.json', DATA/'holdout_release.json'
    if path.exists():
        value, release = read_json(path), read_json(release_path)
        if release['freeze_digest'] != digest(manifest) or release['suite_digest'] != digest(value):
            raise RuntimeError('Released holdout does not match freeze')
        return value
    if release_authorization is None:
        raise RuntimeError('Root must explicitly authorize holdout release after freeze and validation selection')
    authorization = read_json(release_authorization)
    if authorization != {'action': 'release_new_holdout', 'freeze_digest': digest(manifest)}:
        raise RuntimeError('Release authorization must exactly bind this frozen manifest')
    seed = bytes.fromhex(read_json(DATA/'private_seed.json')['seed_hex'])
    if hashlib.sha256(seed).hexdigest() != read_json(DATA/'seed_commitment.json')['sha256']:
        raise RuntimeError('Private seed commitment mismatch')
    rules, counter = public_rules()['holdout'], 0
    def choose(values):
        nonlocal counter
        value = int.from_bytes(hashlib.sha256(seed+counter.to_bytes(8,'big')).digest(),'big')
        counter += 1
        return values[value % len(values)]
    def singular():
        return q('abs_power', choose(rules['negative_amplitudes']), exponent=choose(rules['negative_exponents']))
    items = [case('H01','spectrum',q(amplitude='0',offset=choose(rules['constant_offsets'])),'1/10000000000',False),
             case('H02','spectrum',q(amplitude='0',offset=choose(rules['constant_offsets'])),'1/10000000000'),
             case('H03','spectrum',q(amplitude=choose(rules['positive_step_amplitudes'])),'1/1000'),
             case('H04','approximation',q(amplitude=choose(rules['positive_step_amplitudes'])),'1/10000000000')]
    for tolerance in rules['tolerances']:
        for scope in (True,False):
            items.append(case('H'+str(len(items)+1).zfill(2),'spectrum',singular(),tolerance,scope,'negative_power_spectrum'))
    for tolerance in rules['tolerances']:
        items.append(case('H'+str(len(items)+1).zfill(2),'approximation',singular(),tolerance,family='negative_power_approximation'))
    items.append(case('H14','spectrum',q(amplitude=choose(rules['negative_amplitudes'])),'1/100000000',family='step_spectrum_hard'))
    for item in items:
        item['classification'] = 'new_within_family_parameters_only'
    write_json(path, items)
    write_json(release_path, {'freeze_digest': digest(manifest), 'suite_digest': digest(items), 'released_utc': now(),
                             'authorization_digest': digest(authorization), 'seed_commitment_verified': True,
                             'seed_disclosed': False, 'one_shot_interpretation': True})
    return items


def run_evaluation(phase):
    initialize()
    selection = read_json(DATA/'validation_selection.json')
    if selection['protocol_digest'] != digest(public_rules()):
        raise RuntimeError('Selection mismatch')
    cases = holdout_cases() if phase=='holdout' else validation_cases() if phase=='validation' else regression_cases()
    if not cases:
        raise RuntimeError('No cases')
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    destination = DATA/(phase+'_'+stamp)
    sessions = []
    for repetition in range(REPETITIONS):
        for method in METHODS[repetition:]+METHODS[:repetition]:
            index = selection['selected_config_index'] if method=='validation_selected_fixed' else None
            session = run_session(method,cases,repetition,destination,index)
            sessions.append(session)
            write_json(destination/'progress.json', {'completed_sessions':len(sessions),'total_sessions':len(METHODS)*REPETITIONS,
                                                   'latest_method':method,'latest_repetition':repetition})
    if phase=='holdout':
        checked_manifest()
    report = {'version': VERSION, 'phase': phase, 'created_utc': now(), 'protocol_digest': digest(public_rules()),
              'cases': cases, 'session_files': [s['result_file'] for s in sessions],
              'results': summarize(sessions,cases,selection),
              'environment': {'python': sys.version, 'executable':sys.executable,'platform':platform.platform(),
                  'logical_cpus':os.cpu_count(),'concurrent_other_activity':'unknown, shared host',
                  'external_model_calls_by_evaluator':0,'actual_model_input_output_reasoning_tokens':None,
                  'development_model_usage':None,'financial_cost':None},
              'limitations': ['Only known-family parameter transfer; no new mathematical families.',
                 'Five repetitions quantify machine timing variation, not five independent mathematical tasks.',
                 'Independent evaluation orchestration shares trusted frozen mathematical replay code with solver backends.',
                 'All failures/unsupported cases remain; faster failure is not a benefit.',
                 'Development/model/human costs unknown; preparation, setup, online, replay, exceptions and probes retained.',
                 'Any algorithm change after reveal requires a fresh protocol and untouched seed.']}
    output = destination/'report.json'
    write_json(output,report)
    write_json(DATA/('latest_'+phase+'.json'),{'report_path':str(output),'results':report['results']})
    return {'report_path':str(output),'results':report['results']}


def self_test():
    checks = []
    def check(name, predicate):
        if not predicate:
            raise AssertionError(name)
        checks.append({'name':name,'passed':True})
    def catches_exception():
        try:
            time.sleep(2)
        except Exception:
            return 'swallowed'
    timeout = timed_call(catches_exception,0.03)
    check('BaseException alarm not swallowed by except Exception',timeout['timed_out'] and timeout['wall_seconds']<0.5)
    thread_result = []
    def thread_probe():
        try:
            timed_call(lambda:True,0.1)
        except RuntimeError:
            thread_result.append(True)
    thread = threading.Thread(target=thread_probe)
    thread.start(); thread.join()
    check('Non-main-thread evaluation rejected',thread_result==[True])
    with tempfile.TemporaryDirectory(prefix='geometry-evaluator-selftest-') as tmp:
        adapter, verifier, normalize = load_adapter('frozen_v317',tmp)
        item = validation_cases()[0]
        result = timed_call(lambda:adapter(item['task']),3)
        cert = extract_certificate(result.get('response'))
        check('Real exact constant certificate independently accepted',assess(cert,item,verifier,normalize)['target_met'])
        changed_item = deepcopy(item)
        changed_item['task']['function']['offset']='1/7'
        check('Correct proof of wrong original task rejected',not assess(cert,changed_item,verifier,normalize)['certificate_valid'])
        changed_item = deepcopy(item)
        changed_item['task']['kind']='approximation'
        changed_item['task'].pop('mean_zero')
        check('Wrong conclusion kind rejected',not assess(cert,changed_item,verifier,normalize)['certificate_valid'])
        row = evaluate_row('frozen_v317',item,adapter,verifier,normalize)
        probes = mutation_probes([row],[item],verifier,normalize)
        check('Forged status, wrong original space/function, impossible bound rejected',len(probes)==4 and not any(p['accepted'] for p in probes))
        tiny = deepcopy(item); tiny['budget_seconds']=0.025
        def sleeping(task,index):
            try:
                time.sleep(1)
            except Exception:
                return {'verified':True,'target_met':True}
        row = evaluate_row('release',tiny,sleeping,verifier,normalize)
        check('Timed-out solver never counted successful; partial attempt retained',not row['success'] and row['execution']['timed_out'] and len(row['attempts'])==1)
        def slow_replay(*args,**kwargs):
            time.sleep(1)
            return True
        row = evaluate_row('release',tiny,lambda task,index:cert,slow_replay,normalize)
        check('Independent replay charged to same whole-task deadline',not row['success'] and row['execution']['timed_out'] and row['last_stage']=='independent_replay')
        check('Unknown format cannot use solver success flags',not assess({'format':'fake','function':item['task']['function'],'mean_zero':False,'verified':True,'target_met':True},item,verifier,normalize)['certificate_valid'])
    check('Old private seed cannot be read by helper',_old_seed_read_is_forbidden())
    with tempfile.TemporaryDirectory(prefix='geometry-release-gate-test-') as tmp:
        actual_data, actual_check = globals()['DATA'], globals()['checked_manifest']
        globals()['DATA'] = Path(tmp)
        globals()['checked_manifest'] = lambda:{'test_only':'freeze_stub'}
        try:
            rejected = False
            try:
                holdout_cases()
            except RuntimeError as error:
                rejected = 'explicitly authorize' in str(error)
            check('Missing explicit release authorization rejected before seed read',rejected)
            authorization = Path(tmp)/'wrong_authorization.json'
            write_json(authorization,{'action':'release_new_holdout','freeze_digest':'wrong'})
            rejected = False
            try:
                holdout_cases(authorization)
            except RuntimeError as error:
                rejected = 'exactly bind' in str(error)
            check('Wrong freeze authorization rejected before seed read',rejected)
            check('Gate tests created no holdout',not (Path(tmp)/'holdout_instances.json').exists())
        finally:
            globals()['DATA'], globals()['checked_manifest'] = actual_data, actual_check
    result = {'created_utc':now(),'checks':checks,'passed':len(checks),'timeout_probe':timeout,
              'holdout_generated_or_read':False,'old_private_seed_read':False}
    write_json(DATA/'self_test.json',result)
    return result


def _old_seed_read_is_forbidden():
    try:
        read_json(OLD/'evaluation/private_seed.json')
    except RuntimeError:
        return True
    return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['initialize','self-test','select-validation','run-validation','run-regression','freeze','release-holdout','run-holdout','worker'])
    parser.add_argument('--authorization',type=Path)
    parser.add_argument('--spec',type=Path)
    args = parser.parse_args()
    if args.command=='worker':
        worker(read_json(args.spec)); return
    if args.command=='initialize': result=initialize()
    elif args.command=='self-test': result=self_test()
    elif args.command=='select-validation': result=select_validation()
    elif args.command=='freeze': result=write_freeze_manifest()
    elif args.command=='release-holdout':
        cases=holdout_cases(args.authorization)
        result={'released_cases':len(cases),'suite_digest':digest(cases),'private_seed_disclosed':False}
    else: result=run_evaluation(args.command.removeprefix('run-'))
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
