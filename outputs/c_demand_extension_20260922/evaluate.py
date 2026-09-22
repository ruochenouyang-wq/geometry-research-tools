"""Public C2 extension evaluation: baseline/C1/C2, 24 tasks, 3 cold runs.
Every 10-second deadline covers import, solve, JSON, full replay and raw output.
No old evaluator is executed, no seed is read, and no model API is called.
"""
from copy import deepcopy
from datetime import datetime, timezone
from fractions import Fraction
import argparse
import hashlib
import importlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import selectors
import signal
import statistics
import subprocess
import sys
import time
import traceback
import zipfile

ROOT = Path(__file__).resolve().parent
ORIGINAL = ROOT.parent / 'geometry_parallel_v1_20260921'
PREVIOUS = ROOT.parent / 'parallel_cross_dimension_20260921'
SYSTEMS = ('baseline', 'C1', 'C2')
REPEATS = 3
SECONDS = 10.0
CASE_COUNT = 24
VERSION = 'c_demand_extension_public_v1'


def encode(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('utf-8')


def safe(path):
    path = Path(path)
    if any(word in path.name.lower() for word in ('seed', 'private')):
        raise ValueError('Seed/private files are outside this evaluator')
    return path


def read(path):
    return json.loads(safe(path).read_bytes())


def write(path, value):
    path = safe(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    with temporary.open('wb') as stream:
        stream.write(encode(value) + b'\n')
        stream.flush()
    temporary.replace(path)


def sha(path):
    return hashlib.sha256(safe(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(encode(value)).hexdigest()


def stamp():
    return datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')


def exact(value):
    if type(value) not in (int, str):
        raise ValueError('Require an exact integer or rational string')
    return Fraction(value)


def peak_bytes(usage):
    if sys.platform == 'darwin':
        return int(usage.ru_maxrss)
    if sys.platform.startswith('linux'):
        return int(usage.ru_maxrss) * 1024
    raise RuntimeError('Unknown peak RSS units')


def cpu(usage):
    return usage.ru_utime + usage.ru_stime


def meta(item):
    return {key: deepcopy(item[key]) for key in ('id', 'dimension', 'label', 'classification') if key in item}


def assert_files(files):
    for name, expected in files.items():
        if not Path(name).is_file() or sha(name) != expected:
            raise RuntimeError('Frozen file changed or missing: ' + name)


def assert_snapshot(path, expected):
    if sha(path) != expected:
        raise RuntimeError('Source snapshot archive changed')


def identity():
    pointer_path = PREVIOUS / 'evaluation/LATEST_FREEZE.json'
    pointer = read(pointer_path)
    if sha(pointer['manifest']) != pointer['sha256']:
        raise RuntimeError('Previous manifest changed')
    previous = read(pointer['manifest'])
    old_files = previous['identity']['files']
    if len(old_files) != 604:
        raise RuntimeError('Expected the preserved previous 604-file union')
    assert_files(old_files)
    assert_snapshot(previous['source_archive'], previous['source_archive_sha256'])
    # Explicit paths only: no geometry_* scanner and no expansion of old scope.
    new_files = {str(p.resolve()): sha(p) for p in sorted(ROOT.rglob('*.py'))}
    for name in ('cases.json', 'protocol.json'):
        p = ROOT / name
        new_files[str(p)] = sha(p)
    references = {str(pointer_path): sha(pointer_path), pointer['manifest']: sha(pointer['manifest'])}
    return {'files': {**old_files, **new_files, **references}, 'new_files': new_files,
            'preserved_previous_files': len(old_files), 'previous_freeze': pointer,
            'previous_archive': {'path': previous['source_archive'], 'sha256': previous['source_archive_sha256']}}


def inputs():
    cases, protocol = read(ROOT / 'cases.json'), read(ROOT / 'protocol.json')
    if type(cases) is not list or len(cases) != CASE_COUNT or len({c['id'] for c in cases}) != CASE_COUNT:
        raise ValueError('Require exactly 24 uniquely named fixed cases')
    expected = {'systems': list(SYSTEMS), 'repetitions': REPEATS, 'cold_seconds': 10, 'case_count': CASE_COUNT}
    for key, value in expected.items():
        if key in protocol and protocol[key] != value:
            raise ValueError('Protocol contradicts confirmed value: ' + key)
    for case in cases:
        task = case['task']
        if task['kind'] not in ('spectrum', 'approximation') or exact(task['tolerance']) <= 0:
            raise ValueError('Invalid task or tolerance')
        if task['kind'] == 'spectrum' and type(task.get('mean_zero')) is not bool:
            raise ValueError('Require exact spectrum space bool')
        if case.get('classification') not in ('new_public', 'old_failure_regression', 'expected_unsupported_control'):
            raise ValueError('Each case must identify its predeclared public/regression/control class')
    return cases, protocol


def freeze():
    inputs()
    for name in ('controller.py', 'runtime.py', 'evaluate.py', 'selftest_evaluate.py',
                 'fullspace_extension.py', 'step_extension.py'):
        if not (ROOT / name).is_file():
            raise RuntimeError('Source not ready: ' + name)
    before = identity()
    folder = ROOT / 'evaluation'
    folder.mkdir(exist_ok=True)
    label = stamp()
    archive = folder / ('snapshot_' + label + '.zip')
    with zipfile.ZipFile(archive, 'x', zipfile.ZIP_DEFLATED) as out:
        index = {}
        for number, (name, expected) in enumerate(sorted(before['files'].items())):
            raw = safe(name).read_bytes()
            if hashlib.sha256(raw).hexdigest() != expected:
                raise RuntimeError('Source changed during snapshot')
            member = 'files/' + str(number).zfill(6) + '.blob'
            out.writestr(member, raw)
            index[name] = {'member': member, 'sha256': expected}
        out.writestr('index.json', encode(index))
    if before != identity():
        raise RuntimeError('Identity changed while freezing')
    path = folder / ('freeze_' + label + '.json')
    write(path, {'version': VERSION, 'classification': 'public_development_not_holdout',
                 'identity': before, 'source_archive': str(archive), 'source_archive_sha256': sha(archive)})
    pointer = {'manifest': str(path), 'sha256': sha(path)}
    write(folder / 'LATEST_FREEZE.json', pointer)
    return pointer


def checked_freeze():
    pointer = read(ROOT / 'evaluation/LATEST_FREEZE.json')
    if sha(pointer['manifest']) != pointer['sha256']:
        raise RuntimeError('Freeze manifest changed')
    frozen = read(pointer['manifest'])
    if frozen['identity'] != identity():
        raise RuntimeError('Code, original dependencies or tasks changed')
    assert_snapshot(frozen['source_archive'], frozen['source_archive_sha256'])
    return pointer


def assess(certificate, task, verdict):
    rejected = {'certificate_valid': False, 'target_met': False, 'failure_class': 'invalid_certificate'}
    if type(certificate) is not dict:
        return {**rejected, 'reason': 'missing_certificate', 'failure_class': 'no_certificate'}
    if 'kind' in certificate and certificate['kind'] != task['kind']:
        return {**rejected, 'reason': 'original_kind_mismatch'}
    if certificate.get('function') != task['function']:
        return {**rejected, 'reason': 'original_function_mismatch'}
    if task['kind'] == 'spectrum' and certificate.get('mean_zero') is not task['mean_zero']:
        return {**rejected, 'reason': 'original_space_mismatch'}
    if 'tolerance' in certificate and exact(certificate['tolerance']) != exact(task['tolerance']):
        return {**rejected, 'reason': 'original_tolerance_mismatch'}
    if type(verdict) is not dict or verdict.get('certificate_valid') is not True:
        return {**rejected, 'reason': 'full_replay_rejected'}
    if task['kind'] == 'spectrum':
        metric = exact(certificate['upper']) - exact(certificate['lower'])
        quantity = 'lowest_spectral_interval_width'
    else:
        metric = exact(certificate['error_upper'])
        quantity = 'function_L2_error_upper'
        if 'error_squared' in certificate and not 0 <= exact(certificate['error_squared']) <= metric * metric:
            return {**rejected, 'reason': 'invalid_L2_enclosure'}
    if metric < 0:
        return {**rejected, 'reason': 'negative_metric'}
    reached = metric <= exact(task['tolerance'])
    return {'certificate_valid': True, 'target_met': reached, 'metric': str(metric), 'quantity': quantity,
            'metric_over_tolerance': str(metric / exact(task['tolerance'])), 'tolerance': task['tolerance'],
            'failure_class': None if reached else 'precision_insufficient'}


def cert_of(response):
    return response.get('certificate') if type(response) is dict else None


def provenance(system, response):
    response = response if type(response) is dict else {}
    counters = response.get('counters', {})
    counters = counters if type(counters) is dict else {}
    attempts = response.get('attempts', [])
    if system == 'baseline':
        return {'route': 'baseline_reference', 'solver_reported_route': response.get('route'),
                'fallback_reason': None, 'explanation_complete': True}
    if system == 'C1':
        fallback = bool(counters.get('baseline_fallback_calls') or response.get('baseline_fallback') or response.get('fallback'))
        reason = response.get('fallback_reason')
        inferred = False
        if fallback and not reason:
            failures = [a.get('reason', a.get('error')) for a in attempts if type(a) is dict and (a.get('reason') or a.get('error'))] if type(attempts) is list else []
            reason = '; '.join(str(x) for x in failures) or ('C1 uses its unchanged baseline route for this task' if response.get('route') == 'unchanged_baseline_route' else 'C1 returned through its recorded baseline_fallback_calls route')
            inferred = True
        return {'route': 'baseline_fallback' if fallback else 'native_c1',
                'solver_reported_route': response.get('route'), 'fallback_reason': reason,
                'fallback_reason_inferred_from_record': inferred, 'explanation_complete': not fallback or bool(reason)}
    route = response.get('route')
    recognized = route in ('native_c2', 'reused_c1', 'baseline_fallback')
    reason = response.get('fallback_reason')
    return {'route': route if recognized else 'unclassified', 'solver_reported_route': route,
            'fallback_reason': reason, 'route_recognized': recognized,
            'reuse_reason': response.get('native_kind') if route == 'reused_c1' else None,
            'explanation_complete': recognized and (route != 'baseline_fallback' or bool(reason))}


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class Synthetic:
    def __init__(self, behavior, marker=None):
        self.behavior = behavior
        self.marker = marker
    def solve(self, task, progress=None):
        if progress:
            progress({'stage': 'synthetic_search', 'event': 'attempt_started', 'reason': 'synthetic test'})
        if self.behavior == 'raise':
            raise RuntimeError('synthetic exception')
        child_metrics = None
        if self.behavior == 'child_busy':
            script = "import time,resource,json,os; start=time.process_time(); memory=bytearray(64*1024*1024)\nwhile time.process_time()-start<.12: pass\nu=resource.getrusage(resource.RUSAGE_SELF); print(json.dumps({'cpu':u.ru_utime+u.ru_stime,'maxrss':u.ru_maxrss,'pid':os.getpid()}))"
            completed = subprocess.run([sys.executable,'-B','-c',script],capture_output=True,text=True,check=True)
            child_metrics = json.loads(completed.stdout)
            if progress:
                progress({'stage':'synthetic_child','event':'child_reaped','child_metrics':child_metrics})
        if self.behavior == 'child_group_timeout':
            script = "import time,sys; from pathlib import Path; time.sleep(.6); Path(sys.argv[1]).write_text('descendant escaped process-group kill')"
            child = subprocess.Popen([sys.executable,'-B','-c',script,self.marker])
            if progress:
                progress({'stage':'synthetic_child','event':'child_started','child_pid':child.pid})
            child.wait()
        if self.behavior in ('swallow', 'memory_timeout'):
            if self.behavior == 'memory_timeout':
                self.memory = bytearray(32 * 1024 * 1024)
            def alarm(*_):
                raise BaseException('synthetic alarm')
            signal.signal(signal.SIGALRM, alarm)
            signal.setitimer(signal.ITIMER_REAL, .01, .01)
            while True:
                try:
                    time.sleep(5)
                except BaseException:
                    pass
        if self.behavior == 'unsupported':
            return {'status': 'unsupported_method_domain', 'certificate': None, 'route': 'native_c2', 'failure_stage': 'capability_check'}
        certificate = {'format': 'synthetic', 'function': task['function'], 'mean_zero': task.get('mean_zero'),
                       'lower': '0', 'upper': '1/1000000000', 'tolerance': task['tolerance']}
        return {'certificate': certificate, 'route': 'native_c2', 'synthetic_child_metrics':child_metrics,
                'attempts': [{'stage': 'synthetic_search', 'status': 'done'}]}
    def verify(self, certificate, task):
        if self.behavior == 'verify_timeout':
            time.sleep(30)
        if self.behavior == 'status_only':
            return {'status': 'target_met', 'target_met': True}
        return {'certificate_valid': True, 'target_met': True}


def load_engine(spec):
    sys.dont_write_bytecode = True
    if spec.get('synthetic'):
        module = Synthetic(spec['synthetic'],spec.get('synthetic_marker'))
        return module.solve, module.verify, True
    sys.path.insert(0, str(ROOT))
    if spec['system'] == 'C2':
        module = load_module(ROOT / 'controller.py', 'c2_evaluation_controller')
        runtime = importlib.import_module('runtime')
        if Path(runtime.__file__).resolve() != (ROOT / 'runtime.py').resolve():
            raise RuntimeError('Wrong runtime module loaded')
        return module.solve, runtime.verify, True
    path = ORIGINAL / ('shared_baseline.py' if spec['system'] == 'baseline' else 'c_adaptive/variant.py')
    module = load_module(path, 'c2_evaluation_frozen_' + spec['system'])
    return module.solve, module.verify, False


def worker(spec, deadline):
    event_stream = sys.stdout
    log = open(spec['log_path'], 'w', buffering=1)
    progress_stream = open(spec['progress_path'], 'wb', buffering=0)
    sys.stdout = log
    sys.stderr = log
    row = {'system': spec['system'], 'case': meta(spec['case']), 'task': spec['case']['task'],
           'task_sha256': digest(spec['case']['task']), 'repetition': spec.get('repetition'),
           'success': False, 'status': 'started', 'timings': {}, 'response_json_bytes': 0,
           'certificate_json_bytes': 0, 'failure_stage': 'import'}
    stage = 'import'
    def emit(event, **extra):
        usage = resource.getrusage(resource.RUSAGE_SELF)
        children = resource.getrusage(resource.RUSAGE_CHILDREN)
        message = {'event': event, 'stage': stage, 'monotonic': time.monotonic(),
                   'self_cpu_seconds': cpu(usage), 'reaped_children_cpu_seconds':cpu(children),
                   'self_plus_reaped_children_cpu_seconds':cpu(usage)+cpu(children),
                   'peak_rss_high_water_bytes': peak_bytes(usage),
                   'reaped_children_peak_rss_high_water_bytes':peak_bytes(children), **extra}
        event_stream.write(encode(message).decode('utf-8') + '\n')
        event_stream.flush()
    def progress(message):
        record = {'received_monotonic': time.monotonic(), 'payload': deepcopy(message)}
        progress_stream.write(encode(record) + b'\n')
        progress_stream.flush()
    def measured(name, function):
        nonlocal stage
        stage = name
        row['failure_stage'] = name
        emit('stage_started')
        start, before = time.monotonic(), resource.getrusage(resource.RUSAGE_SELF)
        children_before = resource.getrusage(resource.RUSAGE_CHILDREN)
        try:
            return function()
        finally:
            after = resource.getrusage(resource.RUSAGE_SELF)
            children_after = resource.getrusage(resource.RUSAGE_CHILDREN)
            self_cpu, child_cpu = cpu(after)-cpu(before), cpu(children_after)-cpu(children_before)
            row['timings'][name] = {'wall_seconds': time.monotonic()-start, 'cpu_seconds':self_cpu+child_cpu,
                'self_cpu_seconds':self_cpu,'reaped_children_cpu_seconds':child_cpu,
                'self_user_seconds':after.ru_utime-before.ru_utime,'self_system_seconds':after.ru_stime-before.ru_stime,
                'reaped_children_user_seconds':children_after.ru_utime-children_before.ru_utime,
                'reaped_children_system_seconds':children_after.ru_stime-children_before.ru_stime,
                'process_peak_rss_high_water_before_bytes': peak_bytes(before),
                'process_peak_rss_high_water_after_bytes': peak_bytes(after),
                'reaped_children_peak_rss_high_water_after_bytes':peak_bytes(children_after)}
            emit('stage_finished')
    write(spec['output_path'], row)
    try:
        solve, verify, accepts_progress = measured('import', lambda: load_engine(spec))
        solve_task = deepcopy(spec['case']['task'])
        solve_task['budget'] = {'wall_seconds': max(0.0, min(SECONDS, deadline-time.monotonic()))}
        row['solver_budget_seconds'] = solve_task['budget']['wall_seconds']
        response = measured('solve', lambda: solve(solve_task, progress=progress) if accepts_progress else solve(solve_task))
        def serialize():
            raw = encode(response)
            restored = json.loads(raw)
            certificate = cert_of(restored)
            proof = encode(certificate) if certificate is not None else b''
            for name, data in (('response.json', raw), ('certificate.json', proof)):
                if data:
                    with (Path(spec['output_path']).parent / name).open('wb') as out:
                        out.write(data)
                        out.flush()
            return restored, certificate, len(raw), len(proof)
        restored, certificate, size, proof_size = measured('serialization_and_payload_write', serialize)
        row.update(response=restored, response_json_bytes=size, certificate_json_bytes=proof_size,
                   provenance=provenance(spec['system'], restored))
        write(spec['output_path'], row)
        if certificate is None:
            unsupported = restored.get('status') in ('unsupported','unsupported_method_domain')
            row.update(status='unsupported' if unsupported else 'no_certificate',
                       failure_class='unsupported' if unsupported else 'no_certificate',
                       failure_stage=restored.get('failure_stage') or 'solve',
                       failure_reason=restored.get('reason') or restored.get('fallback_reason') or 'Solver returned no complete certificate')
        else:
            verdict = measured('full_verification', lambda: verify(deepcopy(certificate), deepcopy(spec['case']['task'])))
            row['verifier_result'] = json.loads(encode(verdict))
            assessment = measured('exact_target_check', lambda: assess(certificate, spec['case']['task'], row['verifier_result']))
            row.update(assessment=assessment, success=assessment['target_met'], failure_class=assessment['failure_class'])
            row['status'] = 'target_met' if row['success'] else assessment['failure_class']
            row['failure_stage'] = None if row['success'] else ('exact_target_check' if assessment['certificate_valid'] else 'full_verification')
            row['failure_reason'] = assessment.get('reason') if not row['success'] else None
            if assessment['failure_class'] == 'precision_insufficient':
                row['failure_reason'] = 'Verified bound exceeds the original tolerance'
    except Exception as error:
        row.update(success=False, status='exception', failure_class='execution_failure', failure_stage=stage,
                   failure_reason=str(error), error={'type': type(error).__name__, 'message': str(error), 'traceback': traceback.format_exc()})
    before_write = time.monotonic()
    own_usage, child_usage = resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)
    row['worker_usage_before_final_write'] = {'self_cpu_seconds':cpu(own_usage), 'reaped_children_cpu_seconds':cpu(child_usage),
        'self_plus_reaped_children_cpu_seconds':cpu(own_usage)+cpu(child_usage),
        'self_peak_rss_high_water_bytes':peak_bytes(own_usage),
        'reaped_children_peak_rss_high_water_bytes':peak_bytes(child_usage)}
    write(spec['output_path'], row)
    emit('output_written', output_write_wall_seconds=time.monotonic()-before_write)
    progress_stream.flush()
    log.flush()


def safe_output(path):
    try:
        value = read(path)
        if type(value) is not dict:
            raise ValueError('Expected object')
        return value
    except (OSError, ValueError, TypeError) as error:
        return {'success': False, 'failure_class': 'missing_or_invalid_output',
                'failure_reason': type(error).__name__+': '+str(error)}


def progress_records(path):
    records, errors = [], []
    if not Path(path).exists():
        return records, errors
    for line in Path(path).read_bytes().splitlines():
        try:
            records.append(json.loads(line))
        except (ValueError, TypeError) as error:
            errors.append(str(error))
    return records, errors


def run_process(spec, destination):
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    spec = {**deepcopy(spec), 'log_path': str(destination/'worker.log'), 'progress_path': str(destination/'progress.jsonl'),
            'output_path': str(destination/'output.json')}
    write(destination/'spec.json', spec)
    limit = spec.get('limit', SECONDS)
    events, errors, buffer = [], [], b''
    status = usage = watchdog = None
    start = time.monotonic()
    deadline = start+limit
    stderr = open(destination/'stderr.log', 'wb')
    child = subprocess.Popen([sys.executable, '-B', str(Path(__file__).resolve()), 'worker',
        '--spec', str(destination/'spec.json'), '--deadline', repr(deadline)], cwd=str(ROOT),
        stdout=subprocess.PIPE, stderr=stderr, start_new_session=True,
        env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0'))
    stderr.close()
    selector = selectors.DefaultSelector()
    selector.register(child.stdout, selectors.EVENT_READ)
    with (destination/'events.jsonl').open('wb') as raw_events:
        while True:
            if status is None and (time.monotonic()>=deadline or errors):
                watchdog = 'entire_process_deadline' if not errors else 'invalid_worker_event'
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                _, status, usage = os.wait4(child.pid, 0)
            wait = min(.01, max(0.0, deadline-time.monotonic())) if status is None else 0
            for key, _ in selector.select(wait):
                chunk = os.read(key.fileobj.fileno(), 65536)
                if not chunk:
                    selector.unregister(key.fileobj)
                    continue
                raw_events.write(chunk)
                raw_events.flush()
                buffer += chunk
                while b'\n' in buffer:
                    line, buffer = buffer.split(b'\n', 1)
                    try:
                        event = json.loads(line)
                        if type(event) is not dict or not start <= event['monotonic'] <= time.monotonic():
                            raise ValueError('Invalid timestamp')
                        if event['event'] not in ('stage_started', 'stage_finished', 'output_written'):
                            raise ValueError('Unknown event')
                        events.append(event)
                    except (ValueError, KeyError, TypeError) as error:
                        errors.append(str(error))
            if status is None:
                waited, wait_status, wait_usage = os.wait4(child.pid, os.WNOHANG)
                if waited:
                    status, usage = wait_status, wait_usage
            if status is not None and not selector.get_map():
                break
            if status is not None and time.monotonic()>deadline+.1:
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                errors.append('Worker exited but a descendant retained stdout')
                break
    if buffer:
        errors.append('Incomplete event line')
    selector.close()
    child.stdout.close()
    child.returncode = os.waitstatus_to_exitcode(status)
    elapsed = time.monotonic()-start
    process = {'wall_seconds': elapsed, 'cpu_seconds': cpu(usage), 'user_seconds': usage.ru_utime,
               'system_seconds': usage.ru_stime, 'peak_rss_bytes': peak_bytes(usage), 'ru_maxrss_raw': usage.ru_maxrss,
               'ru_maxrss_raw_unit': 'bytes' if sys.platform=='darwin' else 'KiB', 'platform': sys.platform,
               'returncode': child.returncode, 'watchdog': watchdog, 'protocol_errors': errors,
               'measurement': 'outer wait4 CPU already includes descendants reaped by worker; do not add children CPU again',
               'cpu_coverage': 'worker_and_reaped_descendants; unreaped descendants at forced termination may be missing',
               'forced_termination_cpu_may_be_lower_bound':watchdog is not None,
               'rss_scope':'largest observed worker/reaped-child single-process high water; not simultaneous process-tree RSS sum',
               'outer_wait4_peak_rss_bytes':peak_bytes(usage),
               'forced_termination_rss_descendant_coverage_may_be_incomplete':watchdog is not None}
    row = safe_output(spec['output_path'])
    child_peaks = [e.get('reaped_children_peak_rss_high_water_bytes',0) for e in events]
    child_peaks.append(row.get('worker_usage_before_final_write',{}).get('reaped_children_peak_rss_high_water_bytes',0))
    process['observed_reaped_children_peak_rss_bytes'] = max(child_peaks)
    process['peak_rss_bytes'] = max(process['outer_wait4_peak_rss_bytes'],process['observed_reaped_children_peak_rss_bytes'])
    records, progress_errors = progress_records(spec['progress_path'])
    process_ok = child.returncode==0 and watchdog is None and not errors and elapsed<=limit
    row.update(system=spec['system'], case=meta(spec['case']), task=spec['case']['task'], task_sha256=digest(spec['case']['task']),
               repetition=spec.get('repetition'), process=process, events=events, progress=records,
               progress_parse_errors=progress_errors, last_progress=records[-1] if records else None,
               within_budget=process_ok, success=bool(row.get('success') and process_ok), raw_directory=str(destination))
    last_payload = records[-1].get('payload') if records and type(records[-1]) is dict else None
    row['last_internal_action'] = last_payload.get('action',last_payload.get('stage')) if type(last_payload) is dict else None
    row['last_internal_reason'] = last_payload.get('reason',last_payload.get('reason_code')) if type(last_payload) is dict else None
    row.setdefault('response_json_bytes', 0)
    row.setdefault('certificate_json_bytes', 0)
    if not process_ok:
        row.update(failure_class='timeout' if watchdog=='entire_process_deadline' else 'process_failure',
                   failure_stage=events[-1]['stage'] if events else 'startup_or_import',
                   failure_reason='Parent killed process at fixed deadline' if watchdog=='entire_process_deadline' else 'Worker process failed or violated event protocol')
    row.setdefault('provenance', {'route':'unavailable_before_return', 'fallback_reason':None, 'explanation_complete':False})
    write(destination/'result.json', row)
    return row


def distribution(values):
    return {'raw':values, 'min':min(values), 'median':statistics.median(values), 'max':max(values)} if values else None


def summarize(rows, cases):
    result = {}
    for system in SYSTEMS:
        selected = [r for r in rows if r['system']==system]
        details = []
        for case in cases:
            runs = sorted([r for r in selected if r['case']['id']==case['id']], key=lambda r:r['repetition'])
            stable = len(runs)==REPEATS and {r['repetition'] for r in runs}==set(range(REPEATS)) and all(r['success'] for r in runs)
            details.append({**meta(case), 'kind':case['task']['kind'], 'fixed_repetitions':REPEATS, 'returned_rows':len(runs),
                'stable_3_of_3':stable, 'successful_runs':sum(r['success'] for r in runs),
                'routes':[r['provenance']['route'] for r in runs],
                'fallback_reasons':[r['provenance'].get('fallback_reason') for r in runs],
                'failure_classes':[r.get('failure_class') for r in runs],
                'failure_stages':[r.get('failure_stage') for r in runs],
                'failure_reasons':[r.get('failure_reason') for r in runs],
                'actual_metric_exact':[r.get('assessment',{}).get('metric') for r in runs],
                'metric_over_tolerance':[r.get('assessment',{}).get('metric_over_tolerance') for r in runs],
                'target_tolerance':case['task']['tolerance'],
                **{metric:distribution([r['process'][metric] for r in runs]) for metric in ('wall_seconds','cpu_seconds','peak_rss_bytes')},
                **{metric:distribution([r[metric] for r in runs]) for metric in ('response_json_bytes','certificate_json_bytes')},
                'full_verify_wall_seconds':distribution([r['timings']['full_verification']['wall_seconds'] for r in runs if 'full_verification' in r.get('timings',{})])})
        grouped = {}
        for field in ('classification','dimension','kind'):
            grouped['by_'+field] = {value:{'case_denominator':sum(r.get(field)==value for r in details),
                'run_denominator':REPEATS*sum(r.get(field)==value for r in details),
                'successful_runs':sum(r['successful_runs'] for r in details if r.get(field)==value),
                'stable_success_cases':sum(r.get(field)==value and r['stable_3_of_3'] for r in details)}
                for value in sorted({r.get(field,'unspecified') for r in details})}
        result[system] = {'case_denominator':len(cases), 'run_denominator':len(cases)*REPEATS,
            'returned_rows':len(selected), 'successful_runs':sum(r['success'] for r in selected),
            'stable_success_cases':sum(r['stable_3_of_3'] for r in details), 'per_case':details, **grouped,
            'total_wall_seconds_all_attempts':sum(r['process']['wall_seconds'] for r in selected),
            'total_cpu_seconds_all_attempts':sum(r['process']['cpu_seconds'] for r in selected),
            'total_cpu_all_attempts_may_be_lower_bound':any(r['process'].get('forced_termination_cpu_may_be_lower_bound') for r in selected),
            'cpu_lower_bound_possible_run_count':sum(bool(r['process'].get('forced_termination_cpu_may_be_lower_bound')) for r in selected),
            'route_counts':{route:sum(r['provenance']['route']==route for r in selected) for route in sorted({r['provenance']['route'] for r in selected})},
            'successful_route_counts':{route:sum(r['success'] and r['provenance']['route']==route for r in selected) for route in sorted({r['provenance']['route'] for r in selected})},
            'missing_route_explanations':sum(not r['provenance'].get('explanation_complete') for r in selected)}
        result[system]['expected_unsupported_control_outcomes']={
            'fixed_run_denominator':REPEATS*sum(c.get('classification')=='expected_unsupported_control' for c in cases),
            'reported_unsupported_within_deadline':sum(r['case'].get('classification')=='expected_unsupported_control' and r['within_budget'] and r.get('failure_class')=='unsupported' for r in selected),
            'math_successes':sum(r['case'].get('classification')=='expected_unsupported_control' and r['success'] for r in selected),
            'other_failures_or_timeouts':sum(r['case'].get('classification')=='expected_unsupported_control' and not r['success'] and not(r['within_budget'] and r.get('failure_class')=='unsupported') for r in selected),
            'note':'Expected method rejection is reported separately and never counted as a solved mathematical task.'}
        math_cases=[c for c in cases if c.get('classification')!='expected_unsupported_control']
        result[system]['mathematical_tasks_excluding_method_boundary_controls']={
            'fixed_case_denominator':len(math_cases),'fixed_run_denominator':len(math_cases)*REPEATS,
            'successful_runs':sum(r['case'].get('classification')!='expected_unsupported_control' and r['success'] for r in selected),
            'stable_success_cases':sum(r.get('classification')!='expected_unsupported_control' and r['stable_3_of_3'] for r in details)}
    def stable(system,case_id):
        return next(r for r in result[system]['per_case'] if r['id']==case_id)['stable_3_of_3']
    common = [c['id'] for c in cases if all(stable(s,c['id']) for s in SYSTEMS)]
    for system in SYSTEMS:
        result[system]['all_three_stable_common_subset'] = common
        result[system]['pairwise_comparisons'] = {}
        for reference in ('baseline','C1'):
            subset = [c['id'] for c in cases if stable(reference,c['id']) and stable(system,c['id'])]
            pairs = []
            for case_id in subset:
                ours = sorted([r for r in rows if r['system']==system and r['case']['id']==case_id],key=lambda r:r['repetition'])
                theirs = sorted([r for r in rows if r['system']==reference and r['case']['id']==case_id],key=lambda r:r['repetition'])
                case = next(c for c in cases if c['id']==case_id)
                ratios = {metric+'_ratio':[a['process'][metric]/b['process'][metric] if b['process'][metric] else None for a,b in zip(ours,theirs)] for metric in ('wall_seconds','cpu_seconds','peak_rss_bytes')}
                ratios.update({metric+'_ratio':[a[metric]/b[metric] if b[metric] else None for a,b in zip(ours,theirs)] for metric in ('response_json_bytes','certificate_json_bytes')})
                pairs.append({**meta(case), **ratios, 'routes':[r['provenance']['route'] for r in ours]})
            result[system]['pairwise_comparisons'][reference] = {'stable_same_success_subset':subset,'per_case':pairs,
                'warning':'Different subsets and new/regression strata are not interchangeable; inspect exact achieved precision and route.'}
    return result


def run():
    cases, protocol = inputs()
    pointer = checked_freeze()
    before = identity()
    destination = ROOT/'evaluation'/('final_'+stamp())
    destination.mkdir(parents=True)
    write(destination/'started.json', {'version':VERSION,'freeze':pointer,'identity_before':before,'cases_sha256':digest(cases)})
    rows = []
    for repetition in range(REPEATS):
        order = SYSTEMS[repetition%len(SYSTEMS):]+SYSTEMS[:repetition%len(SYSTEMS)]
        case_order = cases[repetition%len(cases):]+cases[:repetition%len(cases)]
        for case in case_order:
            for system in order:
                row = run_process({'system':system,'case':case,'repetition':repetition},destination/(system+'_'+case['id']+'_r'+str(repetition)))
                rows.append(row)
                entry={'completed_runs':len(rows),'fixed_runs':CASE_COUNT*len(SYSTEMS)*REPEATS,
                    'case_id':case['id'],'system':system,'repetition':repetition,'success':row['success'],
                    'route':row['provenance']['route'],'failure_class':row.get('failure_class'),
                    'wall_seconds':row['process']['wall_seconds'],'cpu_seconds':row['process']['cpu_seconds'],
                    'peak_rss_bytes':row['process']['peak_rss_bytes'],'result_path':str(Path(row['raw_directory'])/'result.json')}
                with (destination/'completed_index.jsonl').open('ab') as index:
                    index.write(encode(entry)+b'\n')
                    index.flush()
                write(destination/'checkpoint.json', {'complete':False,**entry})
                print('['+str(len(rows))+'/'+str(CASE_COUNT*len(SYSTEMS)*REPEATS)+'] '+system+' '+case['id']+
                      ' r='+str(repetition)+' '+('success' if row['success'] else str(row.get('failure_class')))+
                      ' route='+row['provenance']['route']+' wall='+format(row['process']['wall_seconds'],'.3f')+'s',flush=True)
    try:
        after = identity()
        unchanged = before==after and checked_freeze()==pointer
        error = None
    except (OSError,ValueError,RuntimeError) as exception:
        after=None;unchanged=False;error=str(exception)
    unique = {(r['system'],r['case']['id'],r['repetition']) for r in rows}
    complete = len(rows)==CASE_COUNT*len(SYSTEMS)*REPEATS and len(unique)==len(rows)
    report = {'version':VERSION,'classification':'public_development_not_holdout','complete_fixed_denominator':complete,
        'valid_frozen_comparison':complete and unchanged,'freeze':pointer,'identity_before':before,'identity_after':after,
        'identity_error':error,'rows':rows,'summary':summarize(rows,cases),
        'execution_environment':{'python_executable':sys.executable,'python_version':sys.version.split()[0],'platform':sys.platform},
        'actual_model_api_calls':0,'actual_model_tokens':None,'actual_cached_model_tokens':None,
        'limitations':['Three serial interleaved repeats on a shared machine; timing is preliminary.',
            'JSON bytes are serialized volume, not token counts.',
            'New public cases and disclosed old failure regressions must be reported separately.',
            'C2 native, reused C1 and true baseline fallback are distinct routes.',
            'Independent certificate replay shares frozen mathematical dependencies; not a formal proof assistant.']}
    write(destination/'report.json',report)
    write(destination/'checkpoint.json',{'complete':complete,'completed_runs':len(rows),
        'fixed_runs':CASE_COUNT*len(SYSTEMS)*REPEATS,'valid_frozen_comparison':complete and unchanged,
        'report_path':str(destination/'report.json')})
    if not complete or not unchanged:
        raise RuntimeError('Results retained, but freeze identity or denominator is invalid')
    return {'report_path':str(destination/'report.json'),'rows':len(rows),'valid_frozen_comparison':True}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command',choices=['worker','freeze','run','check'])
    parser.add_argument('--spec');parser.add_argument('--deadline',type=float)
    args=parser.parse_args()
    if args.command=='worker':
        worker(read(args.spec),args.deadline);return
    if args.command=='freeze':result=freeze()
    elif args.command=='run':result=run()
    else:
        cases,_=inputs();current=identity()
        result={'cases':len(cases),'preserved_previous_files':current['preserved_previous_files'],'snapshot_union_files':len(current['files']),'all_unchanged':True}
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
