"""Frozen public cross-dimension comparison; no network/model or seed access.

Four systems, five serial interleaved repetitions. Each process has a parent
absolute watchdog, and its own wait4 CPU/peak-RSS accounting, including failures.
"""
from copy import deepcopy
from datetime import datetime, timezone
from fractions import Fraction
import argparse
import hashlib
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
WORKSPACE = ROOT.parents[1]
OLD = ROOT.parent / 'geometry_parallel_v1_20260921'
SYSTEMS = ('baseline', 'A', 'B', 'C')
REPEATS = 5
COLD = 10.0
PREPARE = 30.0
QUERY = 10.0
SESSION = 200.0
REPLAY = 10.0
VERSION = 'cross_dimension_public_development_v1'


def encoded(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False).encode('utf-8')


def allowed(path):
    path = Path(path)
    if 'seed' in path.name.lower() or 'private' in path.name.lower():
        raise ValueError('Seed/private files are outside this evaluator')
    return path


def read(path):
    return json.loads(allowed(path).read_bytes())


def write(path, value):
    path = allowed(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix + '.tmp')
    with temp.open('wb') as stream:
        stream.write(encoded(value) + b'\n')
        stream.flush()
    temp.replace(path)


def sha(path):
    return hashlib.sha256(allowed(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(encoded(value)).hexdigest()


def stamp():
    return datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')


def exact(value):
    if type(value) not in (int, str):
        raise ValueError('Bounds must be exact integers or rational strings')
    return Fraction(value)


def rss_bytes(usage):
    if sys.platform == 'darwin':
        return int(usage.ru_maxrss)
    if sys.platform.startswith('linux'):
        return int(usage.ru_maxrss) * 1024
    raise RuntimeError('Unknown ru_maxrss units on this platform')


def cpu(usage):
    return usage.ru_utime + usage.ru_stime


def metadata(item):
    return {k: deepcopy(item[k]) for k in ('id', 'dimension', 'label') if k in item}


def identity():
    previous_pointer = OLD / 'evaluation/LATEST_FREEZE.json'
    pointer = read(previous_pointer)
    previous_manifest = Path(pointer['manifest'])
    if sha(previous_manifest) != pointer['sha256']:
        raise RuntimeError('Previous frozen manifest identity changed')
    frozen = read(previous_manifest)
    old_files = {**frozen['identity']['sources'], **frozen['identity']['artifacts']}
    for path, expected in old_files.items():
        if not Path(path).is_file() or sha(path) != expected:
            raise RuntimeError('Previous frozen source/artifact changed: ' + path)
    if sha(frozen['source_archive']) != frozen['source_archive_sha256']:
        raise RuntimeError('Previous frozen source archive changed')
    baseline_path = OLD / 'BASELINE.json'
    baseline = read(baseline_path)['files']
    baseline_files = {}
    for name, expected in baseline.items():
        path = (WORKSPACE / name).resolve()
        if not path.is_file() or sha(path) != expected:
            raise RuntimeError('Original baseline changed: ' + name)
        baseline_files[str(path)] = expected
    new_files = {str(p.resolve()): sha(p) for p in sorted(ROOT.rglob('*.py'))}
    for name in ('cases.json', 'families.json', 'protocol.json', 'environment.json', 'notes/ANALYTIC_CONTROLS.md'):
        path = ROOT / name
        new_files[str(path)] = sha(path)
    preserved_metadata = {str(previous_pointer): sha(previous_pointer),
                          str(previous_manifest): sha(previous_manifest),
                          str(baseline_path): sha(baseline_path)}
    return {'files': {**baseline_files, **old_files, **new_files, **preserved_metadata},
            'original_baseline_files': len(baseline),
            'previous_sources': len(frozen['identity']['sources']),
            'previous_artifacts': len(frozen['identity']['artifacts']),
            'previous_manifest': pointer, 'new_files': new_files,
            'previous_archive': {'path': frozen['source_archive'], 'sha256': frozen['source_archive_sha256']}}


def inputs():
    cases, families, protocol = (read(ROOT / name) for name in ('cases.json', 'families.json', 'protocol.json'))
    if type(cases) is not list or len(cases) != 18 or len({c['id'] for c in cases}) != 18:
        raise ValueError('Require exactly 18 uniquely named common cases')
    if sum(c['task']['kind'] == 'spectrum' for c in cases) != 16 or sum(c['task']['kind'] == 'approximation' for c in cases) != 2:
        raise ValueError('Require 16 spectrum and two approximation cases')
    if type(families) is not list or len(families) != 2 or len({f['id'] for f in families}) != 2:
        raise ValueError('Require exactly two parameter families')
    for family in families:
        if len(family['tasks']) != 16 or len({t['id'] for t in family['tasks']}) != 16:
            raise ValueError('Require exactly 16 uniquely named queries per family')
    expected = {'systems': list(SYSTEMS), 'repetitions': 5, 'cold_seconds': 10,
                'prepare_seconds': 30, 'query_seconds': 10, 'session_seconds': 200, 'replay_seconds': 10}
    for key, value in expected.items():
        if key in protocol and protocol[key] != value:
            raise ValueError('Protocol contradicts confirmed evaluation constant: ' + key)
    for item in cases + [q for f in families for q in f['tasks']]:
        task = item['task']
        if exact(task['tolerance']) <= 0:
            raise ValueError('Tolerance must be positive')
        if task['kind'] == 'spectrum' and type(task.get('mean_zero')) is not bool:
            raise ValueError('Spectrum space must be an exact bool')
    return cases, families, protocol


def resolve_family(item, families):
    value = item.get('family')
    if value is None or type(value) is dict:
        return value
    matches = [f['spec'] for f in families if f['id'] == value]
    if len(matches) != 1:
        raise ValueError('Unknown family id: ' + str(value))
    return matches[0]


def freeze():
    inputs()
    for name in ('engines.py', 'a_adapter.py', 'b_adapter.py', 'benchmark.py', 'selftest.py'):
        if not (ROOT / name).is_file():
            raise RuntimeError('Entrypoint not ready: ' + name)
    data = identity()
    folder = ROOT / 'evaluation'
    folder.mkdir(exist_ok=True)
    label = stamp()
    archive = folder / ('source_snapshot_' + label + '.zip')
    with zipfile.ZipFile(archive, 'x', zipfile.ZIP_DEFLATED) as out:
        index = {}
        for i, (name, expected) in enumerate(sorted(data['files'].items())):
            raw = allowed(name).read_bytes()
            if hashlib.sha256(raw).hexdigest() != expected:
                raise RuntimeError('Source changed during snapshot')
            member = 'files/' + str(i).zfill(6) + '.blob'
            out.writestr(member, raw)
            index[name] = {'member': member, 'sha256': expected}
        out.writestr('index.json', encoded(index))
    manifest = {'version': VERSION, 'created_utc': label, 'identity': data,
                'classification': 'frozen_public_development_not_holdout',
                'source_archive': str(archive), 'source_archive_sha256': sha(archive)}
    path = folder / ('freeze_' + label + '.json')
    write(path, manifest)
    pointer = {'manifest': str(path), 'sha256': sha(path)}
    write(folder / 'LATEST_FREEZE.json', pointer)
    return pointer


def checked_freeze():
    pointer = read(ROOT / 'evaluation/LATEST_FREEZE.json')
    if sha(pointer['manifest']) != pointer['sha256']:
        raise RuntimeError('Freeze manifest changed')
    manifest = read(pointer['manifest'])
    if identity() != manifest['identity']:
        raise RuntimeError('Code, tasks or frozen dependencies changed')
    if sha(manifest['source_archive']) != manifest['source_archive_sha256']:
        raise RuntimeError('Source snapshot changed')
    return pointer


def assess(certificate, task, verdict):
    invalid = {'certificate_valid': False, 'target_met': False, 'failure_class': 'invalid_certificate'}
    if type(certificate) is not dict:
        return {**invalid, 'reason': 'missing_certificate', 'failure_class': 'capability_failure'}
    if 'kind' in certificate and certificate['kind'] != task['kind']:
        return {**invalid, 'reason': 'original_kind_mismatch'}
    if certificate.get('function') != task['function']:
        return {**invalid, 'reason': 'original_function_mismatch'}
    if task['kind'] == 'spectrum' and certificate.get('mean_zero') is not task['mean_zero']:
        return {**invalid, 'reason': 'original_space_mismatch'}
    if 'tolerance' in certificate and exact(certificate['tolerance']) != exact(task['tolerance']):
        return {**invalid, 'reason': 'original_tolerance_mismatch'}
    if type(verdict) is not dict or verdict.get('certificate_valid') is not True:
        return {**invalid, 'reason': 'complete_verifier_rejected'}
    if task['kind'] == 'spectrum':
        metric = exact(certificate['upper']) - exact(certificate['lower'])
        quantity = 'lowest_spectral_interval_width'
    elif task['kind'] == 'approximation':
        metric = exact(certificate['error_upper'])
        if 'error_squared' in certificate and not 0 <= exact(certificate['error_squared']) <= metric * metric:
            return {**invalid, 'reason': 'invalid_L2_enclosure'}
        quantity = 'function_L2_error_upper'
    else:
        return {**invalid, 'reason': 'unsupported_task_kind', 'failure_class': 'unsupported'}
    if metric < 0:
        return {**invalid, 'reason': 'negative_metric'}
    met = metric <= exact(task['tolerance'])
    return {'certificate_valid': True, 'target_met': met, 'metric': str(metric),
            'tolerance': task['tolerance'], 'quantity': quantity,
            'failure_class': None if met else 'precision_insufficient'}


def certificate_of(response):
    return response.get('certificate') if type(response) is dict else None


def fallback_of(response):
    if type(response) is not dict:
        return {'fallback': False, 'route': None}
    counters = response.get('counters', {})
    return {'fallback': bool(response.get('fallback') or response.get('baseline_fallback') or
                             (type(counters) is dict and counters.get('baseline_fallback_calls'))),
            'route': response.get('route'), 'reason': response.get('fallback_reason')}


def preparation_payload_sizes(response):
    found = {}
    def scan(obj):
        if type(obj) is dict:
            for key, value in obj.items():
                if key in ('bank', 'certificate', 'family_certificate') and type(value) is dict:
                    found[digest(value)] = len(encoded(value))
                else:
                    scan(value)
        elif type(obj) is list:
            for value in obj:
                scan(value)
    scan(response)
    return {'preparation_certificate_payload_count': len(found),
            'preparation_certificate_json_bytes': sum(found.values()),
            'serialized_library_json_bytes': sum(found.values()),
            'library_bytes_definition': 'sum of unique complete bank/certificate payloads; not heap size',
            'preparation_response_json_bytes': len(encoded(response))}


class SyntheticEngine:
    """Only used by selftest; no mathematical success claims."""
    def __init__(self, behavior):
        self.behavior = behavior
    def prepare(self, seconds):
        return {'ok': self.behavior != 'unsupported', 'status': 'unsupported' if self.behavior == 'unsupported' else 'ready'}
    def solve(self, task):
        if self.behavior == 'raise':
            raise RuntimeError('synthetic solve exception')
        if self.behavior in ('sleep', 'swallow', 'memory_sleep'):
            if self.behavior == 'memory_sleep':
                self.retained_memory = bytearray(32 * 1024 * 1024)
            if self.behavior == 'swallow':
                def alarm(*_):
                    raise BaseException('synthetic alarm')
                signal.signal(signal.SIGALRM, alarm)
                signal.setitimer(signal.ITIMER_REAL, .01, .01)
            while True:
                try:
                    time.sleep(5)
                except BaseException:
                    pass
        cert = {'format': 'synthetic_only', 'function': task['function'], 'mean_zero': task.get('mean_zero'),
                'lower': '0', 'upper': '1/1000000000', 'tolerance': task['tolerance']}
        return {'certificate': cert, 'status': 'synthetic_only', 'route': 'synthetic'}
    def verify(self, cert, task, full=True):
        if self.behavior == 'verify_sleep':
            time.sleep(30)
        if self.behavior == 'status_only':
            return {'status': 'target_met', 'target_met': True}
        return {'certificate_valid': True, 'target_met': True, 'synthetic_only': True}


def load_engine(spec):
    if 'synthetic' in spec:
        return SyntheticEngine(spec['synthetic'])
    sys.dont_write_bytecode = True
    sys.path.insert(0, str(ROOT))
    module_spec = importlib.util.spec_from_file_location('cross_dimension_engines', ROOT / 'engines.py')
    module = importlib.util.module_from_spec(module_spec)
    sys.modules[module_spec.name] = module
    module_spec.loader.exec_module(module)
    return module.make_engine(spec['system'], family_spec=spec.get('family'))


def worker(spec, absolute_deadline):
    event_stream = sys.stdout
    log = open(spec['log_path'], 'w', buffering=1)
    sys.stdout = log
    sys.stderr = log
    def event(kind, name, **extra):
        usage = resource.getrusage(resource.RUSAGE_SELF)
        data = {'event': kind, 'name': name, 'monotonic': time.monotonic(),
                'cpu_seconds': cpu(usage), 'peak_rss_bytes': rss_bytes(usage), **extra}
        event_stream.write(encoded(data).decode('utf-8') + '\n')
        event_stream.flush()
    def timed(function):
        start, before = time.monotonic(), resource.getrusage(resource.RUSAGE_SELF)
        try:
            value = function()
            return value, {'wall_seconds': time.monotonic() - start,
                           'cpu_seconds': cpu(resource.getrusage(resource.RUSAGE_SELF)) - cpu(before),
                           'process_peak_rss_high_water_before_bytes': rss_bytes(before),
                           'process_peak_rss_high_water_after_bytes': rss_bytes(resource.getrusage(resource.RUSAGE_SELF))}
        except Exception as error:
            error._benchmark_timing = {'wall_seconds': time.monotonic() - start,
                                      'cpu_seconds': cpu(resource.getrusage(resource.RUSAGE_SELF)) - cpu(before)}
            raise
    started = time.monotonic()
    engine, import_cost = timed(lambda: load_engine(spec))
    event('ready', 'initialization', import_cost=import_cost)
    def prepare():
        seconds = max(0.0, min(PREPARE, absolute_deadline - time.monotonic()))
        response, elapsed = timed(lambda: engine.prepare(seconds))
        raw = encoded(response)
        response = json.loads(raw)
        row = {'ok': type(response) is dict and response.get('ok') is True, 'response': response,
               'timing': elapsed, **preparation_payload_sizes(response)}
        write(spec['preparation_path'], row)
        return row
    def operation(item, directory, full, replay_certificate=None, query_deadline=None):
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        row = {'system': spec['system'], 'case': metadata(item), 'task': item['task'],
               'task_sha256': digest(item['task']), 'success': False, 'timings': {},
               'verification_full': full, 'status': 'started',
               'response_json_bytes': 0, 'certificate_json_bytes': 0}
        write(directory / 'output.json', row)
        operation_start = time.monotonic()
        try:
            if spec['mode'] == 'replay':
                response = {'certificate': replay_certificate}
            else:
                solve_task = deepcopy(item['task'])
                solve_deadline = min(absolute_deadline, query_deadline) if query_deadline is not None else absolute_deadline
                solve_task['budget'] = {'wall_seconds': max(0.0, min(QUERY, solve_deadline - time.monotonic()))}
                row['solver_budget_seconds'] = solve_task['budget']['wall_seconds']
                response, row['timings']['solve'] = timed(lambda: engine.solve(solve_task))
            def serialize():
                response_bytes = encoded(response)
                restored = json.loads(response_bytes)
                certificate = certificate_of(restored)
                cert_bytes = encoded(certificate) if certificate is not None else b''
                for name, content in (('response.json', response_bytes), ('certificate.json', cert_bytes)):
                    if content:
                        with (directory / name).open('wb') as stream:
                            stream.write(content)
                            stream.flush()
                return restored, certificate, len(response_bytes), len(cert_bytes)
            payload, row['timings']['serialization_and_payload_write'] = timed(serialize)
            restored, certificate, response_size, certificate_size = payload
            row.update(response=restored, response_json_bytes=response_size, certificate_json_bytes=certificate_size,
                       fallback=fallback_of(restored))
            # Preserve returned evidence even if the complete verifier times out.
            write(directory / 'output.json', row)
            if certificate is None:
                row.update(status='unsupported' if restored.get('status') == 'unsupported' or restored.get('unsupported') else 'no_certificate',
                           failure_class='unsupported' if restored.get('status') == 'unsupported' or restored.get('unsupported') else 'capability_failure')
            else:
                verdict, row['timings']['verification'] = timed(lambda: engine.verify(deepcopy(certificate), deepcopy(item['task']), full=full))
                row['verifier_result'] = json.loads(encoded(verdict))
                row['assessment'] = assess(certificate, item['task'], row['verifier_result'])
                row['success'] = row['assessment']['target_met']
                row['failure_class'] = row['assessment']['failure_class']
                row['status'] = 'target_met' if row['success'] else row['failure_class']
        except Exception as error:
            row.update(success=False, status='exception', failure_class='execution_failure',
                       error={'type': type(error).__name__, 'message': str(error), 'traceback': traceback.format_exc(),
                              'timing_until_exception': getattr(error, '_benchmark_timing', None)})
        row['operation_wall_seconds_before_final_write'] = time.monotonic() - operation_start
        usage = resource.getrusage(resource.RUSAGE_SELF)
        row['process_peak_rss_high_water_bytes'] = rss_bytes(usage)
        before_write = time.monotonic()
        write(directory / 'output.json', row)
        event('output_written', item['id'], final_output_write_wall_seconds=time.monotonic() - before_write)
        return row
    if spec['mode'] == 'replay':
        operation(spec['case'], spec['output_directory'], True, spec['certificate'])
        event('finished', 'replay')
        return
    event('phase_started', 'preparation')
    try:
        preparation = prepare()
    except Exception as error:
        preparation = {'ok': False, 'error': {'type': type(error).__name__, 'message': str(error),
                                             'traceback': traceback.format_exc(),
                                             'timing_until_exception': getattr(error, '_benchmark_timing', None)}}
        write(spec['preparation_path'], preparation)
    event('phase_finished', 'preparation', ok=preparation['ok'])
    if not preparation['ok']:
        return
    if spec['mode'] == 'cold':
        operation(spec['case'], spec['output_directory'], True)
        event('finished', 'cold')
    elif spec['mode'] == 'batch':
        for item in spec['cases']:
            query_deadline = min(absolute_deadline, time.monotonic() + QUERY)
            event('phase_started', item['id'])
            row = operation(item, Path(spec['queries_directory']) / item['id'], False, query_deadline=query_deadline)
            event('phase_finished', item['id'], ok=True, target_met=row['success'])
    else:
        raise ValueError('Unknown worker mode')
    log.flush()


def process_metrics(usage, elapsed, status, watchdog, errors, start):
    return {'wall_seconds': elapsed, 'cpu_seconds': cpu(usage),
            'user_seconds': usage.ru_utime, 'system_seconds': usage.ru_stime,
            'peak_rss_bytes': rss_bytes(usage), 'ru_maxrss_raw': usage.ru_maxrss,
            'ru_maxrss_raw_unit': 'bytes' if sys.platform == 'darwin' else 'KiB',
            'platform': sys.platform, 'returncode': os.waitstatus_to_exitcode(status),
            'watchdog': watchdog, 'protocol_errors': errors, 'parent_start_monotonic': start,
            'measurement': 'per-process os.wait4; peak is not a children-max difference'}


def safe_output(path):
    try:
        value = read(path)
        if type(value) is not dict:
            raise ValueError('Output must be an object')
        return value
    except (OSError, ValueError, TypeError) as error:
        return {'success': False, 'status': 'missing_or_invalid_output', 'failure_class': 'execution_failure',
                'output_read_error': type(error).__name__ + ': ' + str(error)}


def run_process(spec, destination):
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    spec = {**deepcopy(spec), 'log_path': str(destination / 'worker.log'),
            'preparation_path': str(destination / 'preparation.json'),
            'output_directory': str(destination), 'queries_directory': str(destination / 'queries')}
    write(destination / 'spec.json', spec)
    mode = spec['mode']
    limit = spec.get('limit', SESSION if mode == 'batch' else REPLAY if mode == 'replay' else COLD)
    start = time.monotonic()
    global_deadline = start + limit
    phase_deadline = min(global_deadline, start + QUERY) if mode == 'batch' else global_deadline
    active_phase = 'initialization' if mode == 'batch' else 'entire_process'
    sequence = ['preparation'] + [item['id'] for item in spec.get('cases', [])]
    next_phase = 0
    events, errors, buffer = [], [], b''
    watchdog = None
    usage = status = None
    stderr_stream = open(destination / 'stderr.log', 'wb')
    child = subprocess.Popen([sys.executable, '-B', str(Path(__file__).resolve()), 'worker',
                              '--spec', str(destination / 'spec.json'), '--deadline', repr(global_deadline)],
                             cwd=str(ROOT), stdout=subprocess.PIPE, stderr=stderr_stream,
                             start_new_session=True,
                             env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1', PYTHONHASHSEED='0'))
    stderr_stream.close()
    selector = selectors.DefaultSelector()
    selector.register(child.stdout, selectors.EVENT_READ)
    with (destination / 'events.jsonl').open('wb') as events_file:
        while True:
            now = time.monotonic()
            if status is None and (now >= min(global_deadline, phase_deadline) or errors):
                watchdog = 'global_deadline' if now >= global_deadline else ('protocol_error' if errors else active_phase)
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                _, status, usage = os.wait4(child.pid, 0)
            select_wait = min(.01, max(0.0, min(global_deadline, phase_deadline)-time.monotonic())) if status is None else 0
            for key, _ in selector.select(select_wait):
                chunk = os.read(key.fileobj.fileno(), 65536)
                if not chunk:
                    selector.unregister(key.fileobj)
                    continue
                events_file.write(chunk)
                events_file.flush()
                buffer += chunk
                while b'\n' in buffer:
                    line, buffer = buffer.split(b'\n', 1)
                    try:
                        msg = json.loads(line)
                        if type(msg) is not dict or not start <= msg['monotonic'] <= time.monotonic():
                            raise ValueError('Invalid event timestamp')
                        events.append(msg)
                        if mode == 'batch' and msg['event'] == 'phase_started':
                            if active_phase not in ('initialization', 'between_phases'):
                                raise ValueError('Nested phase')
                            if next_phase >= len(sequence) or msg['name'] != sequence[next_phase]:
                                raise ValueError('Changed phase order or fixed denominator')
                            next_phase += 1
                            active_phase = msg['name']
                            phase_deadline = min(global_deadline, msg['monotonic'] + (PREPARE if active_phase == 'preparation' else QUERY))
                        elif mode == 'batch' and msg['event'] == 'phase_finished':
                            if msg['name'] != active_phase or msg['monotonic'] > phase_deadline:
                                raise ValueError('Wrong phase or completed after deadline')
                            active_phase = 'between_phases'
                            phase_deadline = min(global_deadline, time.monotonic() + QUERY)
                        elif msg['event'] not in ('ready', 'output_written', 'finished', 'phase_started', 'phase_finished'):
                            raise ValueError('Unknown event')
                    except (ValueError, KeyError, TypeError) as error:
                        errors.append(str(error))
            if status is None:
                waited, wait_status, wait_usage = os.wait4(child.pid, os.WNOHANG)
                if waited:
                    status, usage = wait_status, wait_usage
            if status is not None and not selector.get_map():
                break
            if status is not None and time.monotonic() > global_deadline + .1:
                # A descendant retaining the pipe cannot strand the parent.
                try:
                    os.killpg(child.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
                errors.append('stdout_pipe_held_after_worker_exit')
                break
    if buffer:
        errors.append('unterminated_event_output')
    selector.close()
    child.stdout.close()
    child.returncode = os.waitstatus_to_exitcode(status)  # wait4 already reaped it.
    elapsed = time.monotonic() - start
    measured = process_metrics(usage, elapsed, status, watchdog, errors, start)
    process_ok = child.returncode == 0 and watchdog is None and not errors and elapsed <= limit
    prep = safe_output(spec['preparation_path']) if mode != 'replay' else None
    if prep is not None:
        prep_begins = [e for e in events if e['event'] == 'phase_started' and e['name'] == 'preparation']
        prep_ends = [e for e in events if e['event'] == 'phase_finished' and e['name'] == 'preparation']
        prep['phase_wall_seconds_including_encoding_and_output'] = prep_ends[0]['monotonic']-prep_begins[0]['monotonic'] if prep_begins and prep_ends else None
        prep['phase_cpu_seconds_including_encoding_and_output'] = prep_ends[0]['cpu_seconds']-prep_begins[0]['cpu_seconds'] if prep_begins and prep_ends else None
    if mode == 'batch':
        rows = []
        for item in spec['cases']:
            raw = Path(spec['queries_directory']) / item['id'] / 'output.json'
            row = safe_output(raw)
            row.setdefault('response_json_bytes', 0)
            row.setdefault('certificate_json_bytes', 0)
            begins = [e for e in events if e['event'] == 'phase_started' and e['name'] == item['id']]
            ends = [e for e in events if e['event'] == 'phase_finished' and e['name'] == item['id']]
            completed = len(begins) == len(ends) == 1 and not errors
            if completed:
                completed = ends[0]['monotonic'] - begins[0]['monotonic'] <= QUERY and ends[0]['monotonic'] <= global_deadline
            row.update(case=metadata(item), task=item['task'], task_sha256=digest(item['task']),
                       success=bool(row.get('success') and completed), phase_complete=completed,
                       phase_wall_seconds=ends[0]['monotonic']-begins[0]['monotonic'] if begins and ends else None,
                       phase_cpu_seconds=ends[0]['cpu_seconds']-begins[0]['cpu_seconds'] if begins and ends else None,
                       cumulative_wall_seconds=ends[0]['monotonic']-start if ends else None,
                       cumulative_cpu_seconds=ends[0]['cpu_seconds'] if ends else None,
                       process_peak_rss_high_water_bytes=ends[0]['peak_rss_bytes'] if ends else None,
                       raw_path=str(raw))
            if not completed:
                row['failure_class'] = 'timeout' if watchdog else ('unsupported' if prep and prep.get('response', {}).get('status') == 'unsupported' else 'not_executed_or_failed')
            rows.append(row)
        result = {'system': spec['system'], 'mode': mode, 'family': spec.get('family_metadata'),
                  'repetition': spec.get('repetition'), 'preparation': prep, 'rows': rows,
                  'session_complete': process_ok and next_phase == len(sequence) and all(r['phase_complete'] for r in rows)}
    else:
        result = safe_output(destination / 'output.json')
        result.setdefault('response_json_bytes', 0)
        result.setdefault('certificate_json_bytes', 0)
        result.update(system=spec['system'], mode=mode, case=metadata(spec['case']), task=spec['case']['task'],
                      task_sha256=digest(spec['case']['task']), repetition=spec.get('repetition'),
                      success=bool(result.get('success') and process_ok), within_budget=process_ok,
                      preparation=prep)
        if not process_ok:
            result['failure_class'] = 'timeout' if watchdog else 'process_failure'
        elif mode == 'cold' and not prep.get('ok'):
            result['failure_class'] = 'unsupported' if prep.get('response', {}).get('status') == 'unsupported' else 'preparation_failed'
            result['status'] = result['failure_class']
    result.update(process=measured, events=events, raw_directory=str(destination))
    write(destination / 'result.json', result)
    return result


def distribution(values):
    return {'raw': values, 'min': min(values), 'median': statistics.median(values), 'max': max(values)} if values else None


def common_summary(rows, cases):
    out = {}
    for system in SYSTEMS:
        selected = [r for r in rows if r['system'] == system]
        details = []
        for case in cases:
            runs = [r for r in selected if r['case']['id'] == case['id']]
            good = len(runs) == REPEATS and {r['repetition'] for r in runs} == set(range(REPEATS)) and all(r['success'] for r in runs)
            details.append({**metadata(case), 'kind': case['task']['kind'], 'runs': len(runs),
                'stable_5_of_5': good, 'successful_runs': sum(r['success'] for r in runs),
                'stable_5_of_5_without_fallback': good and not any(r.get('fallback', {}).get('fallback') for r in runs),
                'fallback_runs': sum(r.get('fallback', {}).get('fallback', False) for r in runs),
                'failure_classes': [r.get('failure_class') for r in runs],
                'target_tolerance': case['task']['tolerance'],
                'actual_metric_exact_by_repetition': [{ 'repetition': r['repetition'],
                    'metric':r.get('assessment',{}).get('metric'),
                    'metric_over_tolerance':str(exact(r['assessment']['metric'])/exact(case['task']['tolerance']))
                        if r.get('assessment',{}).get('certificate_valid') else None,
                    'lower':certificate_of(r.get('response')).get('lower') if certificate_of(r.get('response')) else None,
                    'upper':certificate_of(r.get('response')).get('upper') if certificate_of(r.get('response')) else None}
                    for r in runs],
                'wall_seconds': distribution([r['process']['wall_seconds'] for r in runs]),
                'cpu_seconds': distribution([r['process']['cpu_seconds'] for r in runs]),
                'peak_rss_bytes': distribution([r['process']['peak_rss_bytes'] for r in runs]),
                'response_json_bytes': distribution([r['response_json_bytes'] for r in runs if 'response_json_bytes' in r]),
                'certificate_json_bytes': distribution([r['certificate_json_bytes'] for r in runs if 'certificate_json_bytes' in r]),
                'full_verify_wall_seconds': distribution([r['timings']['verification']['wall_seconds'] for r in runs if 'verification' in r.get('timings', {})]),
                'full_verify_cpu_seconds': distribution([r['timings']['verification']['cpu_seconds'] for r in runs if 'verification' in r.get('timings', {})]),
                'routes': [r.get('fallback', {}).get('route') for r in runs]})
        out[system] = {'fixed_case_denominator': len(cases), 'fixed_run_denominator': len(cases)*REPEATS,
                       'returned_rows': len(selected), 'successful_runs': sum(r['success'] for r in selected),
                       'stable_success_cases': sum(r['stable_5_of_5'] for r in details), 'per_case': details,
                       'total_wall_seconds_all_attempts': sum(r['process']['wall_seconds'] for r in selected),
                       'total_cpu_seconds_all_attempts': sum(r['process']['cpu_seconds'] for r in selected),
                       'totals_are_diagnostic_not_mixed_capability_ranking': True,
                       'by_kind': {kind:{'cases':sum(r['kind']==kind for r in details),
                           'run_denominator':REPEATS*sum(r['kind']==kind for r in details),
                           'successful_runs':sum(r['successful_runs'] for r in details if r['kind']==kind),
                           'stable_success_cases':sum(r['kind']==kind and r['stable_5_of_5'] for r in details)}
                           for kind in ('spectrum','approximation')},
                       'by_dimension': {d: {'cases': sum(r.get('dimension') == d for r in details),
                           'run_denominator':REPEATS*sum(r.get('dimension') == d for r in details),
                           'successful_runs':sum(r['successful_runs'] for r in details if r.get('dimension')==d),
                           'stable_success_cases': sum(r.get('dimension') == d and r['stable_5_of_5'] for r in details)}
                           for d in sorted({r.get('dimension', 'unspecified') for r in details})}}
    shared = [case['id'] for case in cases if all(next(r for r in out[s]['per_case'] if r['id'] == case['id'])['stable_5_of_5'] for s in SYSTEMS)]
    for system in SYSTEMS:
        def pairs_for(ids):
            pairs=[]
            for case_id in ids:
                base = sorted([r for r in rows if r['system']=='baseline' and r['case']['id']==case_id],key=lambda r:r['repetition'])
                ours = sorted([r for r in rows if r['system']==system and r['case']['id']==case_id],key=lambda r:r['repetition'])
                ratios={key+'_ratio_by_repetition':[a['process'][key]/b['process'][key] if b['process'][key] else None for a,b in zip(ours,base)]
                        for key in ('wall_seconds','cpu_seconds','peak_rss_bytes')}
                ratios.update({key+'_ratio_by_repetition':[a[key]/b[key] if b[key] else None for a,b in zip(ours,base)]
                               for key in ('response_json_bytes','certificate_json_bytes')})
                original=next(c for c in cases if c['id']==case_id)
                pairs.append({'case_id':case_id,'dimension':original.get('dimension'),'kind':original['task']['kind'],**ratios,
                              'fallback_in_any_run':any(r.get('fallback',{}).get('fallback') for r in ours)})
            return pairs
        out[system]['all_four_systems_5_of_5_common_subset'] = shared
        out[system]['paired_cost_on_common_subset'] = pairs_for(shared)
        pairwise=[c['id'] for c in cases if all(next(r for r in out[s]['per_case'] if r['id']==c['id'])['stable_5_of_5'] for s in ('baseline',system))]
        out[system]['baseline_pairwise_5_of_5_subset']=pairwise
        out[system]['paired_cost_on_baseline_pairwise_subset']=pairs_for(pairwise)
        out[system]['pairwise_subset_warning']='Different pairwise subsets must not be used to rank all methods against each other.'
    return out


def add_replays(session, family, destination):
    replays = []
    for item, row in zip(family['tasks'], session['rows']):
        cert = certificate_of(row.get('response'))
        if cert is None:
            replay = {'success': False, 'not_executed': True, 'reason': 'no_certificate',
                      'process': {'wall_seconds': 0, 'cpu_seconds': 0, 'peak_rss_bytes': None}}
        else:
            replay = run_process({'mode': 'replay', 'system': session['system'], 'case': item,
                                  'certificate': cert, 'repetition': session['repetition']}, destination / item['id'])
        row['full_cold_replay_success'] = replay['success']
        replays.append(replay)
    session['full_cold_replays'] = replays
    session['complete_resource_totals'] = {
        'session_process_wall_seconds':session['process']['wall_seconds'],
        'session_process_cpu_seconds':session['process']['cpu_seconds'],
        'extra_full_cold_replay_wall_seconds':sum(r['process']['wall_seconds'] for r in replays),
        'extra_full_cold_replay_cpu_seconds':sum(r['process']['cpu_seconds'] for r in replays),
        'total_including_session_exit_and_full_cold_replays_wall_seconds':session['process']['wall_seconds']+sum(r['process']['wall_seconds'] for r in replays),
        'total_including_session_exit_and_full_cold_replays_cpu_seconds':session['process']['cpu_seconds']+sum(r['process']['cpu_seconds'] for r in replays),
        'largest_serial_process_peak_rss_bytes':max([session['process']['peak_rss_bytes']]+[r['process']['peak_rss_bytes'] for r in replays if r['process']['peak_rss_bytes'] is not None]),
        'peak_rss_is_max_not_sum':True}
    session['checkpoints'] = {}
    for n in (1, 4, 16):
        last = session['rows'][n-1]
        prefix_complete = last['cumulative_wall_seconds'] is not None
        wall = last['cumulative_wall_seconds'] if prefix_complete else session['process']['wall_seconds']
        cpu_value = last['cumulative_cpu_seconds'] if prefix_complete else session['process']['cpu_seconds']
        replay_wall = sum(r['process']['wall_seconds'] for r in replays[:n])
        replay_cpu = sum(r['process']['cpu_seconds'] for r in replays[:n])
        session['checkpoints'][str(n)] = {'queries': n, 'prefix_complete': prefix_complete,
            'successful_queries': sum(r['success'] and r['full_cold_replay_success'] for r in session['rows'][:n]),
            'including_preparation_wall_seconds': wall, 'including_preparation_cpu_seconds': cpu_value,
            'extra_full_cold_replay_wall_seconds': replay_wall, 'extra_full_cold_replay_cpu_seconds': replay_cpu,
            'total_with_full_cold_replay_wall_seconds': wall+replay_wall,
            'total_with_full_cold_replay_cpu_seconds': cpu_value+replay_cpu,
            'failed_prefix_cost_policy': 'charge all observed session cost if requested prefix was not reached'}
    return session


def batch_summary(sessions, families):
    result = {}
    for family in families:
        by_system = {}
        for system in SYSTEMS:
            runs = [s for s in sessions if s['system'] == system and s['family']['id'] == family['id']]
            data = {'fixed_sessions': REPEATS, 'returned_sessions': len(runs),
                    'fixed_query_denominator': REPEATS*16,
                    'successful_queries_with_cold_replay': sum(r['success'] and r['full_cold_replay_success'] for s in runs for r in s['rows']),
                    'complete_sessions': sum(s['session_complete'] for s in runs),
                    'session_peak_rss_bytes': distribution([s['process']['peak_rss_bytes'] for s in runs]),
                    'total_session_wall_seconds_all_attempts': sum(s['process']['wall_seconds'] for s in runs),
                    'total_session_cpu_seconds_all_attempts': sum(s['process']['cpu_seconds'] for s in runs),
                    'complete_total_wall_seconds_including_exit_and_cold_replays':distribution([s['complete_resource_totals']['total_including_session_exit_and_full_cold_replays_wall_seconds'] for s in runs]),
                    'complete_total_cpu_seconds_including_exit_and_cold_replays':distribution([s['complete_resource_totals']['total_including_session_exit_and_full_cold_replays_cpu_seconds'] for s in runs]),
                    'preparation_response_json_bytes':distribution([s['preparation'].get('preparation_response_json_bytes',0) for s in runs]),
                    'preparation_certificate_json_bytes':distribution([s['preparation'].get('preparation_certificate_json_bytes',0) for s in runs]),
                    'serialized_library_json_bytes':distribution([s['preparation'].get('serialized_library_json_bytes',0) for s in runs]),
                    'preparation_phase_wall_seconds':distribution([s['preparation']['phase_wall_seconds_including_encoding_and_output'] for s in runs if s['preparation']['phase_wall_seconds_including_encoding_and_output'] is not None]),
                    'query_response_json_bytes':distribution([r['response_json_bytes'] for s in runs for r in s['rows']]),
                    'query_certificate_json_bytes':distribution([r['certificate_json_bytes'] for s in runs for r in s['rows']]),
                    'full_cold_replay_peak_rss_bytes':distribution([r['process']['peak_rss_bytes'] for s in runs for r in s['full_cold_replays'] if not r.get('not_executed')]),
                    'checkpoints': {}}
            data['per_query']=[{**metadata(item),
                'fixed_repetitions':REPEATS,
                'successful_repetitions':sum(s['rows'][i]['success'] and s['rows'][i]['full_cold_replay_success'] for s in runs),
                'stable_5_of_5':len(runs)==REPEATS and all(s['rows'][i]['success'] and s['rows'][i]['full_cold_replay_success'] for s in runs),
                'failure_classes':[s['rows'][i].get('failure_class') for s in runs],
                'actual_metric_exact':[s['rows'][i].get('assessment',{}).get('metric') for s in runs],
                'target_tolerance':item['task']['tolerance']}
                for i,item in enumerate(family['tasks'])]
            for n in (1,4,16):
                points = [s['checkpoints'][str(n)] for s in runs]
                data['checkpoints'][str(n)] = {'fixed_query_denominator': REPEATS*n,
                    'successful_queries': sum(p['successful_queries'] for p in points),
                    'all_five_sessions_all_queries_success': len(points) == REPEATS and all(p['successful_queries'] == n for p in points),
                    **{key: distribution([p[key] for p in points]) for key in (
                        'including_preparation_wall_seconds','including_preparation_cpu_seconds',
                        'total_with_full_cold_replay_wall_seconds','total_with_full_cold_replay_cpu_seconds')}}
            for label, indices in (('first_query',range(1)),('later_queries',range(1,16))):
                for metric in ('wall','cpu'):
                    key = 'phase_'+metric+'_seconds'
                    data[label+'_'+metric+'_seconds'] = distribution([s['rows'][i][key] for s in runs for i in indices if s['rows'][i][key] is not None])
            by_system[system] = data
        for system,data in by_system.items():
            data['paired_prefix_cost_vs_baseline']={}
            for n in (1,4,16):
                key=str(n)
                valid=all(by_system[s]['checkpoints'][key]['all_five_sessions_all_queries_success'] for s in ('baseline',system))
                theirs=sorted([s for s in sessions if s['system']=='baseline' and s['family']['id']==family['id']],key=lambda s:s['repetition'])
                ours=sorted([s for s in sessions if s['system']==system and s['family']['id']==family['id']],key=lambda s:s['repetition'])
                data['paired_prefix_cost_vs_baseline'][key]={'eligible':valid,
                    'wall_ratio_including_preparation_and_replays':[a['checkpoints'][key]['total_with_full_cold_replay_wall_seconds']/b['checkpoints'][key]['total_with_full_cold_replay_wall_seconds'] for a,b in zip(ours,theirs)] if valid else [],
                    'cpu_ratio_including_preparation_and_replays':[a['checkpoints'][key]['total_with_full_cold_replay_cpu_seconds']/b['checkpoints'][key]['total_with_full_cold_replay_cpu_seconds'] for a,b in zip(ours,theirs)] if valid else []}
        result[family['id']] = {'family': metadata(family), 'systems': by_system}
    return result


def run():
    cases, families, protocol = inputs()
    pointer = checked_freeze()
    before = identity()
    destination = ROOT / 'evaluation' / ('final_' + stamp())
    destination.mkdir(parents=True)
    write(destination / 'run_started.json', {'version': VERSION, 'freeze': pointer, 'identity_before': before,
                                           'cases_sha256': digest(cases), 'families_sha256': digest(families)})
    common, sessions = [], []
    for repeat in range(REPEATS):
        systems = SYSTEMS[repeat % len(SYSTEMS):] + SYSTEMS[:repeat % len(SYSTEMS)]
        order = cases[repeat % len(cases):] + cases[:repeat % len(cases)]
        for case in order:
            for system in systems:
                row = run_process({'mode':'cold','system':system,'case':case,
                                   'family':resolve_family(case,families),'repetition':repeat},
                                  destination / ('cold_'+system+'_'+case['id']+'_r'+str(repeat)))
                common.append(row)
                write(destination / 'common_checkpoint.json', {'complete':False,'rows':common})
    for repeat in range(REPEATS):
        systems = SYSTEMS[repeat % len(SYSTEMS):] + SYSTEMS[:repeat % len(SYSTEMS)]
        for family in families[repeat % 2:] + families[:repeat % 2]:
            for system in systems:
                session = run_process({'mode':'batch','system':system,'family':family['spec'],
                    'family_metadata':metadata(family),'cases':family['tasks'],'repetition':repeat},
                    destination / ('batch_'+system+'_'+family['id']+'_r'+str(repeat)))
                session = add_replays(session,family,destination / ('replay_'+system+'_'+family['id']+'_r'+str(repeat)))
                sessions.append(session)
                write(destination / 'batch_checkpoint.json', {'complete':False,'sessions':sessions})
    try:
        after = identity()
        unchanged = before == after and checked_freeze() == pointer
        identity_error = None
    except (OSError, ValueError, RuntimeError) as error:
        after = None;unchanged = False;identity_error = str(error)
    complete = len(common)==18*4*5 and len(sessions)==2*4*5 and all(len(s['rows'])==16 and len(s['full_cold_replays'])==16 for s in sessions)
    report = {'version':VERSION,'classification':'public_development_not_holdout','complete_fixed_denominators':complete,
        'valid_frozen_comparison':unchanged and complete,'identity_before':before,'identity_after':after,
        'identity_error':identity_error,'freeze':pointer,'common_rows':common,'common_summary':common_summary(common,cases),
        'batch_sessions':sessions,'batch_summary':batch_summary(sessions,families),
        'actual_model_api_calls':0,'actual_model_tokens':None,'actual_cached_model_tokens':None,
        'metrics_note':'JSON UTF-8 bytes measure serialized volume, never token usage. Batch RSS is cumulative high water, not per-query delta.',
        'timing_note':'Five serial interleaved repetitions on a shared machine; compare only all-system stable-success cases. Background load can affect wall time.',
        'complete_verification_note':'Cold tasks include full replay inside their 10-second whole-process deadline. Batch queries add separately limited, separately charged cold full replays.'}
    write(destination / 'report.json',report)
    if not unchanged or not complete:
        raise RuntimeError('Comparison saved but frozen identity or fixed denominator invalid')
    return {'report_path':str(destination / 'report.json'),'valid_frozen_comparison':True,
            'common_rows':len(common),'batch_sessions':len(sessions),'batch_queries':sum(len(s['rows']) for s in sessions)}


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
        cases,families,_=inputs();current=identity()
        result={'cases':len(cases),'families':len(families),'all_preserved':True,
                'original_baseline_files':current['original_baseline_files'],
                'previous_sources':current['previous_sources'],'previous_artifacts':current['previous_artifacts']}
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
