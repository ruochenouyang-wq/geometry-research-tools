"""Evidence-preserving interaction helpers; no model or network calls.

Executors and verifiers are supplied by the caller. A summary is not a proof;
its evidence_state records whether an explicit verifier actually accepted it.
"""
from copy import deepcopy
from collections import OrderedDict
from contextlib import contextmanager
from fractions import Fraction
import hashlib
import json
import math
import os
from pathlib import Path
import importlib.util
import sys
import time


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False,
                      separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode('utf-8')).hexdigest()


def exact_fraction(value):
    if type(value) not in (int, str):
        raise ValueError('Exact values must be integers or rational strings')
    return Fraction(value)


def outward_decimal(value, places=12, upper=False):
    """Integer-only directed rounding, including negatives and tiny values."""
    if type(places) is not int or not 0 <= places <= 1000:
        raise ValueError('Decimal places must be an integer in 0..1000')
    value = exact_fraction(value)
    scale = 10 ** places
    numerator = value.numerator * scale
    rounded = -((-numerator) // value.denominator) if upper else numerator // value.denominator
    sign = '-' if rounded < 0 else ''
    integer, fractional = divmod(abs(rounded), scale)
    return sign+str(integer)+(('.'+str(fractional).zfill(places)) if places else '')


def exact_bounds(certificate, places=12):
    if 'error_upper' in certificate:
        error = exact_fraction(certificate['error_upper'])
        squared = (exact_fraction(certificate['error_squared'])
                   if 'error_squared' in certificate else None)
        if error < 0 or (squared is not None and (squared < 0 or error*error < squared)):
            raise ValueError('Invalid approximation error bound')
        bounds = {'error_upper': str(error),
                  'decimal_error_upper': outward_decimal(str(error), places, upper=True)}
        if squared is not None:
            bounds['error_squared'] = str(squared)
        return bounds
    lower = exact_fraction(certificate['lower'])
    upper = exact_fraction(certificate['upper'])
    if lower > upper:
        raise ValueError('Reversed spectral interval')
    width = upper-lower
    if 'exact_width' in certificate and exact_fraction(certificate['exact_width']) != width:
        raise ValueError('Width disagrees with exact endpoints')
    return {'lower': str(lower), 'upper': str(upper), 'exact_width': str(width),
            'decimal_lower': outward_decimal(str(lower), places),
            'decimal_upper': outward_decimal(str(upper), places, upper=True)}


def _call(ledger, stage, function, *args):
    if ledger is None:
        return function(*args)
    with ledger.stage(stage):
        return function(*args)


def summarize_result(result, request=None, verifier=None, places=12, ledger=None):
    """Preserve target, scope and verification status, never just a number.

    verifier(certificate) must return the literal boolean True to certify.
    Bind it to the request in the caller when request-specific replay is needed.
    """
    if not isinstance(result, dict):
        raise ValueError('Result must be an object')
    if result.get('ok') is False:
        return {'ok': False, 'status': 'execution_failed',
                'request': deepcopy(request),
                'error': {'type': result.get('error'), 'reason': result.get('reason')}}
    certificate = result.get('certificate')
    if not isinstance(certificate, dict):
        raise ValueError('A mathematical result requires a certificate')
    for field in ('format', 'function', 'scope'):
        if field not in certificate:
            raise ValueError('Certificate missing '+field)
    tolerance = certificate.get('tolerance')
    if isinstance(request, dict) and 'tolerance' in request:
        requested = exact_fraction(request['tolerance'])
        if tolerance is not None and requested != exact_fraction(tolerance):
            raise ValueError('Requested tolerance disagrees with certificate')
        tolerance = str(requested)
    if tolerance is not None:
        tolerance = exact_fraction(tolerance)
        if tolerance <= 0:
            raise ValueError('Tolerance must be positive')
    checked = None if verifier is None else _call(ledger, 'verification', verifier, deepcopy(certificate)) is True
    summary = {
        'ok': checked is not False, 'format': 'interaction_summary_v1',
        'request': deepcopy(request),
        'target': {'tolerance': None if tolerance is None else str(tolerance),
                   'quantity': 'function_L2_error' if 'error_upper' in certificate
                   else 'lowest_spectral_value'},
        'scope': {key: deepcopy(certificate[key]) for key in (
            'function', 'geometry', 'scope', 'mean_zero', 'measure', 'norm',
            'full_infinite_space_covered', 'spectral_transfer_claimed') if key in certificate},
        'evidence_state': 'unchecked' if checked is None else
                          ('verified' if checked else 'rejected'),
        'status': certificate.get('status', 'approximation_bound'),
        'certificate_format': certificate['format'],
        'certificate_sha256': digest(certificate),
        'bounds': exact_bounds(certificate, places),
    }
    if checked is True:
        bound = summary['bounds'].get('exact_width', summary['bounds'].get('error_upper'))
        summary['status'] = ('certified_bound' if tolerance is None else
                            ('certified_met' if exact_fraction(bound) <= tolerance else 'certified_open'))
    if checked is not True:
        summary['status'] = 'unchecked' if checked is None else 'verification_failed'
    return summary


def present_result(result, request=None, verifier=None, store=None,
                   include_certificate=False, places=12, ledger=None):
    """Store certificate once through store(sha256, certificate) -> locator.

    Without a store, include the full certificate so evidence is not lost.
    The store owns durability; a hash alone is not a mathematical verifier.
    """
    summary = summarize_result(result, request, verifier, places, ledger)
    certificate = result.get('certificate')
    if not isinstance(certificate, dict):
        return summary
    if store is not None:
        key = digest(certificate)
        locator = _call(ledger, 'store', store, key, deepcopy(certificate))
        if not isinstance(locator, str) or not locator:
            raise ValueError('Store must return a nonempty evidence locator')
        summary['evidence_ref'] = {'sha256': key, 'locator': locator,
                                   'bytes': len(canonical(certificate).encode('utf-8'))}
    if store is None or include_certificate:
        summary['certificate'] = deepcopy(certificate)
    return summary


def retrieve_evidence(reference, load, verifier=None, require_verified=False, ledger=None):
    """load(locator) -> certificate; reject content substitution before replay."""
    if not isinstance(reference, dict):
        raise ValueError('Evidence reference must be an object')
    key = reference.get('sha256')
    if (not isinstance(key, str) or len(key) != 64 or
            any(character not in '0123456789abcdef' for character in key)):
        raise ValueError('Invalid SHA-256 evidence key')
    if not isinstance(reference.get('locator'), str) or not reference['locator']:
        raise ValueError('Missing evidence locator')
    certificate = _call(ledger, 'load', load, reference['locator'])
    if not isinstance(certificate, dict) or digest(certificate) != key:
        raise ValueError('Evidence content does not match its hash')
    if 'bytes' in reference and reference['bytes'] != len(canonical(certificate).encode('utf-8')):
        raise ValueError('Evidence byte count mismatch')
    checked = None if verifier is None else _call(ledger, 'verification', verifier, deepcopy(certificate)) is True
    if require_verified and checked is not True:
        raise ValueError('Mathematical replay was not accepted')
    return {'certificate': deepcopy(certificate), 'sha256': key,
            'evidence_state': 'unchecked' if checked is None else
                              ('verified' if checked else 'rejected')}


def request_key(request, namespace):
    """Bind ALL supplied semantics; never merge omitted and explicit defaults.

    namespace identifies executor/theorem versions. Exact-input equality is a
    deliberately conservative deduplication rule, not algebraic equivalence.
    """
    if not isinstance(request, dict) or not isinstance(namespace, str) or not namespace:
        raise ValueError('Request object and nonempty executor namespace required')
    def validate(value):
        if isinstance(value, dict):
            if any(type(key) is not str for key in value):
                raise ValueError('JSON object keys must be strings')
            for child in value.values():
                validate(child)
        elif isinstance(value, list):
            for child in value:
                validate(child)
        elif type(value) not in (str, int, float, bool, type(None)):
            raise ValueError('Request must contain JSON values')
    validate(request)
    return digest({'executor_namespace': namespace, 'request': request})


class BudgetExceeded(RuntimeError):
    pass


class UsageLedger:
    """Count every attempted request and every started execution, including errors."""
    def __init__(self, max_requests=None, max_executions=None, max_wall_seconds=None, clock=None):
        for value in (max_requests, max_executions):
            if value is not None and (type(value) is not int or value < 0):
                raise ValueError('Count limits must be nonnegative integers or None')
        self.max_requests = max_requests
        self.max_executions = max_executions
        if max_wall_seconds is not None and (type(max_wall_seconds) not in (int,float) or
                not math.isfinite(max_wall_seconds) or max_wall_seconds < 0):
            raise ValueError('Wall limit must be finite, nonnegative, or None')
        self.max_wall_seconds = max_wall_seconds
        self.clock = clock or time.perf_counter
        self.started_at = self.clock()
        self.stages = {}
        self.counts = dict(received=0, admitted=0, executor_calls=0,
                           successes=0, failures=0, cache_hits=0, budget_rejections=0)

    def admit_request(self):
        self.counts['received'] += 1
        if self.max_wall_seconds is not None and self.clock()-self.started_at >= self.max_wall_seconds:
            self.counts['budget_rejections'] += 1
            raise BudgetExceeded('Elapsed wall-time budget exhausted')
        if self.max_requests is not None and self.counts['admitted'] >= self.max_requests:
            self.counts['budget_rejections'] += 1
            raise BudgetExceeded('Request budget exhausted')
        self.counts['admitted'] += 1

    def start_execution(self):
        if self.max_executions is not None and self.counts['executor_calls'] >= self.max_executions:
            self.counts['budget_rejections'] += 1
            raise BudgetExceeded('Execution budget exhausted')
        self.counts['executor_calls'] += 1

    @contextmanager
    def stage(self, name):
        """Accumulate stage calls and exception time; nested stages overlap."""
        if not isinstance(name, str) or not name:
            raise ValueError('Stage name required')
        start = self.clock()
        failed = False
        try:
            yield
        except BaseException:
            failed = True
            raise
        finally:
            entry = self.stages.setdefault(name, {'calls':0,'exceptions':0,'seconds':0.0})
            entry['calls'] += 1
            entry['exceptions'] += int(failed)
            entry['seconds'] += max(0.0,self.clock()-start)

    def snapshot(self):
        return {'counts': dict(self.counts),
                'limits': {'requests':self.max_requests,'executions':self.max_executions,
                           'wall_seconds':self.max_wall_seconds},
                'elapsed_wall_seconds':max(0.0,self.clock()-self.started_at),
                'stages':deepcopy(self.stages),
                'timing_note':'Elapsed wall includes idle time; stages may overlap and must not be summed as wall time. Admission limits do not interrupt an active call.',
                'actual_api_usage': None, 'model_calls': 0}


class RequestRunner:
    """Deduplicate deterministic calls; disable cache for fresh/stochastic work.

    execute returns a transport envelope whose result is the executor output.
    Cached execution is not mathematical verification. New theorem/code versions
    require a new runner namespace. This small runner is sequential, not thread safe.
    """
    def __init__(self, executor, namespace, cache=True, ledger=None,
                 max_cache_entries=128, max_cache_bytes=16*1024*1024):
        if not callable(executor):
            raise ValueError('Executor must be callable')
        request_key({}, namespace)
        self.executor = executor
        self.namespace = namespace
        self.cache_enabled = bool(cache)
        for value in (max_cache_entries, max_cache_bytes):
            if type(value) is not int or value < 0:
                raise ValueError('Cache limits must be nonnegative integers')
        self.max_cache_entries = max_cache_entries
        self.max_cache_bytes = max_cache_bytes
        self._cache = OrderedDict()
        self.cache_bytes = 0
        self.cache_stats = dict(hits=0, misses=0, writes=0, evictions=0, oversized_skips=0)
        self.ledger = ledger if ledger is not None else UsageLedger()

    def cache_snapshot(self):
        return {**self.cache_stats, 'entries':len(self._cache),
                'serialized_payload_bytes':self.cache_bytes,
                'max_entries':self.max_cache_entries,'max_payload_bytes':self.max_cache_bytes,
                'note':'Payload bounds exclude Python/key overhead. Warm means this runner has cached the exact input, not a fresh-process comparison.'}

    def clear_cache(self):
        self._cache.clear()
        self.cache_bytes = 0

    def _cache_write(self, key, result):
        payload = canonical(result)
        size = len(payload.encode('utf-8'))
        if not self.max_cache_entries or size > self.max_cache_bytes:
            self.cache_stats['oversized_skips'] += 1
            return
        while self._cache and (len(self._cache) >= self.max_cache_entries or
                               self.cache_bytes+size > self.max_cache_bytes):
            _, (_, old_size) = self._cache.popitem(last=False)
            self.cache_bytes -= old_size
            self.cache_stats['evictions'] += 1
        self._cache[key] = (payload, size)
        self.cache_bytes += size
        self.cache_stats['writes'] += 1

    def execute(self, request):
        try:
            self.ledger.admit_request()
            key = _call(self.ledger, 'request_key', request_key, request, self.namespace)
            hit = self.cache_enabled and key in self._cache
            if hit:
                self.cache_stats['hits'] += 1
                self.ledger.counts['cache_hits'] += 1
                self._cache.move_to_end(key)
                result = _call(self.ledger, 'cache_read', json.loads, self._cache[key][0])
            else:
                self.cache_stats['misses'] += 1
                self.ledger.start_execution()
                result = _call(self.ledger, 'execution', self.executor, deepcopy(request))
                if not isinstance(result, dict):
                    raise ValueError('Executor must return an object')
                if self.cache_enabled and result.get('ok') is not False:
                    _call(self.ledger, 'cache_write', self._cache_write, key, result)
            self.ledger.counts['failures' if result.get('ok') is False else 'successes'] += 1
            return {'result': result, 'cache_hit': hit, 'request_key': key,
                    'cache_condition': ('warm_cache_hit' if hit else
                        ('cold_cache_miss' if self.cache_enabled else 'cache_disabled')),
                    'cache': self.cache_snapshot(), 'usage': self.ledger.snapshot()}
        except Exception:
            self.ledger.counts['failures'] += 1
            raise

    def batch(self, requests):
        """Preserve input order; one failure must not cancel later requests."""
        if not isinstance(requests, (list, tuple)):
            raise ValueError('Batch requests must be a finite list or tuple')
        responses = []
        for index, request in enumerate(requests):
            try:
                response = self.execute(request)
            except Exception as error:
                response = {'result': {'ok': False, 'status': 'execution_failed',
                                      'error': type(error).__name__, 'reason': str(error)},
                            'cache_hit': False, 'request_key': None, 'usage':self.ledger.snapshot()}
            responses.append({'index': index, **response})
        return responses


BOUND_PROTOCOL = ('bound_interaction_summary_v1: request_ref is [request_id,SHA256(canonical original request)]. '
    'Inherited scope/target fields come from that visible request; certificate_sha256 comes from evidence_ref. '
    'Restore locator from locator_parts[0]+sha256+locator_parts[1]. '
    'Canonical JSON sorts keys, uses compact separators and UTF-8. Exact fractions and outward decimals remain unchanged.')


def reasonable_summary(summary):
    """Strong standalone baseline: remove repeated request and duplicate digest."""
    result = deepcopy(summary)
    result.pop('request', None)
    if result.get('certificate_sha256') == result.get('evidence_ref', {}).get('sha256'):
        result.pop('certificate_sha256')
    return result


def bind_summary(summary, request, request_id, protocol_known=False):
    """Deduplicate only fields recoverable from the visible, hash-bound request."""
    if not isinstance(request_id, str) or not request_id:
        raise ValueError('A visible request identifier is required')
    if summary.get('format') != 'interaction_summary_v1':
        raise ValueError('Only mathematical summaries can be request-bound')
    if summary.get('request') is not None and canonical(summary['request']) != canonical(request):
        raise ValueError('Summary was made for a different request')
    response = reasonable_summary(summary)
    response['format'] = 'bound_interaction_summary_v1'
    response['request_ref'] = [request_id, digest(request)]
    inherited = []
    for key in ('function', 'mean_zero'):
        if key in request and response['scope'].get(key) == request[key]:
            response['scope'].pop(key)
            inherited.append(key)
    if ('tolerance' in request and response['target'].get('tolerance') is not None and
            exact_fraction(request['tolerance']) == exact_fraction(response['target']['tolerance'])):
        response['target'].pop('tolerance')
        inherited.append('tolerance')
    response['inherit'] = inherited
    reference = response.get('evidence_ref')
    if isinstance(reference,dict):
        locator = reference.get('locator')
        key = reference.get('sha256')
        if isinstance(locator,str) and isinstance(key,str) and len(key) == 64 and locator.count(key) == 1:
            reference['locator_parts'] = locator.split(key)
            reference.pop('locator')
    if not protocol_known:
        response['protocol'] = BOUND_PROTOCOL
    return response


def expand_bound_summary(response, request, request_id):
    """Recover the exact summary semantics; reject a different request/context."""
    if response.get('format') != 'bound_interaction_summary_v1':
        raise ValueError('Not a bound summary')
    if response.get('request_ref') != [request_id, digest(request)]:
        raise ValueError('Request identity or content mismatch')
    inherited = response.get('inherit')
    if (not isinstance(inherited, list) or len(set(inherited)) != len(inherited) or
            not set(inherited) <= {'function','mean_zero','tolerance'}):
        raise ValueError('Unknown inherited fields')
    result = deepcopy(response)
    for key in ('request_ref','inherit','protocol'):
        result.pop(key, None)
    result['format'] = 'interaction_summary_v1'
    result['request'] = deepcopy(request)
    reference = result.get('evidence_ref')
    if isinstance(reference,dict) and 'locator_parts' in reference:
        parts = reference.pop('locator_parts')
        if (not isinstance(parts,list) or len(parts) != 2 or
                not all(isinstance(part,str) for part in parts) or 'locator' in reference):
            raise ValueError('Invalid locator reconstruction')
        reference['locator'] = parts[0]+reference['sha256']+parts[1]
    for key in inherited:
        if key not in request:
            raise ValueError('Inherited field absent from visible request')
        target = result['target'] if key == 'tolerance' else result['scope']
        if key in target:
            raise ValueError('Inherited field also supplied explicitly')
        target[key] = str(exact_fraction(request[key])) if key == 'tolerance' else deepcopy(request[key])
    if 'evidence_ref' in result and 'certificate_sha256' not in result:
        result['certificate_sha256'] = result['evidence_ref']['sha256']
    return result


def offline_text_counter(cache_dir):
    """Load pinned public encodings from verified LOCAL cache only.

    Cache preparation is explicit, outside this function. No downloading or
    character-based fallback is permitted. Returns text -> {encoding: count}.
    """
    cache_dir = Path(cache_dir).resolve()
    expected = {
        '9b5ad71b2ce5302211f9c61530b329a4922fc6a4':
            '223921b76ee99bde995b7ff738513eef100fb51d18c93597a113bcffe865b2a7',
        'fb374d419588a4632f3f557e76b4b70aebbca790':
            '446a9538cb6c348e3516120d7c08b09f57c36495e2acfffe59a5bf8b0cfb1a2d',
    }
    for filename, expected_hash in expected.items():
        path = cache_dir/filename
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != expected_hash:
            raise RuntimeError('Offline encoding cache unavailable or invalid: '+filename)
    root = Path(__file__).resolve().parent
    dependency = root.parents[1]/'work'/'token-deps'
    if dependency.is_dir() and str(dependency) not in sys.path:
        sys.path.insert(0,str(dependency))
    import tiktoken
    if tiktoken.__version__ != '0.11.0':
        raise RuntimeError('Pinned offline token measurement requires tiktoken 0.11.0')
    path = root.parent/'geometry_v36_v90'/'token_meter.py'
    spec = importlib.util.spec_from_file_location('_interaction_frozen_token_meter', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    previous = os.environ.get('TIKTOKEN_CACHE_DIR')
    os.environ['TIKTOKEN_CACHE_DIR'] = str(cache_dir)
    try:
        encoders = {name:module.encoder(name) for name in module.ENCODINGS}
    finally:
        if previous is None:
            os.environ.pop('TIKTOKEN_CACHE_DIR',None)
        else:
            os.environ['TIKTOKEN_CACHE_DIR'] = previous
    return lambda text: {name:len(encoder.encode(text,disallowed_special=()))
                         for name,encoder in encoders.items()}


def measure_conversation(events, counter=None, instructions='', schema=''):
    """Count all explicit items, including requests, failures and detail reads."""
    texts = [instructions, schema]
    for event in events:
        if not isinstance(event,dict) or not isinstance(event.get('text'),str):
            raise ValueError('Every event needs its actual text')
        texts.append(event['text'])
    if any(not isinstance(text,str) for text in texts):
        raise ValueError('Instructions and schema must be strings')
    report = {'items':len(texts),'basis':'Sum of explicit text-item encodings; excludes unknown API framing and hidden reasoning.',
              'model_mapping_asserted':False,'actual_api_usage':None}
    if counter is None:
        return {**report,'available':False,'encodings':None,'reason':'No offline text encoder supplied'}
    totals = None
    for text in texts:
        counts = counter(text)
        if (not isinstance(counts,dict) or not counts or
                any(type(value) is not int or value < 0 for value in counts.values())):
            raise ValueError('Counter must return nonnegative integer encoding counts')
        if totals is None:
            totals = {key:0 for key in counts}
        if totals.keys() != counts.keys():
            raise ValueError('Encoding set changed within one conversation')
        for key,value in counts.items():
            totals[key] += value
    return {**report,'available':True,'encodings':totals}


def choose_response(summary, request, request_id, counter=None, protocol_known=False):
    """Use binding only when no measured encoding becomes longer.

    The original request must already be visible. If the protocol is not already
    in context it is transmitted and counted in the candidate, not assumed free.
    Selection overhead is local computation and should be timed by the caller.
    """
    baseline = reasonable_summary(summary)
    if counter is None:
        return {'strategy':'standalone','response':baseline,'comparison':None,
                'actual_api_usage':None}
    candidate = bind_summary(summary, request, request_id, protocol_known)
    before = measure_conversation([{'text':canonical(baseline)}], counter)['encodings']
    after = measure_conversation([{'text':canonical(candidate)}], counter)['encodings']
    better = all(after[key] <= before[key] for key in before) and any(after[key] < before[key] for key in before)
    return {'strategy':'request_bound' if better else 'standalone',
            'response':candidate if better else baseline,
            'comparison':{'standalone':before,'request_bound':after},
            'actual_api_usage':None}
