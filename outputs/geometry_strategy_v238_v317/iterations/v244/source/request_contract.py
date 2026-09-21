"""Request boundary for the frozen mathematical backend.

This module does not change mathematical algorithms or their certificates.
"""
import json
import math
import re
from copy import deepcopy
import backend

METADATA = {'expected_kind', 'request_id', 'expected_request_digest'}
KINDS = {'spectrum', 'approximation'}
FIELDS = {
    'precise_singular': {'op', 'function', 'max_terms', 'tolerance', 'precision_bits', 'iterations'},
    'spectrum': {'op', 'function', 'mean_zero', 'modes', 'bits', 'max_m', 'tolerance', 'sqrt_bits', 'near_tail'},
    'approximate': {'op', 'function', 'levels', 'degree', 'root_ratio', 'sqrt_bits'},
    'verify': {'op', 'certificate', 'function', 'mean_zero', 'tolerance'},
}
DEFAULTS = {
    'precise_singular': {'max_terms': 16, 'precision_bits': 112, 'iterations': 14, 'tolerance': '1/100000000'},
    'spectrum': {'mean_zero': True, 'modes': 8, 'bits': 32, 'max_m': 4,
                 'tolerance': '1/100000000', 'sqrt_bits': 40, 'near_tail': 0},
    'approximate': {'levels': 4, 'degree': 3, 'root_ratio': '1/2', 'sqrt_bits': 40},
    'verify': {},
}


def semantic_request(request):
    """Normalize explicit mathematical defaults; preserve omitted verify bindings."""
    data = {**DEFAULTS[request['op']], **deepcopy(request)}
    for key in ('request_id', 'expected_request_digest'):
        data.pop(key, None)
    if 'function' in data:
        normalizer = backend.pieces.normalize_function if data['op'] == 'approximate' else backend.direct.normalize
        data['function'] = normalizer(data['function'])
    for key in ('tolerance', 'root_ratio'):
        if key in data:
            data[key] = str(backend.pieces.rational(data[key]))
    return data


def request_digest(request):
    return backend.digest(semantic_request(request))


def certificate_kind(certificate):
    if not isinstance(certificate, dict):
        return None
    return {backend.direct.FULL: 'spectrum', backend.enriched.FORMAT: 'spectrum',
            backend.pieces.FORMAT: 'approximation'}.get(certificate.get('format'))


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('Duplicate JSON key: ' + key)
        result[key] = value
    return result


def parse_request(text):
    return validate_request(json.loads(text, object_pairs_hook=_unique_object,
                                       parse_constant=_reject_constant))


def _reject_constant(value):
    raise ValueError('Non-finite JSON number: ' + value)


def _finite_values(value):
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError('Non-finite request number')
    if isinstance(value, dict):
        for item in value.values():
            _finite_values(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _finite_values(item)


def _tolerance(value):
    value = backend.pieces.rational(value)
    if not backend.F(1, 10**30) <= value <= 1:
        raise ValueError('Tolerance must be in 1e-30..1')
    return value


def _validate_singular(request):
    if 'function' not in request or backend.direct.normalize(request['function']) != backend.enriched.FUNCTION:
        raise ValueError('precise_singular supports the fixed original singular function only')
    for name, default, low, high in (('max_terms', 16, 1, 16),
                                    ('precision_bits', 112, 32, 256),
                                    ('iterations', 14, 0, 64)):
        backend.direct.integer(request.get(name, default), low, high, name)
    _tolerance(request.get('tolerance', '1/100000000'))


def validate_request(request):
    if not isinstance(request, dict):
        raise ValueError('Request must be an object')
    _finite_values(request)
    op = request.get('op')
    if type(op) is not str or op not in FIELDS:
        raise ValueError('Unknown operation')
    if not set(request) <= FIELDS[op] | METADATA:
        raise ValueError('Unknown request field')
    for key in ('function', 'mean_zero', 'tolerance'):
        if key in request and request[key] is None:
            raise ValueError('Explicit null binding: ' + key)
    if 'expected_kind' in request:
        if type(request['expected_kind']) is not str or request['expected_kind'] not in KINDS:
            raise ValueError('expected_kind must be spectrum or approximation')
        produced = {'spectrum': 'spectrum', 'precise_singular': 'spectrum',
                    'approximate': 'approximation'}.get(request.get('op'))
        if produced is not None and produced != request['expected_kind']:
            raise ValueError('Requested operation and expected_kind differ')
    if request.get('op') == 'precise_singular':
        _validate_singular(request)
    elif op in ('spectrum', 'approximate'):
        if 'function' not in request:
            raise ValueError('Missing original function')
        normalizer = backend.pieces.normalize_function if op == 'approximate' else backend.direct.normalize
        normalizer(request['function'])
    elif 'function' in request:
        backend.direct.normalize(request['function'])
    if 'mean_zero' in request and type(request['mean_zero']) is not bool:
        raise ValueError('mean_zero must be boolean')
    if 'tolerance' in request:
        _tolerance(request['tolerance'])
    budgets = ({'modes': (1, 32), 'bits': (8, 96), 'max_m': (0, 12),
                'sqrt_bits': (8, 160), 'near_tail': (0, 32)} if op == 'spectrum' else
               {'levels': (1, 64), 'degree': (0, 32), 'sqrt_bits': (0, 256)} if op == 'approximate' else {})
    for key, (low, high) in budgets.items():
        if key in request:
            backend.direct.integer(request[key], low, high, key)
    if op == 'approximate' and 'root_ratio' in request:
        ratio = backend.pieces.rational(request['root_ratio'])
        if not 0 < ratio < 1 or ratio.denominator > 4096:
            raise ValueError('root_ratio must lie in (0,1) with denominator<=4096')
    if 'request_id' in request:
        identifier = request['request_id']
        if type(identifier) is not str or not 1 <= len(identifier) <= 128 or identifier.strip() != identifier or any(ord(c) < 32 for c in identifier):
            raise ValueError('request_id must be a nonempty identifier of at most 128 characters')
    if 'expected_request_digest' in request:
        value = request['expected_request_digest']
        if type(value) is not str or re.fullmatch('[0-9a-f]{64}', value) is None or value != request_digest(request):
            raise ValueError('Request digest mismatch')
    return request


def execute(request, executor=None):
    request = deepcopy(validate_request(request))
    digest = request_digest(request)
    payload = {key: deepcopy(value) for key, value in request.items() if key not in METADATA}
    result = dict((executor or backend.legacy.handle)(payload))
    if 'expected_kind' in request:
        cert = request.get('certificate') if request.get('op') == 'verify' else result.get('certificate')
        result['kind_matches'] = certificate_kind(cert) == request['expected_kind']
        if request.get('op') == 'verify' and not result['kind_matches']:
            result['verified'] = False
    result['request_digest'] = digest
    if 'request_id' in request:
        result['request_id'] = request['request_id']
    return result


def process_lines(lines, executor=None):
    for line in lines:
        if line.strip():
            yield execute(parse_request(line), executor=executor)
