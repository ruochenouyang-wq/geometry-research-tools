"""Request boundary for the frozen mathematical backend.

This module does not change mathematical algorithms or their certificates.
"""
import json
import math
import backend

METADATA = {'expected_kind'}
KINDS = {'spectrum', 'approximation'}


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


def validate_request(request):
    if not isinstance(request, dict):
        raise ValueError('Request must be an object')
    _finite_values(request)
    for key in ('function', 'mean_zero', 'tolerance'):
        if key in request and request[key] is None:
            raise ValueError('Explicit null binding: ' + key)
    if 'expected_kind' in request:
        if request['expected_kind'] not in KINDS:
            raise ValueError('expected_kind must be spectrum or approximation')
        produced = {'spectrum': 'spectrum', 'precise_singular': 'spectrum',
                    'approximate': 'approximation'}.get(request.get('op'))
        if produced is not None and produced != request['expected_kind']:
            raise ValueError('Requested operation and expected_kind differ')
    return request


def execute(request, executor=None):
    request = validate_request(request)
    payload = {key: value for key, value in request.items() if key not in METADATA}
    result = dict((executor or backend.legacy.handle)(payload))
    if 'expected_kind' in request:
        cert = request.get('certificate') if request.get('op') == 'verify' else result.get('certificate')
        result['kind_matches'] = certificate_kind(cert) == request['expected_kind']
        if request.get('op') == 'verify' and not result['kind_matches']:
            result['verified'] = False
    return result


def process_lines(lines, executor=None):
    for line in lines:
        if line.strip():
            yield execute(parse_request(line), executor=executor)
