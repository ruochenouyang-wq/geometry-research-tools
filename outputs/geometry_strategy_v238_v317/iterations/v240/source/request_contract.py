"""Request boundary for the frozen mathematical backend.

This module does not change mathematical algorithms or their certificates.
"""
import json
import math
import backend


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
    return request


def execute(request, executor=None):
    request = validate_request(request)
    return (executor or backend.legacy.handle)(request)


def process_lines(lines, executor=None):
    for line in lines:
        if line.strip():
            yield execute(parse_request(line), executor=executor)
