"""Request boundary for the frozen mathematical backend.

This module does not change mathematical algorithms or their certificates.
"""
import json
import backend


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('Duplicate JSON key: ' + key)
        result[key] = value
    return result


def parse_request(text):
    return validate_request(json.loads(text, object_pairs_hook=_unique_object))


def validate_request(request):
    return request


def execute(request, executor=None):
    request = validate_request(request)
    return (executor or backend.legacy.handle)(request)


def process_lines(lines, executor=None):
    for line in lines:
        if line.strip():
            yield execute(parse_request(line), executor=executor)
