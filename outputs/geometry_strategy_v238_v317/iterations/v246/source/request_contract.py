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
CERTIFICATE_KINDS = {backend.direct.FULL: 'spectrum', backend.enriched.FORMAT: 'spectrum',
                     backend.pieces.FORMAT: 'approximation'}
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


class ContractError(ValueError):
    """Stable machine-readable rejection; message is explanatory, not an API."""
    def __init__(self, code, message, field=None):
        super().__init__(message)
        self.code = code
        self.field = field


def _integer(value, low, high, name):
    try:
        return backend.direct.integer(value, low, high, name)
    except ValueError as error:
        raise ContractError('invalid_budget', str(error), name) from error


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
    return CERTIFICATE_KINDS.get(certificate.get('format'))


def verify_certificate(cert, expected_function=None, expected_mean_zero=None, expected_tolerance=None):
    request = {'op': 'verify', 'certificate': cert}
    for key, value in (('function', expected_function), ('mean_zero', expected_mean_zero),
                       ('tolerance', expected_tolerance)):
        if value is not None:
            request[key] = value
    try:
        return backend.legacy.handle(request).get('verified') is True
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


def evidence_state(request, certificate, verifier=None):
    """Replay evidence; a solver's ok/verified/status fields are never proof."""
    verifier = verifier or verify_certificate
    op = request['op']
    generated = op != 'verify'
    expected_kind = request.get('expected_kind',
                                {'spectrum': 'spectrum', 'precise_singular': 'spectrum',
                                 'approximate': 'approximation'}.get(op))
    kind = certificate_kind(certificate)
    function = request.get('function')
    mean_zero = request.get('mean_zero', True if op in ('spectrum', 'precise_singular') else None)
    tolerance = request.get('tolerance', '1/100000000' if op in ('spectrum', 'precise_singular') else None)
    checked = verifier(certificate, expected_function=function, expected_mean_zero=mean_zero,
                       expected_tolerance=tolerance) is True
    has_expected = any(value is not None for value in (function, mean_zero, tolerance))
    valid = checked or (has_expected and verifier(certificate) is True)
    kind_matches = expected_kind is None or kind == expected_kind
    bindings_match = checked and kind_matches
    provided = {'function': function is not None, 'kind': expected_kind is not None,
                'mean_zero': mean_zero is not None, 'tolerance': tolerance is not None}
    complete = provided['function'] and provided['kind'] and (
        kind == 'approximation' or (kind == 'spectrum' and provided['mean_zero'] and provided['tolerance']))
    target_met = None
    certificate_target_met = None
    if valid and kind == 'spectrum':
        width = backend.F(certificate['upper']) - backend.F(certificate['lower'])
        if width < 0:
            raise ValueError('Verified certificate has reversed endpoints')
        certificate_target_met = width <= backend.F(certificate['tolerance'])
        if tolerance is not None:
            target_met = bindings_match and width <= _tolerance(tolerance)
    if not valid:
        status = 'invalid_certificate'
    elif not bindings_match:
        status = 'binding_mismatch'
    elif not complete:
        status = 'verified_unbound'
    elif target_met is True:
        status = 'target_met'
    elif target_met is False:
        status = 'certified_open'
    else:
        status = 'certified_no_target'
    if not isinstance(certificate, dict):
        reason = 'missing_certificate' if certificate is None else 'malformed_certificate'
    elif kind is None:
        reason = 'unsupported_certificate_format'
    elif not valid:
        reason = 'invalid_evidence'
    elif not kind_matches:
        reason = 'conclusion_kind_mismatch'
    elif not bindings_match:
        reason = 'binding_mismatch'
    elif not complete:
        reason = 'unbound_verification'
    elif target_met is False:
        reason = 'target_not_met'
    elif target_met is None:
        reason = 'no_target_requested'
    else:
        reason = None
    return {'certificate_valid': valid, 'verified': valid and bindings_match,
            'bindings_match': bindings_match, 'binding_complete': bool(complete),
            'bindings_provided': provided, 'kind_matches': kind_matches,
            'certificate_kind': kind, 'target_met': target_met,
            'certificate_target_met': certificate_target_met, 'status': status,
            'reason_code': reason}


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ContractError('duplicate_key', 'Duplicate JSON key: ' + key, key)
        result[key] = value
    return result


def parse_request(text):
    try:
        parsed = json.loads(text, object_pairs_hook=_unique_object,
                            parse_constant=_reject_constant)
    except json.JSONDecodeError as error:
        raise ContractError('malformed_json', str(error)) from error
    return validate_request(parsed)


def _reject_constant(value):
    raise ContractError('nonfinite_number', 'Non-finite JSON number: ' + value)


def _finite_values(value):
    if isinstance(value, float) and not math.isfinite(value):
        raise ContractError('nonfinite_number', 'Non-finite request number')
    if isinstance(value, dict):
        for item in value.values():
            _finite_values(item)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _finite_values(item)


def _tolerance(value):
    try:
        value = backend.pieces.rational(value)
    except ValueError as error:
        raise ContractError('invalid_rational', str(error), 'tolerance') from error
    if not backend.F(1, 10**30) <= value <= 1:
        raise ContractError('invalid_budget', 'Tolerance must be in 1e-30..1', 'tolerance')
    return value


def _validate_singular(request):
    if 'function' not in request or backend.direct.normalize(request['function']) != backend.enriched.FUNCTION:
        raise ValueError('precise_singular supports the fixed original singular function only')
    for name, default, low, high in (('max_terms', 16, 1, 16),
                                    ('precision_bits', 112, 32, 256),
                                    ('iterations', 14, 0, 64)):
        _integer(request.get(name, default), low, high, name)
    _tolerance(request.get('tolerance', '1/100000000'))


def validate_request(request):
    try:
        return _validate_request(request)
    except ContractError:
        raise
    except (ValueError, TypeError, KeyError, ArithmeticError) as error:
        raise ContractError('invalid_request', str(error)) from error


def _validate_request(request):
    if not isinstance(request, dict):
        raise ContractError('invalid_root', 'Request must be an object')
    _finite_values(request)
    op = request.get('op')
    if type(op) is not str or op not in FIELDS:
        raise ContractError('unknown_operation', 'Unknown operation', 'op')
    if not set(request) <= FIELDS[op] | METADATA:
        raise ContractError('unknown_field', 'Unknown request field')
    for key in ('function', 'mean_zero', 'tolerance'):
        if key in request and request[key] is None:
            raise ContractError('null_binding', 'Explicit null binding: ' + key, key)
    if 'expected_kind' in request:
        if type(request['expected_kind']) is not str or request['expected_kind'] not in KINDS:
            raise ContractError('invalid_kind', 'expected_kind must be spectrum or approximation', 'expected_kind')
        produced = {'spectrum': 'spectrum', 'precise_singular': 'spectrum',
                    'approximate': 'approximation'}.get(request.get('op'))
        if produced is not None and produced != request['expected_kind']:
            raise ContractError('conclusion_kind_mismatch', 'Requested operation and expected_kind differ', 'expected_kind')
    if request.get('op') == 'precise_singular':
        _validate_singular(request)
    elif op in ('spectrum', 'approximate'):
        if 'function' not in request:
            raise ContractError('missing_function', 'Missing original function', 'function')
        normalizer = backend.pieces.normalize_function if op == 'approximate' else backend.direct.normalize
        normalizer(request['function'])
    elif 'function' in request:
        backend.direct.normalize(request['function'])
    if 'mean_zero' in request and type(request['mean_zero']) is not bool:
        raise ContractError('invalid_type', 'mean_zero must be boolean', 'mean_zero')
    if 'tolerance' in request:
        _tolerance(request['tolerance'])
    budgets = ({'modes': (1, 32), 'bits': (8, 96), 'max_m': (0, 12),
                'sqrt_bits': (8, 160), 'near_tail': (0, 32)} if op == 'spectrum' else
               {'levels': (1, 64), 'degree': (0, 32), 'sqrt_bits': (0, 256)} if op == 'approximate' else {})
    for key, (low, high) in budgets.items():
        if key in request:
            _integer(request[key], low, high, key)
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
            raise ContractError('request_digest_mismatch', 'Request digest mismatch', 'expected_request_digest')
    return request


def execute(request, executor=None, verifier=None):
    request = deepcopy(validate_request(request))
    digest = request_digest(request)
    payload = {key: deepcopy(value) for key, value in request.items() if key not in METADATA}
    result = {'ok': True} if request['op'] == 'verify' else dict((executor or backend.legacy.handle)(payload))
    cert = request.get('certificate') if request['op'] == 'verify' else result.get('certificate')
    result['execution_ok'] = result.get('ok') is True
    result.update(evidence_state(request, cert, verifier=verifier))
    result['request_digest'] = digest
    if 'request_id' in request:
        result['request_id'] = request['request_id']
    return result


def process_lines(lines, executor=None, verifier=None):
    for line in lines:
        if line.strip():
            yield execute(parse_request(line), executor=executor, verifier=verifier)
