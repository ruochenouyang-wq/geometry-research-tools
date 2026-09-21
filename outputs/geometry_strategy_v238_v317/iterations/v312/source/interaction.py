"""Evidence-preserving interaction helpers; no model or network calls.

Executors and verifiers are supplied by the caller. A summary is not a proof;
its evidence_state records whether an explicit verifier actually accepted it.
"""
from copy import deepcopy
from fractions import Fraction
import hashlib
import json


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
        squared = exact_fraction(certificate['error_squared'])
        if error < 0 or squared < 0 or error*error < squared:
            raise ValueError('Invalid approximation error bound')
        return {'error_upper': str(error), 'error_squared': str(squared),
                'decimal_error_upper': outward_decimal(str(error), places, upper=True)}
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


def summarize_result(result, request=None, verifier=None, places=12):
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
    checked = None if verifier is None else verifier(deepcopy(certificate)) is True
    summary = {
        'ok': checked is not False, 'format': 'interaction_summary_v1',
        'request': deepcopy(request),
        'target': {'tolerance': certificate.get('tolerance'),
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
    if checked is True and 'lower' in summary['bounds']:
        tolerance = exact_fraction(certificate['tolerance'])
        if tolerance <= 0:
            raise ValueError('Tolerance must be positive')
        summary['status'] = ('certified_met' if exact_fraction(summary['bounds']['exact_width'])
                             <= tolerance else 'certified_open')
    if checked is not True:
        summary['status'] = 'unchecked' if checked is None else 'verification_failed'
    return summary


def present_result(result, request=None, verifier=None, store=None,
                   include_certificate=False, places=12):
    """Store certificate once through store(sha256, certificate) -> locator.

    Without a store, include the full certificate so evidence is not lost.
    The store owns durability; a hash alone is not a mathematical verifier.
    """
    summary = summarize_result(result, request, verifier, places)
    certificate = result.get('certificate')
    if not isinstance(certificate, dict):
        return summary
    if store is not None:
        key = digest(certificate)
        locator = store(key, deepcopy(certificate))
        if not isinstance(locator, str) or not locator:
            raise ValueError('Store must return a nonempty evidence locator')
        summary['evidence_ref'] = {'sha256': key, 'locator': locator,
                                   'bytes': len(canonical(certificate).encode('utf-8'))}
    if store is None or include_certificate:
        summary['certificate'] = deepcopy(certificate)
    return summary


def retrieve_evidence(reference, load, verifier=None, require_verified=False):
    """load(locator) -> certificate; reject content substitution before replay."""
    if not isinstance(reference, dict):
        raise ValueError('Evidence reference must be an object')
    key = reference.get('sha256')
    if (not isinstance(key, str) or len(key) != 64 or
            any(character not in '0123456789abcdef' for character in key)):
        raise ValueError('Invalid SHA-256 evidence key')
    if not isinstance(reference.get('locator'), str) or not reference['locator']:
        raise ValueError('Missing evidence locator')
    certificate = load(reference['locator'])
    if not isinstance(certificate, dict) or digest(certificate) != key:
        raise ValueError('Evidence content does not match its hash')
    if 'bytes' in reference and reference['bytes'] != len(canonical(certificate).encode('utf-8')):
        raise ValueError('Evidence byte count mismatch')
    checked = None if verifier is None else verifier(deepcopy(certificate)) is True
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


class RequestRunner:
    """Deduplicate deterministic calls; disable cache for fresh/stochastic work.

    execute returns a transport envelope whose result is the executor output.
    Cached execution is not mathematical verification. New theorem/code versions
    require a new runner namespace. This small runner is sequential, not thread safe.
    """
    def __init__(self, executor, namespace, cache=True):
        if not callable(executor):
            raise ValueError('Executor must be callable')
        request_key({}, namespace)
        self.executor = executor
        self.namespace = namespace
        self.cache_enabled = bool(cache)
        self._cache = {}

    def execute(self, request):
        key = request_key(request, self.namespace)
        if self.cache_enabled and key in self._cache:
            return {'result': deepcopy(self._cache[key]), 'cache_hit': True, 'request_key': key}
        result = self.executor(deepcopy(request))
        if not isinstance(result, dict):
            raise ValueError('Executor must return an object')
        if self.cache_enabled and result.get('ok') is not False:
            self._cache[key] = deepcopy(result)
        return {'result': result, 'cache_hit': False, 'request_key': key}
