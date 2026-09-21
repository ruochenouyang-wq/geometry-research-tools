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


def summarize_result(result, request=None, verifier=None):
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
        'bounds': {key: deepcopy(certificate[key]) for key in (
            'lower', 'upper', 'exact_width', 'error_squared', 'error_upper')
                   if key in certificate},
    }
    if checked is not True:
        summary['status'] = 'unchecked' if checked is None else 'verification_failed'
    return summary
