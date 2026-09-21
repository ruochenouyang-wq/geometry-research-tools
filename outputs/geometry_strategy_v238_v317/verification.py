"""One explicit registry for replay and original-task binding.

This integration layer does not generate candidates. A valid certificate may
still miss the requested precision; replay and goal assessment stay separate.
"""
from copy import deepcopy
import importlib
import backend

F = backend.F


def conclusion_kind(certificate):
    if type(certificate) is not dict:
        return None
    fmt=certificate.get('format')
    if fmt == backend.pieces.FORMAT:
        return 'approximation'
    if fmt in (backend.direct.FULL, backend.enriched.FORMAT,
               'adaptive_singular_temple_v1'):
        return 'spectrum'
    if fmt == 'geometry_proof_transport_v1' and certificate.get('kind') in ('spectrum','approximation'):
        return certificate['kind']
    return None


def verify_any(certificate, expected_function=None, expected_mean_zero=None,
               expected_tolerance=None, expected_kind=None, _depth=0):
    """Replay exact evidence; tolerance binds an embedded target if present.

Certificates of a mathematical bound need not themselves embed a precision
target. Such certificates may be reused for any target, with goal status
computed separately by assess. Explicit nulls are rejected by request parsing.
"""
    try:
        if _depth > 48 or type(certificate) is not dict:
            return False
        kind=conclusion_kind(certificate)
        if kind is None or (expected_kind is not None and kind != expected_kind):
            return False
        if expected_function is not None:
            if backend.direct.normalize(expected_function) != certificate.get('function'):
                return False
        if expected_mean_zero is not None:
            if type(expected_mean_zero) is not bool or kind != 'spectrum':
                return False
            if certificate.get('mean_zero') is not expected_mean_zero:
                return False
        if expected_tolerance is not None:
            target=backend.pieces.rational(expected_tolerance)
            if not F(1,10**30) <= target <= 1:
                return False
            if 'tolerance' in certificate and backend.pieces.rational(certificate['tolerance']) != target:
                return False
        fmt=certificate.get('format')
        if fmt == backend.pieces.FORMAT:
            return backend.pieces.verify_piecewise(certificate, expected_function)
        if fmt == backend.direct.FULL:
            return backend.direct.verify_full(certificate, expected_function,
                                             expected_mean_zero, expected_tolerance)
        if fmt == backend.enriched.FORMAT:
            return backend.enriched.verify(certificate, expected_tolerance=expected_tolerance)
        if fmt == 'adaptive_singular_temple_v1':
            module=importlib.import_module('adaptive_singular')
            return module.verify(certificate, expected_function=expected_function,
                                 expected_mean_zero=expected_mean_zero,
                                 expected_tolerance=expected_tolerance) is True
        if fmt == 'geometry_proof_transport_v1':
            module=importlib.import_module('proof_transport')
            def child(cert, **bindings):
                return verify_any(cert, _depth=_depth+1, **bindings)
            return module.verify(certificate, verifier=child) is True
        return False
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError,RecursionError):
        return False


def assess(certificate, task):
    """Separate evidence validity from reaching a particular task's target."""
    if type(task) is not dict or task.get('kind') not in ('spectrum','approximation'):
        raise ValueError('Task requires an explicit mathematical conclusion kind')
    target=backend.pieces.rational(task.get('tolerance','1/100000000'))
    if not F(1,10**30) <= target <= 1:
        raise ValueError('Tolerance outside 1e-30..1')
    scope=task.get('mean_zero',True) if task['kind']=='spectrum' else None
    if not verify_any(certificate, expected_function=task['function'],
                      expected_mean_zero=scope, expected_tolerance=str(target),
                      expected_kind=task['kind']):
        return {'certificate_valid':False,'target_met':False,'status':'verification_failed'}
    value=(F(certificate['upper'])-F(certificate['lower']) if task['kind']=='spectrum'
           else F(certificate['error_upper']))
    return {'certificate_valid':True,'target_met':value<=target,
            'status':'target_met' if value<=target else 'certified_open',
            'quantity':'spectral_interval_width' if task['kind']=='spectrum' else 'function_L2_error_upper',
            'value':str(value),'tolerance':str(target),'kind':task['kind'],
            'function':deepcopy(certificate['function'])}


def verify_request(request):
    """Replay a supplied certificate, distinguishing intrinsic and bound validity."""
    import request_contract
    request_contract.validate_request(request)
    certificate=request.get('certificate')
    kind=conclusion_kind(certificate)
    valid=verify_any(certificate)
    bound=valid and verify_any(certificate,
        expected_function=request.get('function'),expected_mean_zero=request.get('mean_zero'),
        expected_tolerance=request.get('tolerance'),expected_kind=request.get('expected_kind'))
    complete=(request.get('function') is not None and request.get('expected_kind') is not None
              and (kind=='approximation' or
                   (kind=='spectrum' and 'mean_zero' in request and 'tolerance' in request)))
    target_met=None
    if bound and 'tolerance' in request:
        quantity=(F(certificate['error_upper']) if kind=='approximation' else
                  F(certificate['upper'])-F(certificate['lower']))
        target_met=quantity<=F(request['tolerance'])
    status=('invalid_certificate' if not valid else 'binding_mismatch' if not bound else
            'verified_unbound' if not complete else 'target_met' if target_met is True else
            'certified_open' if target_met is False else 'certified_no_target')
    result={'ok':True,'execution_ok':True,'certificate_valid':valid,'verified':bound,
            'bindings_match':bound,'binding_complete':bool(complete),'certificate_kind':kind,
            'target_met':target_met,'status':status,
            'request_digest':request_contract.request_digest(request)}
    if 'request_id' in request:
        result['request_id']=request['request_id']
    return result
