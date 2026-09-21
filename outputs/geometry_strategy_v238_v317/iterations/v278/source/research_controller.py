"""Evidence-bound research decisions for two distinct numerical objectives."""
from copy import deepcopy
from time import perf_counter
import backend

F = backend.F


def normalize(request):
    if type(request) is not dict:
        raise ValueError('Request must be an object')
    if request.get('target') not in ('approximation', 'spectrum'):
        raise ValueError('Choose approximation or spectrum explicitly')
    q = backend.direct.normalize(request['function'])
    tolerance = backend.pieces.rational(request.get('tolerance', '1/100000000'))
    if not F(1, 10**30) <= tolerance <= 1:
        raise ValueError('Tolerance must be in 1e-30..1')
    mean_zero = request.get('mean_zero', True)
    if type(mean_zero) is not bool:
        raise ValueError('mean_zero must be boolean')
    return dict(request, function=q, tolerance=str(tolerance), mean_zero=mean_zero)


def default_solvers():
    return {name: backend.legacy.handle for name in ('approximate', 'singular', 'spectrum')}


def verify_certificate(cert, expected_function=None, expected_mean_zero=None,
                       expected_tolerance=None):
    fmt = cert.get('format') if isinstance(cert, dict) else None
    if fmt == backend.pieces.FORMAT:
        return (expected_mean_zero is None and expected_tolerance is None
                and backend.pieces.verify_piecewise(cert, expected_function))
    if fmt == backend.direct.FULL:
        return backend.direct.verify_full(cert, expected_function,
                                          expected_mean_zero, expected_tolerance)
    if fmt == backend.enriched.FORMAT:
        return (expected_function == backend.enriched.FUNCTION
                and expected_mean_zero is True
                and backend.enriched.verify(cert, expected_tolerance=expected_tolerance))
    return False


def _accept(cert, request, verifier):
    if request['target'] == 'approximation':
        if 'error_upper' not in cert or cert.get('spectral_transfer_claimed') is not False:
            raise ValueError('Expected a function approximation certificate')
        valid = verifier(cert, expected_function=request['function'],
                         expected_mean_zero=None, expected_tolerance=None)
        score = F(cert['error_upper'])
    else:
        if not all(k in cert for k in ('lower', 'upper', 'exact_width')):
            raise ValueError('Expected a complete spectral interval')
        valid = verifier(cert, expected_function=request['function'],
                         expected_mean_zero=request['mean_zero'],
                         expected_tolerance=request['tolerance'])
        score = F(cert['upper'])-F(cert['lower'])
        valid = valid and score == F(cert['exact_width'])
    if not valid or score < 0:
        raise ValueError('Certificate replay or original problem binding failed')
    return score


def first_action(request):
    if request['target'] == 'approximation':
        return {'route': 'approximate', 'op': 'approximate',
                'function': request['function'], 'levels': 2, 'degree': 2}
    return {'route': 'spectrum', 'op': 'spectrum', 'function': request['function'],
            'mean_zero': request['mean_zero'], 'tolerance': request['tolerance'],
            'modes': 2, 'bits': 20, 'max_m': 2, 'near_tail': 0}


def diagnose(cert, request):
    return {'kind': 'goal_gap', 'evidence': 'verified_certificate',
            'target': request['target'], 'hypothesis': None}


def solve(request, solvers=None, verifier=None):
    started = perf_counter()
    request = normalize(request)
    solvers = default_solvers() if solvers is None else solvers
    verifier = verify_certificate if verifier is None else verifier
    action = first_action(request)
    raw = solvers[action['route']]({k: deepcopy(v) for k, v in action.items() if k != 'route'})
    cert = raw.get('certificate', raw)
    score = _accept(cert, request, verifier)
    return {'certificate': cert, 'attempts': [{'action': action, 'metric': str(score),
             'verified': True}], 'status': 'target_met' if score <= F(request['tolerance'])
             else 'certified_open', 'diagnosis': diagnose(cert, request),
             'cost': {'elapsed_seconds': perf_counter()-started, 'solver_calls': 1}}
