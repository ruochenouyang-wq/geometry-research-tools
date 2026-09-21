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


def _singular_supported(request, supports=None):
    if request['target'] != 'spectrum' or not request['mean_zero']:
        return False
    return (bool(supports(deepcopy(request))) if supports is not None
            else request['function'] == backend.enriched.FUNCTION)


def first_action(request, singular_supports=None):
    if request['target'] == 'approximation':
        return {'route': 'approximate', 'op': 'approximate',
                'function': request['function'],
                'levels': 1 if request['function']['profile'] == 'step' else 2,
                'degree': 0 if request['function']['profile'] == 'step' else 2}
    route = request.get('route', 'auto')
    if route not in ('auto', 'singular', 'spectrum'):
        raise ValueError('Unknown spectral route')
    if route == 'singular' or (route == 'auto' and _singular_supported(request, singular_supports)):
        return {'route': 'singular', 'op': 'precise_singular', 'function': request['function'],
                'tolerance': request['tolerance'], 'max_terms': 3,
                'precision_bits': 112, 'iterations': 8}
    return {'route': 'spectrum', 'op': 'spectrum', 'function': request['function'],
            'mean_zero': request['mean_zero'], 'tolerance': request['tolerance'],
            'modes': 2, 'bits': 20, 'max_m': 2, 'near_tail': 0}


def legal_actions(request, actions, singular_supports=None):
    """Filter unsupported mathematical operations before executing any solver."""
    request = normalize(request)
    accepted, rejected = [], []
    integer_limits = {'levels': (1, 64), 'degree': (0, 32), 'sqrt_bits': (8, 160),
                      'modes': (1, 32), 'bits': (8, 96), 'max_m': (0, 12),
                      'near_tail': (0, 32), 'max_terms': (1, 16),
                      'precision_bits': (32, 256), 'iterations': (0, 64)}
    fields = {'approximate': {'levels', 'degree', 'sqrt_bits', 'root_ratio'},
              'singular': {'max_terms', 'precision_bits', 'iterations', 'tolerance'},
              'spectrum': {'mean_zero', 'tolerance', 'modes', 'bits', 'max_m',
                           'near_tail', 'sqrt_bits'}}
    for action in actions:
        try:
            route = action['route']
            if route not in fields:
                raise ValueError('Unknown route')
            if not set(action) <= fields[route] | {'route', 'op', 'function'}:
                raise ValueError('Unknown action field')
            if backend.direct.normalize(action['function']) != request['function']:
                raise ValueError('Action changes the original function')
            if (route == 'approximate') != (request['target'] == 'approximation'):
                raise ValueError('Action changes the target kind')
            if action['op'] != {'approximate': 'approximate', 'singular': 'precise_singular',
                                 'spectrum': 'spectrum'}[route]:
                raise ValueError('Route and operation disagree')
            if route == 'approximate':
                backend.pieces.normalize_function(request['function'])
                if not 0 < F(action.get('root_ratio', '1/2')) < 1:
                    raise ValueError('Root ratio must lie in (0,1)')
            else:
                if F(action['tolerance']) != F(request['tolerance']):
                    raise ValueError('Action changes target tolerance')
                if route == 'singular':
                    if not _singular_supported(request, singular_supports):
                        raise ValueError('Enriched route supports its fixed mean-zero problem only')
                elif action['mean_zero'] is not request['mean_zero']:
                    raise ValueError('Action changes the spectral space')
                backend.direct.form_bound(request['function'])
            for key, (lo, hi) in integer_limits.items():
                if key in action:
                    backend.direct.integer(action[key], lo, hi, key)
            accepted.append(deepcopy(action))
        except (ValueError, TypeError, KeyError, ArithmeticError) as exc:
            rejected.append({'action': deepcopy(action), 'reason': str(exc)})
    return accepted, rejected


def diagnose(cert, request):
    """Inspect an already replayed certificate; arithmetic facts vs hypotheses."""
    if request['target'] == 'approximation':
        total = F(cert['error_squared'])
        core = sum((F(c['error_squared_dt']) for c in cert['cells'] if c.get('core')), F(0))
        annuli = total-core
        if annuli < 0:
            raise ValueError('Approximation error partition is inconsistent')
        rounding = F(cert['error_upper'])**2-total
        target = F(request['tolerance'])**2
        bottleneck = 'rounding' if total <= target < total+rounding else (
            'core' if core >= annuli else 'annuli')
        return {'kind': 'approximation_budget', 'evidence': 'exact_error_partition',
                'core_squared': str(core), 'annuli_squared': str(annuli),
                'rounding_squared_slack': str(rounding), 'target_squared': str(target),
                'bottleneck': bottleneck, 'hypothesis': 'Refine the largest controllable contribution'}
    if cert.get('format') == backend.enriched.FORMAT:
        M = [[F(x) for x in row] for row in cert['matrices']['M']]
        H = [[F(x) for x in row] for row in cert['matrices']['H']]
        v = list(map(F, cert['trial']['coefficients']))
        mu = F(cert['statistics']['rayleigh'])
        r = [sum(((H[i][j]-mu*M[i][j])*v[j] for j in range(len(v))), F(0))
             for i in range(len(v))]
        dual = backend.enriched._ldl_solver(M)(r)
        projected = sum((a*b for a, b in zip(r, dual)), F(0))/F(cert['statistics']['mass'])
        complete = F(cert['statistics']['residual_squared'])
        outside = complete-projected
        if projected < 0 or outside < 0:
            raise ArithmeticError('Orthogonal residual decomposition failed')
        return {'kind': 'trial_residual_budget', 'evidence': 'exact_orthogonal_decomposition',
                'projected_squared': str(projected), 'outside_trial_squared': str(outside),
                'complete_squared': str(complete),
                'bottleneck': 'coefficients' if projected > outside else 'representation',
                'hypothesis': 'A larger residual component suggests a corresponding refinement; '
                              'it does not prove spectral interval improvement'}
    if cert.get('format') == backend.direct.FULL:
        sector_low = min(F(s['lower']) for s in cert['sectors'])
        angular = F(cert['angular_tail_lower'])
        limiting = [s['azimuth_m'] for s in cert['sectors']
                    if F(s['lower']) == F(cert['lower'])]
        return {'kind': 'spectral_lower_budget', 'evidence': 'exact_minimum_of_proven_bounds',
                'bottleneck': 'angular_tail' if angular <= sector_low else 'radial_sector',
                'limiting_m': limiting,
                'best_upper_m': min(cert['sectors'], key=lambda s: F(s['upper']))['azimuth_m'],
                'angular_tail_lower': str(angular), 'sector_min_lower': str(sector_low),
                'hypothesis': 'Radial sector width mixes representation, tail inequality and '
                              'endpoint resolution; these causes are not strictly separated'}
    return {'kind': 'goal_gap', 'evidence': 'verified_certificate',
            'target': request['target'], 'hypothesis': None}


def next_action(action, diagnosis):
    action = deepcopy(action)
    if diagnosis['kind'] == 'approximation_budget':
        key, increment = {'core': ('levels', 2), 'annuli': ('degree', 2),
                          'rounding': ('sqrt_bits', 16)}[diagnosis['bottleneck']]
        action[key] = action.get(key, 40)+increment
        return action
    if diagnosis['kind'] == 'trial_residual_budget':
        if diagnosis['bottleneck'] == 'coefficients':
            action['iterations'] = min(64, max(2, 2*action.get('iterations', 14)))
        else:
            counts = [n for n in (1, 3, 6, 10, 16) if n > action.get('max_terms', 3)]
            if not counts:
                return None
            action['max_terms'] = counts[0]
        return action
    if diagnosis['kind'] == 'spectral_lower_budget':
        if diagnosis['bottleneck'] == 'angular_tail':
            action['max_m'] = action.get('max_m', 2)+1
        elif action.get('near_tail', 0) < action.get('modes', 2)*2:
            action['near_tail'] = action.get('near_tail', 0)+4
        else:
            action['modes'] = action.get('modes', 2)+2
        return action
    return None


def solve(request, solvers=None, verifier=None):
    started = perf_counter()
    request = normalize(request)
    solvers = default_solvers() if solvers is None else solvers
    verifier = verify_certificate if verifier is None else verifier
    budget = request.get('budget', {})
    max_attempts = backend.direct.integer(budget.get('max_attempts', 8), 1, 1024, 'max_attempts')
    max_work = backend.direct.integer(budget.get('work_units', 2000000), 1, 10**12, 'work_units')
    supports = solvers.get('singular_supports')
    action = first_action(request, supports)
    attempts, cert, diagnosis, work = [], None, {'kind': 'not_started'}, 0
    reason = 'attempt_budget'
    for _ in range(max_attempts):
        legal, rejected = legal_actions(request, [action], supports)
        if not legal:
            diagnosis = {'kind': 'no_legal_action', 'rejected': rejected}
            reason = 'unsupported'
            break
        action = legal[0]
        charge = work_estimate(action)
        if work+charge > max_work:
            reason = 'work_budget'
            break
        work += charge
        raw = solvers[action['route']]({k: deepcopy(v) for k, v in action.items() if k != 'route'})
        cert = raw.get('certificate', raw)
        score = _accept(cert, request, verifier)
        diagnosis = diagnose(cert, request)
        attempts.append({'action': action, 'metric': str(score), 'verified': True,
                         'work_units_charged': charge, 'diagnosis': diagnosis})
        if score <= F(request['tolerance']):
            reason = 'target_met'
            break
        action = next_action(action, diagnosis)
        if action is None:
            reason = 'no_supported_refinement'
            break
    return {'certificate': cert, 'attempts': attempts,
            'status': 'target_met' if reason == 'target_met' else (
                'certified_open' if cert is not None else reason),
            'stop_reason': reason, 'diagnosis': diagnosis,
            'cost': {'elapsed_seconds': perf_counter()-started, 'solver_calls': len(attempts),
                     'work_units_charged': work, 'work_is_proxy_not_runtime': True}}


def work_estimate(action):
    """Deterministic scheduling proxy; never equated with seconds or token use."""
    if action['route'] == 'approximate':
        return action.get('levels', 2)*(action.get('degree', 2)+1)**3
    if action['route'] == 'singular':
        return 4000+sum(n**3 for n in (1, 3, 6, 10, 16)
                        if n <= action.get('max_terms', 3))*max(1, action.get('iterations', 8))
    n = action.get('modes', 2)
    return n**3*(n+action.get('near_tail', 0))*(action.get('max_m', 2)+1)*action.get('bits', 20)
