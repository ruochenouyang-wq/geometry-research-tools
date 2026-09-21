"""Task-bound entry point for the new certified geometry routes.

The returned success is derived from an independently replayed certificate.
External evaluation supplies the hard deadline; this module also reports its
elapsed time and never treats a deadline as evidence of mathematical success.
"""
from copy import deepcopy
from fractions import Fraction as F
import importlib
import json
from math import isfinite
from pathlib import Path
import re
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parent
VERSION = 'geometry-accuracy80-development'
sys.dont_write_bytecode = True


def rational(value):
    if type(value) is int:
        return F(value)
    if type(value) is str and re.fullmatch(r'-?\d+(?:/[1-9]\d*)?', value):
        return F(value)
    raise ValueError('Expected an exact integer or rational string')


def normalize_function(function):
    if type(function) is not dict or set(function)-{
            'kind', 'profile', 'axis', 'amplitude', 'offset', 'exponent'}:
        raise ValueError('Malformed original function')
    if function.get('kind') != 'axis_profile' or function.get('profile') not in ('step', 'abs_power'):
        raise ValueError('This tool currently supports axis-profile step and power functions')
    axis = function.get('axis', 2)
    if type(axis) is not int or not 0 <= axis <= 2:
        raise ValueError('axis must be 0, 1 or 2')
    q = {'kind': 'axis_profile', 'profile': function['profile'], 'axis': axis,
         'amplitude': str(rational(function.get('amplitude', 1))),
         'offset': str(rational(function.get('offset', 0)))}
    if q['profile'] == 'abs_power':
        q['exponent'] = str(rational(function['exponent']))
    elif 'exponent' in function:
        raise ValueError('A step profile has no exponent')
    return q


def normalize_task(task):
    if type(task) is not dict or set(task)-{
            'kind', 'function', 'tolerance', 'mean_zero', 'id', 'budget'}:
        raise ValueError('Malformed task fields')
    kind = task.get('kind')
    if kind not in ('spectrum', 'approximation'):
        raise ValueError('Choose spectrum or approximation')
    q = normalize_function(task['function'])
    tolerance = rational(task.get('tolerance', '1/100000000'))
    if not F(1, 10**30) <= tolerance <= 1:
        raise ValueError('Tolerance must be between 1e-30 and 1')
    result = {'kind': kind, 'function': q, 'tolerance': str(tolerance)}
    if kind == 'spectrum':
        scope = task.get('mean_zero', True)
        if type(scope) is not bool:
            raise ValueError('mean_zero must be boolean')
        result['mean_zero'] = scope
    elif 'mean_zero' in task:
        raise ValueError('Function approximation has no mean-zero constraint')
    budget = task.get('budget', {})
    if type(budget) is not dict or set(budget)-{'wall_seconds'}:
        raise ValueError('The unified budget accepts only wall_seconds')
    wall = budget.get('wall_seconds', 10)
    if isinstance(wall, bool) or not isfinite(float(wall)) or not 0 <= float(wall) <= 86400:
        raise ValueError('Invalid wall_seconds budget')
    result['budget'] = {'wall_seconds': float(wall)}
    return result


def verify_certificate(certificate, *, expected_function, expected_mean_zero,
                       expected_tolerance, expected_kind):
    """Replay evidence and bind the original task; do not trust status fields."""
    try:
        if type(certificate) is not dict:
            return False
        q = normalize_function(expected_function)
        tol = rational(expected_tolerance)
        if certificate.get('function') != q or not F(1, 10**30) <= tol <= 1:
            return False
        if expected_kind not in ('approximation', 'spectrum'):
            return False
        if 'tolerance' in certificate and rational(certificate['tolerance']) != tol:
            return False
        if expected_kind == 'spectrum':
            if type(expected_mean_zero) is not bool or certificate.get('mean_zero') is not expected_mean_zero:
                return False
            lower, upper = rational(certificate['lower']), rational(certificate['upper'])
            if upper < lower or rational(certificate['exact_width']) != upper-lower:
                return False
        elif expected_mean_zero is not None:
            return False
        fmt = certificate.get('format')
        if fmt == 'general_rational_power_piecewise_L2_v1':
            if expected_kind != 'approximation':
                return False
            import general_approx
            return general_approx.verify(certificate, function=q, tolerance=str(tol)) is True
        if fmt in ('fullspace_m0_fractional_temple_v1', 'step_c1_strong_residual_temple_v1',
                   'singular_complete_gap_temple_v1'):
            if expected_kind != 'spectrum':
                return False
            module = importlib.import_module({
                'fullspace_m0_fractional_temple_v1': 'fullspace_singular',
                'step_c1_strong_residual_temple_v1': 'step_solver',
                'singular_complete_gap_temple_v1': 'singular_solver'}[fmt])
            return module.verify(certificate, expected_function=q,
                                 expected_mean_zero=expected_mean_zero,
                                 expected_tolerance=str(tol)) is True
        import compat_controller
        return compat_controller.verification.verify_any(
            certificate, expected_function=q, expected_mean_zero=expected_mean_zero,
            expected_tolerance=str(tol), expected_kind=expected_kind) is True
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError, RecursionError):
        return False


def assess(certificate, task):
    kind = task['kind']
    valid = verify_certificate(certificate, expected_function=task['function'],
                               expected_mean_zero=task.get('mean_zero') if kind == 'spectrum' else None,
                               expected_tolerance=task['tolerance'], expected_kind=kind)
    if not valid:
        return {'certificate_valid': False, 'target_met': False, 'status': 'verification_failed'}
    metric = (rational(certificate['exact_width']) if kind == 'spectrum'
              else rational(certificate['error_upper']))
    if metric < 0:
        return {'certificate_valid': False, 'target_met': False, 'status': 'verification_failed'}
    met = metric <= rational(task['tolerance'])
    return {'certificate_valid': True, 'target_met': met,
            'status': 'target_met' if met else 'certified_open',
            'metric': str(metric), 'tolerance': task['tolerance']}


def _route(task):
    q = task['function']
    if task['kind'] == 'approximation' and q['profile'] == 'abs_power':
        import general_approx
        return 'general_rational_power_piecewise', general_approx.solve(q, task['tolerance'])
    if task['kind'] == 'spectrum' and q['profile'] == 'step':
        import step_solver
        return 'step_c1_candidate', step_solver.solve(task, degree=10)
    if task['kind'] == 'spectrum' and q['profile'] == 'abs_power' and task['mean_zero'] is False:
        import fullspace_singular
        return 'fullspace_m0_fractional', fullspace_singular.solve(
            q, tolerance=task['tolerance'], mean_zero=False,
            wall_seconds=min(8.0, task['budget']['wall_seconds']))
    if task['kind'] == 'spectrum' and q['profile'] == 'abs_power' and task['mean_zero'] is True:
        import singular_solver
        return 'complete_gap_mean_zero_singular', singular_solver.solve(task)
    import compat_controller
    original = {k: v for k, v in task.items() if k != 'budget'}
    if task['kind'] == 'spectrum':
        original['budget'] = task['budget']
    return 'repaired_frozen_controller', compat_controller.solve(original)


def solve(task):
    started = perf_counter()
    task = normalize_task(task)
    limit = task['budget']['wall_seconds']
    if limit == 0:
        return {'version': VERSION, 'certificate': None, 'certificate_valid': False,
                'target_met': False, 'status': 'budget_exceeded', 'within_budget': False,
                'total_elapsed_seconds': perf_counter()-started, 'attempts': []}
    attempts = []
    try:
        route, raw = _route(task)
        certificate = raw.get('certificate', raw)
        # JSON is the public interchange format; replay what a caller receives.
        certificate = json.loads(json.dumps(certificate, allow_nan=False))
        assessment = (assess(certificate, task) if certificate is not None else
                      {'certificate_valid': False, 'target_met': False,
                       'status': 'no_verified_certificate',
                       'solver_status': raw.get('status'), 'stop_reason': raw.get('stop_reason')})
        attempts.append({'route': route, 'assessment': assessment,
                         'details': {k:v for k,v in raw.items() if k != 'certificate'}
                                    if 'certificate' in raw else None})
        result = {'version': VERSION, 'certificate': certificate, **assessment,
                  'attempts': attempts, 'route': route}
    except Exception as exc:
        result = {'version': VERSION, 'certificate': None, 'certificate_valid': False,
                  'target_met': False, 'status': 'execution_failed', 'attempts': attempts,
                  'error': type(exc).__name__+': '+str(exc)}
    elapsed = perf_counter()-started
    result.update(total_elapsed_seconds=elapsed, within_budget=elapsed <= limit,
                  budget_scope='whole_solve_call_cooperative; external evaluator enforces hard deadline')
    if elapsed > limit:
        result.update(status='budget_exceeded', target_met=False)
    return result


if __name__ == '__main__':
    for line in sys.stdin:
        if line.strip():
            try:
                print(json.dumps(solve(json.loads(line)), sort_keys=True, allow_nan=False), flush=True)
            except Exception as exc:
                print(json.dumps({'target_met': False, 'error': type(exc).__name__+': '+str(exc)}), flush=True)
