"""Private copy of the frozen controller with two diagnosed integration repairs.

No source or globals in the frozen controller are modified.  This fallback is
kept separate from the new mathematical routes so its effect can be inspected.
"""
from copy import deepcopy
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'geometry_strategy_v238_v317'
sys.dont_write_bytecode = True
if str(OLD) not in sys.path:
    sys.path.append(str(OLD))
import backend
import research_service
import verification

if Path(backend.__file__).resolve() != OLD / 'backend.py':
    raise RuntimeError('The fallback must use the frozen v317 backend')

_spec = importlib.util.spec_from_file_location(
    '_accuracy80_private_frozen_controller', OLD / 'research_controller.py')
_controller = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_controller)
_original_diagnose = _controller.diagnose
_original_next_action = _controller.next_action


def diagnose(certificate, request):
    if certificate.get('format') != 'adaptive_singular_temple_v1':
        return _original_diagnose(certificate, request)
    import adaptive_singular
    split = adaptive_singular.residual_diagnosis(
        certificate['function'], certificate['powers'], certificate['coefficients'])
    inside = backend.F(split['inside_squared'])
    outside = backend.F(split['outside_squared'])
    return {
        'kind': 'trial_residual_budget',
        'evidence': 'exact_orthogonal_decomposition',
        'projected_squared': str(inside), 'outside_trial_squared': str(outside),
        'complete_squared': split['complete_squared'],
        'bottleneck': 'coefficients' if inside > outside else 'representation',
        'hypothesis': 'Residual components guide refinement; they do not prove the spectral gap.',
    }


def next_action(action, diagnosis):
    if action.get('route') == 'singular' and diagnosis.get('kind') == 'spectral_lower_budget':
        # The singular backend can legitimately return an ordinary source
        # certificate. Refine the requested singular route using its own fields.
        # The controller's independent spectrum continuation remains available.
        counts = [n for n in (1, 3, 6, 10, 16) if n > action.get('max_terms', 3)]
        if not counts:
            return None
        result = deepcopy(action)
        result['max_terms'] = counts[0]
        return result
    return _original_next_action(action, diagnosis)


_controller.diagnose = diagnose
_controller.next_action = next_action


def solve(task):
    request = {k: deepcopy(v) for k, v in task.items() if k not in ('kind', 'id')}
    request['target'] = task['kind']
    if task['kind'] == 'approximation':
        import adaptive_approx
        return adaptive_approx.solve(request['function'], tolerance=request['tolerance'])
    service = research_service.ResearchService
    solvers = {'spectrum': service._spectrum, 'singular': service._singular,
               'approximate': backend.legacy.handle,
               'singular_supports': service._singular_supports}
    return _controller.solve(request, solvers=solvers, verifier=verification.verify_any)


def legal_actions(request, actions):
    return _controller.legal_actions(request, actions,
                                     research_service.ResearchService._singular_supports)
