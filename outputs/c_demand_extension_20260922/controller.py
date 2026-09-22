"""C2: precision-driven search, frozen certificates, explicit fallbacks."""
from fractions import Fraction as F
import json
from time import perf_counter
from runtime import baseline, c1, Deadline, Execution, verify

VERSION = 'c2-demand-extension-v1'


def centered_available(q):
    a, alpha = F(q['amplitude']), F(q['exponent'])
    variance = a*a*alpha*alpha / ((1+2*alpha)*(1+alpha)**2)
    return variance < 1


def solve(raw_task, progress=None):
    started = perf_counter()
    task = baseline.normalize_task(raw_task)
    execution = Execution(started, task['budget']['wall_seconds'], progress)
    certificate, fallback_reason, failure_stage = None, None, None
    unverified_candidates = []
    route, native_kind = 'native_c2', None
    assessment = {'certificate_valid': False, 'target_met': False,
                  'status': 'no_verified_certificate'}
    skip_fallback, fallback_assessment = False, None
    try:
        execution.check()
        q = task['function']
        if task['kind'] == 'approximation':
            fallback_reason = 'approximation_objective_not_extended'
            execution.note(fallback_reason)
        elif q['profile'] == 'step':
            import step_extension
            native_kind = 'adaptive_step_degree'
            execution.reserve_fallback()
            try:
                certificate = execution.run('adaptive_step_search', lambda:
                    step_extension.search(task, execution))
                # Keep any recovered open proof and reserve time for replay.
                # When the full schedule finished, this also avoids repeating
                # the baseline's degree10. An early stop makes no such claim.
                skip_fallback = certificate is not None
                if certificate is None:
                    fallback_reason = 'step_search_no_certificate'
            except (ValueError, ArithmeticError, Deadline) as error:
                fallback_reason = ('adaptive_stage_budget' if isinstance(error, Deadline)
                                   else 'step_search_precondition_or_candidate_failure')
                failure_stage = 'adaptive_step_search'
                execution.note(fallback_reason, error=type(error).__name__, reason=str(error))
            finally:
                execution.release_reserve()
        else:
            alpha, amplitude = F(q['exponent']), F(q['amplitude'])
            if not -F(1,2) < alpha < 0:
                # Strong residual moments used here cease to be integrable at
                # alpha=-1/2. This is a method boundary, not nonexistence proof.
                skip_fallback = True
                assessment['status'] = 'unsupported_method_domain'
                execution.note('strong_residual_exponent_domain', exponent=str(alpha),
                               required='-1/2 < exponent < 0')
                failure_stage = 'capability_check'
            elif not task['mean_zero'] and amplitude >= 0:
                skip_fallback = True
                assessment['status'] = 'unsupported_method_domain'
                execution.note('fullspace_requires_negative_amplitude', amplitude=str(amplitude))
                failure_stage = 'capability_check'
            elif not task['mean_zero'] and not centered_available(q):
                import fullspace_extension
                native_kind = 'certified_coarse_gap'
                extension = execution.run('strong_fullspace_search', lambda:
                    fullspace_extension.solve(task, wall_seconds=execution.remaining()*0.7,
                                              progress=execution.emit))
                certificate = extension.get('certificate')
                if extension.get('best_candidate') is not None:
                    unverified_candidates.append({
                        'origin': 'strong_fullspace_search',
                        'is_certificate': False,
                        'candidate': extension['best_candidate']})
                execution.attempts.extend(extension.get('attempts', []))
                for key, count in extension.get('counters', {}).items():
                    if isinstance(count, (int, float)):
                        execution.counters[key] = execution.counters.get(key, 0)+count
                execution.note('strong_fullspace_result', status_detail=extension.get('status'),
                               stop_reason=extension.get('stop_reason'))
                # Preserve any valid open result. Never label it target-met.
                skip_fallback = certificate is not None
                fallback_reason = None if skip_fallback else extension.get(
                    'stop_reason', 'strong_search_no_certificate')
                if certificate is None:
                    failure_stage = 'strong_fullspace_search'
            else:
                route, native_kind = 'reused_c1', 'existing_demand_threshold'
                execution.reserve_fallback()
                try:
                    certificate = execution.run('c1_demand_search', lambda:
                        c1._mean_zero(task, execution) if task['mean_zero'] else
                        c1._fullspace(task, execution))
                    if certificate is None:
                        fallback_reason = 'c1_candidate_schedule_exhausted'
                except (ValueError, ArithmeticError, Deadline) as error:
                    fallback_reason = ('adaptive_stage_budget' if isinstance(error, Deadline)
                                       else 'c1_search_precondition_or_candidate_failure')
                    failure_stage = 'c1_demand_search'
                    execution.note(fallback_reason, error=type(error).__name__, reason=str(error))
                finally:
                    execution.release_reserve()

        if certificate is None and not skip_fallback:
            execution.release_reserve()
            execution.check()
            route = 'baseline_fallback'
            execution.note('baseline_fallback', fallback_reason=fallback_reason,
                           remaining_seconds=execution.remaining())
            execution.counters['baseline_fallback_calls'] += 1
            remaining_task = {**task, 'budget': {'wall_seconds': execution.remaining()}}
            result = execution.run('baseline_remaining_budget', lambda: baseline.solve(remaining_task))
            certificate = result.get('certificate')
            fallback_assessment = {k: result[k] for k in (
                'certificate_valid', 'target_met', 'status', 'metric', 'tolerance', 'error') if k in result}
            execution.attempts[-1].update(baseline_status=result.get('status'),
                                         baseline_attempts=result.get('attempts'),
                                         baseline_details=result.get('details'))
        if fallback_assessment is not None:
            assessment = fallback_assessment
        elif certificate is not None:
            execution.check()
            certificate = json.loads(json.dumps(certificate, allow_nan=False))
            execution.counters['final_frozen_replays'] += 1
            assessment = execution.run('final_frozen_replay', lambda: verify(certificate, task))
            if not assessment.get('certificate_valid'):
                failure_stage = 'final_frozen_replay'
                execution.note('original_verifier_rejected', assessment=assessment)
            elif not assessment.get('target_met'):
                execution.note('valid_certificate_precision_insufficient', metric=assessment.get('metric'),
                               tolerance=task['tolerance'])
                failure_stage = 'precision_gate'
    except Deadline as error:
        assessment.update(target_met=False, status='budget_exceeded', error=str(error))
        failure_stage = failure_stage or 'whole_solve_budget'
    except Exception as error:
        assessment.update(certificate_valid=False, target_met=False, status='execution_failed',
                          error=type(error).__name__+': '+str(error))
        failure_stage = failure_stage or 'unexpected_execution_exception'
        execution.note('unexpected_execution_exception', error=assessment['error'])
    elapsed = perf_counter()-started
    within = execution.limit > 0 and elapsed <= execution.limit
    if not within:
        assessment.update(target_met=False, status='budget_exceeded')
    return {'version': VERSION, 'route': route, 'native_kind': native_kind,
            'certificate': certificate, **assessment, 'within_budget': within,
            'fallback_reason': fallback_reason, 'failure_stage': failure_stage,
            'unverified_candidates': unverified_candidates,
            'total_elapsed_seconds': elapsed, 'attempts': execution.attempts,
            'counters': execution.counters,
            'counter_scope': 'Search operations; excludes frozen verifier internals.',
            'budget_scope': 'Cooperative whole solve; evaluation enforces process deadline.',
            'actual_model_calls': 0, 'actual_model_tokens': None}


if __name__ == '__main__':
    import sys
    for line in sys.stdin:
        if line.strip():
            print(json.dumps(solve(json.loads(line)), allow_nan=False), flush=True)
