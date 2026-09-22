"""Target-width driven degree search using the unchanged step certificate.

The exact screening bounds guide search only. The selected trial is rebuilt by
the frozen certificate constructor; the controller owns the final task-bound
replay. This module neither runs a baseline fallback nor changes old globals.
"""
from fractions import Fraction as F

from runtime import baseline, Deadline
import step_solver as step


DEGREES = (2, 4, 6, 8, 10)
PRECISION_BITS = 96
ITERATIONS = 12


def _increment(execution, key):
    execution.counters[key] = execution.counters.get(key, 0)+1


def _source(q, mean_zero, tolerance, execution):
    """Identical complete-source settings to frozen step_solver.solve."""
    _increment(execution, 'step_sources_generated')
    source = execution.run('step_complete_source', lambda: step.direct.full_ground(
        q, mean_zero, modes=4, bits=20, max_m=2,
        tolerance=tolerance, near_tail=4))
    m = 1 if mean_zero else 0
    if len(source['sectors']) <= m:
        extra = execution.run('step_explicit_target_sector', lambda:
            step.direct.certify_sector(q, m, mean_zero, 4, 20, 40, 4), azimuth_m=m)
        source = execution.run('step_complete_source_extension', lambda:
            step.direct.full_certificate(q, mean_zero, source['sectors']+[extra], tolerance))
    return source


def _screen(q, mean_zero, source, gap, statistics, tolerance):
    """Exact copy of the frozen certificate's scalar bound composition.

    None of these provisional fields is a certificate. Acceptance still calls
    trial_certificate, which verifies the complete source and recomputes them.
    """
    m = 1 if mean_zero else 0
    beta, mu = F(gap['lower']), F(statistics['rayleigh'])
    variance = F(statistics['residual_squared'])
    result = {'screen_is_certificate': False, 'beta': str(beta),
              'rayleigh': str(mu), 'residual_squared': str(variance),
              'strict_temple_gap': mu < beta, 'tolerance': str(tolerance)}
    if mu >= beta:
        return result
    temple = mu-variance/(beta-mu)
    selected = max(temple, F(source['sectors'][m]['lower']))
    others = [{'azimuth_m': c['azimuth_m'], 'lower': c['lower']}
              for c in source['sectors'] if c['azimuth_m'] != m]
    angular = F(source['angular_tail_lower'])
    pointwise = (2 if mean_zero else 0)+min(F(q['offset']), F(q['offset'])+F(q['amplitude']))
    lower = max(pointwise, min([selected, angular]+[F(c['lower']) for c in others]))
    upper = min(mu, F(source['upper']))
    if lower > upper:
        raise ArithmeticError('Contradictory screened global bounds')
    required_lower = upper-F(tolerance)
    limiting_others = ([dict(c) for c in others if F(c['lower']) < required_lower]
                       if pointwise < required_lower else [])
    angular_limits = pointwise < required_lower and angular < required_lower
    result.update(temple_denominator=str(beta-mu), temple_lower=str(temple),
                  selected_sector_lower=str(selected), other_sector_lowers=others,
                  angular_tail_lower=str(angular), pointwise_global_lower=str(pointwise),
                  lower=str(lower), upper=str(upper), exact_width=str(upper-lower),
                  required_global_lower=str(required_lower),
                  limiting_other_sectors=limiting_others, angular_tail_limits=angular_limits,
                  target_met=upper-lower <= F(tolerance))
    return result


def _assemble(q, task, source, gap, candidate, execution):
    proposal = candidate['proposal']
    def build():
        _increment(execution, 'step_certificates_assembled')
        return step.trial_certificate(q, proposal['basis'], proposal['coefficients'], source,
                                      gap, task['tolerance'], task['mean_zero'])
    certificate = execution.run('step_assemble_frozen_certificate', build, degree=candidate['degree'])
    for key in ('lower', 'upper', 'exact_width'):
        if certificate[key] != candidate['screen'][key]:
            raise ArithmeticError('Frozen certificate disagrees with screening: '+key)
    return certificate


def _assemble_using_reserve_if_needed(q, task, source, gap, candidate, execution):
    try:
        return _assemble(q, task, source, gap, candidate, execution)
    except Deadline as error:
        execution.note('step_assembly_uses_reserve', selected_degree=candidate['degree'],
                       cause=type(error).__name__, reason=str(error),
                       next_action='assemble_existing_candidate_without_more_search')
        execution.release_reserve()
        # A genuinely exhausted whole-call budget still raises Deadline here.
        return _assemble(q, task, source, gap, candidate, execution)


def search(raw_task, execution):
    """Return an original step certificate (possibly open) or None.

    execution provides check/run/remaining/attempts/counters/note and the root's
    release_reserve hook. A stage deadline preserves the best candidate using
    the reserve; without a candidate, or when the whole budget is exhausted,
    Deadline propagates. No task id or measured runtime selects a degree.
    """
    execution.check()
    task = baseline.normalize_task(raw_task)
    if task['kind'] != 'spectrum' or task['function']['profile'] != 'step':
        execution.note('step_unsupported_task', kind=task['kind'],
                       profile=task['function']['profile'])
        return None
    q = step.normalize(task['function'])
    m, tolerance = (1 if task['mean_zero'] else 0), F(task['tolerance'])
    execution.note('step_degree_schedule', degrees=list(DEGREES),
                   precision_bits=PRECISION_BITS, iterations=ITERATIONS,
                   policy='first_exact_global_width_at_or_below_original_tolerance',
                   tolerance=task['tolerance'], azimuth_m=m)
    try:
        source = _source(q, task['mean_zero'], task['tolerance'], execution)
        gap = step._gap(q, m, task['mean_zero'], None)
    except (ValueError, ArithmeticError) as error:
        execution.note('step_source_failed', error=type(error).__name__, reason=str(error))
        return None
    execution.note('step_source_ready', source=source, gap=gap)
    # This is a proved obstruction for this fixed pointwise gap, rather than
    # merely a poor candidate. A Rayleigh value cannot be below its sector floor.
    if F(source['sectors'][m]['lower']) >= F(gap['lower']):
        execution.note('step_pointwise_gap_provably_insufficient',
                       selected_sector_lower=source['sectors'][m]['lower'],
                       beta=gap['lower'], next_action='controller_fallback')
        return None
    best = None
    interrupted = False
    for degree in DEGREES:
        try:
            execution.check()
            _increment(execution, 'candidate_proposals')
            proposal = execution.run('step_candidate', lambda: step.propose_trial(
                q, degree=degree, m=m, precision_bits=PRECISION_BITS, iterations=ITERATIONS),
                degree=degree, azimuth_m=m)
            statistics = execution.run('step_full_strong_residual', lambda:
                step.trial_statistics(q, proposal['basis'], proposal['coefficients'], m), degree=degree)
            _increment(execution, 'step_exact_width_checks')
            screen = _screen(q, task['mean_zero'], source, gap, statistics, tolerance)
        except Deadline as error:
            execution.note('step_stage_budget', interrupted_degree=degree,
                           cause=type(error).__name__, reason=str(error),
                           selected_degree=None if best is None else best['degree'],
                           next_action='controller_fallback' if best is None else 'preserve_best_candidate')
            if best is None:
                raise
            execution.release_reserve()
            interrupted = True
            break
        except (ValueError, ArithmeticError) as error:
            execution.note('step_candidate_failed', degree=degree,
                           error=type(error).__name__, reason=str(error))
            continue
        candidate = {'degree': degree, 'proposal': proposal,
                     'statistics': statistics, 'screen': screen}
        execution.note('step_candidate_screened', **candidate)
        if not screen['strict_temple_gap']:
            execution.note('step_temple_gap_not_strict', degree=degree,
                           rayleigh=screen['rayleigh'], beta=screen['beta'],
                           next_action='increase_degree' if degree != DEGREES[-1] else 'controller_fallback')
            continue
        if screen['limiting_other_sectors'] or screen['angular_tail_limits']:
            execution.note('step_other_sector_limits_candidate', degree=degree,
                           required_global_lower=screen['required_global_lower'],
                           limiting_other_sectors=screen['limiting_other_sectors'],
                           angular_tail_limits=screen['angular_tail_limits'],
                           claim_scope='current_candidate_only_not_impossibility')
        if best is None or F(screen['exact_width']) < F(best['screen']['exact_width']):
            best = candidate
        if screen['target_met']:
            execution.note('step_target_width_reached', degree=degree,
                           exact_width=screen['exact_width'], tolerance=task['tolerance'])
            break
    if best is None:
        execution.note('step_no_temple_eligible_candidate', next_action='controller_fallback')
        return None
    if not best['screen']['target_met'] and not interrupted:
        execution.note('step_degree_schedule_exhausted', selected_degree=best['degree'],
                       exact_width=best['screen']['exact_width'], tolerance=task['tolerance'],
                       next_action='return_original_certified_open_avoid_duplicate_degree10')
    try:
        return _assemble_using_reserve_if_needed(q, task, source, gap, best, execution)
    except (ValueError, ArithmeticError) as error:
        execution.note('step_certificate_assembly_failed', error=type(error).__name__,
                       reason=str(error), next_action='controller_fallback')
        return None
