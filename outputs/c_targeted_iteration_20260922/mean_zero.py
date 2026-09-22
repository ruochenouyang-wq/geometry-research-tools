"""Choose the constrained radial block to refine using certified intervals.

Every returned proof has the unchanged direct_moments.FULL format.  Positive
amplitude is a capability gate, never evidence that m=0 is the lowest block.
The common controller owns the final task-bound frozen replay.
"""
from fractions import Fraction as F

from runtime import baseline, Deadline
import spectral_gap

direct = spectral_gap.direct
SCHEDULE = ((8, 12), (12, 16), (16, 24))
MAX_M = 4
MAX_BISECTIONS = 96


def _increment(execution, name):
    execution.counters[name] = execution.counters.get(name, 0)+1


class _Kernel:
    """Count and budget only this search's exact comparisons."""
    def __init__(self, kernel, execution):
        self.kernel, self.execution = kernel, execution

    def __getattr__(self, name):
        return getattr(self.kernel, name)

    def count(self, value, kind='lower', independent=False):
        self.execution.check()
        matrix = self.kernel.matrix(value, kind)
        _increment(self.execution, 'mean_zero_inertia_checks')
        return (direct.base.inertia(matrix) if independent else
                direct.arithmetic.banded_inertia(matrix))

    def lower_holds(self, value):
        return value < self.beta and self.count(value)[0] == 0

    def upper_holds(self, value):
        negative, zero, _ = self.count(value, 'upper')
        return negative+zero >= 1


def _sector(q, m, modes, near, arithmetic_tolerance, execution):
    def construct():
        _increment(execution, 'mean_zero_kernel_constructions')
        # True is essential: m=0 removes l=0; other blocks remain unchanged.
        return _Kernel(direct.Kernel(q, m, True, modes, 40, near), execution)
    kernel = execution.run('mean_zero_construct_kernel', construct,
                           active_m=m, modes=modes, near_tail=near)
    analytic = direct.tail_bound(kernel.form, kernel.start)
    upper = min(kernel.data['A'][i][i]/kernel.data['mass'][i] for i in range(kernel.n))
    if not kernel.upper_holds(upper):
        raise ArithmeticError('Certified diagonal Ritz bracket failed')
    low, high = analytic-1, upper
    for _ in range(MAX_BISECTIONS):
        if high-low <= arithmetic_tolerance:
            break
        middle = (low+high)/2
        if kernel.upper_holds(middle):
            high = middle
        else:
            low = middle
    upper = high
    distance = F(1)
    for _ in range(64):
        low = analytic-distance
        if kernel.lower_holds(low):
            break
        distance *= 2
    else:
        raise ValueError('Complete lower comparison bracket exhausted')
    high = min(upper, kernel.beta)
    for _ in range(MAX_BISECTIONS):
        if high-low <= arithmetic_tolerance:
            break
        middle = (low+high)/2
        if kernel.lower_holds(middle):
            low = middle
        else:
            high = middle
    # Independently recount both endpoints. No bisection flag is a proof.
    certificate = execution.run('mean_zero_sector_certificate', lambda:
        direct.sector_certificate(kernel, low, upper), active_m=m, modes=modes)
    execution.note('mean_zero_complete_sector', active_m=m, modes=modes,
                   near_tail=near, lower=certificate['lower'], upper=certificate['upper'],
                   exact_width=certificate['exact_width'],
                   projection=certificate['projection'])
    return certificate


def _summary(certificate, tolerance):
    sectors = certificate['sectors']
    upper = F(certificate['upper'])
    angular = F(certificate['angular_tail_lower'])
    active = [s['azimuth_m'] for s in sectors if F(s['lower']) < upper-tolerance]
    dominant = []
    for sector in sectors:
        competitors = [angular]+[F(s['lower']) for s in sectors
                                  if s['azimuth_m'] != sector['azimuth_m']]
        if F(sector['upper']) < min(competitors):
            dominant.append(sector['azimuth_m'])
    return {'lower': certificate['lower'], 'upper': certificate['upper'],
            'exact_width': certificate['exact_width'], 'active_sectors': active,
            'proved_dominant_sectors': dominant,
            'angular_tail_lower': certificate['angular_tail_lower'],
            'angular_tail_limits': angular < upper-tolerance,
            'sector_intervals': [{'m': s['azimuth_m'], 'lower': s['lower'],
                                  'upper': s['upper'], 'modes': s['modes'],
                                  'near_tail': s['near_tail']} for s in sectors]}


def _assemble(q, sectors, tolerance, execution):
    def build():
        _increment(execution, 'mean_zero_full_certificates')
        return direct.full_certificate(q, True, sectors, tolerance)
    return execution.run('mean_zero_assemble_original_full', build)


def _assemble_with_reserve(q, sectors, tolerance, execution):
    try:
        return _assemble(q, sectors, tolerance, execution), False
    except Deadline as error:
        execution.note('mean_zero_assembly_uses_reserve', reason=str(error),
                       next_action='assemble_existing_complete_sectors_without_more_search')
        execution.release_reserve()
        return _assemble(q, sectors, tolerance, execution), True


def search(raw_task, execution):
    """Return a frozen full proof, possibly open, or None for the old route.

    execution implements check/run/note/counters/release_reserve.  No local
    baseline fallback or final assess is called; no case identifier is read.
    """
    execution.check()
    task = baseline.normalize_task(raw_task)
    q = task['function']
    if (task['kind'] != 'spectrum' or task.get('mean_zero') is not True or
            q['profile'] != 'abs_power' or not -F(1, 2) < F(q['exponent']) < 0 or
            F(q['amplitude']) <= 0):
        execution.note('mean_zero_capability_delegated',
                       reason='Requires positive-amplitude singular power in the original mean-zero space')
        return None
    try:
        form = execution.run('mean_zero_centered_form', lambda: direct.form_bound(q, 40))
    except (ValueError, ArithmeticError) as error:
        execution.note('mean_zero_form_delegated', reason=str(error))
        return None
    tolerance = F(task['tolerance'])
    sectors, levels, best = [], {}, None
    try:
        # Coarse complete intervals decide where precision work is useful.
        for m in (0, 1):
            sectors.append(_sector(q, m, *SCHEDULE[0], F(1, 16384), execution))
            levels[m] = 0
        current, used_reserve = _assemble_with_reserve(q, sectors, task['tolerance'], execution)
        best = current
        if used_reserve:
            execution.note('mean_zero_search_stopped_after_reserved_assembly',
                           **_summary(best, tolerance))
            return best
        while True:
            summary = _summary(current, tolerance)
            execution.note('mean_zero_global_interval', **summary)
            if F(current['exact_width']) <= tolerance:
                execution.note('mean_zero_original_target_reached',
                               **summary, final_task_bound_replay_required=True)
                return best
            if summary['angular_tail_limits']:
                m = len(sectors)
                if m > MAX_M:
                    execution.note('mean_zero_angular_schedule_exhausted',
                                   first_omitted_m=m, **summary)
                    return best
                execution.note('mean_zero_extend_angular_coverage', active_m=m,
                               reason='Complete angular tail still limits original global width')
                sectors.append(_sector(q, m, *SCHEDULE[0], F(1, 16384), execution))
                levels[m] = 0
            else:
                # The smallest certified floor controls the global interval.
                # A sign heuristic is never substituted for these endpoints.
                active_m = min(summary['active_sectors'], key=lambda m: F(sectors[m]['lower']))
                next_level = levels[active_m]+1
                if next_level >= len(SCHEDULE):
                    execution.note('mean_zero_active_sector_schedule_exhausted',
                                   active_m=active_m, **summary)
                    return best
                modes, near = SCHEDULE[next_level]
                execution.note('mean_zero_refine_active_sector', active_m=active_m,
                               reason='This certified lower endpoint limits the global interval',
                               next_modes=modes, next_near_tail=near, **summary)
                levels[active_m] = next_level
                _increment(execution, 'mean_zero_sector_refinements')
                try:
                    updated = _sector(q, active_m, modes, near, tolerance/F(32), execution)
                except (ValueError, ArithmeticError) as error:
                    execution.note('mean_zero_sector_refinement_failed', active_m=active_m,
                                   modes=modes, near_tail=near, reason=str(error))
                    continue
                # Preserve all unchanged complete sectors, including m=1.
                sectors[active_m] = updated
            current, used_reserve = _assemble_with_reserve(q, sectors, task['tolerance'], execution)
            if F(current['exact_width']) < F(best['exact_width']):
                best = current
            if used_reserve:
                execution.note('mean_zero_search_stopped_after_reserved_assembly',
                               **_summary(best, tolerance))
                return best
    except Deadline as error:
        execution.note('mean_zero_stage_budget', reason=str(error),
                       retained_complete_certificate=best is not None,
                       retained_width=None if best is None else best['exact_width'])
        if best is None:
            raise
        execution.release_reserve()
        return best
    except (ValueError, ArithmeticError) as error:
        execution.note('mean_zero_complete_search_failed', reason=str(error),
                       retained_complete_certificate=best is not None)
        return best
