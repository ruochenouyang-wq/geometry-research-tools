"""Version C: certify the inequalities needed by the requested precision.

Search changes only. Every returned mathematical certificate uses a frozen
format and is independently replayed by the original task-bound verifier.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import sys
from time import perf_counter

sys.dont_write_bytecode = True
PARENT = Path(__file__).resolve().parent.parent
if str(PARENT) not in sys.path:
    sys.path.insert(0, str(PARENT))
import shared_baseline

baseline = shared_baseline.baseline
VERSION = 'c-proof-demand-threshold-v1.1'


class Deadline(Exception):
    pass


class Execution:
    def __init__(self, started, limit):
        self.started, self.limit = started, limit
        self.attempts = []
        self.counters = {'search_inertia_checks': 0, 'gap_threshold_checks': 0,
                         'kernel_constructions': 0, 'candidate_proposals': 0,
                         'final_frozen_replays': 0, 'baseline_fallback_calls': 0,
                         'adaptive_policy_skips': 0}

    def remaining(self):
        return max(0.0, self.limit-(perf_counter()-self.started))

    def check(self):
        if self.remaining() <= 0:
            raise Deadline('Cooperative whole-solve budget exhausted')

    def run(self, action, function, **details):
        self.check()
        begin = perf_counter()
        row = {'action': action, **details}
        try:
            value = function()
            row['status'] = 'returned'
            return value
        except Exception as error:
            row.update(status='failed', error=type(error).__name__, reason=str(error))
            raise
        finally:
            row['elapsed_seconds'] = perf_counter()-begin
            self.attempts.append(row)


class CountedKernel:
    """Instrument only this execution's explicitly requested comparisons."""
    def __init__(self, kernel, execution):
        self.kernel, self.execution = kernel, execution

    def __getattr__(self, key):
        return getattr(self.kernel, key)

    def count(self, x, kind='lower', independent=False):
        self.execution.check()
        # Tail validation precedes inertia; rejected tails are not counted.
        matrix = self.kernel.matrix(x, kind)
        self.execution.counters['search_inertia_checks'] += 1
        import spectral_gap
        return (spectral_gap.direct.base.inertia(matrix) if independent else
                spectral_gap.direct.arithmetic.banded_inertia(matrix))

    def lower_holds(self, x):
        return x < self.beta and self.count(x)[0] == 0

    def upper_holds(self, x):
        negative, zero, _ = self.count(x, 'upper')
        return negative+zero >= 1


def _ceil_dyadic(value, bits=20):
    value = F(value)
    denominator = 1 << bits
    return F(-((-value.numerator*denominator)//value.denominator), denominator)


def _kernel(q, m, mean_zero, modes, near, form, form_verifier, execution):
    import spectral_gap
    def construct():
        execution.counters['kernel_constructions'] += 1
        return spectral_gap._kernel(q, m, mean_zero, modes, near, 40,
                                    form, form_verifier)
    kernel = execution.run('construct_kernel', construct, m=m, modes=modes, near_tail=near)
    return CountedKernel(kernel, execution)


def _gap_at(kernel, beta, form, execution):
    import spectral_gap
    execution.check()
    beta = F(beta)
    def certify():
        # _certificate performs exactly one independent inertia after matrix().
        # Check its strict tail condition before counting that operation.
        if beta >= kernel.beta:
            raise ValueError('Requested gap reaches the strict radial tail barrier')
        execution.counters['gap_threshold_checks'] += 1
        execution.counters['search_inertia_checks'] += 1
        return spectral_gap._certificate(kernel, beta, form)
    return execution.run('gap_threshold', certify, beta=str(beta),
                         radial_tail_lower=str(kernel.beta))


def _needed_gap(kernel, mu, variance, tolerance, form, execution):
    """A sufficient Temple denominator; no floating-point acceptance."""
    seen = set()
    for portion in (F(15,16), F(1)):
        epsilon = tolerance*portion
        required = mu+max(variance/epsilon, F(1, 1 << 24))
        # A modest dyadic denominator avoids huge fractions in exact inertia.
        beta = _ceil_dyadic(required)
        if beta in seen:
            continue
        seen.add(beta)
        try:
            proof = _gap_at(kernel, beta, form, execution)
            execution.attempts[-1].update(allocated_width=str(epsilon),
                                         required_beta=str(required))
            return proof
        except (ValueError, ArithmeticError):
            execution.attempts[-1].update(allocated_width=str(epsilon),
                                         required_beta=str(required))
    return None


def _fullspace(task, execution):
    import fullspace_singular as full
    import power_forms
    q, tolerance = full.normalize(task['function']), F(task['tolerance'])
    try:
        form = execution.run('centered_form', lambda: full.form_certificate(q, 'centered_L2'))
    except ValueError:
        form = None
    if form is None:
        # Public development exposed a strong-potential failure in the
        # provisional search-beta heuristic. Keep this whole form regime on
        # the original verified-gap-first policy; do not spend its fallback
        # budget on repeated speculative candidates. Mean-zero is separate.
        execution.counters['adaptive_policy_skips'] += 1
        execution.attempts.append({'action':'adaptive_form_regime_gate',
            'status':'baseline_fallback','reason':'centered_form_unavailable',
            'elapsed_seconds':0.0,
            'timing_note':'Use total_elapsed_seconds for the full cost; this row records the decision only.'})
        return None
    modes_schedule = (8,)
    for modes in modes_schedule:
        current_form = form or execution.run('holder_form', lambda: power_forms.select_form(q, modes))
        kernel = _kernel(q, 0, False, modes, 16, current_form, full.verify_form, execution)
        for count in (3, 6, 10, 14, 16):
            execution.check()
            powers = full.generated_powers(q, count)
            execution.counters['candidate_proposals'] += 1
            # kernel.beta is ONLY a candidate-search heuristic here, not a
            # claimed second-eigenvalue bound. _needed_gap supplies the proof.
            proposal = execution.run('fullspace_candidate', lambda: full.proposal(
                q, powers, current_form, kernel.beta, 112, 24), terms=count,
                provisional_search_beta=str(kernel.beta), provisional_beta_is_proof=False)
            stats = execution.run('full_residual', lambda: full.statistics(q, powers, proposal['coefficients']))
            mu, variance = F(stats['rayleigh']), F(stats['residual_squared'])
            gap = _needed_gap(kernel, mu, variance, tolerance, current_form, execution)
            if gap is None:
                continue
            cert = execution.run('assemble_original_certificate', lambda: full.certificate(
                q, powers, proposal['coefficients'], gap, task['tolerance']))
            if F(cert['exact_width']) <= tolerance:
                return cert
    return None


def _coarse_source(kernel, form, execution):
    """Enough source information for the trial shift, without full bisection."""
    import spectral_gap
    direct = spectral_gap.direct
    analytic = direct.tail_bound(kernel.form, kernel.start)
    upper = min(kernel.data['A'][i][i]/kernel.data['mass'][i] for i in range(kernel.n))
    distance = F(1)
    for _ in range(64):
        lower = analytic-distance
        if kernel.lower_holds(lower):
            break
        distance *= 2
    else:
        raise ValueError('Coarse lower bracket exhausted')
    high = min(upper, kernel.beta)
    for _ in range(8):
        middle = (lower+high)/2
        if kernel.lower_holds(middle):
            lower = middle
        else:
            high = middle
    return (spectral_gap._ground_certificate(kernel, lower, upper, form)
            if form is not None else direct.sector_certificate(kernel, lower, upper))


def _ground_at(kernel, lower, form, execution):
    import spectral_gap
    upper = min(kernel.data['A'][i][i]/kernel.data['mass'][i] for i in range(kernel.n))
    def certify():
        if lower > upper:
            raise ValueError('Requested exclusion lower exceeds a certified Ritz upper')
        return (spectral_gap._ground_certificate(kernel, lower, upper, form)
                if form is not None else spectral_gap.direct.sector_certificate(kernel, lower, upper))
    return execution.run('other_sector_threshold', certify, lower=str(lower), m=kernel.m)


def _mean_zero(task, execution):
    import singular_solver as singular
    import power_forms
    q, tolerance = singular.normalize(task['function']), F(task['tolerance'])
    default = singular._default_form_available(q)
    engine = singular._engine(q)
    for modes in ((8,) if default else (12, 16)):
        form = None if default else execution.run('holder_form', lambda: power_forms.select_form(q, modes+1))
        m1_kernel = _kernel(q, 1, True, modes, 12, form, singular._verify_form, execution)
        source = execution.run('coarse_m1_source', lambda: _coarse_source(m1_kernel, form, execution))
        shift = F(source['lower'])-F(1,1024)
        m0_kernel = None
        for count in (3, 6, 10, 16):
            execution.check()
            powers = engine.generated_powers(q, count)
            execution.counters['candidate_proposals'] += 1
            proposal = execution.run('mean_zero_candidate', lambda: engine._inverse_proposal(
                q, powers, 112, 14, shift=shift), terms=count)
            stats = execution.run('full_residual', lambda: engine.statistics(q, powers, proposal['coefficients']))
            mu, variance = F(stats['rayleigh']), F(stats['residual_squared'])
            gap = _needed_gap(m1_kernel, mu, variance, tolerance, form, execution)
            if gap is None:
                continue
            if m0_kernel is None:
                m0_kernel = _kernel(q, 0, True, modes, 12, form, singular._verify_form, execution)
            lower = _ceil_dyadic(mu-tolerance*F(15,16))
            try:
                m0 = _ground_at(m0_kernel, lower, form, execution)
            except (ValueError, ArithmeticError):
                continue
            cert = execution.run('assemble_original_certificate', lambda: singular.certificate(
                q, powers, proposal['coefficients'], gap, m0, source, task['tolerance']))
            if F(cert['exact_width']) <= tolerance:
                return cert
    return None


def verify(certificate, task):
    return baseline.assess(certificate, baseline.normalize_task(task))


def solve(raw_task):
    started = perf_counter()
    task = baseline.normalize_task(raw_task)
    execution = Execution(started, task['budget']['wall_seconds'])
    certificate = None
    route = 'proof_demand_threshold'
    fallback_assessment = None
    assessment = {'certificate_valid': False, 'target_met': False, 'status': 'no_verified_certificate'}
    try:
        execution.check()
        if task['kind'] == 'spectrum' and task['function']['profile'] == 'abs_power':
            try:
                certificate = (_mean_zero(task, execution) if task['mean_zero'] else _fullspace(task, execution))
            except (ValueError, ArithmeticError) as error:
                # A candidate policy failure is not a claim that the task is
                # impossible. Preserve it before the original bounded route.
                execution.attempts.append({'action': 'adaptive_route_failed',
                    'status': 'failed', 'error': type(error).__name__, 'reason': str(error),
                    'elapsed_seconds': 0.0, 'timing_note': 'Only total_elapsed_seconds is the complete cost; action timings may nest and omit preparation.'})
        else:
            route = 'unchanged_baseline_route'
        if certificate is None:
            execution.check()
            remaining_task = {**task, 'budget': {'wall_seconds': execution.remaining()}}
            execution.counters['baseline_fallback_calls'] += 1
            fallback = execution.run('baseline_remaining_budget', lambda: baseline.solve(remaining_task))
            certificate = fallback.get('certificate')
            # The trusted baseline entry point already serializes and performs
            # its task-bound final replay. Preserve that result without adding
            # a duplicate replay to the unchanged step/approximation routes.
            fallback_assessment = {key: fallback[key] for key in (
                'certificate_valid','target_met','status','metric','tolerance','error') if key in fallback}
            execution.attempts[-1].update(baseline_status=fallback.get('status'),
                                         baseline_attempts=fallback.get('attempts'))
        if fallback_assessment is not None:
            assessment = fallback_assessment
        elif certificate is not None:
            execution.check()
            certificate = json.loads(json.dumps(certificate, allow_nan=False))
            execution.counters['final_frozen_replays'] += 1
            assessment = execution.run('final_frozen_replay', lambda: baseline.assess(certificate, task))
    except Deadline:
        assessment = {'certificate_valid': False, 'target_met': False, 'status': 'budget_exceeded'}
    except Exception as error:
        assessment = {'certificate_valid': False, 'target_met': False, 'status': 'execution_failed',
                      'error': type(error).__name__+': '+str(error)}
    elapsed = perf_counter()-started
    within = elapsed <= execution.limit
    if not within or execution.limit == 0:
        assessment.update(target_met=False, status='budget_exceeded')
    return {'version': VERSION, 'route': route, 'certificate': certificate, **assessment,
            'within_budget': within and execution.limit > 0, 'total_elapsed_seconds': elapsed,
            'attempts': execution.attempts, 'counters': execution.counters,
            'counter_scope': 'Explicit search comparisons only; excludes Gram PSD checks and frozen replay internals.',
            'budget_scope': 'whole_solve_call_cooperative; external hard deadline required',
            'actual_model_calls': 0, 'actual_model_tokens': None}


if __name__ == '__main__':
    for line in sys.stdin:
        if line.strip():
            try:
                print(json.dumps(solve(json.loads(line)), allow_nan=False), flush=True)
            except Exception as error:
                print(json.dumps({'target_met': False, 'error': type(error).__name__+': '+str(error)}), flush=True)
