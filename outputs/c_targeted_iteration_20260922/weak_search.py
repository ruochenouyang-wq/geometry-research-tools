"""Keep C1's weak-potential search; reuse exact data for unverified assembly."""
from fractions import Fraction as F

from runtime import baseline, c1


def search(raw_task, execution):
    import fullspace_singular as full
    import certification
    from fast_candidate import prepare_trial

    task = baseline.normalize_task(raw_task)
    if task['kind'] != 'spectrum' or task['mean_zero'] is not False:
        raise ValueError('Weak fullspace search preserves the full-space spectrum task')
    q, tolerance = full.normalize(task['function']), F(task['tolerance'])
    try:
        form = execution.run('centered_form', lambda: full.form_certificate(q, 'centered_L2'))
    except ValueError:
        execution.note('weak_route_requires_centered_form', next_action='baseline_fallback')
        return None
    kernel = c1._kernel(q, 0, False, 8, 16, form, full.verify_form, execution)
    for count in (3, 6, 10, 14, 16):
        execution.check()
        powers = full.generated_powers(q, count)
        execution.counters['candidate_proposals'] += 1
        proposal = execution.run('fullspace_candidate', lambda: full.proposal(
            q, powers, form, kernel.beta, 112, 24), terms=count,
            provisional_search_beta=str(kernel.beta), provisional_beta_is_proof=False)
        prepared = execution.run('prepare_exact_trial_once', lambda: prepare_trial(q, powers, full=full))
        stats = execution.run('full_residual', lambda: full.statistics(
            q, powers, proposal['coefficients'], (prepared.M, prepared.H, prepared.R)))
        gap = c1._needed_gap(kernel, F(stats['rayleigh']), F(stats['residual_squared']),
                             tolerance, form, execution)
        if gap is None:
            continue
        def assemble():
            execution.counters['deferred_assemblies'] = execution.counters.get('deferred_assemblies', 0)+1
            return certification.prepare_fullspace(task, gap).assemble(
                powers, proposal['coefficients'], prepared=prepared, statistics=stats)
        pending = execution.run('assemble_pending_certificate', assemble)
        certificate = pending['certificate']
        if F(certificate['exact_width']) <= tolerance:
            return certificate
    return None
