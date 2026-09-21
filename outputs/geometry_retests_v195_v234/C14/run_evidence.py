"""Re-run original C14 and all ten substantive V205--V214 capabilities."""
from pathlib import Path
import io
import json
import sys
import time
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from support import ROOT, F, save, original_case, old_spectrum
import general_spectrum as g


def main():
    began = time.perf_counter()
    original = original_case('C14')
    save('C14/original_case.json', original)
    baseline = []
    for mean_zero in (True, False):
        try:
            previous = old_spectrum.adaptive_spectrum(original['input']['q'], k=9, mean_zero=mean_zero)['certificate']
            count = old_spectrum.count_below(previous['assembly'])
            trace = old_spectrum.negative_trace(previous['assembly'])
            baseline.append({'mean_zero': mean_zero, 'status': 'accepted',
                'count': [count['count_lower'], count['count_upper']],
                'trace': [trace['lower'], trace['upper']]})
        except (ValueError, TypeError) as exc:
            baseline.append({'mean_zero': mean_zero, 'status': 'rejected',
                             'exception': type(exc).__name__, 'message': str(exc)})
    save('C14/baseline.json', {'original_case': original, 'actual_retest': baseline,
         'diagnosis': 'The mean-zero half is already solved; the original unprojected half is rejected. Nonconstant q requires genuinely different m=0 radial compression.'})
    old = old_spectrum.adaptive_spectrum([-20, 0, 1], k=1, modes=6, max_modes=8,
        max_m=4, max_radial=6, bits=32, tolerance='1/1000000')['certificate']
    old_trace = old_spectrum.negative_trace(old['assembly'], '1/1000000')
    save('C14/scheduler_baseline.json', {'q': ['-20', '0', '1'], 'mean_zero': True, 'k': 1,
         'old_ordered_status': old['status'], 'old_trace_status': old_trace['status'],
         'old_trace_width': old_trace['width'],
         'diagnosis': 'Closing only the first eigenvalue is insufficient for a requested full negative trace.'})
    stages = []
    saved_proofs = []

    def proof(path, cert):
        if not g.verify(cert):
            raise AssertionError('Certificate rejected before saving: '+path)
        save(path, cert); saved_proofs.append(path)
        return path

    def stage(version, capability, function, cert, outcome, extras=()):
        path = proof('C14/v%d.json' % version, cert)
        stages.append({'version': version, 'capability': capability,
                       'callable': 'general_spectrum.'+function,
                       'evidence': [path]+list(extras), 'outcome': outcome})

    stage(205, 'Global analytic ordered bounds include l=0 with multiplicity one in the unprojected space.',
          'analytic_bounds', g.analytic_bounds([0, 1], 4, False),
          'q=t: first eigenvalue in [-1,1], next three in [1,3], before any sector solve.')
    full_sector = g.enumerate_sector([0, 1], m=0, radial_count=2, mean_zero=False, modes=6, bits=32)
    projected_path = proof('C14/v206_projected.json', g.enumerate_sector([0, 1], m=0,
        radial_count=2, mean_zero=True, modes=6, bits=32))
    stage(206, 'Distinct radial enumeration for actual full and mean-zero operators.',
          'enumerate_sector', full_sector,
          'For q=t the full m=0 basis starts at l=0 and has a negative ground state; the projected basis starts at l=1.', [projected_path])
    linear = {mean: g.adaptive_spectrum([0, 1], k=6 if not mean else 5,
        mean_zero=mean, tolerance='1/1000000', bits=32)['certificate'] for mean in (False, True)}
    stage(207, 'Full kth eigenvalue merge with uncomputed angular and radial min-max intervals.',
          'ordered_merge', linear[False]['answers']['ordered'],
          'All six requested full-space eigenvalues for q=t certified to tolerance 1e-6.')
    full_model, mean_model = g.assemble([-7], [], False), g.assemble([-7], [], True)
    stage(208, 'Strict threshold counting with the constant mode and exact threshold ties.',
          'strict_count', g.strict_count(full_model), 'Original C14 unprojected negative count exactly 9.')
    stage(209, 'Complete negative trace including every potentially negative omitted mode.',
          'negative_trace', g.negative_trace(full_model),
          'Original C14 unprojected trace exactly 27; l>=3 has lower bound 5 and cannot be negative.')
    trace_aware = g.adaptive_spectrum([-20, 0, 1], k=1, mean_zero=True,
        quantities=['ordered', 'trace'], tolerance='1/1000000', max_m=4,
        max_radial=6, bits=32)
    schedule = save('C14/trace_schedule.json', {
        'search_history': trace_aware['search_history'], 'sector_solves': trace_aware['sector_solves'],
        'budget_or_precision_exhausted': trace_aware['budget_or_precision_exhausted']})
    stage(210, 'Trace-aware adaptive scheduling continues after requested ordered eigenvalues close.',
          'adaptive_spectrum', trace_aware['certificate'],
          'k=1 ordering closes at step 2; full negative trace closes at step 9.', [schedule])
    low = g.adaptive_spectrum([-20, 0, 1], k=1, mean_zero=True,
        quantities=['ordered', 'trace'], tolerance='1/1000000', bits=32, max_steps=2)
    proof('C14/budget_exhaustion.json', low['certificate'])
    stage(211, 'Projection-preserving threshold shift computes both count and sum(T-lambda)+.',
          'spectral_shift', g.spectral_shift(g.assemble([0], [], False), 6),
          'Full q=0 at T=6: strict count 4, shifted negative trace 18; five ties are excluded.')
    tiny = g.negative_trace(g.assemble(['-1/100'], [], False))
    stage(212, 'Scope-bound fixed-V trace quotient with standard-area factor 4 and constant-mode obstruction.',
          'fixed_potential_quotient', g.fixed_potential_quotient(tiny, ['1/100']),
          'V=1/100 gives pi-scaled full trace quotient 25; the 1/(4 epsilon) divergence is explicit.')
    stage(213, 'Verified codimension-one interlacing for the same nonconstant potential.',
          'projection_interlacing', g.projection_interlacing(linear[False]['answers']['ordered'],
              linear[True]['answers']['ordered']),
          'q=t: five interlacing relations bind the actual full and projected operators.')
    cache = {}
    first = g.research_driver(original['input']['q'], spaces=original['input']['spaces'], cache=cache)
    second = g.research_driver(original['input']['q'], spaces=original['input']['spaces'], cache=cache)
    cache_path = save('C14/cache_execution.json', {'first': first['execution'], 'second': second['execution'],
         'cache_binding': 'q, projection, k, quantities, threshold, tolerance and all resource bounds',
         'note': 'Exact constant-potential calibration needs no radial Schur solve.'})
    stage(214, 'One driver gives separate requested-quantity decisions and verified input-bound cache reuse.',
          'research_driver', second['certificate'],
          'Original C14: mean-zero count 8 / trace 20; unprojected count 9 / trace 27; second run reuses both verified results.', [cache_path])
    save('C14/STAGES.json', stages)
    mapping = {'original_input': original['input'], 'driver_input': {'q': original['input']['q'],
        'spaces': original['input']['spaces'], 'k': 9,
        'quantities': ['ordered', 'count', 'trace'], 'threshold': '0'},
        'projection_mapping': {'full_mean_zero': True, 'full_unprojected': False},
        'unchanged_original_potential_and_spaces': True,
        'negative_trace_convention': 'sum_j max(0,-lambda_j), including real angular multiplicity'}
    save('C14/original_input_mapping.json', mapping)
    got = {}
    for result in second['certificate']['results']:
        c, t = result['answers']['count'], result['answers']['trace']
        got[result['request']['space']] = {'count_lower': c['count_lower'], 'count_upper': c['count_upper'],
            'trace_lower': t['lower'], 'trace_upper': t['upper'],
            'unlisted_eigenvalue_floor': t['tails']['unlisted_floor']}
    assert got['full_mean_zero']['count_lower'] == got['full_mean_zero']['count_upper'] == 8
    assert got['full_mean_zero']['trace_lower'] == got['full_mean_zero']['trace_upper'] == '20'
    assert got['full_unprojected']['count_lower'] == got['full_unprojected']['count_upper'] == 9
    assert got['full_unprojected']['trace_lower'] == got['full_unprojected']['trace_upper'] == '27'
    replays = [g.verify(json.loads((ROOT/path).read_text())) for path in saved_proofs]
    suite = unittest.defaultTestLoader.loadTestsFromName('test_general_spectrum')
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    (ROOT/'C14/unit_tests.log').write_text(stream.getvalue())
    save('C14/unit_tests.json', {'all_passed': result.wasSuccessful(), 'tests_run': result.testsRun,
         'failures': len(result.failures), 'errors': len(result.errors), 'skipped': len(result.skipped)})
    save('C14/RESULTS.json', {'case': 'C14', 'versions': [205, 214], 'original_case_completed': True,
         'original_input_mapping': 'C14/original_input_mapping.json', 'results': got,
         'all_passed': result.wasSuccessful() and all(replays), 'tests_run': result.testsRun,
         'saved_certificate_replays': len(replays), 'certificate_replays': replays,
         'certificates': saved_proofs, 'wall_seconds': time.perf_counter()-began,
         'scope_limit': 'Fixed axisymmetric polynomial potentials on unit S2; no universal unprojected L2 trace constant is claimed.',
         'model_speed_or_token_claim': None})
    print({'original_C14_completed': True, 'tests': result.testsRun,
           'certificate_replays': len(replays), 'all_passed': result.wasSuccessful() and all(replays)})


if __name__ == '__main__':
    main()
