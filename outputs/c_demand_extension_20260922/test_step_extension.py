"""Small functional checks of the new search, never timing acceptance tests."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import sys
from time import perf_counter
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from runtime import baseline, Deadline
import step_extension as extension


class Execution:
    def __init__(self, seconds=20):
        self.started, self.seconds = perf_counter(), seconds
        self.attempts, self.counters = [], {}
        self.reserve_released = False

    def remaining(self):
        return max(0, self.seconds-(perf_counter()-self.started))

    def check(self):
        if self.remaining() <= 0:
            raise Deadline('Functional check budget exhausted')

    def run(self, action, function, **details):
        self.check()
        row = {'action': action, **details}
        try:
            out = function()
            row['status'] = 'returned'
            return out
        except Exception as error:
            row.update(status='failed', reason=str(error))
            raise
        finally:
            self.attempts.append(row)

    def note(self, reason_code, **details):
        self.attempts.append({'reason_code': reason_code, **details})

    def release_reserve(self):
        self.reserve_released = True


def task(mean_zero=False, tolerance='1/10000'):
    return {'kind': 'spectrum', 'mean_zero': mean_zero,
            'function': {'kind': 'axis_profile', 'profile': 'step', 'axis': 1,
                         'amplitude': '-7/12', 'offset': '2/9'},
            'tolerance': tolerance, 'budget': {'wall_seconds': 20}}


class StepExtensionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = {}
        for mean_zero in (False, True):
            execution = Execution()
            cert = extension.search(task(mean_zero), execution)
            if cert is None:
                raise AssertionError(execution.attempts)
            cls.results[mean_zero] = (cert, execution)

    def test_both_spaces_use_original_certificate_and_stop_early(self):
        for mean_zero, (cert, execution) in self.results.items():
            judgment = baseline.assess(cert, baseline.normalize_task(task(mean_zero)))
            self.assertTrue(judgment['certificate_valid'])
            self.assertTrue(judgment['target_met'])
            self.assertEqual(cert['format'], extension.step.FORMAT)
            self.assertEqual(cert['trial']['azimuth_m'], int(mean_zero))
            self.assertLess(max(k for _, k in cert['trial']['basis']), 10)
            self.assertEqual(execution.counters['step_sources_generated'], 1)
            self.assertEqual(execution.counters['step_certificates_assembled'], 1)

    def test_screen_rebuild_matches_original_at_every_recorded_candidate(self):
        for mean_zero, (cert, execution) in self.results.items():
            for row in execution.attempts:
                if row.get('reason_code') != 'step_candidate_screened' or not row['screen']['strict_temple_gap']:
                    continue
                p = row['proposal']
                rebuilt = extension.step.trial_certificate(cert['function'], p['basis'], p['coefficients'],
                    cert['source'], cert['gap'], cert['tolerance'], mean_zero)
                for key in ('lower', 'upper', 'exact_width'):
                    self.assertEqual(row['screen'][key], rebuilt[key])

    def test_tight_target_keeps_open_and_exact_original_tolerance(self):
        t = task(False, '1/1000000000000000000000000000000')
        execution = Execution()
        with patch.object(extension, 'DEGREES', (2,)):
            cert = extension.search(t, execution)
        self.assertIsNotNone(cert)
        judgment = baseline.assess(cert, baseline.normalize_task(t))
        self.assertTrue(judgment['certificate_valid'])
        self.assertFalse(judgment['target_met'])
        self.assertEqual(cert['tolerance'], t['tolerance'])
        self.assertTrue(any(row.get('reason_code') == 'step_degree_schedule_exhausted'
                            for row in execution.attempts))

    def test_tighter_precision_really_advances_the_degree(self):
        t = task(False, '1/10000000000')
        execution = Execution()
        cert = extension.search(t, execution)
        self.assertIsNotNone(cert)
        judgment = baseline.assess(cert, baseline.normalize_task(t))
        self.assertTrue(judgment['certificate_valid'])
        self.assertTrue(judgment['target_met'])
        self.assertGreater(execution.counters['candidate_proposals'],
                           self.results[False][1].counters['candidate_proposals'])
        self.assertLessEqual(execution.counters['candidate_proposals'], len(extension.DEGREES))

    def test_task_binding_and_tamper_rejected_by_frozen_verifier(self):
        cert = self.results[False][0]
        for field, value in [('mean_zero', True), ('tolerance', '1/100000')]:
            changed = dict(task(), **{field: value})
            self.assertFalse(baseline.assess(cert, baseline.normalize_task(changed))['certificate_valid'])
        changed = task()
        changed['function']['axis'] = 2
        self.assertFalse(baseline.assess(cert, baseline.normalize_task(changed))['certificate_valid'])
        for field in ('lower', 'exact_width'):
            bad = deepcopy(cert)
            bad[field] = '0'
            self.assertFalse(baseline.assess(bad, baseline.normalize_task(task()))['certificate_valid'])

    def test_strict_gap_and_other_sector_diagnostics_are_exact(self):
        cert = self.results[False][0]
        bad_stats = dict(cert['statistics'], rayleigh=cert['gap']['lower'])
        screen = extension._screen(cert['function'], False, cert['source'], cert['gap'], bad_stats, F(1, 10000))
        self.assertFalse(screen['strict_temple_gap'])
        source = deepcopy(cert['source'])
        source['angular_tail_lower'] = '-100'
        screen = extension._screen(cert['function'], False, source, cert['gap'], cert['statistics'], F(1, 10**20))
        self.assertTrue(screen['angular_tail_limits'])
        self.assertFalse(screen['target_met'])

    def test_failed_candidate_preserved_and_never_assessed_locally(self):
        execution = Execution()
        with patch.object(extension.step, 'propose_trial', side_effect=ValueError('synthetic proposal failure')), \
             patch.object(baseline, 'assess', side_effect=AssertionError('controller owns replay')):
            self.assertIsNone(extension.search(task(), execution))
        rows = [r for r in execution.attempts if r.get('reason_code') == 'step_candidate_failed']
        self.assertEqual(len(rows), len(extension.DEGREES))
        self.assertTrue(all(r['reason'] == 'synthetic proposal failure' for r in rows))

    def test_temple_failure_has_explicit_reason_and_retains_residual(self):
        def inadmissible_gap_candidate(q, degree, **kwargs):
            terms = extension.step.basis(degree)
            # The valid trial r=z has a Rayleigh value above the pointwise beta.
            return {'basis': list(map(list, terms)),
                    'coefficients': ['0', '1']+['0']*(len(terms)-2)}
        execution = Execution()
        with patch.object(extension.step, 'propose_trial', side_effect=inadmissible_gap_candidate):
            self.assertIsNone(extension.search(task(), execution))
        reasons = [r for r in execution.attempts if r.get('reason_code') == 'step_temple_gap_not_strict']
        self.assertEqual(len(reasons), len(extension.DEGREES))
        screened = [r for r in execution.attempts if r.get('reason_code') == 'step_candidate_screened']
        self.assertEqual(len(screened), len(reasons))
        self.assertTrue(all('operator_norm_squared' in r['statistics'] for r in screened))

    def test_stage_deadline_preserves_best_without_more_search(self):
        class StageDeadline(Execution):
            def check(self):
                super().check()
                if not self.reserve_released and any(
                        r.get('reason_code') == 'step_candidate_screened' for r in self.attempts):
                    raise Deadline('deterministic stage boundary')
        t = task(False, '1/10000000000')
        execution = StageDeadline()
        cert = extension.search(t, execution)
        self.assertIsNotNone(cert)
        self.assertTrue(execution.reserve_released)
        self.assertEqual(execution.counters['candidate_proposals'], 1)
        self.assertEqual(execution.counters['step_certificates_assembled'], 1)
        judgment = baseline.assess(cert, baseline.normalize_task(t))
        self.assertTrue(judgment['certificate_valid'])
        self.assertFalse(judgment['target_met'])
        row = next(r for r in execution.attempts if r.get('reason_code') == 'step_stage_budget')
        self.assertEqual(row['reason'], 'deterministic stage boundary')
        self.assertEqual(row['selected_degree'], 2)
        self.assertFalse(any(r.get('reason_code') == 'step_degree_schedule_exhausted'
                             for r in execution.attempts))

    def test_deadline_at_assembly_uses_reserve_once(self):
        class AssemblyDeadline(Execution):
            def run(self, action, function, **details):
                if action == 'step_assemble_frozen_certificate' and not self.reserve_released:
                    raise Deadline('deterministic assembly boundary')
                return super().run(action, function, **details)
        execution = AssemblyDeadline()
        cert = extension.search(task(), execution)
        self.assertTrue(execution.reserve_released)
        self.assertEqual(execution.counters['step_certificates_assembled'], 1)
        self.assertTrue(baseline.assess(cert, baseline.normalize_task(task()))['certificate_valid'])
        self.assertTrue(any(r.get('reason_code') == 'step_assembly_uses_reserve'
                            for r in execution.attempts))

    def test_unsupported_and_deadline_propagate_without_fallback(self):
        t = task()
        t['kind'] = 'approximation'
        del t['mean_zero']
        execution = Execution()
        with patch.object(baseline, 'solve', side_effect=AssertionError('no local fallback')):
            self.assertIsNone(extension.search(t, execution))
            self.assertEqual(execution.attempts[-1]['reason_code'], 'step_unsupported_task')
            with self.assertRaises(Deadline):
                extension.search(task(), Execution(0))


if __name__ == '__main__':
    unittest.main()
