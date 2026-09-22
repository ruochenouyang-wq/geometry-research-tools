"""C3 acceptance, fallback, and public-task integration controls.

Functional checks may use a larger budget during concurrent development;
the frozen evaluator separately imposes the original ten-second deadline.
"""
import json
from pathlib import Path
import unittest
from unittest.mock import patch

import controller
import fast_candidate
import strong_search
import weak_search

OLD_CASES = Path(__file__).resolve().parent.parent/'c_demand_extension_20260922'/'cases.json'


def task(case_id, budget=30):
    data = json.loads(OLD_CASES.read_text())
    cases = data['cases'] if isinstance(data, dict) else data
    found = next(row for row in cases if row['id'] == case_id)
    value = json.loads(json.dumps(found['task']))
    value['budget'] = {'wall_seconds':budget}
    return value


class AcceptanceControls(unittest.TestCase):
    def test_unverified_assembly_cannot_claim_success(self):
        raw = task('E09')
        forged = {'format':'fullspace_m0_fractional_temple_v1', 'status':'target_met',
                  'certificate_valid':True, 'target_met':True, 'exact_width':'0'}
        with patch.object(weak_search,'search',return_value=forged):
            result = controller.solve(raw)
        self.assertFalse(result['certificate_valid'])
        self.assertFalse(result['target_met'])
        self.assertEqual(result['failure_stage'],'final_frozen_replay')
        self.assertEqual(result['counters']['final_frozen_replays'],1)

    def test_fast_rejection_preserves_rational_fallback_and_reason(self):
        raw = task('E19')
        run = strong_search.Run(30)
        full = strong_search._modules()[1]
        powers = full.generated_powers(raw['function'],6)
        failure = {'accepted':False,'reason':'numpy_unavailable','candidate':None}
        marker = {'candidate_kind':'old_rational'}
        with patch.object(fast_candidate,'try_candidate',return_value=failure), \
             patch.object(strong_search,'_run_rational_candidate',return_value=marker) as fallback:
            result = strong_search._run_candidate(raw['function'],powers,{}, {},'1/100000000',run)
        self.assertIs(result,marker)
        fallback.assert_called_once()
        self.assertEqual(run.counts['rational_after_fast'],1)
        self.assertTrue(any(row.get('reason')=='numpy_unavailable' for row in run.attempts))

    def test_sixteen_powers_go_directly_to_rational_path(self):
        raw = task('E19')
        powers = strong_search._modules()[1].generated_powers(raw['function'],16)
        with patch.object(fast_candidate,'try_candidate',side_effect=AssertionError('Unexpected float call')), \
             patch.object(strong_search,'_run_rational_candidate',return_value={}) as fallback:
            strong_search._run_candidate(raw['function'],powers,{}, {},'1/100000000',strong_search.Run(30))
        fallback.assert_called_once()

    def test_method_boundary_remains_explicit(self):
        result = controller.solve(task('E23'))
        self.assertFalse(result['target_met'])
        self.assertIsNone(result['certificate'])
        self.assertEqual(result['status'],'unsupported_method_domain')


class PublicIntegration(unittest.TestCase):
    def check_original(self, case_id):
        raw = task(case_id)
        result = controller.solve(raw)
        self.assertTrue(result['certificate_valid'],result.get('error',result.get('status')))
        self.assertTrue(result['target_met'],result.get('status'))
        assessment = controller.verify(result['certificate'],raw)
        self.assertTrue(assessment['certificate_valid'])
        self.assertTrue(assessment['target_met'])
        return result

    def test_weak_original_target_and_no_numpy_import(self):
        with patch.object(fast_candidate,'_numpy',side_effect=AssertionError('Weak route imported NumPy')):
            result = self.check_original('E09')
        self.assertEqual(result['native_kind'],'weak_fullspace_deferred_assembly')
        self.assertGreater(result['counters']['deferred_assemblies'],0)

    def test_strong_original_target(self):
        result = self.check_original('E19')
        self.assertEqual(result['native_kind'],'hybrid_fullspace_deferred_assembly')
        self.assertGreater(result['counters']['fast_candidate_attempts'],0)

    def test_mean_zero_active_block_original_target(self):
        result = self.check_original('E18')
        self.assertEqual(result['native_kind'],'active_mean_zero_sector')
        self.assertIs(result['certificate']['mean_zero'],True)
        selected = [a['active_m'] for a in result['attempts']
                    if a.get('reason_code')=='mean_zero_refine_active_sector']
        self.assertIn(0,selected)


if __name__ == '__main__':
    unittest.main()
