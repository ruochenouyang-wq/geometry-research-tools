"""Control-flow regression tests; acceptance is always owned by verifier."""
import copy
from pathlib import Path
import json
import sys
import types
import unittest
from unittest.mock import patch
import controller

CASES = json.loads((Path(__file__).resolve().parent/'cases.json').read_text())
OK = {'certificate_valid':True,'target_met':True,'status':'certified_target_met',
      'metric':'1/1000000000000'}


class ControllerTests(unittest.TestCase):
    def task(self,index=0):
        return copy.deepcopy(CASES[index]['task'])

    def test_unsupported_is_not_mathematical_success(self):
        for index in (22,23):
            with patch.object(controller.baseline,'solve') as fallback:
                result=controller.solve(self.task(index))
            self.assertFalse(result['target_met'])
            self.assertEqual(result['status'],'unsupported_method_domain')
            fallback.assert_not_called()

    def test_no_budget_cannot_succeed(self):
        task=self.task();task['budget']['wall_seconds']=0
        with patch.object(controller.baseline,'solve') as fallback:
            result=controller.solve(task)
        fallback.assert_not_called()
        self.assertEqual(result['status'],'budget_exceeded')

    def test_candidate_claim_does_not_override_final_replay(self):
        fake=types.SimpleNamespace(search=lambda task,ex: {'untrusted':'candidate'})
        rejected={'certificate_valid':False,'target_met':False,'status':'invalid_certificate'}
        events=[]
        with patch.dict(sys.modules,{'step_extension':fake}), patch.object(controller,'verify',return_value=rejected):
            result=controller.solve(self.task(),progress=events.append)
        self.assertFalse(result['target_met'])
        self.assertEqual(result['failure_stage'],'final_frozen_replay')
        self.assertTrue(any(x.get('status')=='started' for x in events))

    def test_pending_fullspace_candidate_gets_real_verification(self):
        fake=types.SimpleNamespace(solve=lambda task,**kw:{'certificate':{'candidate':1},
            'target_met':False,'verification_pending':True,'attempts':[]})
        with patch.dict(sys.modules,{'fullspace_extension':fake}), patch.object(controller,'verify',return_value=OK):
            result=controller.solve(self.task(18))
        self.assertTrue(result['target_met'])
        self.assertEqual(result['route'],'native_c2')

    def test_recoverable_candidate_survives_fallback_without_becoming_proof(self):
        candidate={'powers':['0','8/5'],'coefficients':['1','1/4']}
        fake=types.SimpleNamespace(solve=lambda task,**kw:{'certificate':None,
            'best_candidate':candidate,'stop_reason':'protected_budget_stop','attempts':[]})
        failed={'certificate':None,'certificate_valid':False,'target_met':False,'status':'no_certificate'}
        with patch.dict(sys.modules,{'fullspace_extension':fake}), patch.object(controller.baseline,'solve',return_value=failed):
            result=controller.solve(self.task(18))
        self.assertFalse(result['target_met'])
        self.assertIsNone(result['certificate'])
        self.assertFalse(result['unverified_candidates'][0]['is_certificate'])
        self.assertEqual(result['unverified_candidates'][0]['candidate'],candidate)

    def test_fallback_has_remaining_budget_and_reason(self):
        captured=[]
        fake=types.SimpleNamespace(search=lambda task,ex:None)
        def fallback(task):
            captured.append(task)
            return {**OK,'certificate':{'candidate':2}}
        with patch.dict(sys.modules,{'step_extension':fake}), patch.object(controller.baseline,'solve',side_effect=fallback):
            result=controller.solve(self.task())
        self.assertEqual(result['route'],'baseline_fallback')
        self.assertEqual(result['fallback_reason'],'step_search_no_certificate')
        self.assertLessEqual(captured[0]['budget']['wall_seconds'],10)
        self.assertEqual(result['counters']['baseline_fallback_calls'],1)

    def test_valid_open_certificate_kept_without_duplicate_baseline(self):
        fake=types.SimpleNamespace(search=lambda task,ex:{'candidate':3})
        opened={**OK,'target_met':False,'metric':'1/10','status':'certified_open'}
        with patch.dict(sys.modules,{'step_extension':fake}), patch.object(controller,'verify',return_value=opened), patch.object(controller.baseline,'solve') as fallback:
            result=controller.solve(self.task())
        fallback.assert_not_called()
        self.assertTrue(result['certificate_valid'])
        self.assertFalse(result['target_met'])
        self.assertIsNotNone(result['certificate'])
        self.assertEqual(result['failure_stage'],'precision_gate')

    def test_unexpected_implementation_error_visible(self):
        def fail(task,ex):raise TypeError('synthetic implementation defect')
        fake=types.SimpleNamespace(search=fail)
        with patch.dict(sys.modules,{'step_extension':fake}):
            result=controller.solve(self.task())
        self.assertEqual(result['status'],'execution_failed')
        self.assertIn('synthetic implementation defect',result['error'])
        self.assertEqual(result['failure_stage'],'unexpected_execution_exception')

    def test_approximation_marked_intentional_baseline(self):
        with patch.object(controller.baseline,'solve',return_value={**OK,'certificate':{}}):
            result=controller.solve(self.task(20))
        self.assertEqual(result['route'],'baseline_fallback')
        self.assertEqual(result['fallback_reason'],'approximation_objective_not_extended')


if __name__=='__main__':unittest.main(verbosity=2)
