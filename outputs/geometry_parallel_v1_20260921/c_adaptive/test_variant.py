"""Public, deterministic correctness checks; no held-out seeds or inputs."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import variant as c


def task(mean_zero=False, amplitude='-1', exponent='-1/4', tolerance='1/100000000'):
    return {'kind':'spectrum', 'function':{'kind':'axis_profile','profile':'abs_power',
             'axis':2,'amplitude':amplitude,'offset':'0','exponent':exponent},
            'mean_zero':mean_zero,'tolerance':tolerance,'budget':{'wall_seconds':10}}


class VariantTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.results = {mean_zero:c.solve(task(mean_zero)) for mean_zero in (False,True)}

    def test_two_spaces_frozen_replay_and_target(self):
        for mean_zero,result in self.results.items():
            self.assertTrue(result['target_met'], result)
            self.assertTrue(c.verify(result['certificate'],task(mean_zero))['target_met'])
            self.assertEqual(result['counters']['baseline_fallback_calls'],0)
            self.assertEqual(result['counters']['final_frozen_replays'],1)
            self.assertLessEqual(F(result['metric']),F(task()['tolerance']))

    def test_full_residual_exact_threshold_and_tail(self):
        for mean_zero,result in self.results.items():
            cert=result['certificate'];stats=cert['statistics'];gap=cert['gap']
            mu,variance,beta=map(F,(stats['rayleigh'],stats['residual_squared'],gap['beta']))
            self.assertGreater(beta,mu)
            self.assertLess(beta,F(gap['radial_tail_lower']))
            self.assertLessEqual(variance/(beta-mu),F(task()['tolerance']))
            if mean_zero:
                self.assertGreaterEqual(F(cert['m0_source']['lower']),mu-F(task()['tolerance']))

    def test_binding_and_tampering_rejected(self):
        for mz,result in self.results.items():
            cert=result['certificate'];original=task(mz)
            for change in ('space','function','tolerance'):
                wrong=deepcopy(original)
                if change=='space':wrong['mean_zero']=not mz
                elif change=='function':wrong['function']['amplitude']='-2'
                else:wrong['tolerance']='1/1000000'
                self.assertFalse(c.verify(cert,wrong)['certificate_valid'])
            for key in ('beta','lower','radial_tail_lower'):
                wrong=deepcopy(cert);wrong['gap'][key]='999'
                self.assertFalse(c.verify(wrong,original)['certificate_valid'])
            wrong=deepcopy(cert);wrong['statistics']['residual_squared']='0'
            self.assertFalse(c.verify(wrong,original)['certificate_valid'])

    def test_zero_budget_starts_no_route(self):
        original=task();original['budget']['wall_seconds']=0
        with patch.object(c,'_fullspace',side_effect=AssertionError('must not start')):
            result=c.solve(original)
        self.assertEqual(result['status'],'budget_exceeded')
        self.assertFalse(result['target_met']);self.assertEqual(result['attempts'],[])
        self.assertEqual(result['counters']['candidate_proposals'],0)

    def test_expired_cooperative_budget_has_no_fallback(self):
        original=task();original['budget']['wall_seconds']=0.001
        times=iter((0.0,0.0,0.002,0.003))
        with patch.object(c,'perf_counter',side_effect=lambda:next(times)), \
             patch.object(c,'_fullspace',side_effect=c.Deadline('expired')):
            result=c.solve(original)
        self.assertFalse(result['target_met'])
        self.assertEqual(result['counters']['baseline_fallback_calls'],0)

    def test_fallback_receives_remaining_budget(self):
        original=task();seen=[]
        def fallback(t):
            seen.append(t)
            return {'certificate':None,'status':'no_verified_certificate','attempts':[]}
        with patch.object(c,'_fullspace',return_value=None), \
             patch.object(c.baseline,'solve',side_effect=fallback):
            result=c.solve(original)
        self.assertEqual(len(seen),1)
        self.assertLess(seen[0]['budget']['wall_seconds'],10)
        self.assertEqual(seen[0]['function'],original['function'])
        self.assertEqual(seen[0]['tolerance'],original['tolerance'])
        self.assertEqual(result['counters']['baseline_fallback_calls'],1)

    def test_zero_residual_strict_gap_and_failed_check_count(self):
        import spectral_gap
        original=task();q=c.baseline.normalize_task(original)['function']
        from time import perf_counter
        execution=c.Execution(perf_counter(),10)
        k=c._kernel(q,0,False,4,4,None,None,execution)
        proof=c._needed_gap(k,F(-3),F(0),F(1,1000),None,execution)
        self.assertIsNotNone(proof)
        self.assertGreater(F(proof['beta']),F(-3))
        before=execution.counters['search_inertia_checks']
        with self.assertRaises(ValueError):c._gap_at(k,k.beta,None,execution)
        self.assertEqual(execution.counters['search_inertia_checks'],before)
        self.assertEqual(execution.attempts[-1]['status'],'failed')
        self.assertTrue(spectral_gap.verify_gap(proof,q,0,False))

    def test_unchanged_step_and_approximation_routes(self):
        examples=[{'kind':'approximation','function':task()['function'],'tolerance':'1/1000000'},
                  {'kind':'spectrum','function':{'kind':'axis_profile','profile':'step',
                    'amplitude':'-1','offset':'0'},'mean_zero':True,'tolerance':'1/1000000'}]
        for original in examples:
            result=c.solve(original)
            self.assertTrue(result['target_met'],result)
            self.assertEqual(result['route'],'unchanged_baseline_route')
            self.assertTrue(c.verify(result['certificate'],original)['target_met'])

    def test_strong_fullspace_gate_preserves_baseline_budget(self):
        import fullspace_singular
        original=task(False,amplitude='-3',exponent='-2/5')
        with patch.object(fullspace_singular,'proposal',side_effect=AssertionError('C must skip proposals')), \
             patch.object(c.baseline,'solve',return_value={
                 'certificate':None,'certificate_valid':False,'target_met':False,
                 'status':'no_verified_certificate','attempts':[]}):
            result=c.solve(original)
        self.assertEqual(result['counters']['candidate_proposals'],0)
        self.assertEqual(result['counters']['kernel_constructions'],0)
        self.assertEqual(result['counters']['adaptive_policy_skips'],1)
        self.assertEqual(result['counters']['baseline_fallback_calls'],1)


if __name__=='__main__':unittest.main()
