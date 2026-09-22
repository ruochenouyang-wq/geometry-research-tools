"""Bounded public functional and failure tests; not performance evaluation."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent))
import fullspace_extension as e


def task():
    return {'kind':'spectrum','function':{'kind':'axis_profile','profile':'abs_power',
        'axis':1,'amplitude':'-7/3','offset':'1/5','exponent':'-3/8'},
        'mean_zero':False,'tolerance':'1/100000000','budget':{'wall_seconds':10}}


class FullspaceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.events=[]
        cls.result=e.solve(task(),wall_seconds=7,progress=cls.events.append)
        cls.cert=cls.result['certificate']

    def test_frozen_replay_and_complete_residual(self):
        self.assertEqual(self.result['status'],'candidate_ready',self.result)
        self.assertIsNone(self.result['certificate_valid'])
        self.assertFalse(self.result['target_met'])
        self.assertTrue(self.result['candidate_target_met'])
        self.assertTrue(e.verify(self.cert,task())['target_met'])
        mu=F(self.cert['statistics']['rayleigh']);v=F(self.cert['statistics']['residual_squared'])
        beta=F(self.cert['gap']['beta'])
        self.assertLess(mu,beta)
        self.assertLess(beta,F(self.cert['gap']['radial_tail_lower']))
        self.assertEqual(F(self.cert['exact_width']),v/(beta-mu))
        self.assertLessEqual(F(self.cert['exact_width']),F(task()['tolerance']))

    def test_binding_and_tampering(self):
        for field,value in (('mean_zero',True),('tolerance','1/1000000')):
            wrong=task();wrong[field]=value
            self.assertFalse(e.verify(self.cert,wrong)['certificate_valid'])
        wrong=task();wrong['function']['axis']=0
        self.assertFalse(e.verify(self.cert,wrong)['certificate_valid'])
        for section,field in (('gap','beta'),('gap','radial_tail_lower'),
                              ('statistics','residual_squared')):
            bad=deepcopy(self.cert);bad[section][field]='999'
            self.assertFalse(e.verify(bad,task())['certificate_valid'])

    def test_progress_has_begin_end_and_failures(self):
        self.assertEqual(self.events[0]['event'],'begin')
        self.assertEqual(self.events[-1]['action'],'component_complete')
        self.assertTrue(any(x['event']=='error' and x['action']=='demand_gap_threshold' for x in self.events))
        self.assertTrue(any(x['event']=='end' and x['action']=='bounded_candidate' for x in self.events))
        self.assertTrue(all(x['component']=='fullspace_extension' for x in self.events))

    def test_zero_budget_starts_no_math(self):
        with patch.object(e,'_modules',side_effect=AssertionError('must not import math')):
            result=e.solve(task(),wall_seconds=0)
        self.assertFalse(result['target_met'])
        self.assertEqual(result['counters']['candidate_processes'],0)
        self.assertTrue(result['fallback_recommended'])
        self.assertIn('allocation exhausted',result['fallback_reason'])

    def test_worker_timeout_is_hard_and_recorded(self):
        run=e.Run(4)
        with patch.object(e.subprocess,'run',side_effect=subprocess.TimeoutExpired('candidate',2)) as child:
            with self.assertRaises(e.SearchBudget):
                run.action('bounded_candidate',lambda:e._run_candidate({},[F(0)],{}, {},F(1,100),run))
        self.assertEqual(run.counts['candidate_timeouts'],1)
        self.assertEqual(run.attempts[-1]['status'],'failed')
        self.assertGreater(run.remaining(),3)
        self.assertLessEqual(child.call_args.kwargs['timeout'],2)

    def test_worker_rejects_unproved_gap(self):
        payload={'function':task()['function'],'powers':['0','13/8'],
                 'gap':deepcopy(self.cert['gap']),'shift_source':{},'tolerance':task()['tolerance']}
        payload['gap']['beta']='999'
        with self.assertRaisesRegex(ValueError,'proved second-eigenvalue'):
            e._candidate_worker(payload)

    def test_no_hidden_baseline_fallback(self):
        baseline=e._modules()[0]
        with patch.object(baseline,'solve',side_effect=AssertionError('component must not fallback')):
            wrong=task();wrong['mean_zero']=True
            result=e.solve(wrong,wall_seconds=1)
        self.assertEqual(result['status'],'extension_failed')
        self.assertTrue(result['fallback_recommended'])
        self.assertIn('mean_zero=false',result['fallback_reason'])

    def test_budget_not_reserved_twice_and_invalid_budget(self):
        run=e.Run(10)
        self.assertEqual(run.reserve,0)
        self.assertGreater(run.available(),9.9)
        for bad in (True,-1,float('inf')):
            with self.assertRaises(ValueError):e.solve(task(),wall_seconds=bad)

    def test_budget_after_best_retains_complete_open_certificate(self):
        original=task();original['tolerance']='1/1000000000000000000000000000000'
        candidate={'powers':self.cert['trial']['powers'],
            'coefficients':self.cert['trial']['coefficients'],
            'statistics':self.cert['statistics'],'hard_timeout_seconds':0.5,'iterations':8}
        events=[]
        with patch.object(e,'_run_candidate',side_effect=[candidate,e.SearchBudget('second candidate allocation expired')]):
            result=e.solve(original,wall_seconds=7,progress=events.append)
        self.assertIsNotNone(result['certificate'])
        self.assertEqual(result['search_stop_reason'],'protected_budget_stop')
        self.assertEqual(result['counters']['budget_recovery_assemblies'],1)
        self.assertTrue(e.verify(result['certificate'],original)['certificate_valid'])
        self.assertFalse(e.verify(result['certificate'],original)['target_met'])
        self.assertTrue(any(x['action']=='recover_best_certificate' and
            x.get('detail',{}).get('permit_expired_allocation') for x in events))

    def test_failed_recovery_retains_rebuildable_trial(self):
        original=task();original['tolerance']='1/1000000000000000000000000000000'
        candidate={'powers':self.cert['trial']['powers'],
            'coefficients':self.cert['trial']['coefficients'],
            'statistics':self.cert['statistics'],'hard_timeout_seconds':0.5,'iterations':8}
        with patch.object(e,'_run_candidate',side_effect=[candidate,e.SearchBudget('candidate allocation expired')]), \
             patch.object(e,'_assemble_best',side_effect=ValueError('synthetic construction failure')):
            result=e.solve(original,wall_seconds=7)
        self.assertIsNone(result['certificate'])
        retained=result['best_candidate']
        self.assertEqual(retained['powers'],candidate['powers'])
        self.assertEqual(retained['coefficients'],candidate['coefficients'])
        self.assertEqual(retained['function'],original['function'])
        self.assertEqual(retained['tolerance'],original['tolerance'])
        self.assertTrue(retained['requires_original_certificate_construction_and_verification'])

    def test_budget_inside_threshold_keeps_current_open_trial(self):
        original=task();original['tolerance']='1/1000000000000000000000000000000'
        candidate={'powers':self.cert['trial']['powers'],
            'coefficients':self.cert['trial']['coefficients'],
            'statistics':self.cert['statistics'],'hard_timeout_seconds':0.5,'iterations':8}
        original_gap=e._gap_at;calls=[]
        def interrupt_threshold(*args):
            calls.append(1)
            if len(calls)>1:raise e.SearchBudget('budget reached during demand threshold')
            return original_gap(*args)
        with patch.object(e,'_run_candidate',return_value=candidate), \
             patch.object(e,'_gap_at',side_effect=interrupt_threshold):
            result=e.solve(original,wall_seconds=7)
        self.assertEqual(result['search_stop_reason'],'protected_budget_stop')
        self.assertIsNotNone(result['certificate'])
        self.assertTrue(e.verify(result['certificate'],original)['certificate_valid'])


class SchedulingTests(unittest.TestCase):
    """Pure mocks: never import a mathematics backend or start a child."""
    def _schedule(self,mu,threshold_succeeds=False):
        from types import SimpleNamespace
        seen=[]
        fake_baseline=SimpleNamespace(normalize_task=lambda value:deepcopy(value))
        fake_full=SimpleNamespace(normalize=lambda value:value,verify_form=lambda *a,**k:True,
            generated_powers=lambda q,n:list(range(n)))
        fake_sg=SimpleNamespace(_kernel=lambda q,m,mz,n,*args:SimpleNamespace(n=n,beta=F(100)))
        fake_forms=SimpleNamespace(select_form=lambda *args:{'mock_form':True})
        def coarse(kernel,form,run):
            return {'beta':'2','mock_modes':kernel.n},{'lower':'0'}
        def candidate(q,powers,gap,ground,tolerance,run,*,terminal=False):
            seen.append((gap['mock_modes'],len(powers)))
            return {'powers':['0'],'coefficients':['1'],'statistics':{
                'rayleigh':str(mu),'residual_squared':'1'},'hard_timeout_seconds':0.5,'iterations':4}
        def threshold(kernel,requested,form,run):
            if threshold_succeeds:return {'beta':str(requested),'mock_modes':kernel.n}
            raise ValueError('Synthetic threshold rejection')
        def assemble(best,q,tolerance,run,permit_expired=False):
            return {'format':'mock-only'},{'certificate_valid':None,'target_met':False,
                'candidate_target_met':threshold_succeeds,'verification_pending':True}
        with patch.object(e,'_modules',return_value=(fake_baseline,fake_full,fake_sg,fake_forms)), \
             patch.object(e,'_coarse_gap_and_shift',side_effect=coarse), \
             patch.object(e,'_run_candidate',side_effect=candidate), \
             patch.object(e,'_gap_at',side_effect=threshold), \
             patch.object(e,'_assemble_best',side_effect=assemble), \
             patch.object(e.subprocess,'run',side_effect=AssertionError('No mathematics child allowed')):
            result=e.solve(task(),wall_seconds=10)
        return result,seen

    def test_above_gap_prioritizes_larger_comparison(self):
        result,seen=self._schedule(mu=3)
        self.assertEqual(seen,[(12,6),(16,16)])
        self.assertEqual(result['counters']['gap_priority_upgrades'],1)
        decisions=[a for a in result['attempts'] if a['action']=='prioritize_stronger_gap']
        self.assertEqual(decisions[0]['reason'],'candidate_above_current_proved_gap_prioritize_stronger_gap')
        self.assertTrue(decisions[0]['not_a_proof_that_other_candidates_are_impossible'])

    def test_below_gap_preserves_basis_schedule(self):
        result,seen=self._schedule(mu=1)
        self.assertEqual(seen,[(12,6),(12,10),(12,14),(12,16),(16,16)])
        self.assertEqual(result['counters']['gap_priority_upgrades'],0)

    def test_successful_demand_check_precedes_schedule_decision(self):
        result,seen=self._schedule(mu=3,threshold_succeeds=True)
        self.assertEqual(seen,[(12,6)])
        self.assertEqual(result['counters']['gap_priority_upgrades'],0)
        self.assertTrue(result['candidate_target_met'])

    def test_final_candidate_does_not_reserve_for_nonexistent_later_trials(self):
        self.assertAlmostEqual(e._candidate_seconds(2.73),1.365)
        self.assertEqual(e._candidate_seconds(2.73,terminal=True),2.0)
        self.assertAlmostEqual(e._candidate_seconds(1,terminal=True),0.85)
        self.assertEqual(e._candidate_seconds(0.1,terminal=True),0)


if __name__=='__main__':unittest.main()
