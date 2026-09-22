"""Small semantic tests; no timing comparison or private evaluation inputs."""
from copy import deepcopy
from fractions import Fraction as F
import os
from pathlib import Path
import unittest
from unittest.mock import patch
import a_adapter as adapter


def task(offset='-6',mean_zero=False,tolerance='1/1000000',wall=10):
    return {'kind':'spectrum','function':{'kind':'axis_profile','profile':'step','axis':2,
            'amplitude':'0','offset':offset},'mean_zero':mean_zero,'tolerance':tolerance,
            'budget':{'wall_seconds':wall}}


class AdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine=adapter.Engine()
        cls.t=task()
        with patch.object(cls.engine.source.shared,'solve',side_effect=AssertionError('No baseline.solve')):
            cls.result=cls.engine.solve(cls.t)
        cls.cert=cls.result['certificate']

    def test_native_ground_projection_and_replay(self):
        self.assertTrue(self.result['certificate_valid']);self.assertTrue(self.result['target_met'])
        self.assertLessEqual(F(self.cert['lower']),-6);self.assertGreaterEqual(F(self.cert['upper']),-6)
        self.assertFalse(self.result['fallback']);self.assertFalse(self.cert['baseline_fallback_used'])
        self.assertTrue(self.engine.verify(self.cert,self.t)['target_met'])
        self.assertLessEqual(self.result['source_run_task']['budget']['wall_seconds'],10)
        self.assertIsNotNone(self.result['source_response'])
        self.assertTrue(self.result['metadata']['adapter_solve_extra_full_source_replay'])

    def test_wrong_task_and_projection_mutations(self):
        for key,value in (('mean_zero',True),('tolerance','1/1000'),('kind','approximation')):
            t=deepcopy(self.t);t[key]=value
            self.assertFalse(self.engine.verify(self.cert,t)['certificate_valid'])
        t=deepcopy(self.t);t['function']['axis']=1
        self.assertFalse(self.engine.verify(self.cert,t)['certificate_valid'])
        for edit in (lambda c:c.update(lower='-7'),
                     lambda c:c.update(eligible_sectors=[0]),
                     lambda c:c['source_certificate']['sectors'].pop(0),
                     lambda c:c['source_task'].update(tolerance='1/1000'),
                     lambda c:c['source_certificate'].update(mean_zero=True)):
            bad=deepcopy(self.cert);edit(bad)
            self.assertFalse(self.engine.verify(bad,self.t)['certificate_valid'])

    def test_source_sum_target_false_can_still_prove_ground_target(self):
        source=deepcopy(self.cert['source_certificate'])
        ground_width=F(self.cert['exact_width']);sum_width=F(source['sum_width'])
        self.assertLess(ground_width,sum_width)
        tol=(ground_width+sum_width)/2
        t=task(tolerance=str(tol));source.update(tolerance=str(tol),sum_target_met=False,target_met=False)
        native={'certificate':source,'target_met':False,'status':'count_certified_sum_open',
                'cost_note':'retained native failure metadata'}
        with patch.object(self.engine.source,'solve',return_value=native):
            r=self.engine.solve(t)
        self.assertTrue(r['certificate_valid']);self.assertTrue(r['target_met'])
        self.assertFalse(r['source_sum_target_met']);self.assertEqual(r['source_response'],native)
        self.assertTrue(self.engine.verify(r['certificate'],t)['certificate_valid'])

    def test_zero_negative_count_does_not_claim_zero_ground(self):
        r=self.engine.solve(task(offset='1'))
        self.assertIsNone(r['certificate']);self.assertFalse(r['target_met'])
        self.assertEqual(r['source_response']['negative_count'],0)
        self.assertIn('no_certified_strictly_negative',r['reason'])
        self.assertEqual(r['status'],'unsupported')

    def test_unresolved_count_is_preserved_and_rejected(self):
        t=task(offset='-3/2',mean_zero=True);t['function']['amplitude']='-1'
        source_task=self.engine._source_task(self.engine._task(t))
        native=self.engine.source.negative_spectrum(source_task,modes=(1,),near_tail=0,bits=8)
        self.assertEqual((native['count_lower'],native['count_upper']),(0,3))
        with patch.object(self.engine.source,'solve',return_value=native):r=self.engine.solve(t)
        self.assertIsNone(r['certificate']);self.assertFalse(r['target_met'])
        self.assertEqual(r['reason'],'negative_count_unresolved')
        self.assertEqual(r['status'],'no_verified_ground_certificate')
        self.assertEqual(r['source_response'],native)

    def test_count_only_source_cannot_be_projected(self):
        src=deepcopy(self.cert['source_certificate']);rows=src['sectors']
        for row in rows:
            if row['count_lower']:row['eigenvalue_intervals']=None
        goal=self.cert['source_task']
        src=self.engine.source.assemble(goal,src['form_bound'],src['form_certificate'],rows)
        with patch.object(self.engine.source,'solve',return_value={'certificate':src}):r=self.engine.solve(self.t)
        self.assertIsNone(r['certificate']);self.assertIn('intervals_missing',r['reason'])
        self.assertEqual(r['status'],'no_verified_ground_certificate')

    def test_missing_source_is_inconclusive_not_unsupported(self):
        native={'certificate':None,'status':'stage_failed','attempts':[{'reason':'search failed'}]}
        with patch.object(self.engine.source,'solve',return_value=native):r=self.engine.solve(self.t)
        self.assertEqual(r['status'],'no_verified_ground_certificate')
        self.assertEqual(r['source_response'],native)

    def test_source_replay_failure_is_explicit(self):
        source=deepcopy(self.cert['source_certificate']);source['negative_count']=999
        native={'certificate':source,'target_met':True}
        with patch.object(self.engine.source,'solve',return_value=native):r=self.engine.solve(self.t)
        self.assertEqual(r['status'],'verification_failed')
        self.assertEqual(r['reason'],'source_certificate_failed_full_replay')
        self.assertFalse(r['source_assessment']['certificate_valid'])
        self.assertEqual(r['source_response'],native)

    def test_mean_zero_has_its_own_ground(self):
        t=task(mean_zero=True);r=self.engine.solve(t)
        self.assertTrue(r['certificate_valid'])
        self.assertLessEqual(F(r['lower']),-4);self.assertGreaterEqual(F(r['upper']),-4)
        self.assertFalse(self.engine.verify(r['certificate'],self.t)['certificate_valid'])

    def test_every_verification_replays_even_full_false(self):
        original=self.engine.source.verify
        with patch.object(self.engine.source,'verify',wraps=original) as replay,\
             patch.object(self.engine.source,'solve',side_effect=AssertionError('no search')):
            self.assertTrue(self.engine.verify(self.cert,self.t,full=True)['certificate_valid'])
            self.assertTrue(self.engine.verify(self.cert,self.t,full=False)['certificate_valid'])
            self.assertEqual(replay.call_count,2)

    def test_reduced_execution_budget_is_not_a_different_math_task(self):
        cert=deepcopy(self.cert);cert['source_task']['budget']['wall_seconds']=9.0
        self.assertTrue(self.engine.verify(cert,self.t)['certificate_valid'])
        cert['source_task']['budget']['wall_seconds']=11.0
        self.assertFalse(self.engine.verify(cert,self.t)['certificate_valid'])

    def test_zero_budget_and_non_spectrum_never_call_solver(self):
        with patch.object(self.engine.source,'solve',side_effect=AssertionError('No computation')):
            r=self.engine.solve(task(wall=0))
            self.assertFalse(r['target_met']);self.assertEqual(r['status'],'budget_exceeded')
            t=task();t['kind']='approximation'
            self.assertEqual(self.engine.solve(t)['status'],'unsupported')

    def test_cwd_independent_loading_and_no_preparation(self):
        old=os.getcwd()
        try:
            os.chdir('/')
            engine=adapter.Engine({'ignored_family':True})
            self.assertTrue(Path(engine.source_path).is_absolute())
            self.assertEqual(engine.prepare(2),{'ok':True,'no_preparation':True})
        finally:os.chdir(old)


if __name__=='__main__':unittest.main()
