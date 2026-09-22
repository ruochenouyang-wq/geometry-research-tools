"""Mock fallback boundaries plus exactly two predeclared public examples."""
from dataclasses import FrozenInstanceError
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('c3_fast_under_test',ROOT/'fast_candidate.py')
fast=importlib.util.module_from_spec(spec);spec.loader.exec_module(fast)
Q={'kind':'axis_profile','profile':'abs_power','axis':0,'amplitude':'-2/3','offset':'0','exponent':'-2/7'}
SS=tuple(map(F,(0,2,3,4,5,6)))


class FakeFull:
    canonical=staticmethod(lambda q:json.dumps(q,sort_keys=True,separators=(',',':'),allow_nan=False))
    normalize=staticmethod(lambda q:dict(q))
    powers=staticmethod(lambda ss:tuple(map(F,ss)))
    matvec=staticmethod(lambda a,x:[sum(v*z for v,z in zip(row,x)) for row in a])
    statistics=staticmethod(lambda *args:{'rayleigh':'0','residual_squared':'1/100'})
    @staticmethod
    def trial_matrices(q,ss):
        matrix=tuple(tuple(F(i==j) for j in range(len(ss))) for i in range(len(ss)))
        return matrix,matrix,matrix
    @staticmethod
    def cancellation_transform(q,ss):return FakeFull.trial_matrices(q,ss)[0]
    @staticmethod
    def congruence(matrix,transform):return matrix


def gap(q=None,beta='1'):
    return {'function':dict(Q if q is None else q),'beta':beta,'azimuth_m':0,'mean_zero':False}


def fake_candidate(*args):
    return {'powers':list(map(str,SS)),'coefficients':['1']*6,'precision_bits':53}


class FastCandidateControls(unittest.TestCase):
    def mock_try(self,**kwargs):
        np=SimpleNamespace(float64='float64',array=lambda x,**kw:x)
        with patch.object(fast,'_numpy',return_value=np),patch.object(fast,'_float_candidate',side_effect=fake_candidate):
            return fast.try_candidate(Q,SS,gap(),kwargs.pop('tolerance','1/10'),full=FakeFull,**kwargs)

    def test_00_two_fixed_public_examples(self):
        # First actual NumPy load of this test process is measured in E09/10.
        # Reuse public previously certified gaps; no new spectral search/tests.
        data=json.loads((ROOT.parent/'c_external_comparison_20260922'/'numerical'/'RESULTS.json').read_text())
        full=fast._resolve_full();diagnostics=[]
        for ident,count,accepted in (('E09',10,True),('E19',6,False)):
            task=next(x['task'] for x in data['cases'] if x['id']==ident)
            context=next(x['context'] for x in data['common_contexts'] if x['case_id']==ident)
            ss=full.generated_powers(task['function'],count)
            result=fast.try_candidate(task['function'],ss,context['gap'],task['tolerance'],full=full)
            self.assertEqual(result['accepted'],accepted,result)
            self.assertEqual(result['fallback_required'],not accepted)
            candidate=result['candidate'];self.assertIsNotNone(candidate)
            rebuilt=full.statistics(task['function'],ss,candidate['coefficients'])
            self.assertEqual(candidate['statistics'],rebuilt)
            cert=full.certificate(task['function'],ss,candidate['coefficients'],context['gap'],task['tolerance'])
            self.assertTrue(full.verify(cert,task['function'],False,task['tolerance']))
            self.assertEqual(candidate['screening_width'],cert['exact_width'])
            self.assertNotIn('certificate_valid',result)
            diagnostics.append({'case_id':ident,'powers_count':count,'accepted':result['accepted'],
                'reason':result['reason'],'exact_stats_match_rebuild':True,'old_full_verify':True,
                'exact_width':cert['exact_width'],'numpy_preloaded_at_entry':result['numpy_preloaded_at_entry'],
                'phases':result['phases'],'elapsed_seconds':result['elapsed_seconds'],'cpu_seconds':result['cpu_seconds']})
        print('PUBLIC_FUNCTIONAL_CHECKS='+json.dumps(diagnostics,allow_nan=False))

    def test_prepare_is_exact_only_and_immutable(self):
        with patch.object(fast,'_numpy',side_effect=AssertionError('must not import')):
            p=fast.prepare_trial(Q,SS,full=FakeFull)
        self.assertEqual(p.function_key,FakeFull.canonical(Q))
        with self.assertRaises(FrozenInstanceError):p.powers=()
        self.assertIsInstance(p.M,tuple);self.assertIsInstance(p.M[0],tuple)
        self.assertIsInstance(p.M[0][0],F)

    def test_prepared_function_binding(self):
        wrong=dict(Q,axis=1);p=fast.prepare_trial(wrong,SS,full=FakeFull)
        result=self.mock_try(prepared=p)
        self.assertFalse(result['accepted']);self.assertEqual(result['reason'],'exact_prepare_failed')
        self.assertNotIn('numpy_import',result['phases'])

    def test_prepared_powers_binding(self):
        p=fast.prepare_trial(Q,(F(0),F(2),F(3),F(4),F(5),F(7)),full=FakeFull)
        self.assertEqual(self.mock_try(prepared=p)['reason'],'exact_prepare_failed')

    def test_gap_context_binding(self):
        for bad in (gap(dict(Q,axis=1)),dict(gap(),mean_zero=True),dict(gap(),azimuth_m=1)):
            result=fast.try_candidate(Q,SS,bad,'1/10',full=FakeFull)
            self.assertFalse(result['accepted']);self.assertNotIn('numpy_import',result['phases'])

    def test_sixteen_terms_skip_without_numpy_or_full_load(self):
        with patch.object(fast,'_resolve_full',side_effect=AssertionError('skip must be cheap')):
            result=fast.try_candidate(Q,list(range(16)),gap(),'1/10')
        self.assertEqual(result['reason'],'unsupported_power_count')

    def test_import_failure_requests_fallback(self):
        with patch.object(fast,'_numpy',side_effect=ImportError('NumPy unavailable')):
            result=fast.try_candidate(Q,SS,gap(),'1/10',full=FakeFull)
        self.assertEqual(result['reason'],'numpy_unavailable');self.assertTrue(result['fallback_required'])
        self.assertEqual(result['phases']['numpy_import']['status'],'failed')

    def test_cholesky_failure_requests_fallback(self):
        def failed(np,p,arrays,state,full):
            state['substage']='cholesky';raise ArithmeticError('not positive definite')
        np=SimpleNamespace(float64='float64',array=lambda x,**kw:x)
        with patch.object(fast,'_numpy',return_value=np),patch.object(fast,'_float_candidate',side_effect=failed):
            result=fast.try_candidate(Q,SS,gap(),'1/10',full=FakeFull)
        self.assertEqual(result['reason'],'cholesky_failed');self.assertTrue(result['fallback_required'])

    def test_precise_tolerance_failure_retains_trial(self):
        result=self.mock_try(tolerance='1/1000')
        self.assertFalse(result['accepted']);self.assertTrue(result['fallback_required'])
        self.assertEqual(result['reason'],'exact_width_exceeds_tolerance')
        self.assertEqual(result['candidate']['statistics']['residual_squared'],'1/100')
        self.assertEqual(result['candidate']['screening_width'],'1/100')

    def test_mu_equals_beta_never_accepted_even_zero_residual(self):
        with patch.object(FakeFull,'statistics',return_value={'rayleigh':'1','residual_squared':'0'}):
            result=self.mock_try()
        self.assertFalse(result['accepted']);self.assertEqual(result['reason'],'mu_not_below_screening_beta')

    def test_negative_complete_residual_failure(self):
        with patch.object(FakeFull,'statistics',return_value={'rayleigh':'0','residual_squared':'-1'}):
            result=self.mock_try()
        self.assertFalse(result['accepted']);self.assertEqual(result['reason'],'exact_statistics_failed')

    def test_zero_budget_does_not_import(self):
        with patch.object(fast,'_numpy',side_effect=AssertionError('budget already exhausted')):
            result=fast.try_candidate(Q,SS,gap(),'1/10',full=FakeFull,deadline=0.0)
        self.assertEqual(result['reason'],'deadline_exceeded');self.assertEqual(result['phases'],{})

    def test_deadline_after_stats_retains_candidate(self):
        checks=[]
        def check(deadline):
            checks.append(None)
            if len(checks)==8:raise fast.CandidateDeadline('mock elapsed after stats')
        with patch.object(fast,'_check_deadline',side_effect=check):result=self.mock_try()
        self.assertEqual(result['reason'],'deadline_exceeded')
        self.assertIn('statistics',result['candidate']);self.assertFalse(result['accepted'])

    def test_no_false_final_validity_claim(self):
        result=self.mock_try();self.assertTrue(result['accepted'])
        self.assertTrue(result['screening_only']);self.assertFalse(result['gap_proof_checked_here'])
        self.assertNotIn('certificate_valid',result)

    def test_float_matrix_inputs_rejected(self):
        raw=[[[float(i==j) for j in range(6)] for i in range(6)]]*3
        with self.assertRaises(TypeError):fast.prepare_exact(Q,SS,full=FakeFull,matrices=raw)

    def test_nonfinite_or_zero_float_vectors_rejected(self):
        import numpy as np
        p=fast.prepare_trial(Q,SS,full=FakeFull);arrays=(np.eye(6),np.eye(6))
        for values,vectors in ((np.full(6,np.nan),np.eye(6)),(np.ones(6),np.zeros((6,6)))):
            with patch.object(np.linalg,'eigh',return_value=(values,vectors)):
                with self.assertRaises(fast.InvalidFloatCandidate):
                    fast._float_candidate(np,p,arrays,{},FakeFull)


if __name__=='__main__':unittest.main()
