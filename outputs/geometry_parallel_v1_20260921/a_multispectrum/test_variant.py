"""Public exact calibrations and certificate attacks, without private data."""
from copy import deepcopy
from fractions import Fraction as F
from unittest.mock import patch
import unittest
import variant as v


def task(c='6',mean_zero=False,amplitude='0',tolerance='1/1000000'):
    return {'kind':'negative_spectrum','function':{'kind':'axis_profile','axis':2,
            'profile':'step','amplitude':amplitude,'offset':str(-F(c))},
            'mean_zero':mean_zero,'tolerance':tolerance,'budget':{'wall_seconds':10}}


def oracle(c,mean_zero):
    values=[(2*l+1,F(l*(l+1))-F(c)) for l in range(1 if mean_zero else 0,20)]
    return sum(d for d,x in values if x<0),sum((d*x for d,x in values if x<0),F(0))


class MultispectrumTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.t=task();cls.result=v.solve(cls.t);cls.cert=cls.result['certificate']

    def test_constant_thresholds_via_generic_kernel(self):
        for c in ('599/100','6','601/100','7'):
            for mean in (False,True):
                with self.subTest(c=c,mean_zero=mean):
                    t=task(c,mean);r=v.solve(t);cert=r['certificate'];expected,s=oracle(c,mean)
                    self.assertTrue(r['certificate_valid']);self.assertTrue(r['target_met'])
                    self.assertEqual(r['negative_count'],expected)
                    self.assertLessEqual(F(r['sum_lower']),s);self.assertGreaterEqual(F(r['sum_upper']),s)
                    self.assertEqual(F(r['abs_negative_sum_lower']),-F(r['sum_upper']))
                    self.assertFalse(cert['constant_shortcut_used'])
                    self.assertTrue(all(row['kernel_evidence'] for row in cert['sectors']))

    def test_zero_modes_and_equal_angular_tail(self):
        self.assertEqual(self.cert['negative_count'],4)
        self.assertEqual(self.cert['angular_tail_lower'],'0')
        self.assertEqual(self.cert['angular_tail_m_start'],2)
        self.assertTrue(any(row['ritz_at_zero_inertia'][1]>0 for row in self.cert['sectors']))

    def test_nonconstant_step_independent_count_oracle(self):
        # -5 <= q <= -4: precisely l=0,1 can be negative, with multiplicity.
        for mean,count in ((False,4),(True,3)):
            t=task('4',mean,'-1');r=v.solve(t)
            self.assertTrue(r['certificate_valid']);self.assertTrue(r['count_certified'])
            self.assertEqual(r['negative_count'],count)
            self.assertTrue(v.verify(r['certificate'],t)['certificate_valid'])

    def test_nonconstant_negative_power(self):
        for mean,count in ((False,9),(True,8)):
            t=task('6',mean);t['function'].update(profile='abs_power',exponent='-1/4',amplitude='-1/5')
            r=v.solve(t)
            self.assertTrue(r['certificate_valid']);self.assertTrue(r['count_certified'])
            self.assertEqual(r['negative_count'],count)
            self.assertTrue(v.verify(r['certificate'],t)['certificate_valid'])
            self.assertEqual(r['certificate']['function'],t['function'])

    def test_extreme_precision_is_not_falsely_met(self):
        t=task('6',False,tolerance='1/100000000000000000000')
        r=v.negative_spectrum(t,modes=(8,),bits=8)
        self.assertTrue(r['certificate_valid']);self.assertTrue(r['count_certified'])
        self.assertTrue(r['sum_certified']);self.assertFalse(r['sum_target_met']);self.assertFalse(r['target_met'])

    def test_underresolved_radial_tail_is_failure(self):
        r=v.negative_spectrum(task('601/100'),modes=(1,),near_tail=0,bits=8)
        self.assertFalse(r['target_met']);self.assertFalse(r['certificate_valid'])
        self.assertEqual(r['attempts'][0]['status'],'stage_failed')
        self.assertIn('Strictly positive radial tail',r['attempts'][0]['error'])

    def test_zero_ritz_value_retains_unresolved_count_interval(self):
        # The one-mode Ritz value is exactly zero; a coupling correction can
        # allow a negative eigenvalue. Neither zero nor three is certified.
        t=task('3/2',True,'-1')
        r=v.negative_spectrum(t,modes=(1,),near_tail=0,bits=8)
        self.assertTrue(r['certificate_valid']);self.assertFalse(r['count_certified'])
        self.assertEqual((r['count_lower'],r['count_upper']),(0,3))
        self.assertIsNone(r['negative_count']);self.assertIsNone(r['sum_lower'])
        self.assertFalse(r['sum_target_met']);self.assertFalse(r['target_met'])

    def test_replay_without_proposals_or_search(self):
        with patch.object(v,'solve',side_effect=AssertionError('search forbidden')),\
             patch.object(v,'eigenvalue_interval',side_effect=AssertionError('search forbidden')),\
             patch.object(v,'form_for',side_effect=AssertionError('proposal forbidden')):
            self.assertTrue(v.verify(self.cert,self.t)['target_met'])

    def test_certificate_mutations_rejected(self):
        mutations=[
            lambda c:c.update(quantity='absolute_negative_sum'),
            lambda c:c.update(mean_zero=True),
            lambda c:c.update(negative_threshold='1'),
            lambda c:c.update(negative_count=9),
            lambda c:c.update(sum_lower='-100'),
            lambda c:c.update(abs_negative_sum_lower='100'),
            lambda c:c.update(angular_tail_lower='100'),
            lambda c:c['function'].update(axis=0),
            lambda c:c['form_bound'].update(energy_factor='2'),
            lambda c:c['sectors'][1].update(multiplicity=1),
            lambda c:c['sectors'][0].update(projection='remove_l0'),
            lambda c:c['sectors'][0].update(ritz_at_zero_inertia=[0,0,8]),
            lambda c:c['sectors'][0]['eigenvalue_intervals'][0].update(index=2),
            lambda c:c['sectors'][0]['eigenvalue_intervals'][0].update(lower='-5'),
            lambda c:c['sectors'][0]['kernel_evidence']['A'][0].__setitem__(0,'0'),
            lambda c:c['sectors'].pop(0),
            lambda c:c['sectors'].append(deepcopy(c['sectors'][-1])),
        ]
        for i,edit in enumerate(mutations):
            with self.subTest(mutation=i):
                bad=deepcopy(self.cert);edit(bad)
                self.assertFalse(v.verify(bad,self.t)['certificate_valid'])

    def test_task_binding_and_input_rejections(self):
        for key,value in (('mean_zero',True),('tolerance','1/100'),('geometry','torus'),('domain','interval')):
            t=deepcopy(self.t);t[key]=value
            self.assertFalse(v.verify(self.cert,t)['certificate_valid'])
        for value in (0,1,None,'False'):
            t=deepcopy(self.t);t['mean_zero']=value
            self.assertFalse(v.solve(t)['certificate_valid'])
        t=deepcopy(self.t);t['function']['offset']=-6.0
        self.assertFalse(v.solve(t)['certificate_valid'])
        t=deepcopy(self.t);t['function']['axis']=True
        self.assertFalse(v.solve(t)['certificate_valid'])

    def test_count_only_certificate_does_not_claim_sum(self):
        t=v.normalize_task(self.t);form,source=v.form_for(t['function'],8)
        rows=[v.radial_evidence(v.kernel(t['function'],m,False,8,16,form,source))
              for m in range(v.angular_cutoff(form))]
        cert=v.assemble(t,form,source,rows);r=v.verify(cert,t)
        self.assertTrue(r['certificate_valid']);self.assertTrue(r['count_certified'])
        self.assertFalse(r['sum_certified']);self.assertFalse(r['target_met'])

    def test_empty_negative_spectrum_sum_is_zero(self):
        r=v.solve(task('0'))
        self.assertTrue(r['target_met']);self.assertEqual(r['negative_count'],0)
        self.assertEqual(r['sum_lower'],'0');self.assertEqual(r['sum_upper'],'0')

    def test_zero_budget_has_no_false_certificate(self):
        t=task();t['budget']['wall_seconds']=0;r=v.solve(t)
        self.assertFalse(r['certificate_valid']);self.assertFalse(r['target_met'])

    def test_baseline_delegation_unchanged(self):
        t={'kind':'spectrum','function':task()['function'],'mean_zero':False,'tolerance':'1/1000'}
        with patch.object(v.shared,'solve',return_value={'certificate':{'old':True},'target_met':False}) as mock:
            r=v.solve(t);mock.assert_called_once_with(t)
            self.assertTrue(r['fallback']);self.assertEqual(r['certificate'],{'old':True})


if __name__=='__main__':unittest.main()
