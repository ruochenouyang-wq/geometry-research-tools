import copy
import json
import unittest
from unittest.mock import patch
from dependencies import F, wide, constraints
import function_models as f
import spectral_transfer as s


class SpectralTransferTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.step_function = {'kind':'axis_profile','profile':'step','axis':2}
        cls.step = f.l2_model(cls.step_function, degree=0)
        cls.source = s.polynomial_source(cls.step['polynomial'], bits=24, modes=3, max_modes=3)
        cls.bound = s.transfer_l2(cls.step, cls.source)

    def test_exact_constant_linf_has_no_extra_width(self):
        model = f.analytic_model({'kind':'analytic_sum','polynomial':{'0,0,0':'3/2'}})
        for mz, expected in [(True,F(7,2)),(False,F(3,2))]:
            source = s.polynomial_source(model['polynomial'],mean_zero=mz,modes=2,max_modes=2,bits=12)
            c = s.transfer_linf(model,source)
            self.assertEqual((F(c['lower']),F(c['upper'])),(expected,expected))
            self.assertTrue(s.verify(c,expected_mean_zero=mz))

    def test_exact_exp_zero_is_one(self):
        model = f.analytic_model({'kind':'analytic_sum','terms':[{'function':'exp','argument':{}}]},order=4)
        c = s.transfer_linf(model,s.polynomial_source(model['polynomial'],modes=2,max_modes=2,bits=12))
        self.assertEqual((c['lower'],c['upper']),('3','3'))

    def test_real_exponential_full_original_function(self):
        function = {'kind':'analytic_sum','terms':[{'function':'exp','argument':{'0,0,1':'1'}}]}
        c = s.solve(function,order=12,source_options={'modes':8,'max_modes':12,'bits':40})
        self.assertTrue(s.verify(c,expected_function=function))
        self.assertEqual(c['status'],'target_met')
        self.assertLessEqual(F(c['exact_width']),F(1,10**8))
        self.assertEqual(F(c['error_budget']['representation_widening']),2*F(c['model']['error_upper']))

    def test_multiaxis_source_is_complete_cartesian_not_single_sector(self):
        p = {'1,0,0':'1/10','0,1,0':'1/20','0,0,1':'1/30'}
        source = s.polynomial_source(p,L=2,bits=16)
        self.assertEqual(source['backend'],'cartesian_constraints')
        self.assertTrue(s.verify_source(source,expected_polynomial=p))
        self.assertTrue(source['full_infinite_space_covered'])

    def test_axis_isometry_preserves_complete_spaces(self):
        a = s.polynomial_source({'1,0,0':'1/10'},modes=3,max_modes=3,bits=16)
        b = s.polynomial_source({'0,0,1':'1/10'},modes=3,max_modes=3,bits=16)
        self.assertEqual((a['lower'],a['upper']),(b['lower'],b['upper']))
        self.assertFalse(s.verify_source(a,expected_polynomial=b['polynomial']))

    def test_step_probability_l2_exact_calibration(self):
        self.assertEqual(F(self.step['error_squared']),F(1,4))
        self.assertEqual(F(self.step['error_upper']),F(1,2))
        self.assertEqual((self.bound['lower'],self.bound['upper']),('1','4'))
        self.assertEqual(self.bound['status'],'certified_open')
        self.assertTrue(s.verify(self.bound))

    def test_zero_l2_error_recovers_source(self):
        model = f.l2_model({'kind':'axis_profile','profile':'abs_power','exponent':'2'},degree=2)
        self.assertEqual(F(model['error_upper']),0)
        source = s.polynomial_source(model['polynomial'],modes=4,max_modes=4,bits=20)
        result = s.transfer_l2(model,source)
        self.assertEqual((result['lower'],result['upper']),(source['lower'],source['upper']))

    def test_unbounded_negative_profile_fixed_reference(self):
        function = {'kind':'axis_profile','profile':'abs_power','exponent':'-1/4','amplitude':'-1'}
        c = s.solve(function,method='fixed',degree=4,reference={'0,0,0':'-4/3'},
                    source_options={'modes':4,'max_modes':4,'bits':16,'tolerance':'1/1000'})
        self.assertTrue(s.verify(c,expected_function=function))
        self.assertEqual(F(c['coercivity']['shift']),F(7,3))
        self.assertEqual(F(c['coercivity']['distance_squared'])+F(c['model']['error_squared']),F(2,9))
        self.assertLess(F(c['error_evidence']['relative_form_error']),1)

    def test_same_reference_reduces_to_direct_l2_formula(self):
        c = s.transfer_fixed(self.step,self.source,self.source['polynomial'])
        self.assertEqual(c['coercivity']['distance_squared'],'0')
        self.assertEqual((c['lower'],c['upper']),(self.bound['lower'],self.bound['upper']))

    def test_reference_distance_one_rejected(self):
        with self.assertRaises(ValueError):
            s.fixed_reference({'0,0,0':'1'},{})

    def test_relative_error_one_rejected(self):
        with self.assertRaises(ValueError):
            s.transfer_fixed(self.step,self.source,{'0,0,0':'0'})

    def test_sqrt_upper_is_outward_and_exact_at_squares(self):
        for q in [F(0),F(1,4),F(2,9),F(10**20,3),F(1,10**30)]:
            root=s.sqrt_upper(q,bits=40)
            self.assertGreaterEqual(root*root,q)
            if root: self.assertLess((root-F(1,2**40))**2,q)
        self.assertEqual(s.sqrt_upper(F(1,4)),F(1,2))

    def test_source_input_scope_and_axis_tampering_rejected(self):
        for key, value in [('polynomial',{}),('mean_zero',False),('eigenvalue_index',2),
                           ('scope','single_axisymmetric_sector'),('lower','99')]:
            c=copy.deepcopy(self.source);c[key]=value
            self.assertFalse(s.verify_source(c),key)
        c=copy.deepcopy(self.source);c['conversion']['axis']=0
        self.assertFalse(s.verify_source(c))

    def test_single_sector_is_not_full_source(self):
        sector=wide.certify_sector(['1/2'],m=0,modes=3,bits=12)
        with self.assertRaises(ValueError):
            s._source_record({'0,0,0':'1/2'},True,'axis_isometry_wide',sector)

    def test_model_mismatch_and_norm_substitution_rejected(self):
        other=s.polynomial_source({'0,0,0':'2'},modes=2,max_modes=2,bits=12)
        with self.assertRaises(ValueError):s.transfer_l2(self.step,other)
        with self.assertRaises(ValueError):s.transfer_linf(self.step,self.source)
        model=f.analytic_model({'kind':'analytic_sum'})
        with self.assertRaises(ValueError):s.transfer_l2(model,self.source)

    def test_error_original_function_and_target_binding(self):
        self.assertFalse(s.verify(self.bound,expected_mean_zero=False))
        self.assertFalse(s.verify(self.bound,expected_width='1/100'))
        self.assertFalse(s.verify(self.bound,expected_function={'kind':'axis_profile','profile':'step','axis':0}))
        for mutate in [lambda c:c['model'].__setitem__('error_upper','0'),
                       lambda c:c.__setitem__('upper',c['lower']),
                       lambda c:c.__setitem__('measure','d_sigma'),
                       lambda c:c.__setitem__('status','target_met'),
                       lambda c:c['error_budget'].__setitem__('representation_widening','0')]:
            c=copy.deepcopy(self.bound);mutate(c);self.assertFalse(s.verify(c))

    def test_embedding_payload_has_no_shared_mutable_state(self):
        c=s.transfer_l2(self.step,self.source)
        c['error_evidence']['embedding']['energy_coefficient']='0'
        self.assertFalse(s.verify(c))
        self.assertTrue(s.verify(self.bound))
        self.assertEqual(s.EMBEDDING['energy_coefficient'],'1')

    def test_fixed_reference_tampering_rejected(self):
        c=s.transfer_fixed(self.step,self.source,self.source['polynomial'])
        for k,v in [('distance_squared','1/9'),('shift','0'),('coercivity_factor','2')]:
            bad=copy.deepcopy(c);bad['coercivity'][k]=v;self.assertFalse(s.verify(bad))

    def test_verification_replays_without_optimization_search(self):
        c=s.transfer_fixed(self.step,self.source,self.source['polynomial'])
        with patch.object(wide,'full_ground',side_effect=AssertionError('search forbidden')), \
             patch.object(constraints,'certify',side_effect=AssertionError('search forbidden')), \
             patch.object(wide,'potential_range',side_effect=AssertionError('range search forbidden')):
            self.assertTrue(s.verify(json.loads(json.dumps(c))))

    def test_exact_inputs_only_and_unknown_options_rejected(self):
        for bad in [True,0.1]:
            with self.assertRaises((ValueError,TypeError)):s.poly({'0,0,0':bad})
        with self.assertRaises(ValueError):s.solve(self.step_function,method='l2',source_options={'ignored':1})
        with self.assertRaises(ValueError):s.solve(self.step_function,method='unknown')


if __name__=='__main__':unittest.main()
