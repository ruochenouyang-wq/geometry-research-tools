"""Bounded exact tests; every fixture uses at most eight retained modes."""
import copy
from fractions import Fraction as F
from math import comb
import unittest
from unittest.mock import patch
import wide_potential as w
from common import projected, sectors

Q12 = [0]*12+[1]
ASYM = [F(comb(12,k),2**12) for k in range(13)]

class WideTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.range = w.sturm_range(Q12)
        cls.sector = w.certify_sector(Q12,modes=8,bits=20)
        cls.full = w.full_ground(Q12,modes=8,max_modes=8,bits=20,max_m=2)
        cls.exp = w.exponential_ground(1,degree=12,modes=8,max_modes=8,bits=20)

    def mutate(self, cert, verifier, path, value):
        cert=copy.deepcopy(cert); cursor=cert
        for key in path[:-1]: cursor=cursor[key]
        cursor[path[-1]]=value
        self.assertFalse(verifier(cert))

    def test_cheap_degree24_actual_sphere_bound(self):
        cert=w.cheap_ground_enclosure([0]*24+[1])
        self.assertTrue(w.verify_cheap_ground(cert,expected_q=[0]*24+[1],expected_mean_zero=True))
        self.assertEqual((F(cert['lower']),F(cert['upper'])),(F(2),F(451,225)))
        self.assertEqual([x['sphere_function'] for x in cert['trials']],['z','x'])
    def test_cheap_unconstrained_uses_constant_trial(self):
        cert=w.cheap_ground_enclosure([0]*24+[1],mean_zero=False)
        self.assertTrue(w.verify_cheap_ground(cert,expected_mean_zero=False))
        self.assertEqual((F(cert['lower']),F(cert['upper'])),(F(0),F(1,25)))
        self.assertEqual(cert['trials'][0]['sphere_function'],'1')
    def test_cheap_constant_potential_is_exact(self):
        for mean_zero,answer in [(False,-3),(True,-1)]:
            cert=w.cheap_ground_enclosure([-3],mean_zero)
            self.assertEqual((F(cert['lower']),F(cert['upper'])),(answer,answer))
            self.assertEqual(cert['status'],'certified_exact')
    def test_cheap_bounds_enclose_refined_full_bound(self):
        cert=w.cheap_ground_enclosure(Q12)
        self.assertLessEqual(F(cert['lower']),F(self.full['lower']))
        self.assertGreaterEqual(F(cert['upper']),F(self.full['upper']))
    def test_cheap_wrong_original_potential_and_scope_rejected(self):
        cert=w.cheap_ground_enclosure(Q12)
        self.assertFalse(w.verify_cheap_ground(cert,expected_q=[0,0,1]))
        self.assertFalse(w.verify_cheap_ground(cert,expected_mean_zero=False))
        self.mutate(cert,w.verify_cheap_ground,['scope'],'single_azimuth_sector_only')
        self.mutate(cert,w.verify_cheap_ground,['mean_zero'],False)
    def test_cheap_trial_moment_and_mass_tamper_rejected(self):
        cert=w.cheap_ground_enclosure(Q12)
        self.mutate(cert,w.verify_cheap_ground,['trials',1,'potential_integral'],'0')
        self.mutate(cert,w.verify_cheap_ground,['trials',1,'radial_mass'],'1')
        self.mutate(cert,w.verify_cheap_ground,['trials',1,'sphere_function'],'1')
    def test_cheap_range_and_endpoint_tamper_rejected(self):
        cert=w.cheap_ground_enclosure(Q12)
        self.mutate(cert,w.verify_cheap_ground,['range_proof','q_coefficients'],['0'])
        self.mutate(cert,w.verify_cheap_ground,['lower'],'3')
        self.mutate(cert,w.verify_cheap_ground,['upper'],'2')
    def test_cheap_replay_without_generator_or_spectral_search(self):
        cert=w.cheap_ground_enclosure(Q12,range_proof=self.range)
        with patch.object(w,'cheap_ground_enclosure',side_effect=AssertionError('producer')), \
             patch.object(w,'certify_sector',side_effect=AssertionError('spectrum')), \
             patch.object(w.exact_extrema,'maximize',side_effect=AssertionError('range search')):
            self.assertTrue(w.verify_cheap_ground(cert))
    def test_cheap_projection_requires_boolean(self):
        with self.assertRaises(ValueError): w.cheap_ground_enclosure(Q12,1)
    def test_degree24_accepted(self):
        self.assertEqual(len(w.potential([0]*24+[1])),25)
    def test_old_degree12_rejected(self):
        with self.assertRaises(ValueError): projected.certify_sector(Q12,modes=2,bits=8)
    def test_degree25_rejected(self):
        with self.assertRaises(ValueError): w.potential([0]*25+[1])
    def test_inexact_and_boolean_coefficients_rejected(self):
        for coefficient in [0.5,True,None,'1.5','__import__("os").system("false")']:
            with self.subTest(coefficient=coefficient), self.assertRaises(ValueError):
                w.potential([coefficient])
    def test_coefficient_budgets(self):
        self.assertEqual(w.potential([F(1,2**100)]),[F(1,2**100)])
        for q in [[10**6+1],[F(1,2**2049)]]:
            with self.assertRaises(ValueError): w.potential(q)
    def test_multiplication_not_projected_between_powers(self):
        # t^2 P_1 = 3/5 P_1 + 2/5 P_3, including the intermediate l=0 path.
        self.assertEqual(w.q_times_basis([0,0,1],0,1),{1:F(3,5),3:F(2,5)})
    def test_degree24_entire_multiplication(self):
        expanded=w.q_times_basis([0]*24+[1],0,1)
        self.assertEqual(max(expanded),25)
        self.assertEqual(sum(expanded.values()),1) # P_l(1)=1
    def test_assembly_independent_moment(self):
        # P_1=t, so <P_1,t^12 P_1>= integral t^14=2/15.
        a=w.assemble(Q12,modes=1)
        self.assertEqual(a['M'],[F(2,3)])
        self.assertEqual(a['A'][0][0],F(4,3)+F(2,15))
    def test_associated_mass_and_legacy_low_degree_matrix(self):
        self.assertEqual(w.basis_mass(2,2),F(48,5))
        self.assertEqual(w.assemble([1,-2,3],2,2,4),sectors.assemble([1,-2,3],2,2,4))
    def test_high_degree_associated_analytic_moment(self):
        a=w.assemble([0]*24+[1],m=1,start=1,modes=1)
        # (P_1^1)^2=1-t^2: q moment is 2/25-2/27.
        self.assertEqual(a['A'][0][0],2*a['M'][0]+F(2,25)-F(2,27))
    def test_full_searches_range_only_once(self):
        with patch.object(w,'potential_range',wraps=w.potential_range) as search:
            cert=w.full_ground(Q12,modes=3,max_modes=3,bits=8,max_m=2)
            self.assertTrue(w.verify_full(cert))
            self.assertEqual(search.call_count,1)
    def test_sturm_missing_partition_rejected(self):
        self.mutate(self.range,w.verify_range,['proof','positive','partition'],[])
    def test_asymmetric_wide_matrix_is_symmetric(self):
        a=w.assemble(ASYM,m=1,start=1,modes=4)['A']
        self.assertEqual(a,[list(row) for row in zip(*a)])
        self.assertNotEqual(a[0][1],0)
    def test_sturm_whole_domain_exact_bounds(self):
        self.assertTrue(w.verify_range(self.range,Q12))
        self.assertEqual((self.range['lower'],self.range['upper']),('0','1'))
    def test_bernstein_asymmetric_bounds(self):
        cert=w.bernstein_range(ASYM)
        self.assertTrue(w.verify_range(cert,ASYM))
        self.assertEqual((cert['lower'],cert['upper']),('0','1'))
    def test_auto_cost_dispatch(self):
        self.assertEqual(w.potential_range(Q12)['method'],'bernstein')
        self.assertEqual(w.potential_range([0,0,1])['method'],'sturm')
        with self.assertRaises(ValueError): w.potential_range(Q12,'untrusted_backend')
    def test_range_scope_tamper(self):
        self.mutate(self.range,w.verify_range,['interval'],['0','1'])
    def test_sturm_evidence_tamper(self):
        self.mutate(self.range,w.verify_range,['proof','positive','maximum_upper'],'0')
    def test_bernstein_partition_tamper(self):
        self.mutate(w.bernstein_range(Q12),w.verify_range,['proof','cells'],[])
    def test_range_reuse_does_no_extrema_search(self):
        with patch.object(w,'potential_range',side_effect=AssertionError('range search')), \
             patch.object(w.exact_extrema,'maximize',side_effect=AssertionError('root search')):
            cert=w.certify_sector(Q12,modes=4,bits=8,range_proof=self.range)
            self.assertTrue(w.verify_sector(cert))
    def test_reused_range_wrong_q_rejected(self):
        with self.assertRaises(ValueError): w.reuse_range([0,0,1],self.range)
    def test_reused_range_is_deep_copy(self):
        replay=w.reuse_range(Q12,self.range); replay['proof']['positive']['maximum_upper']='100'
        self.assertTrue(w.verify_range(self.range))
    def test_kernel_all_degree12_tail_columns(self):
        cert=w.kernel_certificate(Q12,modes=4)
        self.assertEqual([x['degree'] for x in cert['tail_couplings']],list(range(5,17)))
        self.assertTrue(w.verify_kernel(cert))
    def test_kernel_tail_column_omission_rejected(self):
        cert=w.kernel_certificate(Q12,modes=4); cert['tail_couplings'].pop()
        self.assertFalse(w.verify_kernel(cert))
    def test_kernel_tail_mass_tamper(self):
        self.mutate(w.kernel_certificate(Q12,modes=4),w.verify_kernel,['tail_couplings',0,'mass'],'1')
    def test_sector_verifies_in_original_high_degree_problem(self):
        self.assertTrue(w.verify_sector(self.sector,expected_q=Q12,expected_m=0,expected_mean_zero=True,expected_k=1))
        self.assertFalse(w.verify_sector(self.sector,expected_q=[0,0,1]))
        self.assertFalse(w.verify_sector(self.sector,expected_m=1))
        self.assertFalse(w.verify_sector(self.sector,expected_mean_zero=False))
        self.assertFalse(w.verify_sector(self.sector,expected_k=2))
    def test_sector_target_bound_tamper(self):
        self.mutate(self.sector,w.verify_sector,['lower'],'100')
    def test_sector_full_scope_confusion(self):
        self.mutate(self.sector,w.verify_sector,['scope'],w.scope(True))
        self.assertFalse(w.verify_full(self.sector))
    def test_sector_mass_tamper(self):
        self.mutate(self.sector,w.verify_sector,['kernel_evidence','mass',0],'1')
    def test_full_reconstructs_all_sphere_angular_tail(self):
        self.assertTrue(w.verify_full(self.full,expected_q=Q12,expected_mean_zero=True))
        self.assertEqual(self.full['scope'],w.scope(True))
        self.assertGreaterEqual(F(self.full['angular_tail_lower']),F(self.full['upper']))
    def test_full_missing_sector_rejected(self):
        self.mutate(self.full,w.verify_full,['sectors'],self.full['sectors'][1:])
    def test_full_angular_tail_tamper(self):
        self.mutate(self.full,w.verify_full,['angular_tail_lower'],'1000000')
    def test_finite_budget_returns_valid_open(self):
        cert=w.full_ground(Q12,modes=2,max_modes=2,bits=8,max_m=0,tolerance=F(1,10**18))
        self.assertTrue(w.verify_full(cert))
        self.assertEqual(cert['status'],'certified_bound_open_gap')
        self.assertFalse(cert['angular_tail_cannot_improve_best_upper'])
    def test_low_degree_cross_check(self):
        q=[0,0,2]
        old=projected.full_ground(q,modes=6,max_modes=6,bits=20,max_m=2)
        new=w.full_ground(q,modes=6,max_modes=6,bits=20,max_m=2)
        self.assertTrue(w.verify_full(new))
        self.assertLessEqual(max(F(old['lower']),F(new['lower'])),min(F(old['upper']),F(new['upper'])))
    def test_constant_spectrum_exact(self):
        cert=w.full_ground([3],modes=2,max_modes=2,bits=8)
        self.assertEqual((cert['lower'],cert['upper']),('5','5'))
    def test_taylor_exact_remainder(self):
        cert=w.exponential_approximation(1,12)
        self.assertEqual(F(cert['absolute_error_bound']),F(3,6227020800))
        self.assertTrue(w.verify_approximation(cert,1))
    def test_taylor_invalid_budget(self):
        for a,n in [(4,12),(True,12),(1,25),(1,True),(1,-1)]:
            with self.assertRaises(ValueError): w.exponential_approximation(a,n)
    def test_taylor_rational_parameter_at_max_degree(self):
        cert=w.exponential_approximation(F(1,10),24)
        self.assertTrue(w.verify_approximation(cert,F(1,10)))
        self.assertEqual(len(cert['polynomial_coefficients']),25)
    def test_taylor_zero_exact(self):
        cert=w.exponential_approximation(0,24)
        self.assertEqual(cert['absolute_error_bound'],'0')
        self.assertTrue(w.verify_approximation(cert))
    def test_taylor_delta_tamper(self):
        self.mutate(w.exponential_approximation(1),w.verify_approximation,['absolute_error_bound'],'0')
    def test_taylor_polynomial_tamper(self):
        self.mutate(w.exponential_approximation(1),w.verify_approximation,['polynomial_coefficients',2],'0')
    def test_exponential_original_binding_and_enlargement(self):
        cert=self.exp; delta=F(cert['perturbation_bound']); inner=cert['polynomial_spectrum']
        self.assertTrue(w.verify_exponential(cert,1))
        self.assertFalse(w.verify_exponential(cert,-1))
        self.assertEqual(F(cert['lower']),F(inner['lower'])-delta)
        self.assertEqual(F(cert['upper']),F(inner['upper'])+delta)
    def test_exponential_delta_tamper(self):
        self.mutate(self.exp,w.verify_exponential,['perturbation_bound'],'0')
        self.mutate(self.exp,w.verify_exponential,['approximation','absolute_error_bound'],'0')
    def test_exponential_scope_tamper(self):
        self.mutate(self.exp,w.verify_exponential,['mean_zero'],False)
    def test_exponential_reflection_spectra_agree(self):
        reflected=w.exponential_ground(-1,degree=12,modes=8,max_modes=8,bits=20)
        self.assertTrue(w.verify_exponential(reflected,-1))
        self.assertEqual((self.exp['lower'],self.exp['upper']),(reflected['lower'],reflected['upper']))
    def test_perturbation_rejects_wrong_polynomial(self):
        with self.assertRaises(ValueError):
            w.perturb_certificate(w.exponential_approximation(1),self.full)
    def test_verify_does_no_spectral_or_range_search(self):
        with patch.object(w,'certify_sector',side_effect=AssertionError('spectral search')), \
             patch.object(w,'potential_range',side_effect=AssertionError('range search')), \
             patch.object(w.exact_extrema,'maximize',side_effect=AssertionError('root search')):
            self.assertTrue(w.verify_full(self.full))
            self.assertTrue(w.verify_exponential(self.exp))

if __name__=='__main__': unittest.main()
