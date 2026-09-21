"""Regression and adversarial checks for radius mathematics, not schema echoes."""
import unittest
from copy import deepcopy
from common import F, projected
import sphere_scale as s
import ordered_spectrum as ordered

FAST=dict(modes=2,max_modes=4,max_m=3,bits=20)


class RadiusTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.zero=s.backend_ground([0],2,**FAST)
        cls.quadratic=s.backend_ground([0,0,1],2,**FAST)
        cls.search=s.adaptive_minimum([0,0,1],[1,2],max_splits=1,tolerance='1/100',**FAST)

    def test_laplace_R2_is_half(self):
        self.assertEqual(self.zero['lower'],'1/2');self.assertEqual(self.zero['upper'],'1/2')
        self.assertTrue(s.verify(self.zero,expected_radius=2))

    def test_R1_reproduces_unit(self):
        c=s.backend_ground([0,1,1],1,**FAST)
        self.assertEqual(c['lower'],c['unit_certificate']['lower'])
        self.assertEqual(c['upper'],c['unit_certificate']['upper'])

    def test_negative_constant_shift(self):
        c=s.backend_ground([-6],2,**FAST)
        self.assertEqual(c['lower'],'-11/2');self.assertEqual(c['upper'],'-11/2')

    def test_physical_power_pullback(self):
        p=s.physical_pullback([1,-2,3],2)
        self.assertEqual(p['pulled_back_potential'],['1','-4','12'])
        self.assertEqual(p['unit_operator_potential'],['4','-16','48'])
        self.assertTrue(s.verify(p))

    def test_coordinate_energies_differ(self):
        norm=s.scale_trial([0,0,1],2,[1])
        physical=s.scale_trial([0,0,1],2,[1],coordinates='ambient_x3')
        self.assertEqual(norm['rayleigh_quotient'],'11/10')
        self.assertEqual(physical['rayleigh_quotient'],'29/10')
        self.assertTrue(s.verify(norm));self.assertTrue(s.verify(physical))

    def test_area_and_dirichlet_scaling(self):
        a=s.scale_trial([0],1,[1]);b=s.scale_trial([0],2,[1])
        self.assertEqual(F(b['mass_divided_by_azimuth_factor']),4*F(a['mass_divided_by_azimuth_factor']))
        self.assertEqual(b['dirichlet_energy_divided_by_azimuth_factor'],a['dirichlet_energy_divided_by_azimuth_factor'])

    def test_wide_trial_exact_integral(self):
        t=s.scale_trial([0]*12+[1],2,[1])
        # u=t: potential/mass = E[t^14]/E[t^2] = 1/5.
        self.assertEqual(t['rayleigh_quotient'],'7/10');self.assertTrue(s.verify(t))

    def test_wide_trial_m1(self):
        t=s.scale_trial([0,0,2],1,[1,'-1/36'],m=1,degrees=[1,3])
        self.assertEqual(t['rayleigh_quotient'],'722/303')

    def test_ordered_negative_trace_scaling(self):
        u=ordered.adaptive_spectrum([-6],k=8,modes=2,max_modes=3,max_radial=3,bits=20)['certificate']
        c=s.scale_ordered(['-3/2'],2,u)
        self.assertEqual(c['negative_trace_lower'],'3');self.assertEqual(c['negative_trace_upper'],'3')
        self.assertEqual(c['negative_count_lower'],3);self.assertEqual(c['negative_count_upper'],3)
        self.assertEqual(c['ordered_intervals'][3]['lower'],'0');self.assertTrue(s.verify(c))

    def test_wrong_unit_potential_rejected(self):
        with self.assertRaises(ValueError):s.scale_ground([1],2,self.zero['unit_certificate'])

    def test_non_mean_zero_scope_rejected(self):
        unit=projected.full_ground([0],mean_zero=False,modes=2,max_modes=2,bits=20)
        with self.assertRaises(ValueError):s.scale_ground([0],2,unit)

    def test_radius_validation(self):
        for r in [0,-1,True,0.5,'0']:
            with self.subTest(r=r),self.assertRaises((ValueError,TypeError)):
                s.backend_ground([0],r,**FAST)

    def test_invalid_coordinate_validation(self):
        with self.assertRaises(ValueError):s.physical_pullback([0],0)
        with self.assertRaises(ValueError):s.scale_trial([0],1,[1],coordinates='physical?')

    def test_mean_zero_trial_rejected(self):
        with self.assertRaises(ValueError):s.scale_trial([0],2,[1],degrees=[0])

    def test_fixed_radius_expected_binding(self):
        self.assertFalse(s.verify(self.zero,expected_radius=1))
        self.assertFalse(s.verify(self.zero,expected_coordinates='ambient_x3'))
        self.assertFalse(s.verify(self.zero,expected_q=[1]))

    def test_scope_tamper(self):
        c=deepcopy(self.zero);c['problem']['mean_zero']=False
        self.assertFalse(s.verify(c))

    def test_physical_coordinate_tamper(self):
        c=deepcopy(self.quadratic);c['problem']['coordinates']='ambient_x3'
        self.assertFalse(s.verify(c))

    def test_radius_interval_laplace(self):
        c=s.radius_cell([0],[1,2],s.backend_ground([0],'3/2',**FAST))
        self.assertEqual(c['lower'],'1/2');self.assertEqual(c['upper'],'2')
        self.assertTrue(s.verify(c))

    def test_sign_aware_interval(self):
        c=s.radius_cell([-6],[1,2],s.backend_ground([-6],'3/2',**FAST))
        self.assertEqual(c['lower'],'-11/2');self.assertEqual(c['upper'],'-4')
        products=list(map(F,c['four_signed_products']))
        unitlo,unithi=map(F,c['unit_eigenvalue_interval'])
        # Negative lower uses largest inverse-square, not the smallest.
        self.assertEqual(min(products),unitlo)
        self.assertNotEqual(min(products),unitlo/4)
        self.assertTrue(s.verify(c))

    def test_invalid_interval_and_anchor(self):
        for interval in [[2,1],[1,1],[0,2]]:
            with self.subTest(interval=interval),self.assertRaises(ValueError):
                s.radius_cell([0],interval,self.zero)
        with self.assertRaises(ValueError):s.radius_cell([0],[3,4],self.zero)

    def test_exact_global_laplace_minmax(self):
        g=s.adaptive_minimum([0],[1,2],max_splits=0,**FAST)['certificate']
        self.assertEqual((g['lower'],g['upper']),('1/2','1/2'))
        self.assertEqual((g['maximum_lower'],g['maximum_upper']),('2','2'))
        self.assertEqual(g['best_feasible_radius'],'2')

    def test_budget_open_not_fake_precision(self):
        g=s.adaptive_minimum([0,0,1],[1,2],max_splits=0,tolerance='1/1000000',**FAST)['certificate']
        self.assertEqual(g['status'],'certified_global_bound_open_gap');self.assertTrue(s.verify(g))

    def test_nonconstant_split_improves(self):
        h=self.search['search_history']
        self.assertLess(F(h[-1]['width']),F(h[0]['width']))
        self.assertTrue(s.verify(self.search['certificate']))

    def test_partition_gap_rejected(self):
        c=deepcopy(self.search['certificate']);c['cells'].pop()
        self.assertFalse(s.verify(c))

    def test_partition_overlap_rejected(self):
        c=deepcopy(self.search['certificate']);c['cells'].append(deepcopy(c['cells'][0]))
        self.assertFalse(s.verify(c))

    def test_feasible_upper_binding(self):
        c=deepcopy(self.search['certificate'])
        best=next(p for p in c['feasible_points'] if p['radius']==c['best_feasible_radius'])
        self.assertEqual(c['upper'],best['upper'])
        c['upper']=c['maximum_upper'];self.assertFalse(s.verify(c))

    def test_cached_proof_wrong_q(self):
        with self.assertRaises(ValueError):s.AnchorCache([1],entries=[self.zero])
        with self.assertRaises(ValueError):s.AnchorCache([0],'ambient_x3',entries=[self.zero])

    def test_resume_reuses_and_is_monotone(self):
        old=self.search['certificate'];cp=s.checkpoint(old)
        new=s.resume_minimum(cp,max_splits=1,**FAST)
        self.assertGreaterEqual(F(new['certificate']['lower']),F(old['lower']))
        self.assertLessEqual(F(new['certificate']['upper']),F(old['upper']))
        self.assertEqual(new['work']['new_spectral_anchor_computations'],2)
        self.assertGreaterEqual(new['work']['verified_anchor_cache_hits'],3)
        self.assertTrue(s.verify(new['certificate']))

    def test_checkpoint_domain_and_digest(self):
        cp=s.checkpoint(self.search['certificate'])
        self.assertFalse(s.verify_checkpoint(cp,expected_domain=[1,3]))
        cp['certificate_digest']='0'*64;self.assertFalse(s.verify_checkpoint(cp))

    def test_explicit_radius_violation(self):
        t=s.scale_trial([0,0,'1/2'],2,[1,'-1/36'],m=1,degrees=[1,3])
        c=s.transfer_counterexample(t,'3/5')
        self.assertEqual(c['status'],'explicit_counterexample')
        self.assertEqual(c['rayleigh_quotient'],'361/606')
        self.assertTrue(s.verify(c,expected_threshold='3/5'))
        self.assertFalse(s.verify(c,expected_threshold='1/2'))

    def test_threshold_tamper(self):
        c=s.transfer_counterexample(s.scale_trial([0],2,[1]),1)
        c['threshold']='1/2';self.assertFalse(s.verify(c))

    def test_backends_by_accuracy_and_degree(self):
        p=s.backend_ground([0],2,tolerance='1/1000000000000000000',**FAST)
        self.assertEqual(p['unit_backend'],'precision_v114');self.assertTrue(s.verify(p))
        w=s.backend_ground([0]*12+[1],1,**FAST)
        self.assertEqual(w['unit_backend'],'wide_v124');self.assertTrue(s.verify(w))

    def test_job_all_three_outcomes(self):
        proved=s.research_job([0],'1/2',domain=[1,2],max_splits=0,**FAST)['certificate']
        violated=s.research_job([0],1,radius=2,trial={'coefficients':[1]},**FAST)['certificate']
        undecided=s.research_job([0,0,1],'3/5',domain=[1,2],max_splits=0,**FAST)['certificate']
        self.assertEqual(proved['status'],'proved');self.assertEqual(violated['status'],'refuted_by_explicit_function')
        self.assertEqual(undecided['status'],'undetermined')
        for c in [proved,violated,undecided]:self.assertTrue(s.verify(c))

    def test_job_proof_nodes_and_scope(self):
        c=s.research_job([0],1,radius=2,**FAST)['certificate']
        bad=deepcopy(c);bad['radius_request']['radius']='1';self.assertFalse(s.verify(bad))
        bad=deepcopy(c);bad['threshold']='1/2';self.assertFalse(s.verify(bad))
        bad=deepcopy(c);key=bad['proof_roots']['spectral'];bad['proof_nodes'][key]['radius']='1'
        self.assertFalse(s.verify(bad))

    def test_job_conflicting_trial_radius(self):
        with self.assertRaises(ValueError):
            s.research_job([0],1,radius=2,trial={'coefficients':[1],'radius':1},**FAST)


if __name__=='__main__':unittest.main(verbosity=2)
