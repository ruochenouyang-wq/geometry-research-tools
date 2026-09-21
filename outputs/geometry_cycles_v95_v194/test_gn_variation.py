import copy
from fractions import Fraction as F
import unittest
from unittest.mock import patch
import gn_variation as g


U=[F(0),F(2,5),F(0),F(1)]  # P1+(2/5)P3


class GNVariationTests(unittest.TestCase):
    def test_degree36_moments(self):
        u=[F(0)]*36+[F(1)]
        c=g.moments(u)
        self.assertEqual(F(c['mean']),F(1,37))
        self.assertEqual(F(c['mass']),F(1,73))
        self.assertEqual(F(c['quartic']),F(1,145))
        self.assertEqual(F(c['energy']),36**2*(F(1,71)-F(1,73)))
        self.assertTrue(g.verify(c,u))

    def test_reject_degree37_and_float(self):
        for u in [[0]*37+[1],[0,1.0],[False,1]]:
            with self.assertRaises(ValueError): g.moments(u)

    def test_center_and_signed_scale(self):
        a=g.normalize([7,0,2]);b=g.normalize([-2,0,-6])
        self.assertEqual(a['normalized_u'],['-1/3','0','1'])
        self.assertEqual(a['normalized_u'],b['normalized_u'])
        self.assertTrue(g.verify(a,[7,0,2]))
        self.assertFalse(g.verify(a,[-2,0,-6]))

    def test_zero_and_constant_rejected(self):
        for u in [[0],[3]]:
            with self.assertRaises(ValueError): g.normalize(u)
            with self.assertRaises(ValueError): g.ratio(u)
        with self.assertRaises(ValueError): g.ratio([1,1])

    def test_ratio_scale_invariant(self):
        a=g.ratio(U);b=g.ratio(g.scale(U,F(-7,3)))
        self.assertEqual(a['pi_K_lower'],b['pi_K_lower'])
        self.assertEqual(F(g.ratio([0,1])['pi_K_lower']),F(9,40))
        self.assertFalse(g.verify(a,g.scale(U,2)))

    def test_complete_high_order_residual(self):
        r=g.residual(U)
        self.assertEqual(r['full_degree'],9)
        self.assertNotEqual(F(r['residual'][-1]),0)
        self.assertEqual(F(r['inner_u']),0)
        self.assertEqual(g.mean([F(x) for x in r['residual']]),0)
        self.assertTrue(g.verify(r,U))

    def test_degree36_residual_degree108(self):
        u=g.centered_cap(36)
        r=g.residual(u)
        self.assertEqual(r['full_degree'],108)
        self.assertTrue(g.verify(r,u))

    def test_direction_sign_and_identity(self):
        c=g.direction(U)
        self.assertGreater(F(c['derivative']),0)
        self.assertEqual(c['derivative'],c['identity_derivative'])
        reversed_direction=g.direction(U,g.scale([F(x) for x in c['v']],-1))
        self.assertLess(F(reversed_direction['derivative']),0)
        self.assertTrue(g.verify(c,U))
        c['derivative']=str(-F(c['derivative']))
        self.assertFalse(g.verify(c,U))

    def test_derivative_independent_coefficient_identity(self):
        v=g.centered_cap(4);c=g.direction(U,v)
        m=g.chart_moments(U,v);n,d=m['numerator'],m['denominator']
        expected=(n[1]*d[0]-n[0]*d[1])/d[0]**2
        self.assertEqual(F(c['derivative']),expected)

    def test_projected_zero_is_not_full_stationarity(self):
        c=g.project_residual([0,1],1)
        self.assertTrue(c['projected_stationary'])
        self.assertFalse(c['full_stationary'])
        self.assertEqual(c['outside_support'],[3])
        self.assertTrue(g.verify(c,[0,1]))
        c['full_stationary']=True
        self.assertFalse(g.verify(c,[0,1]))

    def test_projection_pythagoras(self):
        c=g.project_residual(U,5)
        self.assertEqual(F(c['cross_inner']),0)
        self.assertEqual(F(c['full']['norm_squared']),F(c['inside_norm_squared'])+F(c['outside_norm_squared']))
        self.assertTrue(set(c['outside_support'])=={7,9})

    def test_two_chart_global_certificate(self):
        f,gp=[0,1],[0,F(-3,2),0,F(5,2)]
        c=g.plane([f,gp])
        self.assertTrue(g.verify_plane(c,[f,gp]))
        for x,y in [(1,0),(0,1),(1,7),(-3,1),(5,-8),(-2,-9)]:
            val=F(g.ratio(g.add(g.scale(f,x),g.scale(gp,y)))['pi_K_lower'])
            self.assertLessEqual(val,F(c['upper']))
        c['charts'].pop()
        self.assertFalse(g.verify_plane(c))

    def test_plane_high_degree_not_old_adapter(self):
        basis=[g.centered_cap(8),g.centered_cap(12)]
        c=g.plane(basis,F(1,10**6))
        self.assertTrue(g.verify_plane(c,basis))
        self.assertTrue(c['target_met'])

    def test_plane_binding_and_sign_tamper(self):
        c=g.plane([[0,1],U])
        self.assertFalse(g.verify_plane(c,[[0,2],U]))
        c['upper']=c['lower']
        self.assertFalse(g.verify_plane(c))

    def test_dependent_plane_rejected(self):
        with self.assertRaises(ValueError): g.plane([[0,1],[0,2]])

    def test_monotone_update_improves_fixed_seed(self):
        c=g.refine(U,iterations=1,degree_budget=9,tolerance=F(1,10**7))
        self.assertTrue(g.verify(c,U))
        self.assertTrue(c['steps'][0]['accepted'])
        self.assertGreater(F(c['final']['pi_K_lower']),F(g.ratio(U)['pi_K_lower']))
        self.assertFalse(c['full_function_space_global_optimum_claimed'])

    def test_zero_iteration_and_projection_stop(self):
        c=g.refine(U,iterations=0,degree_budget=3)
        self.assertTrue(g.verify(c,U))
        self.assertEqual(F(c['final']['pi_K_lower']),F(g.ratio(U)['pi_K_lower']))
        d=g.refine([0,1],iterations=2,degree_budget=1)
        self.assertEqual(d['stop'],'projected_stationary_with_nonzero_full_residual')
        self.assertTrue(g.verify(d,[0,1]))

    def test_iteration_tamper_and_binding(self):
        c=g.refine(U,iterations=0)
        self.assertFalse(g.verify(c,[0,1]))
        c['final']['pi_K_lower']='100'
        self.assertFalse(g.verify(c,U))

    def test_failed_plane_preserves_old_trial(self):
        with patch.object(g,'plane',side_effect=ArithmeticError('test bounded search exhausted')):
            c=g.refine(U,iterations=2,degree_budget=9)
        self.assertTrue(g.verify(c,U))
        self.assertEqual(c['stop'],'plane_certificate_budget_failed')
        self.assertFalse(c['steps'][0]['accepted'])
        self.assertEqual(c['final']['pi_K_lower'],g.ratio(U)['pi_K_lower'])

    def test_stop_metadata_cannot_claim_stationary(self):
        c=g.refine(U,iterations=0)
        for status in ('full_stationary','no_strict_improvement_found','plane_certificate_budget_failed'):
            altered=copy.deepcopy(c);altered['stop']=status
            self.assertFalse(g.verify(altered,U))

    def test_cap_exact_moments(self):
        for n in (2,4,8,12):
            u=g.centered_cap(n);c=g.moments(u)
            self.assertEqual(F(c['mean']),0)
            self.assertEqual(F(c['mass']),F(1,2*n+1)-F(1,(n+1)**2))
            self.assertEqual(F(c['energy']),F(n,2*(2*n+1)))
            expected=sum(F(comb4,1)*(-F(1,n+1))**(4-j)*F(1,n*j+1) for j,comb4 in enumerate((1,4,6,4,1)))
            self.assertEqual(F(c['quartic']),expected)

    def test_cap_screening_budget_and_target(self):
        c=g.screen_caps(budget=4)
        self.assertEqual([x['degree'] for x in c['trials']],[2,4,8,12])
        self.assertTrue(g.verify(c))
        d=g.screen_caps(budget=2)
        self.assertEqual(d['untried'],[8,12])
        self.assertEqual(d['stop'],'screening_budget_exhausted')
        e=g.screen_caps(target=F(1,100))
        self.assertEqual(len(e['trials']),1)

    def test_degree24_compilation_identity(self):
        u=g.centered_cap(12);c=g.compile_potential(u)
        self.assertEqual(c['q_degree'],24)
        self.assertEqual(c['trial_rayleigh'],c['required_identity'])
        self.assertEqual(F(c['potential_energy']),-2*F(c['ratio']['moments']['energy']))
        self.assertTrue(g.verify(c,u))
        c['q_coefficients'][0]='0'
        self.assertFalse(g.verify(c,u))

    def test_spectral_diagnostic_replay_scope(self):
        c=g.spectral_diagnostic(U,modes=4,max_modes=4,bits=18,max_m=3,tolerance=F(1,1000))
        self.assertIn('spectral',c)
        self.assertTrue(g.verify(c,U))
        self.assertFalse(c['universal_GN_claim'])
        c['spectral']['mean_zero']=False
        self.assertFalse(g.verify(c,U))

    def test_compile_degree_limit_and_diagnostic_rejection(self):
        with self.assertRaises(ValueError): g.compile_potential(g.centered_cap(13))
        c=g.spectral_diagnostic(U,max_m=65)
        self.assertEqual(c['status'],'rejected_or_budget_failed')
        self.assertNotIn('spectral',c)
        self.assertTrue(g.verify(c,U))
        self.assertFalse(g.verify(c,[0,1]))


if __name__=='__main__':
    unittest.main()
