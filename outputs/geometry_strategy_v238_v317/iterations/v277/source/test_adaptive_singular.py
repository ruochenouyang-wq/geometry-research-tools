import json
from copy import deepcopy
import unittest
from unittest.mock import patch

import adaptive_singular as a


class AdaptiveSingularTests(unittest.TestCase):
    def test_generated_powers_recover_structure_not_literal_table(self):
        self.assertEqual(a.generated_powers(a.ORIGINAL, 16), a.enriched.DEFAULT_EXPONENTS)
        with self.assertRaises(ValueError):
            a.generated_powers(a.ORIGINAL, 17)

    def test_original_source_and_global_three_term_proof(self):
        r = a.solve(a.ORIGINAL, max_terms=3, iterations=8)
        c = r['certificate']
        self.assertTrue(a.verify(c, a.ORIGINAL, True, '1/100000000'))
        self.assertEqual(r['status'], 'certified_open')
        self.assertLess(a.F(c['exact_width']), a.F(r['attempts'][0]['exact_width']))
        bad = deepcopy(c);bad['function']['amplitude'] = '1'
        self.assertFalse(a.verify(bad))

    def test_certificate_replay_does_not_search(self):
        c = a.solve(a.ORIGINAL, max_terms=3, iterations=6)['certificate']
        with patch.object(a, 'proposal', side_effect=AssertionError('search')):
            self.assertTrue(a.verify(json.loads(json.dumps(c))))
        self.assertFalse(a.verify(c, expected_mean_zero=False))

    def test_residual_frontier_finds_the_missing_fractional_power(self):
        frontier = a.residual_frontier(a.ORIGINAL, [0], [1])
        self.assertEqual(frontier[0], a.F(7,4))
        terms = a.residual_terms(a.ORIGINAL, [0], [1])
        integral = sum((x*y*a.enriched.weighted_moment(r+s)
                        for r,x in terms.items() for s,y in terms.items()), a.F(0))
        stats = a.statistics(a.ORIGINAL, [0], [1])
        self.assertEqual(integral/a.F(stats['mass']), a.F(stats['residual_squared']))

    def test_greedy_budget_and_original_problem_certificate(self):
        r = a.solve(a.ORIGINAL, max_terms=3, iterations=6,
                    basis_policy='residual', candidate_budget=2)
        self.assertEqual([x['terms'] for x in r['attempts']], [1,2,3])
        self.assertTrue(a.verify(r['certificate']))
        self.assertTrue(all(len(x['frontier_trials']) <= 2 for x in r['attempts']))

    def test_leading_cancellation_is_exact_after_every_rounding_step(self):
        ss = a.generated_powers(a.ORIGINAL, 4)
        p = a.proposal(a.ORIGINAL, ss, 48, 3, leading_cancellation=True)
        v = list(map(a.F,p['coefficients']))
        self.assertEqual(v[1], -a.F(16,21)*v[0])
        self.assertNotIn(-a.F(1,4), a.residual_terms(a.ORIGINAL,ss,v))
        self.assertTrue(p['leading_cancellation_applied'])
        r = a.solve(a.ORIGINAL,max_terms=3,iterations=4,leading_cancellation=True)
        self.assertTrue(a.verify(r['certificate']))

    def test_cancellation_requires_actual_basis_support(self):
        with self.assertRaises(ValueError):
            a.cancellation_transform(a.ORIGINAL, [0,2])
        self.assertFalse(a.proposal(a.ORIGINAL,[0],iterations=1,
                                    leading_cancellation=True)['leading_cancellation_applied'])

    def test_complete_residual_splits_with_mass_inverse(self):
        ss = a.generated_powers(a.ORIGINAL,3)
        raw = a.residual_diagnosis(a.ORIGINAL,ss,[1,0,0])
        self.assertGreater(a.F(raw['inside_squared']),0)
        p = a.proposal(a.ORIGINAL,ss,iterations=8)
        good = p['diagnosis']
        self.assertLess(a.F(good['inside_squared']), a.F(raw['inside_squared']))
        self.assertGreater(a.F(good['outside_squared']),0)
        for d in (raw,good):
            self.assertEqual(a.F(d['complete_squared']), a.F(d['inside_squared'])+a.F(d['outside_squared']))

    def test_adaptive_iteration_stops_using_complete_residual_split(self):
        p = a.proposal(a.ORIGINAL,a.generated_powers(a.ORIGINAL,3),iterations=14,
                       adaptive_iterations=True)
        self.assertLess(p['iterations'],14)
        self.assertEqual(p['diagnosis']['action'],'basis_or_gap')
        r = a.solve(a.ORIGINAL,max_terms=3,adaptive_iterations=True)
        self.assertTrue(a.verify(r['certificate']))

    def test_certified_shift_improves_early_convergence(self):
        ss = a.generated_powers(a.ORIGINAL,6)
        plain = a.proposal(a.ORIGINAL,ss,iterations=1,shift_policy='none')
        shifted = a.proposal(a.ORIGINAL,ss,iterations=1,shift_policy='certified')
        self.assertLess(a.F(shifted['diagnosis']['inside_squared']),
                        a.F(plain['diagnosis']['inside_squared']))
        source = a.source_certificate(a.ORIGINAL)
        self.assertLess(a.F(shifted['inverse_shift']), a.F(source['sectors'][1]['lower']))
        bad = deepcopy(source);bad['lower'] = '100'
        with self.assertRaises(ValueError):
            a.proposal(a.ORIGINAL,ss,iterations=1,shift_policy='certified',source=bad)

    def test_inverse_shift_does_not_change_original_operator_proof(self):
        r = a.solve(a.ORIGINAL,max_terms=3,iterations=2,shift_policy='certified')
        c = r['certificate']
        self.assertTrue(a.verify(c))
        self.assertEqual(c['function'],a.ORIGINAL)

    def test_mass_orthogonalization_exactly_preserves_full_forms(self):
        ss = a.generated_powers(a.ORIGINAL,5)
        C = a.mass_orthogonal_transform(a.ORIGINAL,ss)
        M,H,R = a.trial_matrices(a.ORIGINAL,ss)
        transformed = a._congruence(M,C)
        self.assertTrue(all(transformed[i][j] == 0 for i in range(5) for j in range(5) if i != j))
        w = [a.F(1),a.F(2),a.F(-1),a.F(1,3),a.F(0)]
        v = a._matvec(C,w)
        for form in (M,H,R):
            lhs = sum((x*y for x,y in zip(v,a._matvec(form,v))),a.F(0))
            rhs = sum((x*y for x,y in zip(w,a._matvec(a._congruence(form,C),w))),a.F(0))
            self.assertEqual(lhs,rhs)

    def test_orthogonalization_can_preserve_cancellation_constraint(self):
        ss = a.generated_powers(a.ORIGINAL,4)
        p = a.proposal(a.ORIGINAL,ss,precision_bits=32,iterations=2,
                       leading_cancellation=True,orthogonalize=True,shift_policy='certified')
        v = list(map(a.F,p['coefficients']))
        self.assertEqual(v[1], -a.F(16,21)*v[0])
        c = a.certificate(a.ORIGINAL,ss,v,a.source_certificate(a.ORIGINAL))
        self.assertTrue(a.verify(c))

    def test_full_residual_refinement_preserves_or_improves_temple_width(self):
        ss = a.generated_powers(a.ORIGINAL,4)
        p = a.proposal(a.ORIGINAL,ss,iterations=3)
        r = a.residual_refinement(a.ORIGINAL,ss,p['coefficients'],steps=2)
        source = a.source_certificate(a.ORIGINAL)
        before = a.certificate(a.ORIGINAL,ss,p['coefficients'],source)
        after = a.certificate(a.ORIGINAL,ss,r['coefficients'],source)
        self.assertTrue(a.verify(after))
        self.assertLessEqual(a.F(after['exact_width']),a.F(before['exact_width']))
        self.assertTrue(r['trace'])

    def test_residual_refinement_does_not_drop_constraint_or_budget(self):
        ss = a.generated_powers(a.ORIGINAL,4)
        p = a.proposal(a.ORIGINAL,ss,iterations=2,leading_cancellation=True,
                       orthogonalize=True,residual_steps=1)
        v = list(map(a.F,p['coefficients']))
        self.assertEqual(v[1],-a.F(16,21)*v[0])
        self.assertEqual(len(p['residual_refinement']),1)
        with self.assertRaises(ValueError):
            a.residual_refinement(a.ORIGINAL,ss,v,steps=5)

    def test_even_gap_keeps_the_complete_odd_component(self):
        form = a.direct.form_bound(a.ORIGINAL,40)
        gap = a.gap_bound(a.ORIGINAL,'even_odd')
        self.assertEqual(a.F(gap['beta']),12*a.F(form['energy_factor'])+a.F(form['mass_offset']))
        self.assertEqual(a.F(gap['odd_component_lower']),6*a.F(form['energy_factor'])+a.F(form['mass_offset']))
        ss = a.generated_powers(a.ORIGINAL,3)
        p = a.proposal(a.ORIGINAL,ss,iterations=6)
        source = a.source_certificate(a.ORIGINAL)
        full = a.certificate(a.ORIGINAL,ss,p['coefficients'],source)
        split = a.certificate(a.ORIGINAL,ss,p['coefficients'],source,gap_policy='even_odd')
        self.assertLess(a.F(split['exact_width']),a.F(full['exact_width']))
        self.assertTrue(a.verify(split))
        bad = deepcopy(split);bad['gap']['odd_component_lower']='100'
        self.assertFalse(a.verify(bad))

    def test_parity_scope_and_second_eigenvalue_are_bound(self):
        c = a.solve(a.ORIGINAL,max_terms=3,iterations=4)['certificate']
        for key,value in [('free_second_eigenvalue','6'),('space','cos_plus_sin'),
                          ('even_and_odd_cover_full_radial_m1',False)]:
            bad = deepcopy(c);bad['gap'][key]=value
            self.assertFalse(a.verify(bad))

    def test_general_amplitude_and_offset_bind_all_three_forms(self):
        q = dict(a.ORIGINAL,amplitude='-1/2',offset='-3')
        ss = a.generated_powers(q,3)
        M,A,R = a.trial_matrices(q,ss)
        base = dict(q,offset='0')
        M0,A0,R0 = a.trial_matrices(base,ss)
        self.assertEqual(M,M0)
        for i in range(3):
            for j in range(3):
                self.assertEqual(A[i][j],A0[i][j]-3*M[i][j])
                self.assertEqual(R[i][j],R0[i][j]-6*A0[i][j]+9*M[i][j])
        r = a.solve(q,max_terms=3,iterations=3)
        self.assertTrue(a.verify(r['certificate'],q,True))
        self.assertFalse(a.verify(r['certificate'],a.ORIGINAL,True))

    def test_constant_shift_reuses_bounds_without_new_search(self):
        base = dict(a.ORIGINAL,amplitude='-1/2')
        source = a.source_certificate(base)
        q = dict(base,offset='-7/3')
        with patch.object(a.direct,'certify_sector',side_effect=AssertionError('unexpected bisection')):
            shifted = a.source_certificate(q)
        self.assertEqual(a.F(shifted['lower']),a.F(source['lower'])-a.F(7,3))
        self.assertEqual(a.F(shifted['upper']),a.F(source['upper'])-a.F(7,3))
        self.assertTrue(a.direct.verify_full(shifted,q,True))

    def test_constant_potential_and_unsafe_amplitude(self):
        q = dict(a.ORIGINAL,amplitude='0',offset='-3')
        r = a.solve(q,max_terms=1,iterations=1)
        self.assertTrue(a.verify(r['certificate'],q,True))
        # The m1 witness is exact, but the independently retained m0 source
        # still has a small bisection width. It cannot be silently discarded.
        self.assertLessEqual(a.F(r['certificate']['lower']),-1)
        self.assertEqual(a.F(r['certificate']['upper']),-1)
        self.assertLess(a.F(r['certificate']['exact_width']),a.F(1,10**8))
        with self.assertRaises(ValueError):
            a.solve(dict(a.ORIGINAL,amplitude='-3'),max_terms=1)

    def test_new_exponents_generate_their_own_first_singular_correction(self):
        for alpha,amp in [('-1/8','-1'),('-1/3','-1/2'),('-3/8','-1/2')]:
            q = dict(a.ORIGINAL,exponent=alpha,amplitude=amp)
            ss = a.generated_powers(q,4)
            self.assertEqual(ss[1],2+a.F(alpha))
            self.assertGreater(ss[1],a.F(3,2))
            ratio = a.F(amp)/(ss[1]*(ss[1]-1))
            terms = a.residual_terms(q,ss,[1,ratio,0,0])
            self.assertNotIn(a.F(alpha),terms)
            p = a.proposal(q,ss,iterations=3)
            c = a.certificate(q,ss,p['coefficients'],a.source_certificate(q),gap_policy='even_odd')
            self.assertTrue(a.verify(c,q,True))

    def test_exponent_amplitude_and_axis_are_bound(self):
        q = dict(a.ORIGINAL,exponent='-1/3',amplitude='-1/2',axis=0)
        r = a.solve(q,max_terms=3,iterations=3)
        self.assertTrue(a.verify(r['certificate'],q,True))
        for update in ({'exponent':'-1/5'},{'amplitude':'1/2'},{'axis':2}):
            self.assertFalse(a.verify(r['certificate'],dict(q,**update),True))
        for alpha in ('-1/2','0','1/4'):
            with self.assertRaises(ValueError):
                a.solve(dict(q,exponent=alpha),max_terms=1)
        with self.assertRaises(ValueError):
            a.solve(q,max_terms=1,candidate_budget=0)

    def test_zero_potential_excited_trial_cannot_enter_temple_gap(self):
        q = dict(a.ORIGINAL,amplitude='0')
        # f(z)=5z^2-1 is the l=3,m=1 even radial free eigenfunction.
        stats = a.statistics(q,[0,2],[-1,5])
        self.assertEqual(a.F(stats['rayleigh']),12)
        self.assertEqual(a.F(stats['residual_squared']),0)
        with self.assertRaisesRegex(ValueError,'Temple'):
            a.certificate(q,[0,2],[-1,5],a.source_certificate(q),gap_policy='even_odd')

    def test_certificate_scope_and_matrix_digest_mutation(self):
        c = a.solve(a.ORIGINAL,max_terms=3,iterations=3)['certificate']
        for key,value in [('matrix_digest','wrong'),('full_residual_not_projected_residual',False),
                          ('mean_zero',False),('scope','axisymmetric_only')]:
            bad = deepcopy(c);bad[key]=value
            self.assertFalse(a.verify(bad))

    def test_finite_L2_form_does_not_imply_a_sufficient_temple_gap(self):
        q = dict(a.ORIGINAL,exponent='-1/3')
        self.assertLess(a.F(a.direct.form_bound(q)['eta_upper']),1)
        for policy in ('nested','residual'):
            r = a.solve(q,max_terms=3,iterations=2,basis_policy=policy,candidate_budget=1)
            self.assertTrue(a.verify(r['certificate'],q,True))
            self.assertEqual(r['status'],'certified_open')
            self.assertTrue(any(t['status']=='temple_gap_unavailable' for t in r['attempts']))


if __name__ == '__main__':
    unittest.main()
