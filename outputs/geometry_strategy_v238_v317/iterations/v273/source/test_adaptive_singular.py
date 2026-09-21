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
        plain = a.proposal(a.ORIGINAL,ss,iterations=1)
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


if __name__ == '__main__':
    unittest.main()
