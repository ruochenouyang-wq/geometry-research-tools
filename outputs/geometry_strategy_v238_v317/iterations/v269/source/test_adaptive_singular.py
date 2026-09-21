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


if __name__ == '__main__':
    unittest.main()
