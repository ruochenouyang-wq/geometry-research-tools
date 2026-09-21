from fractions import Fraction as F
import unittest
import adaptive_approx as adaptive

SINGULAR = {'kind':'axis_profile','profile':'abs_power','exponent':'-1/4','amplitude':'-1'}
STEP = {'kind':'axis_profile','profile':'step','amplitude':'3','offset':'-2'}


class ApproximationTests(unittest.TestCase):
    def test_decomposition_matches_full_integrated_certificate(self):
        for ratio in ['1/2', '2/3', '3/4']:
            terms = adaptive.error_terms(3, 5, ratio, '-3/2')
            function = dict(SINGULAR, amplitude='-3/2', offset='9/7', axis=1)
            cert = adaptive.pieces.piecewise_model(function, 5, 3, ratio)
            self.assertEqual(F(terms['error_squared']), F(cert['error_squared']))
            self.assertEqual(F(terms['core_error_squared'])+F(terms['annulus_error_squared']),
                             F(terms['error_squared']))

    def test_step_is_exact_and_bound_to_input(self):
        result = adaptive.solve(STEP)
        self.assertEqual(result['status'], 'met')
        self.assertEqual(result['certificate']['error_squared'], '0')
        self.assertTrue(adaptive.pieces.verify_piecewise(result['certificate'], STEP))

    def test_invalid_budget_and_target_are_rejected(self):
        for kwargs in [{'tolerance':'0'}, {'max_levels':0}, {'max_degree':True}, {'root_ratio':'1'}]:
            with self.assertRaises(ValueError):
                adaptive.solve(SINGULAR, **kwargs)

    def test_level_selection_is_minimal_not_rounded_logarithm(self):
        result = adaptive.solve(SINGULAR, tolerance='3/100')
        cert = result['certificate']
        self.assertEqual(result['status'], 'met')
        self.assertGreater(cert['levels'], 4)
        previous = adaptive.error_terms(cert['degree'], cert['levels']-1, cert['root_ratio'])
        self.assertGreater(F(previous['error_squared']), F(3,100)**2)
        self.assertLessEqual(result['cost']['level_comparisons'], 7)

    def test_unreachable_fixed_degree_is_proved_before_level_search(self):
        result = adaptive.solve(SINGULAR, max_degree=3)
        terms = result['diagnostics']['error_decomposition']
        self.assertEqual(result['status'], 'open')
        self.assertEqual(terms['limiting_component'], 'fixed_degree_floor')
        self.assertGreater(F(terms['infinite_levels_error_floor_squared']), F(1,10**8)**2)
        self.assertEqual(result['cost']['level_comparisons'], 0)

    def test_degree_increases_only_until_exact_target_is_attainable(self):
        result = adaptive.solve(SINGULAR, tolerance='1/10000')
        self.assertEqual(result['status'], 'met')
        self.assertGreater(result['certificate']['degree'], 3)
        self.assertLess(result['certificate']['degree'], 24)
        self.assertEqual(result['cost']['certificate_builds'], 1)
        self.assertTrue(all(a['limiting_component'] == 'fixed_degree_floor'
                            for a in result['attempts'][:-1]))


if __name__ == '__main__':
    unittest.main()
