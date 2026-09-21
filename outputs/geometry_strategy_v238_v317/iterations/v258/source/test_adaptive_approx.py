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


if __name__ == '__main__':
    unittest.main()
