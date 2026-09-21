import copy
import json
import unittest
from fractions import Fraction as F
import witness as w
from common import ROOT


def source():
    return json.loads((ROOT/'round01/baseline.json').read_text())['old_certificate']['evidence']


class WitnessStages(unittest.TestCase):
    def setUp(self):
        self.one = w.rayleigh_certificate([0, 0, 2], [1], m=1, threshold=F(12, 5))
        self.simple = w.rayleigh_certificate([0, 0, 2], [1, F(-1, 36)], m=1,
                                            degrees=[1, 3], threshold=F(12, 5))

    def test_v95_exact_trial(self):
        self.assertTrue(w.verify(self.simple))
        self.assertEqual(self.simple['rayleigh_quotient'], '722/303')
        self.assertEqual(self.simple['strict_margin'], '26/1515')
        self.assertEqual(self.simple['norm_squared_divided_by_azimuth_factor'], '505/378')
        self.assertEqual(self.simple['energy_divided_by_azimuth_factor'], '1805/567')
        tiny = w.rayleigh_certificate([0], [F(1, 10**300)], m=3)
        self.assertEqual(tiny['rayleigh_quotient'], '12')
        tinier = w.rayleigh_certificate([0], [F(1, 10**3000)], m=1)
        self.assertTrue(w.verify(tinier))
        self.assertEqual(tinier['rayleigh_quotient'], '2')

    def test_v96_sector_selection(self):
        selection = w.select_sector(source(), F(12, 5))
        self.assertTrue(w.verify(selection))
        self.assertEqual(selection['selected_m'], 1)
        self.assertTrue(selection['spectral_existence_below_threshold'])

    def test_v97_mass_normalized_proposal(self):
        proposal = w.mass_normalized_proposal([0, 0, 2], m=1, modes=4)
        self.assertAlmostEqual(sum(x*x for x in proposal['mass_normalized_eigenvector']), 1)
        self.assertLess(proposal['numerical_rayleigh'], 2.4)
        self.assertEqual(proposal['degrees'], [1, 2, 3, 4])

    def test_v98_quantization_keeps_best(self):
        result = w.quantize_proposal(w.mass_normalized_proposal([0, 0, 2], m=1, modes=4), threshold=F(12, 5))
        self.assertTrue(w.verify(result['certificate']))
        best = F(result['certificate']['rayleigh_quotient'])
        for row in result['trace']:
            self.assertTrue(w.verify(row['certificate']))
            self.assertLessEqual(best, F(row['certificate']['rayleigh_quotient']))
        self.assertEqual(result['trace'][1]['status'], 'retained_old_trial')

    def test_v99_sparse_acceptance_and_rejection(self):
        noisy = w.rayleigh_certificate([0, 0, 2], [1, F(-1, 36), F(1, 1000)],
                                       m=1, degrees=[1, 3, 7], threshold=F(12, 5))
        result = w.sparsify(noisy, (1, 2))
        self.assertEqual(result['trace'][0]['status'], 'rejected_no_gain')
        self.assertEqual(result['trace'][1]['status'], 'accepted')
        self.assertLess(F(result['certificate']['rayleigh_quotient']), F(noisy['rayleigh_quotient']))

    def test_v100_complete_projected_residual(self):
        candidate = w.rayleigh_certificate([0, 1], [1], m=0)
        residual = w.residual_certificate(candidate)
        self.assertTrue(w.verify(residual))
        self.assertEqual(residual['removed_l0_coefficient'], '1/3')
        self.assertEqual(residual['residual_coefficients'], [{'degree': 2, 'coefficient': '2/3'}])
        self.assertEqual(residual['normalized_residual_squared'], '4/15')
        self.assertEqual(residual['outside_support_residual_squared'], '4/15')
        self.assertEqual(residual['orthogonality_to_trial'], '0')
        unprojected = w.residual_certificate(w.rayleigh_certificate([0, 1], [1], m=0,
                                                                   mean_zero=False, degrees=[1]))
        self.assertEqual(unprojected['normalized_residual_squared'], '3/5')
        high = w.residual_certificate(w.rayleigh_certificate([0]*6+[1], [1], m=1))
        self.assertEqual(max(row['degree'] for row in high['residual_coefficients']), 7)

    def test_v101_residual_driven_support_growth(self):
        result = w.residual_enrich(self.one)
        self.assertEqual(result['added_degrees'], [3])
        self.assertEqual(result['status'], 'accepted_improvement')
        self.assertLess(F(result['certificate']['rayleigh_quotient']), F(self.one['rayleigh_quotient']))
        exact = w.rayleigh_certificate([0], [1], m=1)
        self.assertEqual(w.residual_enrich(exact)['status'], 'support_did_not_grow')

    def test_v102_two_trial_ritz_scope_and_actual_upper(self):
        second = w.rayleigh_certificate([0, 0, 2], [1], m=1, degrees=[3], threshold=F(12, 5))
        result = w.two_trial_ritz(self.one, second)
        self.assertTrue(w.verify(result))
        self.assertTrue(w.verify(result['certificate']))
        self.assertLess(F(result['ritz_gap']), F(1, 10**8))
        self.assertLess(F(result['actual_trial_upper']), F(12, 5))
        bad = copy.deepcopy(result)
        bad['scope'] = 'full_sphere_ground_lower'
        self.assertFalse(w.verify(bad))
        with self.assertRaises(ValueError):
            w.two_trial_ritz(self.one, self.one)

    def test_v103_negative_residual_line_search(self):
        result = w.residual_line_search(self.one)
        self.assertEqual(result['status'], 'accepted_improvement')
        self.assertLess(F(result['certificate']['rayleigh_quotient']), F(self.one['rayleigh_quotient']))
        self.assertLess(F(result['derivative_numerator_divided_by_two'][0]), 0)
        self.assertTrue(all(w.verify(row['certificate']) for row in result['trace']))
        zero = w.residual_line_search(w.rayleigh_certificate([0], [1], m=1))
        self.assertEqual(zero['status'], 'zero_residual_or_degree_budget')

    def test_v104_explicit_refutation_and_honest_not_found(self):
        result = w.find_counterexample([0, 0, 2], F(12, 5), full_certificate=source())
        self.assertTrue(w.verify(result, expected_q=[0, 0, 2], expected_m=1,
                                  expected_mean_zero=True, expected_threshold=F(12, 5)))
        self.assertEqual(result['status'], 'explicit_counterexample_verified')
        self.assertTrue(result['explicit_real_harmonic_expansion'])
        failed = w.find_counterexample([0], F(1), modes=2, max_steps=0)
        self.assertTrue(w.verify(failed))
        self.assertEqual(failed['status'], 'not_found')

    def test_binding_scope_and_tampering(self):
        for field in ['scope', 'norm_squared_divided_by_azimuth_factor', 'rayleigh_quotient',
                      'full_space_spectral_upper', 'status', 'strict_margin', 'azimuth_integral_factor']:
            bad = copy.deepcopy(self.simple)
            bad[field] = 'wrong'
            self.assertFalse(w.verify(bad), field)
        for options in [dict(expected_q=[0, 0, 3]), dict(expected_m=0),
                        dict(expected_mean_zero=False), dict(expected_coefficients=[1, 2]),
                        dict(expected_threshold=F(5, 2)), dict(expected_m=True)]:
            self.assertFalse(w.verify(self.simple, **options))
        bad = w.residual_certificate(self.simple)
        bad['residual_coefficients'] = []
        self.assertFalse(w.verify(bad))

    def test_input_rejection(self):
        for kwargs in [dict(coefficients=[0]), dict(coefficients=[1.0]),
                       dict(coefficients=[1], m=True), dict(coefficients=[1], mean_zero=1),
                       dict(coefficients=[1], degrees=[0]), dict(coefficients=[1], degrees=[251])]:
            with self.assertRaises((ValueError, TypeError)):
                w.rayleigh_certificate([0], **kwargs)
        with self.assertRaises(ValueError):
            w.rayleigh_certificate([0]*7+[1], [1], m=1)
        with self.assertRaises(ValueError):
            w.find_counterexample([0, 0, 3], F(12, 5), full_certificate=source())


if __name__ == '__main__':
    unittest.main()
