"""Independent GN review: analytic calibrations, scope and complete residuals."""
import copy
import json
from pathlib import Path
from fractions import Fraction as F
import unittest
import gn_variation as g
import wide_potential as wide

HERE = Path(__file__).resolve().parent


def certificate_nodes(value):
    if isinstance(value, dict):
        if str(value.get('format', '')).startswith('gn_'):
            yield value
        for child in value.values():
            yield from certificate_nodes(child)
    elif isinstance(value, list):
        for child in value:
            yield from certificate_nodes(child)


class IndependentGNReview(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plane = g.plane([[0, 1], [0, F(-3, 2), 0, F(5, 2)]])
        cls.diagnostic = g.spectral_diagnostic([0, 1], modes=4, max_modes=4,
                            bits=24, max_m=4, tolerance=F(1, 1000))

    def test_area_measure_factor_four_independent_integrals(self):
        # For t=x3: area mass=4pi/3, area energy=8pi/3,
        # area quartic=4pi/5. Hence pi*K=(4/5)/((4/3)*(8/3)).
        expected = F(4, 5)/(F(4, 3)*F(8, 3))
        c = g.ratio([0, 1])
        self.assertEqual(expected, F(9, 40))
        self.assertEqual(F(c['pi_K_lower']), expected)
        bad = copy.deepcopy(c)
        bad['pi_K_lower'] = str(4*expected)
        self.assertFalse(g.verify(bad))

    def test_euler_residual_closed_form_for_coordinate(self):
        c = g.residual([0, 1])
        self.assertEqual(list(map(F, c['residual'])), [0, 4, 0, F(-20, 3)])
        self.assertEqual(F(c['norm_squared']), F(64, 63))
        self.assertEqual(F(c['inner_u']), 0)
        self.assertFalse(c['full_stationary'])
        bad = copy.deepcopy(c)
        bad['residual'] = ['0', '4']
        self.assertFalse(g.verify(bad))

    def test_complete_degree15_euler_residual_for_t5(self):
        c = g.residual([0]*5+[1])
        expected = [F(0)]*16
        expected[3], expected[5], expected[15] = F(-20), F(320, 9), F(-700, 33)
        self.assertEqual(list(map(F, c['residual'])), expected)
        self.assertEqual(c['full_degree'], 15)
        self.assertTrue(g.verify(c))

    def test_projected_stationarity_not_stationarity(self):
        c = g.project_residual([0, 1], 1)
        self.assertTrue(c['projected_stationary'])
        self.assertFalse(c['full_stationary'])
        self.assertEqual(list(map(F, c['legendre_coefficients'])), [0, 0, 0, F(-8, 3)])
        self.assertEqual(F(c['outside_norm_squared']), F(64, 63))
        result = g.refine([0, 1], iterations=2, degree_budget=1)
        self.assertEqual(result['stop'], 'projected_stationary_with_nonzero_full_residual')
        bad = copy.deepcopy(result)
        bad['stop'] = 'full_stationary'
        self.assertFalse(g.verify(bad))

    def test_local_derivative_does_not_certify_arbitrary_step(self):
        c = g.direction([0, 1])
        self.assertEqual(F(c['derivative']), F(24, 35))
        v = list(map(F, c['v']))
        large_step = g.add([0, 1], g.scale(v, 100))
        self.assertLess(F(g.ratio(large_step)['pi_K_lower']), F(9, 40))
        self.assertIn('not_a_finite_step_guarantee', c['claim'])

    def test_plane_ceiling_is_not_full_H1_ceiling(self):
        self.assertTrue(g.verify_plane(self.plane))
        outside = g.ratio(g.centered_cap(12))
        self.assertGreater(F(outside['pi_K_lower']), F(self.plane['upper']))
        bad = copy.deepcopy(self.plane)
        bad['universal_GN_upper_bound'] = True
        self.assertFalse(g.verify(bad))

    def test_complete_projective_chart_and_sign_bindings(self):
        for coefficients in [(1, 10**6), (10**6, 1), (-1, 10**6), (0, -1)]:
            f, h = [list(map(F, b)) for b in self.plane['basis']]
            trial = g.add(g.scale(f, coefficients[0]), g.scale(h, coefficients[1]))
            self.assertLessEqual(F(g.ratio(trial)['pi_K_lower']), F(self.plane['upper']))
        bad = copy.deepcopy(self.plane)
        bad['charts'][1] = copy.deepcopy(bad['charts'][0])
        self.assertFalse(g.verify(bad))
        bad = copy.deepcopy(self.plane)
        bad['charts'][0]['upper_nonpositive']['positive_multiplier'] = '-1'
        self.assertFalse(g.verify(bad))

    def test_compile_coordinate_exact_operator_identity(self):
        c = g.compile_potential([0, 1])
        self.assertEqual(list(map(F, c['q_coefficients'])), [0, 0, F(-20, 3)])
        self.assertEqual(F(c['trial_rayleigh']), -2)
        self.assertEqual(F(c['potential_energy']), F(-4, 3))
        self.assertTrue(g.verify(c, expected_u=[0, 1]))
        self.assertFalse(g.verify(c, expected_u=[0, 2]))

    def test_compiled_potential_and_wide_spectral_proof_bound_together(self):
        c = self.diagnostic
        self.assertIn('spectral', c)
        self.assertTrue(g.verify(c, expected_u=[0, 1]))
        wrong = wide.full_ground([0], mean_zero=True, modes=2, max_modes=2,
                                 bits=16, max_m=2, tolerance=F(1, 100))
        self.assertTrue(wide.verify_full(wrong))
        bad = copy.deepcopy(c)
        bad['spectral'] = wrong
        self.assertFalse(g.verify(bad))
        bad = copy.deepcopy(c)
        bad['spectral']['mean_zero'] = False
        self.assertFalse(g.verify(bad))
        bad = copy.deepcopy(c)
        bad['universal_GN_claim'] = True
        self.assertFalse(g.verify(bad))

    def test_budget_diagnostic_makes_no_mathematical_success_claim(self):
        c = g.spectral_diagnostic([0, 1], modes=0)
        self.assertEqual(c['status'], 'rejected_or_budget_failed')
        self.assertNotIn('spectral', c)
        self.assertFalse(c['universal_GN_claim'])
        self.assertTrue(g.verify(c))
        bad = copy.deepcopy(c)
        bad['status'] = 'lower_spectral_ground_separated'
        self.assertFalse(g.verify(bad))

    def test_saved_certificate_roundtrips(self):
        found = 0
        formats = set()
        for path in sorted((HERE/'round08').rglob('*.json')):
            if path.name == 'INDEPENDENT_TESTS.json':
                continue
            for c in certificate_nodes(json.loads(path.read_text())):
                self.assertTrue(g.verify(c), str(path)+':'+c['format'])
                found += 1
                formats.add(c['format'])
        self.assertGreaterEqual(found, 8, 'Need actual saved mathematical certificates')
        self.assertGreaterEqual(len(formats), 5)


if __name__ == '__main__':
    unittest.main()
