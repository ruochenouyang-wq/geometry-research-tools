"""Independent Temple, projection, continuation and full-sphere scope review."""
from copy import deepcopy
from fractions import Fraction as F
import unittest
import precision as p
import witness as w
from common import projected


class PrecisionIndependentReview(unittest.TestCase):
    def test_temple_matches_analytic_projected_variance(self):
        trial = w.rayleigh_certificate([0, 1], [1], degrees=[1])
        separator = projected.certify_sector([0, 1], m=0, mean_zero=True,
                                             k=2, modes=4, bits=20)
        cert = p.temple_certificate(trial, separator)
        g = F(separator['lower'])
        self.assertGreater(g, 2)
        self.assertEqual(F(cert['residual_variance']), F(4, 15))
        self.assertEqual(F(cert['lower']), 2-F(4, 15)/(g-2))
        self.assertEqual(F(cert['upper']), 2)
        self.assertEqual(cert['complete_residual']['residual_coefficients'],
                         [{'degree': 2, 'coefficient': '2/3'}])
        self.assertTrue(p.verify_temple(cert))

    def test_zero_and_negative_gap_never_produce_temple_lower(self):
        separator = projected.certify_sector([0], m=0, k=2, modes=2, bits=8)
        for l in (2, 3):
            trial = w.rayleigh_certificate([0], [1], degrees=[l])
            cert = p.temple_certificate(trial, separator)
            self.assertTrue(p.verify_temple(cert))
            self.assertFalse(cert['gap_accepted'])
            self.assertIsNone(cert['lower'])
            self.assertEqual(cert['status'], 'gap_failed_rayleigh_only')
            self.assertEqual(F(cert['upper']), l*(l+1))

    def test_same_operator_separator_binding(self):
        trial = w.rayleigh_certificate([0], [1], m=0, degrees=[1])
        wrong = [projected.certify_sector([0], m=1, k=2, modes=2, bits=8),
                 projected.certify_sector([0], m=0, mean_zero=False, k=2, modes=2, bits=8),
                 projected.certify_sector([1], m=0, k=2, modes=2, bits=8),
                 projected.certify_sector([0], m=0, k=1, modes=2, bits=8)]
        for separator in wrong:
            with self.assertRaises(ValueError):
                p.temple_certificate(trial, separator)

    def test_gap_recovery_actually_strengthens_separator(self):
        trial = w.rayleigh_certificate([0], [1], degrees=[1])
        kernel = projected.Kernel([0], modes=2)
        weak = projected.sector_certificate(kernel, 2, 0, 6, 'upper_ritz')
        self.assertFalse(p.temple_certificate(trial, weak)['gap_accepted'])
        recovered = p.recover_gap(trial, modes=2, bits=8, previous_gap=weak)
        self.assertTrue(recovered['gap_accepted'])
        self.assertEqual(F(recovered['lower']), 2)

    def test_temple_scope_and_residual_tampering(self):
        trial = w.rayleigh_certificate([0, 1], [1], degrees=[1])
        cert = p.recover_gap(trial, modes=3, bits=16)
        for key, value in [('scope', 'all_mean_zero_H1_functions_on_unit_S2'),
                           ('lower', '2'), ('residual_variance', '0'),
                           ('next_eigenvalue_index', 1), ('mean_zero', False)]:
            bad = deepcopy(cert)
            bad[key] = value
            self.assertFalse(p.verify_temple(bad), key)
        bad = deepcopy(cert)
        bad['complete_residual']['residual_coefficients'] = []
        self.assertFalse(p.verify_temple(bad))

    def test_residual_window_is_not_sparse_support(self):
        trial = w.rayleigh_certificate([0, 0, 2], [1], m=1, degrees=[1])
        sparse = w.residual_certificate(trial)
        self.assertGreater(F(sparse['outside_support_residual_squared']), 0)
        inside = p.residual_schedule(trial, modes=3)
        outside = p.residual_schedule(trial, modes=1)
        self.assertEqual(F(inside['window_tail_variance']), 0)
        self.assertEqual(F(outside['window_finite_variance']), 0)
        self.assertEqual(inside['action'], 'increase_rational_coefficient_bits')
        self.assertEqual(outside['action'], 'increase_radial_modes')

    def test_budget_exit_keeps_omitted_angular_lower(self):
        cert = p.full_ground([0, 1, 1], modes=2, max_modes=2, max_m=0, max_steps=1)
        self.assertTrue(p.verify_full(cert, expected_q=[0, 1, 1], expected_mean_zero=True))
        self.assertEqual(cert['status'], 'certified_bound_open_gap')
        self.assertEqual(F(cert['lower']), F(7, 4))
        self.assertLess(F(cert['lower']), F(cert['sectors'][0]['lower']))
        self.assertLessEqual(F(cert['upper']), F(13, 5))

    def test_full_aggregation_cannot_promote_trial_scope(self):
        cert = p.full_ground([0, 0, 2], modes=2, max_modes=2, max_steps=2)
        self.assertTrue(p.verify_full(cert))
        trial = w.rayleigh_certificate([0, 0, 2], [1], m=1)
        self.assertFalse(p.verify_sector(trial))
        bad = deepcopy(cert)
        bad['sectors'][1] = trial
        self.assertFalse(p.verify_full(bad))
        for key, value in [('scope', 'single_azimuth_sector_only'), ('lower', '100'),
                           ('angular_tail_lower', '1000')]:
            bad = deepcopy(cert)
            bad[key] = value
            self.assertFalse(p.verify_full(bad), key)
        bad = deepcopy(cert)
        bad['sectors'] = bad['sectors'][1:]
        self.assertFalse(p.verify_full(bad))

    def test_checkpoint_replays_binding_and_retains_enclosure(self):
        old = p.full_ground([0, 1, 1], modes=2, max_modes=3, max_steps=1)
        saved = p.checkpoint(old)
        self.assertTrue(p.verify_checkpoint(saved, [0, 1, 1], True))
        resumed = p.full_ground([0, 1, 1], modes=2, max_modes=3, max_steps=2, saved=saved)
        self.assertTrue(p.verify_full(resumed))
        self.assertGreaterEqual(F(resumed['lower']), F(old['lower']))
        self.assertLessEqual(F(resumed['upper']), F(old['upper']))
        for q, projection in [([0, 1], True), ([0, 1, 1], False)]:
            with self.assertRaises(ValueError):
                p.full_ground(q, mean_zero=projection, saved=saved)
        for key, value in [('scope', 'single_azimuth_sector_only'), ('lower', '100'),
                           ('certificate_digest', '0'*64)]:
            bad = deepcopy(saved)
            bad[key] = value
            self.assertFalse(p.verify_checkpoint(bad), key)

    def test_refinement_and_intersection_preserve_actual_bounds(self):
        rough = p.refine_coefficients([0, 0, 2], m=1, modes=4, coefficient_bits=8, iterations=0)
        refined = p.refine_coefficients([0, 0, 2], m=1, modes=4, coefficient_bits=64, iterations=4)
        self.assertTrue(w.verify(refined))
        self.assertLess(F(w.residual_certificate(refined)['normalized_residual_squared']),
                        F(w.residual_certificate(rough)['normalized_residual_squared']))
        schur = projected.certify_sector([0, 0, 2], m=1, modes=4, bits=20)
        combined = p.intersect_bounds(schur, p.recover_gap(refined, modes=4, bits=24))
        self.assertTrue(p.verify_sector(combined))
        self.assertGreaterEqual(F(combined['lower']), F(schur['lower']))
        self.assertLessEqual(F(combined['upper']), F(schur['upper']))


if __name__ == '__main__':
    unittest.main()
