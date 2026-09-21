from fractions import Fraction as F
import copy
import unittest

import precision as p
import witness
from common import projected


class PrecisionTests(unittest.TestCase):
    def test_absolute_bits(self):
        self.assertGreater(p.choose_bits([0, 0, 2], F(1, 10**18)), 44)
        self.assertGreater(p.choose_bits([0, 0, 100], F(1, 10**18)),
                           p.choose_bits([0, 0, 2], F(1, 10**18)))
        self.assertGreater(p.choose_bits([0, 0, 2], F(1, 10**24)),
                           p.choose_bits([0, 0, 2], F(1, 10**18)))

    def test_diagnostic_bit_then_radial(self):
        cert = projected.certify_sector([0, 0, 2], m=1, modes=4, bits=44)
        self.assertEqual(p.diagnose_precision(cert, F(1, 10**18), 44)['action'], 'increase_bisection_bits')
        self.assertEqual(p.diagnose_precision(cert, F(1, 10**18), 100)['action'], 'increase_radial_modes')

    def test_exact_temple(self):
        trial = witness.rayleigh_certificate([0], [1], m=1)
        gap = projected.certify_sector([0], m=1, k=2, modes=2)
        cert = p.temple_certificate(trial, gap)
        self.assertTrue(p.verify_temple(cert))
        self.assertEqual(F(cert['lower']), 2)
        self.assertEqual(F(cert['upper']), 2)

    def test_gap_failure_has_no_lower(self):
        trial = witness.rayleigh_certificate([0], [1], m=1, degrees=[3])
        cert = p.recover_gap(trial, 4, 64)
        self.assertTrue(p.verify_temple(cert))
        self.assertFalse(cert['gap_accepted'])
        self.assertIsNone(cert['lower'])

    def test_wrong_gap_binding_rejected(self):
        trial = witness.rayleigh_certificate([0], [1], m=1)
        for gap in (projected.certify_sector([0], m=0, k=2, modes=2),
                    projected.certify_sector([0], m=1, k=1, modes=2),
                    projected.certify_sector([1], m=1, k=2, modes=2)):
            with self.assertRaises(ValueError):
                p.temple_certificate(trial, gap)

    def test_complete_projected_residual_binding(self):
        trial = witness.rayleigh_certificate([0, 1], [1], m=0, mean_zero=True)
        cert = p.recover_gap(trial, 4, 64)
        self.assertEqual(F(cert['complete_residual']['removed_l0_coefficient']), F(1, 3))
        self.assertEqual(F(cert['residual_variance']), F(4, 15))
        cert['complete_residual']['normalized_residual_squared'] = '0'
        self.assertFalse(p.verify_temple(cert))

    def test_radial_vs_finite_residual(self):
        tail = witness.rayleigh_certificate([0, 0, 2], [1], m=1)
        self.assertEqual(p.residual_schedule(tail, 1)['action'], 'increase_radial_modes')
        finite = witness.rayleigh_certificate([0], [1, F(1, 10)], m=1, degrees=[1, 3])
        decision = p.residual_schedule(finite, 3)
        self.assertEqual(decision['action'], 'increase_rational_coefficient_bits')
        self.assertEqual(F(decision['window_tail_variance']), 0)

    def test_coefficient_refinement(self):
        rough = p.refine_coefficients([0, 0, 2], 1, True, 8, 8, 0)
        fine = p.refine_coefficients([0, 0, 2], 1, True, 8, 80, 8)
        rough_var = F(witness.residual_certificate(rough)['normalized_residual_squared'])
        fine_var = F(witness.residual_certificate(fine)['normalized_residual_squared'])
        self.assertLess(fine_var, rough_var)
        self.assertTrue(witness.verify(fine))

    def test_intersection_never_widens(self):
        schur = projected.certify_sector([0, 0, 2], m=1, modes=6, bits=64)
        trial = p.refine_coefficients([0, 0, 2], 1, True, 6, 80)
        temple = p.recover_gap(trial, 6, 64)
        cert = p.intersect_bounds(schur, temple)
        self.assertTrue(p.verify_sector(cert))
        self.assertGreaterEqual(F(cert['lower']), F(schur['lower']))
        self.assertLessEqual(F(cert['upper']), F(schur['upper']))

    def test_precision_target_and_selective_refinement(self):
        cert = p.full_ground([0, 0, 2])
        self.assertTrue(p.verify_full(cert, [0, 0, 2], True))
        self.assertEqual(cert['status'], 'certified_target_met')
        self.assertLessEqual(F(cert['exact_width']), F(1, 10**18))
        self.assertEqual(cert['sectors'][0]['schur']['modes'], 4)
        self.assertLess(cert['sectors'][1]['schur']['modes'], 32)

    def test_zero_budget_open_is_valid(self):
        cert = p.full_ground([0, 0, 2], max_steps=0)
        self.assertTrue(p.verify_full(cert))
        self.assertEqual(cert['status'], 'certified_bound_open_gap')

    def test_checkpoint_monotonic_and_binding(self):
        old = p.full_ground([0, 1, 1], max_steps=1)
        saved = p.checkpoint(old)
        self.assertTrue(p.verify_checkpoint(saved, [0, 1, 1], True))
        new = p.full_ground([0, 1, 1], saved=saved)
        self.assertGreaterEqual(F(new['lower']), F(old['lower']))
        self.assertLessEqual(F(new['upper']), F(old['upper']))
        for q, mean_zero in (([0, 0, 2], True), ([0, 1, 1], False)):
            with self.assertRaises(ValueError):
                p.full_ground(q, mean_zero=mean_zero, saved=saved)
        tampered = copy.deepcopy(saved)
        tampered['lower'] = '100'
        self.assertFalse(p.verify_checkpoint(tampered))

    def test_sector_gap_and_status_tampering(self):
        cert = p.full_ground([0, 0, 2], max_steps=2)
        for mutate in ('omit', 'wrong_m', 'lower', 'status'):
            bad = copy.deepcopy(cert)
            if mutate == 'omit':
                bad['sectors'] = bad['sectors'][1:]
            elif mutate == 'wrong_m':
                bad['sectors'][0]['azimuth_m'] = 1
            elif mutate == 'lower':
                bad['lower'] = '100'
            else:
                bad['status'] = 'certified_target_met'
            self.assertFalse(p.verify_full(bad))

    def test_parameter_validation(self):
        for kwargs in ({'max_steps': True}, {'max_modes': 2}, {'max_m': -1},
                       {'tolerance': 0}, {'use_temple': 'yes'}):
            with self.assertRaises(ValueError):
                p.full_ground([0, 0, 2], **kwargs)


if __name__ == '__main__':
    unittest.main()
