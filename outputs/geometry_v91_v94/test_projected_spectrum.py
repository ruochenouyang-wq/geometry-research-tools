"""Mathematical invariance, tail, and exact-limit regression checks."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import projected_spectrum as p
import v03_tail


class ProjectedSpectrumTests(unittest.TestCase):
    def test_constant_potential_exact_spectra_and_mean_zero_indexing(self):
        for m, zero, k in [(0, False, 1), (0, True, 1), (0, True, 3),
                           (1, True, 3), (2, False, 2)]:
            with self.subTest(m=m, zero=zero, k=k):
                cert = p.certify_sector([F(-7, 3)], m=m, mean_zero=zero, k=k, modes=4)
                start = 1 if m == 0 and zero else m
                l = start+k-1
                expected = l*(l+1)-F(7, 3)
                self.assertEqual(F(cert['lower']), expected)
                self.assertEqual(F(cert['upper']), expected)
                self.assertTrue(p.verify_sector(cert))
        self.assertEqual(p.full_ground([F(-7, 3)])['lower'], '-1/3')

    def test_additive_potential_shift_preserves_all_full_spectral_endpoints(self):
        a = p.full_ground([0, 1, 1])
        b = p.full_ground([F(7, 3), 1, 1])
        for key in ('lower', 'upper', 'angular_tail_lower'):
            self.assertEqual(F(b[key])-F(a[key]), F(7, 3))
        self.assertEqual(a['exact_width'], b['exact_width'])

    def test_equator_reflection_is_unitary_even_with_projection(self):
        a = p.full_ground([0, 1, 1])
        b = p.full_ground([0, -1, 1])
        self.assertEqual((a['lower'], a['upper']), (b['lower'], b['upper']))

    def test_degree_six_sector_agrees_with_frozen_unconstrained_core(self):
        q = [F(1, 7), F(1, 3), 0, F(-1, 5), 0, 0, F(1, 11)]
        cert = p.certify_sector(q, m=0, mean_zero=False, k=2, modes=12)
        old = v03_tail.certify(q, k=2, modes=12, bits=44)
        self.assertLessEqual(max(F(cert['lower']), F(old['lower'])),
                             min(F(cert['upper']), F(old['upper'])))
        self.assertLess(F(cert['exact_width']), F(1, 10**9))

    def test_global_projection_cannot_be_relabelled_as_unconstrained(self):
        cert = p.full_ground([0, 1])
        self.assertTrue(p.verify_full(cert, expected_q=[0, 1], expected_mean_zero=True))
        self.assertFalse(p.verify_full(cert, expected_mean_zero=False))
        bad = deepcopy(cert)
        bad['mean_zero'] = False
        bad['scope'] = p.scope(False)
        self.assertFalse(p.verify_full(bad))

    def test_angular_budget_exhaustion_still_encloses_full_ground(self):
        limited = p.full_ground([0, 0, 2], max_m=0)
        complete = p.full_ground([0, 0, 2])
        self.assertTrue(p.verify_full(limited))
        self.assertFalse(limited['angular_tail_cannot_improve_best_upper'])
        self.assertEqual(limited['status'], 'certified_bound_open_gap')
        self.assertLessEqual(F(limited['lower']), F(complete['lower']))
        self.assertGreaterEqual(F(limited['upper']), F(complete['upper']))
        self.assertEqual(limited['lower'], '2')

    def test_radial_budget_does_not_claim_more_accuracy_than_achieved(self):
        limited = p.full_ground([0, 0, 100], modes=1, max_modes=1, max_m=1)
        complete = p.full_ground([0, 0, 100])
        self.assertEqual(limited['status'], 'certified_bound_open_gap')
        self.assertTrue(p.verify_full(limited))
        self.assertLessEqual(F(limited['lower']), F(complete['lower']))
        self.assertGreaterEqual(F(limited['upper']), F(complete['upper']))

    def test_potential_monotonicity_on_constrained_full_space(self):
        no_potential = p.full_ground([0])
        positive = p.full_ground([0, 0, 2])
        self.assertGreater(F(positive['lower']), F(no_potential['upper']))
        # Multiplication by 2*t^2 lies between 0 and 2.
        self.assertLess(F(positive['upper']), F(no_potential['upper'])+2)

    def test_verification_replays_evidence_without_search(self):
        cert = p.weighted_poincare([0, 0, 2], '119/50')
        with patch.object(p, 'certify_sector', side_effect=AssertionError('search called')), \
             patch.object(p, 'full_ground', side_effect=AssertionError('search called')):
            self.assertTrue(p.verify_inequality(cert, expected_q=[0, 0, 2], expected_threshold='119/50'))

    def test_weighted_inequality_proof_refutation_and_undetermined_are_distinct(self):
        proved = p.weighted_poincare([0, 0, 2], '119/50')
        refuted = p.weighted_poincare([0, 0, 2], '12/5')
        open_case = p.weighted_poincare([0, 0, 2], '12/5', max_m=0)
        self.assertEqual(proved['status'], 'proved')
        self.assertEqual(refuted['status'], 'refuted_by_spectral_existence')
        self.assertEqual(open_case['status'], 'undetermined')
        for cert in (proved, refuted, open_case):
            self.assertTrue(p.verify_inequality(cert))
            self.assertTrue(cert['not_a_universal_Gagliardo_Nirenberg_constant_certificate'])
        self.assertTrue(refuted['refutation_does_not_include_an_explicit_function'])

    def test_reject_different_external_problem_or_threshold(self):
        cert = p.weighted_poincare([0, 0, 2], '119/50')
        self.assertFalse(p.verify_inequality(cert, expected_q=[0, 0, 4]))
        self.assertFalse(p.verify_inequality(cert, expected_threshold='12/5'))
        sector = p.certify_sector([0, 1], m=0, mean_zero=True)
        self.assertFalse(p.verify_sector(sector, expected_k=2))
        self.assertFalse(p.verify_sector(sector, expected_m=1))

    def test_sparse_or_reordered_sector_coverage_rejected(self):
        cert = p.full_ground([0, 0, 2])
        self.assertEqual(len(cert['sectors']), 2)
        for sectors in [cert['sectors'][1:], list(reversed(cert['sectors'])),
                        [cert['sectors'][0], cert['sectors'][0]]]:
            bad = deepcopy(cert)
            bad['sectors'] = sectors
            self.assertFalse(p.verify_full(bad))

    def test_strict_tail_positivity_and_evidence_mutations(self):
        kernel = p.Kernel([0, 0, 100], modes=1)
        with self.assertRaises(ValueError):
            kernel.matrix(kernel.beta, 'lower')
        cert = p.full_ground([0, 0, 2])
        mutations = {'lower': '99', 'angular_tail_lower': '99', 'omitted_azimuth_m_start': 100,
                     'mean_zero': 1, 'status': 'solved_open_problem', 'formal_assistant_checked': True}
        for key, value in mutations.items():
            with self.subTest(field=key):
                bad = deepcopy(cert)
                bad[key] = value
                self.assertFalse(p.verify_full(bad))

    def test_invalid_geometry_controls_and_nonrational_potentials(self):
        for kwargs in [dict(m=-1), dict(m=True), dict(mean_zero=1), dict(k=0),
                       dict(modes=0), dict(bits=7)]:
            with self.assertRaises(ValueError):
                p.certify_sector([0, 1], **kwargs)
        for q in [[0, .5], [False, 1], [0]*8]:
            with self.assertRaises(ValueError):
                p.full_ground(q)


if __name__ == '__main__':
    unittest.main()
