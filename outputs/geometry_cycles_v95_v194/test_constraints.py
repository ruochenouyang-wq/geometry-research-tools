"""Acceptance and adversarial checks for V145--V154."""
from copy import deepcopy
from fractions import Fraction as F
import json
import unittest
import constraints as c


class ConstraintTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.full = c.full_exclusion_ground([0, 0, 2], cutoff=1, modes=8, bits=40)
        cls.moment = c.certify_moment_sector([0, 0, 2], [[1, 1, 1]], modes=8, bits=40)
        cls.q0 = c.full_exclusion_ground([0], cutoff=1)

    def test_analytic_exclusion_exact_values(self):
        for cutoff, value in [(-1, 0), (0, 2), (1, 6), (2, 12)]:
            cert = c.analytic_exclusion_bound([0], cutoff)
            self.assertEqual(F(cert['lower']), value)
            self.assertEqual(F(cert['upper']), value)
            self.assertTrue(c.verify_analytic_exclusion(cert, [0], cutoff))

    def test_analytic_exclusion_actual_ritz_function(self):
        cert = c.analytic_exclusion_bound([0, 0, 2], 1)
        self.assertEqual(F(cert['lower']), 6)
        self.assertEqual(F(cert['upper']), F(148, 21))
        self.assertEqual(cert['ritz_witness']['harmonic_degree'], 2)
        self.assertTrue(c.verify(cert))

    def test_analytic_exclusion_bound_input_binding(self):
        cert = c.analytic_exclusion_bound([0], 1)
        self.assertFalse(c.verify_analytic_exclusion(cert, expected_cutoff=2))
        self.assertFalse(c.verify_analytic_exclusion(cert, expected_q=[1]))

    def test_analytic_exclusion_witness_tamper(self):
        cert = c.analytic_exclusion_bound([0], 1)
        cert['ritz_witness']['harmonic_degree'] = 1
        self.assertFalse(c.verify(cert))

    def test_analytic_moment_variational_bracket(self):
        cert = c.analytic_moment_bound([0, 0, 2], [[1, 1, 1]], modes=6)
        self.assertEqual(F(cert['lower']), 0)
        self.assertEqual(F(cert['upper']), 4)
        self.assertEqual(cert['variational_trial_dimension'], 2)
        self.assertTrue(c.verify_analytic_moment(cert, [0, 0, 2], [[1, 1, 1]], 0))

    def test_analytic_moment_redundancy_keeps_safe_bound(self):
        cert = c.analytic_moment_bound([0], [[1], [2], [0]], modes=4)
        self.assertEqual(F(cert['upper']), 12)
        self.assertEqual(cert['constraint_count_upper_bound'], 3)
        self.assertTrue(c.verify(cert))

    def test_analytic_moment_dimension_and_scope_tamper(self):
        cert = c.analytic_moment_bound([0], [[1, 0, 1]])
        bad = deepcopy(cert)
        bad['variational_trial_dimension'] = 1
        self.assertFalse(c.verify(bad))
        bad = deepcopy(cert)
        bad['scope'] = 'all_H1_functions_on_unit_S2'
        self.assertFalse(c.verify(bad))
        self.assertFalse(c.verify_analytic_moment(cert, expected_rows=[[1]]))

    def test_general_moments_bind_real_angular_component(self):
        cert = c.certify_moment_sector([0], [[1]], m=1, modes=4, bits=24)
        self.assertEqual(cert['scope'], 'single_real_cosine_azimuth_component_only')
        self.assertIn('cos(m*phi)', cert['constraint_spec']['real_angular_component'])
        bad = deepcopy(cert)
        bad['constraint_spec']['real_angular_component'] = 'both sine and cosine together'
        self.assertFalse(c.verify(bad))

    def test_complete_harmonic_exclusion_q0(self):
        for cutoff in (-1, 0, 1, 2, 5):
            cert = c.full_exclusion_ground([0], cutoff)
            value = (cutoff+1)*(cutoff+2)
            self.assertEqual(F(cert['lower']), value)
            self.assertEqual(F(cert['upper']), value)
            self.assertTrue(c.verify_full(cert, expected_cutoff=cutoff))

    def test_shift_covariance(self):
        cert = c.full_exclusion_ground(['-7/3'], cutoff=2)
        self.assertEqual(F(cert['lower']), F(29, 3))
        self.assertEqual(cert['lower'], cert['upper'])

    def test_sector_degree_start(self):
        self.assertEqual(c.sector_assembly([0], 0, 2, 3)['degrees'], [3, 4, 5])
        self.assertEqual(c.sector_assembly([0], 5, 2, 3)['degrees'], [5, 6, 7])

    def test_constrained_angular_floor(self):
        floors = c.tail_floors([0], cutoff=2, m=0, modes=4, next_m=1)
        self.assertEqual(F(floors['angular_floor']), 12)
        self.assertEqual(F(floors['radial_floor']), 56)

    def test_full_nonzero_potential_and_tail(self):
        self.assertTrue(c.verify_full(self.full, [0, 0, 2], 1))
        self.assertLess(F(self.full['exact_width']), F(1, 10**9))
        self.assertGreater(F(self.full['lower']), 6)
        self.assertLess(F(self.full['upper']), F(13, 2))
        self.assertGreater(len(self.full['sectors']), 1)
        self.assertTrue(self.full['angular_tail_cannot_improve_best_upper'])

    def test_open_angular_gap_is_not_hidden(self):
        cert = c.full_exclusion_ground([0, 0, 2], cutoff=1, max_m=0)
        self.assertFalse(cert['angular_tail_cannot_improve_best_upper'])
        self.assertEqual(cert['status'], 'certified_bound_open_gap')
        self.assertTrue(c.verify_full(cert))

    def test_moment_weighting_known_eigenvalue(self):
        # g=P0+P2: <u,g>=2*a0+(2/5)*a2=0. The mixed eigenvector
        # P2-P0/5 has quotient 5, while P1 has quotient 2.
        cert = c.certify_moment_sector([0], [[1, 0, 1]], k=2, modes=6, bits=40)
        self.assertLessEqual(F(cert['lower']), 5)
        self.assertGreaterEqual(F(cert['upper']), 5)
        self.assertLess(F(cert['exact_width']), F(1, 10**9))
        null = cert['forms']['nullspace']
        self.assertEqual(null['weighted_rows'][0][:3], ['2', '0', '2/5'])
        self.assertEqual(null['T'][0][1], '-1/5')

    def test_constant_alone_violates_moment(self):
        null = c.mass_nullspace([[1, 0, 1]], [F(2), F(2, 3), F(2, 5)])
        self.assertNotEqual(F(null['weighted_rows'][0][0]), 0)
        self.assertEqual(null['rank'], 1)

    def test_non_diagonal_reduced_mass(self):
        forms = c.constrained_forms([0], 0, 5, [[1, 1, 1]])
        self.assertFalse(forms['mass_is_diagonal'])
        self.assertEqual(F(forms['reduced_mass'][0][1]), F(2, 15))
        kernel = c.ConstraintKernel([0], 0, -1, 5, [[1, 1, 1]])
        self.assertEqual(kernel.matrix(F(1), 'upper_ritz')[0][1], -F(2, 15))

    def test_non_diagonal_known_quadratic_root(self):
        # Independent 2x2 hand calculation gives 23*x^2-156*x+180=0;
        # the remaining allowed pure modes start at eigenvalue 12.
        cert = c.certify_moment_sector([0], [[1, 1, 1]], modes=5, bits=40)
        lo, hi = F(cert['lower']), F(cert['upper'])
        polynomial = lambda x: 23*x*x-156*x+180
        self.assertGreaterEqual(polynomial(lo), 0)
        self.assertLessEqual(polynomial(hi), 0)
        self.assertLess(hi, F(78, 23))

    def test_constraint_scaling_invariance(self):
        a = c.certify_moment_sector([0, 0, 2], [[1, 1, 1]], modes=6, bits=32)
        b = c.certify_moment_sector([0, 0, 2], [[-3, -3, -3]], modes=6, bits=32)
        self.assertEqual((a['lower'], a['upper']), (b['lower'], b['upper']))
        self.assertEqual(a['forms']['reduced_mass'], b['forms']['reduced_mass'])

    def test_duplicate_constraint_invariance(self):
        a = c.certify_moment_sector([0], [[1, 0, 1]], modes=6, bits=32)
        b = c.certify_moment_sector([0], [[1, 0, 1], [2, 0, 2], [0]], modes=6, bits=32)
        self.assertEqual((a['lower'], a['upper']), (b['lower'], b['upper']))
        self.assertEqual(b['forms']['nullspace']['rank'], 1)

    def test_two_moments_remove_low_modes(self):
        a = c.certify_moment_sector([0], [[1], [0, 1]], modes=6, bits=40)
        self.assertLessEqual(F(a['lower']), 6)
        self.assertGreaterEqual(F(a['upper']), 6)
        self.assertEqual(a['forms']['reduced_dimension'], 4)

    def test_nonzero_m_has_correct_mass(self):
        cert = c.certify_moment_sector([0], [[1]], m=1, modes=5, bits=36)
        self.assertEqual(cert['forms']['head_mass'][0], '4/3')
        self.assertLessEqual(F(cert['lower']), 6)
        self.assertGreaterEqual(F(cert['upper']), 6)

    def test_general_moments_entire_radial_tail(self):
        self.assertTrue(c.verify_sector(self.moment, expected_rows=[[1, 1, 1]], expected_k=1))
        self.assertEqual(len(self.moment['forms']['couplings']), 2)
        self.assertGreater(F(self.moment['radial_tail_lower']), F(self.moment['upper']))
        self.assertLess(F(self.moment['exact_width']), F(1, 10**8))

    def test_expected_inputs_bind(self):
        self.assertFalse(c.verify_sector(self.moment, expected_rows=[[1, 0, 1]]))
        self.assertFalse(c.verify_sector(self.moment, expected_cutoff=-1))
        self.assertFalse(c.verify_sector(self.moment, expected_m=1))
        self.assertFalse(c.verify_sector(self.moment, expected_k=2))
        self.assertFalse(c.verify_full(self.full, expected_cutoff=0))

    def test_row_tamper(self):
        bad = deepcopy(self.moment)
        bad['constraint_spec']['rows'][0][0] = '2'
        self.assertFalse(c.verify(bad))

    def test_rank_tamper(self):
        bad = deepcopy(self.moment)
        bad['forms']['nullspace']['rank'] = 0
        self.assertFalse(c.verify(bad))

    def test_nullspace_tamper(self):
        bad = deepcopy(self.moment)
        bad['forms']['nullspace']['T'][0][0] = '0'
        self.assertFalse(c.verify(bad))

    def test_mass_tamper(self):
        bad = deepcopy(self.moment)
        bad['forms']['reduced_mass'][0][1] = '0'
        self.assertFalse(c.verify(bad))

    def test_tail_tamper(self):
        bad = deepcopy(self.moment)
        bad['radial_tail_lower'] = '1000000'
        self.assertFalse(c.verify(bad))

    def test_lower_comparison_cannot_certify_upper_endpoint(self):
        bad = deepcopy(self.q0['sectors'][0])
        bad['upper_kind'] = 'lower'
        self.assertFalse(c.verify(bad))
        kernel = c.ConstraintKernel([0], cutoff=1)
        with self.assertRaises(ValueError):
            c._sector_evidence(kernel, 1, F(6), F(6), 'lower')

    def test_index_tamper(self):
        bad = deepcopy(self.moment)
        bad['eigenvalue_index_in_constrained_sector'] = 2
        self.assertFalse(c.verify(bad))

    def test_scope_tamper(self):
        bad = deepcopy(self.moment)
        bad['scope'] = 'all_H1_functions_on_unit_S2'
        self.assertFalse(c.verify(bad))

    def test_cutoff_tamper(self):
        bad = deepcopy(self.full)
        bad['constraint_spec']['cutoff'] = 2
        self.assertFalse(c.verify(bad))

    def test_angular_tail_tamper(self):
        bad = deepcopy(self.full)
        bad['angular_minimum_degree'] += 1
        self.assertFalse(c.verify(bad))

    def test_nonconsecutive_sector_tamper(self):
        bad = deepcopy(self.full)
        bad['sectors'] = bad['sectors'][1:]
        self.assertFalse(c.verify(bad))

    def test_boolean_rejected(self):
        with self.assertRaises(ValueError):
            c.canonical_exclusion(True)
        with self.assertRaises(ValueError):
            c.canonical_moments([[True]])
        bad = deepcopy(self.moment)
        bad['eigenvalue_index_in_constrained_sector'] = True
        self.assertFalse(c.verify(bad))

    def test_floats_rejected(self):
        with self.assertRaises(ValueError):
            c.canonical_moments([[1.0]])
        with self.assertRaises(ValueError):
            c.full_exclusion_ground([0.0])

    def test_support_not_silently_truncated(self):
        with self.assertRaises(ValueError):
            c.canonical_moments([[1, 0, 0, 1]], modes=3)
        with self.assertRaises(ValueError):
            c.mass_nullspace([[1, 0], [0, 1]], [1, 1])

    def test_inequality_three_way_status(self):
        self.assertEqual(c.inequality_certificate(self.q0, 6)['status'], 'proved')
        self.assertEqual(c.inequality_certificate(self.q0, 7)['status'], 'refuted_by_spectral_existence')
        midpoint = (F(self.full['lower'])+F(self.full['upper']))/2
        self.assertEqual(c.inequality_certificate(self.full, midpoint)['status'], 'undetermined')

    def test_single_sector_inequality_scope(self):
        cert = c.inequality_certificate(self.moment, 1)
        self.assertIn('single_azimuth', cert['quantifier'])
        self.assertTrue(c.verify_inequality(cert))
        bad = deepcopy(cert)
        bad['quantifier'] = 'all_H1_functions_on_unit_S2'
        self.assertFalse(c.verify_inequality(bad))

    def test_full_cannot_accept_general_moment_sector(self):
        bad = deepcopy(self.q0)
        bad['sectors'] = [self.moment]
        self.assertFalse(c.verify_full(bad))

    def test_higher_index_cannot_prove_ground_inequality(self):
        cert = c.certify_moment_sector([0], [[1]], k=2, modes=4, bits=24)
        with self.assertRaises(ValueError):
            c.inequality_certificate(cert, 5)

    def test_service_moment_scope_explicit(self):
        with self.assertRaises(ValueError):
            c.constrained_inequality([0], 1, rows=[[1]], modes=4)
        cert = c.constrained_inequality([0], 1, cutoff=-1, rows=[[1]], modes=4, bits=24)
        self.assertTrue(c.verify(cert))

    def test_json_roundtrip(self):
        for cert in (self.full, self.moment, self.q0):
            self.assertTrue(c.verify(json.loads(json.dumps(cert))))


if __name__ == '__main__':
    unittest.main()
