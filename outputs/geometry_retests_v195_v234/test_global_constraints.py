from copy import deepcopy
from fractions import Fraction as F
import json
import unittest
import global_constraints as g

C = {'0,0,0': 1}
X = {'1,0,0': 1}
Y = {'0,1,0': 1}
Z = {'0,0,1': 1}


class GlobalConstraintTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c12 = g.certify({}, [C, Z], L=2)
        cls.cross = g.certify({'1,0,0': 1, '0,0,1': 1}, [C, Z], L=2, bits=28)

    def test_original_c12_exact_full_ground(self):
        self.assertEqual((self.c12['lower'], self.c12['upper']), ('2', '2'))
        self.assertEqual(self.c12['scope'], g.SCOPE)

    def test_original_x_and_y_are_admissible(self):
        for trial in (X, Y):
            cert = g.trial_enclosure({}, [C, Z], trial, support_L=1)
            self.assertEqual(cert['upper'], '2')
            self.assertEqual(cert['moment_inner_products'], ['0', '0'])
            self.assertTrue(g.verify(cert))

    def test_z_is_not_admissible(self):
        with self.assertRaises(ValueError):
            g.trial_enclosure({}, [C, Z], Z)

    def test_all_linear_moments_raise_ground(self):
        cert = g.certify({}, [C, X, Y, Z], L=2)
        self.assertEqual((cert['lower'], cert['upper']), ('6', '6'))

    def test_cross_m_linear_constraint_retains_ground_two(self):
        cert = g.certify({}, [C, {'1,0,0': 1, '0,0,1': 1}], L=2)
        self.assertEqual((cert['lower'], cert['upper']), ('2', '2'))

    def test_empty_constraints_include_constant(self):
        cert = g.certify({}, [], L=0)
        self.assertEqual((cert['lower'], cert['upper']), ('0', '0'))

    def test_constant_shift(self):
        cert = g.certify({'0,0,0': '-7/3'}, [C, Z], L=2)
        self.assertEqual((cert['lower'], cert['upper']), ('-1/3', '-1/3'))

    def test_all_components_basis(self):
        data = g.full_harmonics(2)
        self.assertEqual(data['dimension'], 9)
        self.assertEqual(data['basis'][0]['polynomial'], {'0,0,0': '1'})
        self.assertTrue(any(b['degree'] == 1 and b['m'] == 1 and b['part'] == 'imag' for b in data['basis']))
        self.assertEqual(data['gram'][0][0], '1')

    def test_constant_term_is_not_constant_rowspace(self):
        rows = [{'0,0,0': 1, '1,0,0': 1}]
        cert = g.trial_enclosure({}, rows, {'1,0,0': 1, '0,0,0': '-1/3'}, support_L=1)
        self.assertFalse(cert['constant_in_row_span'])
        self.assertEqual(cert['lower'], '0')
        self.assertEqual(cert['upper'], '3/2')
        spectral = g.certify({}, rows, L=2, bits=36)
        self.assertLessEqual(F(spectral['lower']), F(3, 2))
        self.assertGreaterEqual(F(spectral['upper']), F(3, 2))

    def test_combined_rows_imply_mean_zero(self):
        rows = [{'0,0,0': 1, '1,0,0': 1}, {'0,0,0': 1, '1,0,0': -1}]
        system = g.constraint_nullspace(rows, 1)
        self.assertTrue(system['constant_in_row_span'])
        self.assertEqual(system['constant_span_coefficients'], ['1/2', '1/2'])
        cert = g.certify({}, rows, L=2)
        self.assertEqual(cert['lower'], '2')

    def test_three_component_mixed_known_value(self):
        rows = [{'0,0,0': 1, '1,0,0': 1, '0,1,0': 2, '0,0,1': 3}]
        cert = g.certify({}, rows, L=2, bits=36)
        self.assertLessEqual(F(cert['lower']), F(6, 17))
        self.assertGreaterEqual(F(cert['upper']), F(6, 17))
        self.assertFalse(cert['forms']['constraint_system']['constant_in_row_span'])

    def test_dense_mass_is_retained(self):
        data = g.reduced_forms({}, [{'0,0,0': 1, '1,0,0': 1, '0,0,1': 1}], L=1)
        self.assertFalse(data['mass_is_diagonal'])
        self.assertTrue(any(F(row[j]) for i, row in enumerate(data['reduced_mass']) for j in range(i)))

    def test_tail_contains_all_real_modes(self):
        tail = g.transformed_tail({'1,0,0': 1}, [C, Z], L=2)
        self.assertEqual(len(tail['columns']), 7)
        self.assertEqual({x['degree'] for x in tail['columns']}, {3})
        self.assertTrue(any(x['part'] == 'imag' for x in tail['columns']))
        self.assertTrue(any(F(v) for x in tail['columns'] for v in x['entries']))

    def test_nonconstant_full_sphere_certificate(self):
        self.assertTrue(g.verify(self.cross, expected_q={'1,0,0': 1, '0,0,1': 1}, expected_constraints=[C, Z]))
        self.assertLess(F(self.cross['upper']), 2)
        self.assertGreater(F(self.cross['lower']), 1)

    def test_snap_requires_exact_inertia(self):
        self.assertIsNone(g.exact_level({}, [{'0,0,0': 1, '1,0,0': 1}], L=2))
        self.assertFalse(self.cross['exact_level_verified_by_inertia'])

    def test_equivalent_scaled_duplicate_rows(self):
        target = [{'0,0,0': -3}, {'0,0,1': 2}, {'0,0,1': 4}, {}]
        transfer = g.equivalent_constraints(self.c12, target)
        self.assertEqual(transfer['lower'], '2')
        self.assertTrue(g.verify(transfer, expected_constraints=target))

    def test_nonequivalent_transfer_rejected(self):
        with self.assertRaises(ValueError):
            g.equivalent_constraints(self.c12, [C, X, Y, Z])

    def test_support_cannot_be_truncated(self):
        with self.assertRaises(ValueError):
            g.constraint_nullspace([{'0,0,2': 1}], L=1)

    def test_sphere_identity_constraint_is_redundant(self):
        identity = {'2,0,0': 1, '0,2,0': 1, '0,0,2': 1, '0,0,0': -1}
        system = g.constraint_nullspace([identity], L=0)
        self.assertEqual(system['rank'], 0)

    def test_adaptive_keeps_all_constraints(self):
        result = g.adaptive({}, [C, X, Y, Z], start_L=1, max_L=2)
        self.assertEqual(result['attempts'][0]['status'], 'constraints_eliminate_head')
        self.assertEqual(result['final_certificate']['lower'], '6')
        self.assertTrue(g.verify(result))

    def test_budget_failure_remains_visible(self):
        result = g.adaptive({}, [C, X, Y, Z], start_L=1, max_L=1)
        self.assertEqual(result['status'], 'unsupported_within_budget')
        self.assertIsNone(result['final_certificate'])
        self.assertTrue(g.verify(result))

    def test_gap_open_remains_visible(self):
        result = g.adaptive({'1,0,0': 1, '0,0,1': 1}, [C, Z], tolerance=F(1, 10**12), start_L=1, max_L=1, bits=16)
        self.assertEqual(result['status'], 'certified_bound_open_gap')
        self.assertTrue(g.verify(result))

    def test_inequality_proves_and_explicitly_refutes(self):
        proved = g.inequality({}, [C, Z], 2, L=2)
        refuted = g.inequality({}, [C, Z], 3, L=2)
        self.assertEqual(proved['status'], 'proved')
        self.assertEqual(refuted['status'], 'refuted_with_trial')
        self.assertLess(F(refuted['trial_certificate']['upper']), 3)
        self.assertTrue(g.verify(refuted))

    def test_undetermined_status(self):
        midpoint = (F(self.cross['lower'])+F(self.cross['upper']))/2
        result = g.inequality_record(self.cross, midpoint)
        self.assertEqual(result['status'], 'undetermined')
        self.assertTrue(g.verify(result))

    def test_exact_negative_direction(self):
        for matrix in ([[F(0), F(1)], [F(1), F(0)]], [[F(2), F(3)], [F(3), F(2)]]):
            v = g._negative_direction(matrix)
            value = sum(v[i]*matrix[i][j]*v[j] for i in range(2) for j in range(2))
            self.assertLess(value, 0)

    def test_explicit_refutation_can_lift_schur_tail(self):
        q = {'1,0,0': 1, '0,0,1': 1}
        threshold = F(39, 20)
        kernel = g.Kernel(q, [C, Z], L=1)
        self.assertIsNone(g._negative_direction(kernel.matrix(threshold, 'upper_ritz')))
        trial = g.violating_trial(q, [C, Z], threshold, L=1)
        self.assertIsNotNone(trial)
        self.assertTrue(any(sum(map(int, e.split(','))) == 2 for e in trial['trial']))
        self.assertLess(F(trial['upper']), threshold)
        self.assertTrue(g.verify(trial))

    def test_wrong_scope_rejected(self):
        bad = deepcopy(self.c12); bad['scope'] = 'single_azimuth_sector_only'
        self.assertFalse(g.verify(bad))

    def test_wrong_constant_claim_rejected(self):
        cert = g.trial_enclosure({}, [], C, support_L=1)
        cert['constant_in_row_span'] = True
        self.assertFalse(g.verify(cert))

    def test_tail_tamper_rejected(self):
        bad = deepcopy(self.cross); bad['tail']['columns'].pop()
        self.assertFalse(g.verify(bad))

    def test_tail_floor_tamper_rejected(self):
        bad = deepcopy(self.cross); bad['infinite_tail_lower'] = '100000'
        self.assertFalse(g.verify(bad))

    def test_constraint_tamper_rejected(self):
        bad = deepcopy(self.c12); bad['constraints'][1] = X
        self.assertFalse(g.verify(bad))

    def test_mass_tamper_rejected(self):
        bad = deepcopy(self.cross); bad['forms']['reduced_mass'][0][0] = '10'
        self.assertFalse(g.verify(bad))

    def test_lower_comparison_cannot_prove_upper(self):
        bad = deepcopy(self.c12); bad['upper_kind'] = 'lower'
        self.assertFalse(g.verify(bad))

    def test_false_exact_level_rejected(self):
        bad = deepcopy(self.cross); bad['exact_level_verified_by_inertia'] = True
        self.assertFalse(g.verify(bad))

    def test_adaptive_skipped_support_failure_rejected(self):
        bad = g.adaptive({}, [C, X, Y, Z], start_L=1, max_L=2)
        bad['attempts'].pop(0)
        self.assertFalse(g.verify(bad))

    def test_input_boolean_and_float_rejected(self):
        for raw in (True, 1.0):
            with self.assertRaises(ValueError):
                g.certify({'0,0,0': raw}, [C])
            with self.assertRaises(ValueError):
                g.constraint_nullspace([C], L=raw)

    def test_json_replay(self):
        self.assertTrue(g.verify(json.loads(json.dumps(self.cross))))


if __name__ == '__main__':
    unittest.main()
