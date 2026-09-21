"""Independent review: weighted forms, polynomial integrals and full tails."""
import copy
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
import unittest

import constraints as c


def legendre(l):
    # Rodrigues' closed coefficient formula, independent of recurrence assembly.
    result = [F(0)]*(l+1)
    for k in range(l//2+1):
        result[l-2*k] = F((-1)**k*factorial(2*l-2*k),
                          2**l*factorial(k)*factorial(l-k)*factorial(l-2*k))
    return result


def multiply(a, b):
    result = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return result


def integral(poly):
    return sum((x*F(2, i+1) for i, x in enumerate(poly) if i % 2 == 0), F(0))


def weighted(q, i, j):
    return integral(multiply(q, multiply(legendre(i), legendre(j))))


class ConstraintsIndependentReview(unittest.TestCase):
    def test_weighted_nullspace_is_not_unweighted(self):
        masses = [F(2), F(2, 3), F(2, 5)]
        null = c.mass_nullspace([[1, 1, 1]], masses)
        t = [[F(x) for x in row] for row in null['T']]
        self.assertEqual(t, [[F(-1, 3), F(-1, 5)], [F(1), F(0)], [F(0), F(1)]])
        for j in range(2):
            self.assertEqual(sum((masses[i]*t[i][j] for i in range(3)), F(0)), 0)
        self.assertNotEqual(sum(row[0] for row in t), 0)

    def test_full_non_diagonal_mass(self):
        forms = c.constrained_forms([0], 0, 3, rows=[[1, 1, 1]])
        self.assertEqual([[F(x) for x in row] for row in forms['reduced_mass']],
                         [[F(8, 9), F(2, 15)], [F(2, 15), F(12, 25)]])
        self.assertEqual([[F(x) for x in row] for row in forms['reduced_A']],
                         [[F(4, 3), F(0)], [F(0), F(12, 5)]])
        self.assertFalse(forms['mass_is_diagonal'])

    def test_generalized_shift_uses_off_diagonal_mass(self):
        kernel = c.ConstraintKernel([0], 0, -1, 3, rows=[[1, 1, 1]])
        matrix = kernel.matrix(F(3, 2), 'upper_ritz')
        self.assertEqual(matrix[0][1], -F(1, 5))
        self.assertEqual(matrix[1][0], -F(1, 5))
        # The independent determinant is proportional to 23*x^2-156*x+180.
        x = F(3, 2)
        determinant = matrix[0][0]*matrix[1][1]-matrix[0][1]**2
        self.assertEqual(determinant, F(4, 225)*(23*x*x-156*x+180))

    def test_reduced_potential_forms_by_rodrigues_integrals(self):
        q, modes = [F(1, 7), F(-2, 3), F(3, 5), F(1, 9)], 5
        forms = c.constrained_forms(q, 0, modes, rows=[[1, 1, 1], [0, 1, 0, 2]])
        t = [[F(x) for x in row] for row in forms['nullspace']['T']]
        dimension = len(t[0])
        for a in range(dimension):
            for b in range(dimension):
                expected = F(0)
                for i in range(modes):
                    for j in range(modes):
                        entry = weighted(q, i, j)
                        if i == j:
                            entry += i*(i+1)*weighted([F(1)], i, i)
                        expected += t[i][a]*entry*t[j][b]
                self.assertEqual(F(forms['reduced_A'][a][b]), expected)

    def test_transformed_tail_couplings_by_independent_integrals(self):
        q, modes = [F(0), F(1, 3), F(0), F(1)], 5
        forms = c.constrained_forms(q, 0, modes, rows=[[1, 1, 1]])
        t = [[F(x) for x in row] for row in forms['nullspace']['T']]
        self.assertEqual([entry['degree'] for entry in forms['couplings']], [5, 6, 7])
        for entry in forms['couplings']:
            l = entry['degree']
            self.assertEqual(F(entry['mass']), weighted([F(1)], l, l))
            for a, actual in enumerate(entry['column']):
                expected = sum((t[i][a]*weighted(q, i, l) for i in range(modes)), F(0))
                self.assertEqual(F(actual), expected)

    def test_complete_multiplication_before_exclusion(self):
        data = c.sector_assembly([0, 0, 1], m=0, cutoff=1, modes=3)
        self.assertEqual(data['degrees'], [2, 3, 4])
        self.assertEqual(data['A'][0][0]/data['M'][0]-6, F(11, 21))
        self.assertEqual(data['A'][0][0]-6*data['M'][0], weighted([0, 0, 1], 2, 2))

    def test_general_moment_exact_algebraic_ground(self):
        cert = c.certify_moment_sector([0], [[1, 1, 1]], modes=3, bits=48)
        lo, hi = F(cert['lower']), F(cert['upper'])
        polynomial = lambda x: 23*x*x-156*x+180
        self.assertTrue(1 < lo < hi < 2)
        self.assertGreaterEqual(polynomial(lo), 0)
        self.assertLessEqual(polynomial(hi), 0)
        self.assertTrue(c.verify_sector(cert, expected_rows=[[1, 1, 1]], expected_k=1))
        # A different moment leaves P1 at 2 and shifts the P0/P2 combination to 5.
        second = c.certify_moment_sector([0], [[1, 0, 1]], k=2, modes=3, bits=48)
        self.assertLessEqual(F(second['lower']), 5)
        self.assertGreaterEqual(F(second['upper']), 5)

    def test_redundant_constraints_preserve_reduced_forms(self):
        first = c.constrained_forms([0, 1, 2], 1, 5, rows=[[1, 2, 3]])
        duplicate = c.constrained_forms([0, 1, 2], 1, 5,
                                        rows=[[2, 4, 6], [-1, -2, -3], [0]])
        self.assertEqual(first['reduced_A'], duplicate['reduced_A'])
        self.assertEqual(first['reduced_mass'], duplicate['reduced_mass'])
        self.assertEqual(first['couplings'], duplicate['couplings'])

    def test_constraint_support_not_silently_truncated(self):
        with self.assertRaises(ValueError):
            c.certify_moment_sector([0], [[0, 0, 0, 1]], modes=3)
        with self.assertRaises(ValueError):
            c.constrained_forms([0], 0, 3, rows=[[1]], cutoff=1)
        with self.assertRaises(ValueError):
            c.certify_moment_sector([0], [[1], [0, 1]], modes=2)

    def test_complete_cutoff_raises_all_angular_floors(self):
        for cutoff in (0, 1, 2, 8):
            cert = c.full_exclusion_ground([0], cutoff=cutoff, modes=1, max_m=0)
            expected = (cutoff+1)*(cutoff+2)
            self.assertEqual(F(cert['lower']), expected)
            self.assertEqual(F(cert['upper']), expected)
            self.assertEqual(F(cert['angular_tail_lower']), expected)
            self.assertTrue(c.verify_full(cert, expected_cutoff=cutoff))

    def test_incomplete_angular_budget_remains_open(self):
        cert = c.full_exclusion_ground([0, 0, 2], cutoff=1, modes=5, max_m=0)
        self.assertEqual(F(cert['lower']), 6)
        self.assertGreater(F(cert['upper']), 7)
        self.assertEqual(cert['status'], 'certified_bound_open_gap')
        self.assertTrue(c.verify_full(cert))

    def test_single_sector_scope_and_tamper_rejection(self):
        cert = c.certify_moment_sector([0, 1], [[1, 1, 1]], modes=4, bits=32)
        inequality = c.inequality_certificate(cert, 0)
        self.assertIn('single_azimuth_cosine_component', inequality['quantifier'])
        self.assertEqual(cert['scope'], 'single_real_cosine_azimuth_component_only')
        self.assertFalse(c.verify_full(cert))
        for key in ('reduced_mass', 'couplings', 'nullspace'):
            bad = copy.deepcopy(cert)
            if key == 'reduced_mass':
                bad['forms'][key][0][1] = '0'
            elif key == 'couplings':
                bad['forms'][key][0]['column'][0] = '100'
            else:
                bad['forms'][key]['T'][0][0] = '0'
            self.assertFalse(c.verify_sector(bad))
        inequality['quantifier'] = 'all_H1_functions_on_unit_S2'
        self.assertFalse(c.verify_inequality(inequality))

    def test_replay_all_saved_round_six_certificates(self):
        paths = sorted((Path(__file__).parent/'round06'/'certificates').glob('*.json'))
        self.assertGreaterEqual(len(paths), 17)
        for path in paths:
            with self.subTest(path=path.name):
                self.assertTrue(c.verify(json.loads(path.read_text())))

    def test_analytic_exclusion_ritz_by_explicit_integral(self):
        cert = c.analytic_exclusion_bound([0, 0, 1], cutoff=1)
        self.assertEqual(F(cert['lower']), 6)
        self.assertEqual(F(cert['upper']), F(137, 21))
        expected = 6+weighted([0, 0, 1], 2, 2)/weighted([1], 2, 2)
        self.assertEqual(F(cert['upper']), expected)
        self.assertTrue(c.verify_analytic_exclusion(cert, [0, 0, 1], 1))
        self.assertFalse(c.verify_analytic_exclusion(cert, expected_cutoff=0))

    def test_analytic_moment_dimension_argument_can_use_tail(self):
        # Four redundant conditions supported at l0; head alone is exhausted,
        # but the infinite radial space survives and the dimension bound is valid.
        cert = c.analytic_moment_bound([0], [[1], [2], [3], [4]], modes=1)
        self.assertEqual(F(cert['lower']), 0)
        self.assertEqual(F(cert['upper']), 20)
        self.assertEqual(cert['variational_trial_dimension'], 5)
        self.assertEqual(cert['trial_highest_harmonic_degree'], 4)
        self.assertTrue(c.verify_analytic_moment(cert, expected_rows=[[1], [2], [3], [4]]))
        self.assertFalse(c.verify_analytic_moment(cert, expected_rows=[[0], [2], [3], [4]]))
        cert['scope'] = 'all_H1_functions_on_unit_S2'
        self.assertFalse(c.verify(cert))


if __name__ == '__main__':
    unittest.main()
