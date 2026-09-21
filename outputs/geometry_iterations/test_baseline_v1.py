from fractions import Fraction as F
from math import comb, factorial
import random
import unittest

import spectral_certifier as s


class ExactAlgebraTests(unittest.TestCase):
    def test_legendre_assembly_against_rodrigues_and_integration(self):
        def p(n):
            out = [F(0)] * (n + 1)
            for j in range(n + 1):
                if 2 * j >= n:
                    out[2 * j - n] = F(comb(n, j) * (-1)**(n - j) * factorial(2 * j),
                                        factorial(2 * j - n) * 2**n * factorial(n))
            return out
        def product(a, b):
            out = [F(0)] * (len(a) + len(b) - 1)
            for i, x in enumerate(a):
                for j, y in enumerate(b):
                    out[i + j] += x * y
            return out
        q = [F(3, 7), F(-2, 5), F(7, 4), F(1, 3)]
        data = s.assemble(q, 5)
        for i in range(5):
            for j in range(5):
                integrand = product(product(p(i), p(j)), q)
                integral = sum((F(2) * c / (k + 1) for k, c in enumerate(integrand) if k % 2 == 0), F(0))
                expected = integral + (F(2 * i * (i + 1), 2 * i + 1) if i == j else 0)
                self.assertEqual(data["A"][i][j], expected)

    def test_two_by_two_zero_diagonal_pivot(self):
        self.assertEqual(s.inertia([[0, 3], [3, 0]]), [1, 0, 1])

    def test_zero_and_indefinite_matrices(self):
        self.assertEqual(s.inertia([[0, 0], [0, 0]]), [0, 2, 0])
        self.assertEqual(s.inertia([[1, 2], [2, 1]]), [1, 0, 1])
        self.assertEqual(s.inertia([[0, 0, 0], [0, -F(1, 10**90), 0], [0, 0, 1]]), [1, 1, 1])

    def test_inertia_preserved_under_exact_congruences(self):
        rng = random.Random(192)
        for _ in range(20):
            n = 5
            diagonal = [rng.choice([-3, -1, 0, 2, 4]) for _ in range(n)]
            t = [[F(1) if i == j else F(rng.randint(-3, 3), rng.randint(1, 5)) if i < j else F(0)
                  for j in range(n)] for i in range(n)]
            a = [[sum((t[k][i] * diagonal[k] * t[k][j] for k in range(n)), F(0)) for j in range(n)] for i in range(n)]
            self.assertEqual(s.inertia(a), [sum(d < 0 for d in diagonal), sum(d == 0 for d in diagonal), sum(d > 0 for d in diagonal)])

    def test_linear_potential_tail_correction_by_hand(self):
        data = s.assemble([F(0), F(1)], 2)
        self.assertEqual(data["beta"], 5)
        self.assertEqual(data["G"], [[F(0), F(0)], [F(0), F(8, 45)]])

    def test_decimal_rounding_is_outward_for_both_signs(self):
        for value in [F(1, 3), F(-1, 3), F(2), F(-1, 10**15)]:
            lower = F(s.outward_decimal(value))
            upper = F(s.outward_decimal(value, upper=True))
            self.assertLessEqual(lower, value)
            self.assertGreaterEqual(upper, value)


class SpectralProofTests(unittest.TestCase):
    def test_known_sphere_eigenvalues(self):
        for k in (1, 2, 4):
            certificate = s.certify(["0"], k, 6, 24)
            exact = (k - 1) * k
            self.assertLessEqual(F(certificate["lower"]), exact)
            self.assertGreaterEqual(F(certificate["upper"]), exact)
            self.assertTrue(s.verify(certificate))

    def test_constant_potential_shifts_spectrum(self):
        for c in (F(-7, 3), F(11, 5)):
            certificate = s.certify([str(c)], 3, 6, 24)
            self.assertLessEqual(F(certificate["lower"]), 6 + c)
            self.assertGreaterEqual(F(certificate["upper"]), 6 + c)

    def test_omitted_tail_prevents_false_lower_bound(self):
        data = s.assemble([F(0), F(1)], 2)
        false_lower = F(-155, 1000)
        self.assertEqual(s.inertia(s.shifted(data, false_lower))[0], 0)
        self.assertFalse(s.lower_holds(data, 1, false_lower))
        fine = s.certify(["0", "1"], 1, 8, 28)
        self.assertLess(F(fine["upper"]), false_lower)

    def test_more_modes_improves_certified_gap(self):
        coarse = s.certify(["0", "1"], 1, 2, 28)
        fine = s.certify(["0", "1"], 1, 8, 28)
        self.assertLess(F(fine["exact_width"]), F(coarse["exact_width"]) / 1000)
        self.assertTrue(s.verify(fine))

    def test_reflection_invariance(self):
        a = s.certify(["1", "2", "3"], 1, 8, 20)
        b = s.certify(["1", "-2", "3"], 1, 8, 20)
        self.assertEqual((a["lower"], a["upper"]), (b["lower"], b["upper"]))

    def test_strong_coupling_establishes_lower_bracket(self):
        data = s.assemble([F(0), F(1000)], 2)
        self.assertFalse(s.lower_holds(data, 1, data["qmin"] - 1))
        cert = s.certify(["0", "1000"], 1, 2, 24)
        self.assertTrue(s.lower_holds(data, 1, F(cert["lower"])))
        self.assertTrue(s.verify(cert))
        self.assertLessEqual(F(cert["lower"]), F(cert["upper"]))

    def test_wrong_eigenvalue_index_rejected(self):
        cert = s.certify(["0"], 2, 6, 20)
        cert["eigenvalue_index"] = 1
        self.assertFalse(s.verify(cert))

    def test_tampered_tail_bound_rejected(self):
        cert = s.certify(["0", "1"], 1, 4, 20)
        cert["evidence"]["tail_lower"] = "100000"
        self.assertFalse(s.verify(cert))

    def test_wrong_operator_rejected(self):
        cert = s.certify(["0", "1"], 1, 4, 20)
        cert["q_coefficients"] = ["10"]
        self.assertFalse(s.verify(cert))

    def test_false_lower_and_upper_rejected(self):
        for key, value in [("lower", "1"), ("upper", "-1")]:
            cert = s.certify(["0", "1"], 1, 4, 20)
            cert[key] = value
            self.assertFalse(s.verify(cert))

    def test_scope_cannot_be_promoted_to_arbitrary_geometry(self):
        cert = s.certify(["0"], 1, 4, 20)
        cert["scope"] = "all minimal hypersurfaces"
        self.assertFalse(s.verify(cert))

    def test_presented_decimal_cannot_disagree_with_certificate(self):
        cert = s.certify(["0", "1"], 1, 4, 20)
        cert["decimal_enclosure"] = ["0", "0"]
        self.assertFalse(s.verify(cert))

    def test_invalid_sizes_and_inputs(self):
        for args in [(0, 8), (2, 1), (1, 100), (True, 2)]:
            with self.assertRaises(ValueError):
                s.check_sizes(*args)
        for raw in [[0.1], [True], ["1/0"], [], ["1"] * 8]:
            with self.assertRaises(ValueError):
                s.potential(raw)


if __name__ == "__main__":
    unittest.main()
