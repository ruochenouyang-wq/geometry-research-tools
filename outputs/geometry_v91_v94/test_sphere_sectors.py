"""Independent low-degree weighted integrals and frozen m=0 cross-checks."""
from fractions import Fraction as F
import unittest

import sphere_sectors as sectors


# Explicit ordinary Legendre polynomials, independently differentiated below.
LEGENDRE = [
    [F(1)], [F(0), F(1)], [F(-1, 2), F(0), F(3, 2)],
    [F(0), F(-3, 2), F(0), F(5, 2)],
    [F(3, 8), F(0), F(-15, 4), F(0), F(35, 8)],
    [F(0), F(15, 8), F(0), F(-35, 4), F(0), F(63, 8)],
]


def multiply(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def derivative(p, order):
    for _ in range(order):
        p = [(i + 1) * p[i + 1] for i in range(len(p) - 1)]
    return p


def weighted_integral(q, m, i, j):
    # P_i^m P_j^m = (1-t^2)^m (d^m P_i)(d^m P_j): the signs cancel.
    p = multiply(multiply(q, derivative(LEGENDRE[i], m)), derivative(LEGENDRE[j], m))
    for _ in range(m):
        p = multiply(p, [F(1), F(0), F(-1)])
    return sum((coefficient * F(2, degree + 1)
                for degree, coefficient in enumerate(p) if degree % 2 == 0), F(0))


class SectorAssemblyTests(unittest.TestCase):
    def test_m_zero_matches_frozen_forms(self):
        for q in ([0], ["2/3"], [0, 1], [1, "-3/7", 2, 0, 0, "1/13", "2/5"]):
            for modes in (1, 2, 8):
                legacy = sectors._legacy.assemble(sectors.potential(q), modes)
                data = sectors.assemble(q, 0, 0, modes)
                self.assertEqual(data["A"], legacy["A"])
                self.assertEqual(data["M"], legacy["M"])
                self.assertEqual(data["degrees"], list(range(modes)))

    def test_full_m_zero_multiplication_matches_frozen(self):
        q = sectors.potential([1, "2/7", 3, "-2/3", 0, 2, 1])
        for l in (0, 1, 5, 17):
            self.assertEqual(sectors.q_times_basis(q, 0, l), sectors._legacy.q_times_p(q, l))

    def test_explicit_m_one_two_integrals(self):
        for m in (1, 2):
            for q in ([F(0)], [F(7, 5)], [F(1, 3), F(-2), F(3, 7)],
                      [F(0), F(0), F(0), F(0), F(0), F(0), F(2, 9)]):
                data = sectors.assemble(q, m, m, 6 - m)
                for a, i in enumerate(data["degrees"]):
                    self.assertEqual(data["M"][a], weighted_integral([F(1)], m, i, i))
                    for b, j in enumerate(data["degrees"]):
                        expected = weighted_integral(q, m, i, j)
                        if i == j:
                            expected += i * (i + 1) * data["M"][a]
                        self.assertEqual(data["A"][a][b], expected)

    def test_lowest_degree_boundary(self):
        for m in (0, 1, 2, 64):
            self.assertEqual(sectors.q_times_basis([0, 1], m, m), {m + 1: F(1, 2 * m + 1)})
            self.assertNotIn(m - 1, sectors.q_times_basis([0, 0, 0, 1], m, m))

    def test_zero_mean_projects_after_whole_polynomial(self):
        q = [F(0), F(0), F(1)]
        data = sectors.assemble(q, 0, 1, 5)
        full = sectors._legacy.assemble(q, 6)
        self.assertEqual(data["A"], [row[1:] for row in full["A"][1:]])
        self.assertEqual(data["M"], full["M"][1:])
        # t^2 P_1 has P_1 coefficient 3/5; the path through P_0 is essential.
        self.assertEqual(data["A"][0][0] / data["M"][0] - 2, F(3, 5))
        self.assertEqual(sectors.q_times_basis([0, 1], 0, 1)[0], F(1, 3))

    def test_arbitrary_window_is_full_form_submatrix(self):
        for m in (0, 1, 2, 5):
            q = [1, "-1/3", 0, 2, 0, 0, "3/5"]
            full = sectors.assemble(q, m, m, 11)
            data = sectors.assemble(q, m, m + 3, 5)
            self.assertEqual(data["A"], [row[3:8] for row in full["A"][3:8]])
            self.assertEqual(data["M"], full["M"][3:8])

    def test_constant_potential_and_exact_types(self):
        for m, start, modes in ((0, 1, 4), (2, 2, 4), (64, 128, 64)):
            data = sectors.assemble(["-7/11"], m, start, modes)
            for i, l in enumerate(data["degrees"]):
                self.assertIsInstance(data["M"][i], F)
                for j in range(modes):
                    self.assertIsInstance(data["A"][i][j], F)
                    self.assertEqual(data["A"][i][j],
                                     (l * (l + 1) - F(7, 11)) * data["M"][i] if i == j else 0)

    def test_symmetry_bandwidth_and_parity(self):
        for m in (0, 1, 2, 17, 64):
            data = sectors.assemble([1, 0, "-3/2", 0, 0, 0, "2/5"], m, m, 12)
            for i in range(12):
                for j in range(12):
                    self.assertEqual(data["A"][i][j], data["A"][j][i])
                    if abs(i - j) > 6 or (i - j) % 2:
                        self.assertEqual(data["A"][i][j], 0)

    def test_helper_limit_covers_tail_degrees(self):
        series = sectors.q_times_basis([0, 0, 0, 0, 0, 0, 1], 64, 250)
        self.assertEqual(max(series), 256)
        self.assertTrue(all(sectors.basis_mass(64, l) > 0 for l in series))

    def test_invalid_parameters(self):
        for m, start, modes in ((-1, 0, 1), (65, 65, 1), (True, 1, 1),
                                (2, 1, 1), (0, 129, 1), (0, 0.0, 1),
                                (0, 0, 0), (0, 0, 65), (0, 0, True)):
            with self.subTest(m=m, start=start, modes=modes), self.assertRaises(ValueError):
                sectors.assemble([0], m, start, modes)
        for q in ([], [0] * 8, [1.0], [True], ["nan"], [1001], ["1/1000000001"]):
            with self.subTest(q=q), self.assertRaises(ValueError):
                sectors.assemble(q, 0, 0, 1)
        with self.assertRaises(ValueError):
            sectors.basis_mass(1, 0)
        with self.assertRaises(ValueError):
            sectors.q_times_basis([0], 0, 251)


if __name__ == "__main__":
    unittest.main()
