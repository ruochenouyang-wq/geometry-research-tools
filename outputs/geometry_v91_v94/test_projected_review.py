"""Independent mathematical regression checks for the projected spectrum adapter.

Form matrices below are assembled by polynomial differentiation/integration,
not by the adapter's associated-Legendre multiplication recurrence.
"""
from copy import deepcopy
from fractions import Fraction as F
from math import factorial
import unittest

import projected_spectrum as p


def add(a, b):
    c = [F(0)] * max(len(a), len(b))
    for i, v in enumerate(a):
        c[i] += v
    for i, v in enumerate(b):
        c[i] += v
    return c


def scale(a, c):
    return [F(c) * v for v in a]


def mul(a, b):
    c = [F(0)] * (len(a) + len(b) - 1)
    for i, u in enumerate(a):
        for j, v in enumerate(b):
            c[i+j] += u*v
    return c


def derivative(a):
    return [i*a[i] for i in range(1, len(a))] or [F(0)]


def integral(a):
    return sum((2*v/F(i+1) for i, v in enumerate(a) if i % 2 == 0), F(0))


def legendre(l):
    # Rodrigues' formula, independent of the multiplication recurrence.
    a = [F(1)]
    for _ in range(l):
        a = mul(a, [-1, 0, 1])
    for _ in range(l):
        a = derivative(a)
    return scale(a, F(1, 2**l * factorial(l)))


def direct_forms(q, m, start, n):
    last = start+n+len(q)-2
    basis = {}
    for l in range(start, last+1):
        a = legendre(l)
        for _ in range(m):
            a = derivative(a)
        basis[l] = a
    w = [F(1)]
    for _ in range(m):
        w = mul(w, [1, 0, -1])
    degrees = list(range(start, start+n))
    def inner(i, j, potential):
        return integral(mul(mul(w, mul(basis[i], basis[j])), potential))
    mass = [inner(l, l, [1]) for l in degrees]
    a = [[inner(i, j, q) + (i*(i+1)*mass[ii] if i == j else 0)
          for j in degrees] for ii, i in enumerate(degrees)]
    tails = [(l, inner(l, l, [1]), [inner(i, l, q) for i in degrees])
             for l in range(start+n, last+1)]
    return a, mass, tails


class ProjectedIndependentReview(unittest.TestCase):
    def test_direct_rodrigues_forms_and_tail_mass(self):
        for m, start, q in [(0, 1, [0, 1, 1]), (1, 1, [1, -2, 3]),
                            (2, 2, [1, -2, 0, 3, 0, -1, 1])]:
            with self.subTest(m=m):
                kernel = p.Kernel(q, m=m, mean_zero=True, modes=4)
                a, mass, tails = direct_forms(q, m, start, 4)
                self.assertEqual(kernel.a, a)
                self.assertEqual(kernel.mass, mass)
                expected = [(l, h, [(i, b) for i, b in enumerate(column) if b])
                            for l, h, column in tails]
                self.assertEqual([(l, h, sorted(column)) for l, h, column in kernel.couplings], expected)
                x = kernel.beta-1
                for kind, bound in [('lower', kernel.qlo), ('upper_tail', kernel.qhi)]:
                    expected_matrix = [row[:] for row in a]
                    for i in range(4):
                        expected_matrix[i][i] -= x*mass[i]
                    for l, h, column in tails:
                        for i in range(4):
                            for j in range(4):
                                expected_matrix[i][j] -= column[i]*column[j]/(h*(l*(l+1)+bound-x))
                    self.assertEqual(kernel.matrix(x, kind), expected_matrix)

    def test_polynomial_applied_before_projection(self):
        kernel = p.Kernel([0, 0, 1], modes=2)
        # Integral t^4 = 2/5 includes the l=1 -> 0 -> 1 path.
        self.assertEqual(kernel.a[0][0]-2*kernel.mass[0], F(2, 5))

    def test_independently_calculated_coarse_endpoints(self):
        cases = [
            ([0, 1], 0, False, 2, F(548077, 262144), F(68511, 32768)),
            ([0, 1], 0, True, 1, F(506981, 262144), F(31687, 16384)),
            ([0, 1], 1, True, 1, F(127815, 65536), F(511271, 262144)),
            ([0, 1, 1], 0, False, 2, F(347469, 131072), F(694949, 262144)),
            ([0, 1, 1], 0, True, 1, F(663027, 262144), F(331519, 131072)),
            ([0, 1, 1], 1, True, 1, F(70469, 32768), F(563763, 262144)),
        ]
        for q, m, projected, k, lo, hi in cases:
            with self.subTest(q=q, m=m, projected=projected):
                kernel = p.Kernel(q, m=m, mean_zero=projected, modes=6)
                cert = p.sector_certificate(kernel, k, lo, hi, 'upper_tail', independent=True)
                self.assertTrue(p.verify_sector(cert, expected_q=q, expected_m=m,
                                               expected_mean_zero=projected, expected_k=k))

    def test_tail_positivity_is_strict(self):
        kernel = p.Kernel([0, 1], modes=2)
        for kind in ('lower', 'upper_tail'):
            with self.assertRaises(ValueError):
                kernel.matrix(kernel.beta, kind)
        kernel.matrix(kernel.beta, 'upper_ritz')

    def test_full_constant_and_degenerate_angular_equality(self):
        for c in [0, F(-7, 3), F(5, 2)]:
            cert = p.full_ground([c], modes=2, max_modes=2, bits=8, max_m=0)
            self.assertTrue(p.verify_full(cert, expected_q=[c], expected_mean_zero=True))
            self.assertEqual(F(cert['lower']), 2+c)
            self.assertEqual(F(cert['upper']), 2+c)
            # Equality does not exclude the omitted m=1 eigenvalue tie.
            self.assertEqual(F(cert['angular_tail_lower']), 2+c)

    def test_resource_limit_retains_angular_floor(self):
        cert = p.full_ground([0, 1, 1], modes=3, max_modes=3, bits=12, max_m=0)
        self.assertTrue(p.verify_full(cert))
        self.assertFalse(cert['angular_tail_cannot_improve_best_upper'])
        self.assertEqual(cert['status'], 'certified_bound_open_gap')
        self.assertEqual(F(cert['lower']), F(cert['angular_tail_lower']))
        self.assertLess(F(cert['lower']), F(cert['sectors'][0]['lower']))

    def test_full_ground_can_be_nonaxisymmetric(self):
        cert = p.full_ground([0, 1, 1], modes=6, max_modes=6, bits=18,
                             max_m=1, tolerance=F(1, 1000))
        self.assertTrue(p.verify_full(cert))
        self.assertEqual(len(cert['sectors']), 2)
        self.assertTrue(cert['angular_tail_cannot_improve_best_upper'])
        self.assertLess(F(cert['sectors'][1]['upper']), F(cert['sectors'][0]['lower']))
        self.assertGreater(F(cert['lower']), F(215, 100))
        self.assertLess(F(cert['upper']), F(216, 100))

    def test_scope_input_and_endpoint_tampering(self):
        cert = p.certify_sector([0, 1], modes=4, bits=12)
        self.assertFalse(p.verify_sector(cert, expected_q=[0, 1, 1]))
        self.assertFalse(p.verify_sector(cert, expected_mean_zero=False))
        self.assertFalse(p.verify_sector(cert, expected_m=1))
        for key, value in [('projection', 'none'), ('degree_start', 0),
                           ('scope', 'all_mean_zero_H1_functions_on_unit_S2'),
                           ('lower', '3'), ('upper', '1'), ('radial_tail_lower', '1000')]:
            bad = deepcopy(cert)
            bad[key] = value
            self.assertFalse(p.verify_sector(bad), key)

    def test_missing_and_reordered_sectors_rejected(self):
        cert = p.full_ground([0, 1, 1], modes=3, max_modes=3, bits=10, max_m=1)
        for sectors in [cert['sectors'][1:], list(reversed(cert['sectors']))]:
            bad = deepcopy(cert)
            bad['sectors'] = sectors
            self.assertFalse(p.verify_full(bad))
        for key, value in [('angular_tail_lower', '1000'), ('mean_zero', False),
                           ('scope', 'single_azimuth_sector_only'), ('lower', '3')]:
            bad = deepcopy(cert)
            bad[key] = value
            self.assertFalse(p.verify_full(bad), key)

    def test_fixed_potential_inequality_semantics(self):
        for threshold, status in [(F(2), 'proved'), (F(21, 10), 'refuted_by_spectral_existence')]:
            cert = p.weighted_poincare([0], threshold, modes=2, max_modes=2, bits=8)
            self.assertEqual(cert['status'], status)
            self.assertTrue(cert['not_a_universal_Gagliardo_Nirenberg_constant_certificate'])
            self.assertTrue(p.verify_inequality(cert, expected_q=[0], expected_threshold=threshold))
            self.assertFalse(p.verify_inequality(cert, expected_threshold=threshold+1))
            self.assertFalse(p.verify_inequality(cert, expected_q=[1]))
        cert = p.weighted_poincare([0, 1, 1], 2, modes=3, max_modes=3, bits=10, max_m=0)
        self.assertEqual(cert['status'], 'undetermined')


if __name__ == '__main__':
    unittest.main()
