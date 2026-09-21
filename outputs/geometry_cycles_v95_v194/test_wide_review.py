"""Independent degree-24 tail and original exponential-potential review."""
from copy import deepcopy
from fractions import Fraction as F
from math import factorial
from unittest.mock import patch
import unittest
import wide_potential as w
from round01.test_witness_review import associated_polynomial, mul, integral


def direct_kernel(q, m, start, n):
    end = start+n+len(q)-2
    bases = {l: associated_polynomial(l, m) for l in range(start, end+1)}
    weight = [F(1)]
    for _ in range(m):
        weight = mul(weight, [1, 0, -1])
    def inner(i, j, polynomial):
        return integral(mul(mul(weight, mul(bases[i], bases[j])), polynomial))
    degrees = list(range(start, start+n))
    mass = [inner(l, l, [1]) for l in degrees]
    form = [[inner(i, j, q)+(i*(i+1)*mass[ii] if i == j else 0)
             for j in degrees] for ii, i in enumerate(degrees)]
    tails = [(l, inner(l, l, [1]), [(i, inner(k, l, q)) for i, k in enumerate(degrees)
                                  if inner(k, l, q)]) for l in range(start+n, end+1)]
    return form, mass, tails


class WideIndependentReview(unittest.TestCase):
    def test_v115_is_a_standalone_full_sphere_bound(self):
        q = [0]*24+[1]
        with patch.object(w, 'assemble', side_effect=AssertionError('matrix assembly forbidden')):
            cert = w.cheap_ground_enclosure(q)
            self.assertTrue(w.verify_cheap_ground(cert, expected_q=q, expected_mean_zero=True))
        self.assertEqual(F(cert['lower']), 2)
        self.assertEqual(F(cert['upper']), F(451, 225))
        x = next(row for row in cert['trials'] if row['sphere_function'] == 'x')
        self.assertEqual(F(x['potential_integral']), F(4, 675))
        self.assertEqual(F(x['radial_mass']), F(4, 3))
        ordinary = w.cheap_ground_enclosure([-5], mean_zero=False)
        self.assertEqual(F(ordinary['lower']), -5)
        self.assertEqual(F(ordinary['upper']), -5)
        self.assertTrue(w.verify_cheap_ground(ordinary, expected_mean_zero=False))
        self.assertFalse(w.verify_cheap_ground(ordinary, expected_mean_zero=True))
        for key, value in [('scope', 'single_azimuth_sector_only'), ('upper', '2'), ('lower', '3')]:
            bad = deepcopy(cert)
            bad[key] = value
            self.assertFalse(w.verify_cheap_ground(bad), key)

    def test_degree24_all_tail_columns_with_independent_integrals(self):
        q = [0]*24+[1]
        for m in (0, 1, 2):
            with self.subTest(m=m):
                k = w.Kernel(q, m=m, modes=2)
                form, mass, tail = direct_kernel(q, m, max(1, m), 2)
                self.assertEqual(k.a, form)
                self.assertEqual(k.mass, mass)
                self.assertEqual([(l, h, sorted(column)) for l, h, column in k.couplings], tail)
                self.assertEqual(len(tail), 24)
                self.assertEqual(tail[-1][0], k.tail_start+23)
                self.assertTrue(tail[-1][2])
                x = F(-1)
                expected = [row[:] for row in form]
                for i in range(2):
                    expected[i][i] -= x*mass[i]
                for l, h, column in tail:
                    for i, a in column:
                        for j, b in column:
                            expected[i][j] -= a*b/(h*(l*(l+1)+k.qlo-x))
                self.assertEqual(k.matrix(x, 'lower'), expected)

    def test_kernel_missing_last_column_and_mass_tampering(self):
        cert = w.kernel_certificate([0]*24+[1], m=1, modes=2)
        self.assertTrue(w.verify_kernel(cert, expected_q=[0]*24+[1]))
        bad = deepcopy(cert)
        bad['tail_couplings'] = bad['tail_couplings'][:-1]
        self.assertFalse(w.verify_kernel(bad))
        bad = deepcopy(cert)
        bad['tail_couplings'][-1]['mass'] = '1'
        self.assertFalse(w.verify_kernel(bad))
        self.assertFalse(w.verify_kernel(cert, expected_q=[0]*24+[2]))

    def test_ranges_cover_whole_interval_and_bind_input(self):
        q = [0]*24+[1]
        cert = w.bernstein_range(q, depth=1)
        self.assertEqual(F(cert['lower']), 0)
        self.assertEqual(F(cert['upper']), 1)
        self.assertTrue(w.verify_range(cert, expected_q=q))
        bad = deepcopy(cert)
        bad['proof']['cells'].pop()
        self.assertFalse(w.verify_range(bad))
        bad = deepcopy(cert)
        bad['lower'] = '1/2'
        self.assertFalse(w.verify_range(bad))
        self.assertFalse(w.verify_range(cert, expected_q=[0]*24+[-1]))
        with self.assertRaises(ValueError):
            w.Kernel([0]*24+[-1], modes=2, range_proof=cert)

    def test_reused_sturm_range_does_not_run_search(self):
        q = [0]*8+[1]
        cert = w.sturm_range(q)
        with patch.object(w.exact_extrema, 'maximize', side_effect=AssertionError('search forbidden')):
            replayed = w.reuse_range(q, cert)
            kernel = w.Kernel(q, modes=2, range_proof=replayed)
            self.assertEqual(kernel.qlo, 0)
            self.assertEqual(kernel.qhi, 1)

    def test_full_angular_coverage_and_open_limit(self):
        q = [0]*12+[1]
        cert = w.full_ground(q, modes=2, max_modes=2, bits=12, max_m=1, tolerance=F(1, 1000))
        self.assertTrue(w.verify_full(cert, expected_q=q, expected_mean_zero=True))
        self.assertEqual(len(cert['sectors']), 2)
        self.assertTrue(cert['angular_tail_cannot_improve_best_upper'])
        self.assertLess(F(cert['sectors'][1]['upper']), F(cert['sectors'][0]['lower']))
        limited = w.full_ground(q, modes=2, max_modes=2, bits=12, max_m=0)
        self.assertTrue(w.verify_full(limited))
        self.assertEqual(limited['status'], 'certified_bound_open_gap')
        self.assertEqual(F(limited['lower']), F(limited['angular_tail_lower']))
        bad = deepcopy(cert)
        bad['sectors'] = bad['sectors'][1:]
        self.assertFalse(w.verify_full(bad))
        bad = deepcopy(cert)
        bad['scope'] = 'single_azimuth_sector_only'
        self.assertFalse(w.verify_full(bad))

    def test_taylor_formula_uniform_absolute_remainder(self):
        for a, degree in [(F(0), 0), (F(-3), 4), (F(1, 2), 12), (F(1, 10), 24)]:
            with self.subTest(a=a, degree=degree):
                cert = w.exponential_approximation(a, degree)
                self.assertTrue(w.verify_approximation(cert, expected_a=a))
                polynomial = list(map(F, cert['polynomial_coefficients']))
                self.assertEqual(polynomial, [a**k/F(factorial(k)) for k in range(degree+1)] if a else [F(1)])
                ceil = (abs(a).numerator+abs(a).denominator-1)//abs(a).denominator
                self.assertEqual(F(cert['absolute_error_bound']), F(3**ceil)*abs(a)**(degree+1)/factorial(degree+1))
                self.assertEqual(cert['interval'], ['-1', '1'])
                self.assertFalse(cert['formal_assistant_checked'])

    def test_taylor_source_and_error_tampering(self):
        cert = w.exponential_approximation(F(1, 2), 4)
        self.assertFalse(w.verify_approximation(cert, expected_a=F(-1, 2)))
        for key, value in [('absolute_error_bound', '0'), ('interval', ['0', '1']),
                           ('exponential_envelope', '1'), ('degree', 5),
                           ('original_potential', {'kind': 'exp_linear', 'a': '-1/2'})]:
            bad = deepcopy(cert)
            bad[key] = value
            self.assertFalse(w.verify_approximation(bad), key)

    def test_original_exponential_operator_transfer_and_source_binding(self):
        approximation = w.exponential_approximation(F(1, 2), 4)
        q = approximation['polynomial_coefficients']
        spectrum = w.full_ground(q, modes=2, max_modes=2, bits=12, max_m=1)
        cert = w.perturb_certificate(approximation, spectrum)
        self.assertTrue(w.verify_exponential(cert, expected_a=F(1, 2)))
        delta = F(approximation['absolute_error_bound'])
        self.assertEqual(F(cert['lower']), F(spectrum['lower'])-delta)
        self.assertEqual(F(cert['upper']), F(spectrum['upper'])+delta)
        self.assertEqual(F(cert['exact_width']), F(spectrum['exact_width'])+2*delta)
        self.assertFalse(w.verify_exponential(cert, expected_a=F(-1, 2)))
        for key, value in [('perturbation_bound', '0'), ('scope', 'single_azimuth_sector_only'),
                           ('original_potential', {'kind': 'exp_linear', 'a': '1'}), ('lower', '100')]:
            bad = deepcopy(cert)
            bad[key] = value
            self.assertFalse(w.verify_exponential(bad), key)
        other_q = w.full_ground([1], modes=2, max_modes=2, bits=8)
        with self.assertRaises(ValueError):
            w.perturb_certificate(approximation, other_q)
        ordinary = w.full_ground(q, mean_zero=False, modes=2, max_modes=2, bits=10)
        with self.assertRaises(ValueError):
            w.perturb_certificate(approximation, ordinary)

    def test_exponential_zero_exact_calibration_and_unsupported_inputs(self):
        cert = w.exponential_ground(0, degree=0, modes=2, max_modes=2, bits=8)
        self.assertTrue(w.verify_exponential(cert, expected_a=0))
        self.assertEqual(F(cert['lower']), 3)
        self.assertEqual(F(cert['upper']), 3)
        for q in ([0]*25+[1], [0.5], [True]):
            with self.assertRaises(ValueError):
                w.potential(q)
        for a, degree in [(4, 4), (1, 25), (True, 4)]:
            with self.assertRaises(ValueError):
                w.exponential_approximation(a, degree)


if __name__ == '__main__':
    unittest.main()
