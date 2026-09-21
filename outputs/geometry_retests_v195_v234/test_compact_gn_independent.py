"""C16 independent review via ordinary t-polynomial integration and identities."""
import copy
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import unittest

import compact_gn as c


def multiply(a, b):
    out = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return out


def probability_mean(p):
    return sum((x/F(i+1) for i, x in enumerate(p) if i % 2 == 0), F(0))


def t_moments(p):
    derivative = [i*p[i] for i in range(1, len(p))]
    square = multiply(p, p)
    return {'mean': probability_mean(p), 'mass': probability_mean(square),
            'energy': probability_mean(multiply([1, 0, -1], multiply(derivative, derivative))),
            'quartic': probability_mean(multiply(square, square))}


def cap_t(n):
    p = [F(comb(n, j)) for j in range(n+1)]
    p[0] -= F(2**n, n+1)
    return p


def y_to_t(terms):
    degree = max((n for n, _ in terms), default=0)
    p = [F(0)]*(degree+1)
    for n, coefficient in terms:
        coefficient = F(coefficient)
        for j in range(n+1):
            p[j] += coefficient*F(comb(n, j), 2**n)
    return p


def centered_quartic_kernel(indices):
    result = F(0)
    for mask in range(16):
        degree, coefficient = 0, F(1)
        for i, n in enumerate(indices):
            if mask & (1 << i):
                degree += n
            else:
                coefficient *= -F(1, n+1)
        result += coefficient/F(degree+1)
    return result


def reference_Q(n):
    return (24*n**4+30*n*n+18*n)/(12*n**4+31*n**3+27*n*n+9*n+1)


class CompactGNIndependentReview(unittest.TestCase):
    def test_sparse_mean_and_exact_input_validation(self):
        cert = c.mean_certificate([[64, '3/7'], [0, '-2/5'], [64, '4/7']])
        self.assertEqual(cert['terms'], [[0, '-2/5'], [64, '1']])
        self.assertEqual(F(cert['mean']), F(1, 65)-F(2, 5))
        for raw in ([[1, 1.0]], [[True, 1]], [[-1, 1]]):
            with self.assertRaises(ValueError):
                c.mean_certificate(raw)

    def test_centering_and_scale_retain_function_identity(self):
        function = c.centered_function([[1, 2], [3, -1]], scale=7)
        expected_mean_removed = F(2, 2)-F(1, 4)
        self.assertEqual(function['terms'], [[0, str(-7*expected_mean_removed)], [1, '14'], [3, '-7']])
        self.assertEqual(function['mean'], '0')
        self.assertTrue(c.verify(function))

    def test_mass_kernel_by_independent_low_degree_integrals(self):
        for a, b in ((1, 1), (1, 4), (3, 5), (0, 64)):
            expected = F(1, a+b+1)-F(1, (a+1)*(b+1))
            self.assertEqual(c.mass_kernel(a, b), expected)
            if max(a, b) <= 5:
                pa, pb = y_to_t([[a, 1], [0, -F(1, a+1)]]), y_to_t([[b, 1], [0, -F(1, b+1)]])
                self.assertEqual(probability_mean(multiply(pa, pb)), expected)

    def test_energy_kernel_by_t_derivative_integrals(self):
        for a, b in ((1, 1), (2, 5), (3, 4), (0, 0), (0, 64)):
            expected = F(a*b, (a+b)*(a+b+1)) if a and b else F(0)
            self.assertEqual(c.energy_kernel(a, b), expected)
            if a and b and max(a, b) <= 5:
                pa, pb = y_to_t([[a, 1]]), y_to_t([[b, 1]])
                da = [i*pa[i] for i in range(1, len(pa))]
                db = [i*pb[i] for i in range(1, len(pb))]
                self.assertEqual(probability_mean(multiply([1, 0, -1], multiply(da, db))), expected)

    def test_mass_energy_gram_off_diagonal_terms(self):
        weights = [[1, 2], [3, '-1/2'], [7, '3/5']]
        terms = c.centered_function(weights)['terms']
        direct = t_moments(y_to_t(terms))
        self.assertEqual(F(c.mass_form(weights)['value']), direct['mass'])
        self.assertEqual(F(c.energy_form(weights)['value']), direct['energy'])
        self.assertNotEqual(F(c.mass_form(weights)['gram'][0][1]), 0)

    def test_quartic_kernel_all_sixteen_subset_terms(self):
        terms = c.centered_function([[4, 1]])['terms']
        self.assertEqual(F(c.quartic_moment(terms)['quartic']), centered_quartic_kernel([4]*4))
        f, g = c.sparse(c.centered_function([[2, 1]])['terms']), c.sparse(c.centered_function([[5, 1]])['terms'])
        chart = c.chart_moments(f, g)
        for j in range(5):
            expected = comb(4, j)*centered_quartic_kernel([2]*(4-j)+[5]*j)
            self.assertEqual(chart['numerator'][j], expected)

    def test_original_n4_n16_n64_all_moments_by_dense_t_integration(self):
        known = {4: F(6696, 5525), 16: F(1580832, 920465), 64: F(402777216, 209564225)}
        for n in (4, 16, 64):
            with self.subTest(n=n):
                raw = cap_t(n)
                direct = t_moments(raw)
                cert = c.bridge_original(raw, n)
                moments = cert['evaluation']['ratio']['moments']
                for field, expected in direct.items():
                    self.assertEqual(F(moments[field]), expected)
                ratio = cert['evaluation']['ratio']
                self.assertEqual(F(ratio['probability_ratio']), known[n])
                self.assertEqual(F(ratio['pi_times_area_K']), known[n]/4)
                self.assertTrue(c.verify(cert, expected_n=n, expected_original_coefficients=raw))

    def test_original_scale_factors_for_n64(self):
        normalized, original = c.evaluate_cap(64, False), c.evaluate_cap(64, True)
        for field, power in (('mass', 128), ('energy', 128), ('quartic', 256)):
            self.assertEqual(F(original['ratio']['moments'][field]),
                             2**power*F(normalized['ratio']['moments'][field]))
            self.assertEqual(F(original['moment_scale_factors'][field]), 2**power)
        self.assertEqual(original['ratio']['probability_ratio'], normalized['ratio']['probability_ratio'])
        self.assertEqual(original['function']['scale'], str(2**64))

    def test_original_bridge_binds_every_coefficient(self):
        raw = cap_t(64)
        bad = list(raw)
        bad[31] += 1
        with self.assertRaises(ValueError):
            c.bridge_original(bad, 64)
        cert = c.bridge_original(raw, 64)
        self.assertEqual(cert['verified_coefficient_count'], 65)
        self.assertFalse(c.verify(cert, expected_n=16))
        self.assertFalse(c.verify(cert, expected_original_coefficients=bad))
        cert['evaluation']['moment_scale_factors']['quartic'] = '1'
        self.assertFalse(c.verify(cert))

    def test_probability_and_area_pi_factor_four(self):
        cert = c.ratio([[0, '-1/2'], [1, 1]])
        self.assertEqual(F(cert['moments']['mass']), F(1, 12))
        self.assertEqual(F(cert['moments']['energy']), F(1, 6))
        self.assertEqual(F(cert['moments']['quartic']), F(1, 80))
        self.assertEqual(F(cert['probability_ratio']), F(9, 10))
        self.assertEqual(F(cert['pi_times_area_K']), F(9, 40))
        self.assertIn('not_a_full_space_GN_upper_bound', cert['scope'])

    def test_n64_full_residual_explicit_coefficients_and_norm(self):
        n = 64
        terms = [[0, -F(1, n+1)], [n, 1]]
        cert = c.residual(terms)
        m = cert['ratio']['moments']
        mass, energy, quartic = (F(m[field]) for field in ('mass', 'energy', 'quartic'))
        a = F(1, n+1)
        expected = {3*n: -2*energy/quartic, 2*n: 6*a*energy/quartic,
                    n: n*(n+1)+energy/mass-6*a*a*energy/quartic, n-1: F(-n*n)}
        expected[0] = -sum((v/F(power+1) for power, v in expected.items()), F(0))
        self.assertEqual(dict((power, F(value)) for power, value in cert['residual']), expected)
        norm = sum((x*y/F(i+j+1) for i, x in expected.items() for j, y in expected.items()), F(0))
        self.assertEqual(F(cert['residual_norm_squared']), norm)
        self.assertEqual(cert['full_degree'], 192)
        self.assertEqual(cert['residual_mean'], '0')
        self.assertEqual(cert['inner_trial'], '0')
        self.assertFalse(cert['degree_truncation_used'])
        cert['residual'] = [row for row in cert['residual'] if row[0] != 192]
        self.assertFalse(c.verify(cert))

    def test_helpers_keep_powers_beyond_public_input_degree_cap(self):
        n = c.MAX_POWER
        terms = [[0, -F(1, n+1)], [n, 1]]
        cert = c.residual(terms)
        self.assertEqual(cert['full_degree'], 3*n)
        self.assertTrue(any(power == 3*n for power, _ in cert['residual']))
        self.assertGreater(F(cert['residual_norm_squared']), 0)
        self.assertTrue(c.verify(cert))
        direction = c.direction(terms)
        self.assertTrue(direction['strict_local_improvement'])
        self.assertTrue(c.verify(direction))

    def test_direction_identity_has_complete_residual_and_no_finite_step_claim(self):
        terms = c.evaluate_cap(64, False)['function']['terms']
        result = c.direction(terms)
        residual = result['residual_certificate']
        expected = 2*F(residual['ratio']['pi_times_area_K'])*F(residual['residual_norm_squared'])/F(residual['ratio']['moments']['energy'])
        self.assertEqual(F(result['derivative_pi_times_K']), expected)
        self.assertGreater(expected, 0)
        self.assertFalse(result['finite_step_improvement_claimed'])
        self.assertTrue(c.verify(result))

    def test_plane_two_complete_charts_and_scope(self):
        f = c.evaluate_cap(16, False)['function']['terms']
        h = c.evaluate_cap(64, False)['function']['terms']
        plane = c.plane([f, h], tolerance=F(1, 10**6))
        self.assertTrue(c.verify_plane(plane, [f, h]))
        self.assertEqual(len(plane['charts']), 2)
        self.assertFalse(plane['universal_GN_upper_bound'])
        self.assertEqual(plane['scope'], 'global_in_recorded_two_function_plane_only')
        for parameter in (F(-3), F(-1), F(0), F(1, 3), F(1), F(4)):
            terms = list(f)+[[power, parameter*F(value)] for power, value in h]
            value = F(c.ratio(terms)['pi_times_area_K'])
            self.assertLessEqual(value, F(plane['upper_pi_times_K']))
        bad = copy.deepcopy(plane)
        bad['charts'].pop()
        self.assertFalse(c.verify(bad))
        bad = copy.deepcopy(plane)
        bad['universal_GN_upper_bound'] = True
        self.assertFalse(c.verify(bad))

    def test_original_plane_scaling_and_coordinate_pullback_are_bound(self):
        f = c.evaluate_cap(4, True)['function']['terms']
        h = c.evaluate_cap(64, True)['function']['terms']
        plane = c.plane([f, h], F(1, 10**6))
        self.assertTrue(c.verify_plane(plane, [f, h]))
        self.assertEqual(plane['basis'], [f, h])
        self.assertEqual(list(map(F, plane['positive_basis_scales'])), [F(1, 2**4), F(1, 2**64)])
        coefficients = list(map(F, plane['witness']['original_basis_coefficients']))
        recombined = {}
        for terms, factor in zip((f, h), coefficients):
            for power, raw in terms:
                recombined[power] = recombined.get(power, F(0))+factor*F(raw)
        recombined = [[power, str(value)] for power, value in sorted(recombined.items()) if value]
        self.assertEqual(recombined, plane['witness']['terms'])
        for field in ('scale', 'chart_basis', 'pullback'):
            bad = copy.deepcopy(plane)
            if field == 'scale':
                bad['positive_basis_scales'][1] = '1'
            elif field == 'chart_basis':
                bad['chart_basis'][0][0][1] = '0'
            else:
                bad['witness']['original_basis_coefficients'][0] = '1'
            self.assertFalse(c.verify(bad))

    def test_plane_union_budget_rejects_before_returning_unreplayable_proof(self):
        # Regression: each legal input has 33 terms; the actual mixed witness
        # has 65. The first implementation returned a certificate verify rejected.
        f = c.centered_function([[i, 1 if i == 1 else F(1, 1024)] for i in range(1, 33)])['terms']
        h = c.centered_function([[i, 1 if i == 33 else F(1, 1024)] for i in range(33, 65)])['terms']
        self.assertEqual((len(f), len(h)), (33, 33))
        with self.assertRaises(ValueError):
            c.plane([f, h], F(1, 10**5))
        # Centering itself adds a constant term and must respect the same cap.
        with self.assertRaises(ValueError):
            c.centered_function([[i, 1] for i in range(1, c.MAX_TERMS+1)])

    def test_integer_family_law_and_exact_limit_are_trial_only(self):
        previous = F(0)
        for n in (1, 2, 4, 16, 64, 10**12):
            cert = c.family_value(n)
            value = F(cert['probability_ratio'])
            self.assertEqual(value, reference_Q(F(n)))
            self.assertTrue(previous < value < 2)
            previous = value
            self.assertFalse(cert['full_space_sharpness_claimed'])
            self.assertTrue(c.verify(cert, expected_n=n))
        law = c.family_law()
        self.assertEqual(law['n_domain'], 'all_integers_n>=1')
        self.assertEqual(law['limit_probability_ratio'], '2')
        self.assertEqual(law['limit_pi_times_area_K'], '1/2')
        self.assertFalse(law['sharp_constant_claimed'])

    def test_discrete_monotonicity_does_not_claim_continuous_monotonicity(self):
        # Q'(1)=-3/25: real n near 1 decreases, though the integer sequence increases.
        self.assertLess(reference_Q(F(11, 10)), reference_Q(F(1)))
        with self.assertRaises(ValueError):
            c.family_value(F(11, 10))
        law = c.family_law()
        law['n_domain'] = 'all_real_n>=1'
        self.assertFalse(c.verify(law))

    def test_every_saved_c16_certificate_replays(self):
        formats = {'compact_mean_v215', 'compact_centered_v216', 'compact_mass_form_v217',
                   'compact_energy_form_v218', 'compact_quartic_v219', 'compact_moments_v220',
                   'compact_ratio_v220', 'compact_residual_v221', 'compact_direction_v221',
                   'compact_plane_v222', 'compact_family_law_v223', 'compact_family_value_v223',
                   'compact_cap_v224', 'compact_original_bridge_v224'}
        count, original_n = 0, set()
        def walk(value, path):
            nonlocal count
            if isinstance(value, dict):
                if value.get('format') in formats:
                    with self.subTest(path=path):
                        self.assertTrue(c.verify(value))
                    count += 1
                    if value['format'] == 'compact_original_bridge_v224':
                        original_n.add(value['n'])
                    return
                for key, item in value.items():
                    walk(item, path+'/'+str(key))
            elif isinstance(value, list):
                for i, item in enumerate(value):
                    walk(item, path+'/'+str(i))
        for path in sorted((Path(__file__).parent/'C16').rglob('*.json')):
            walk(json.loads(path.read_text()), str(path.name))
        self.assertGreaterEqual(count, 10)
        self.assertTrue({4, 16, 64}.issubset(original_n))


if __name__ == '__main__':
    unittest.main()
