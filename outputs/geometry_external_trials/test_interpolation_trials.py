"""Independent reduction checks and adversarial replay of the two GN records."""
from fractions import Fraction as F
from pathlib import Path
from unittest.mock import patch
import copy
import hashlib
import json
import unittest

import interpolation_trials as adapter


HERE = Path(__file__).resolve().parent
P1 = [F(0), F(1)]
P3 = [F(0), F(-3, 2), F(0), F(5, 2)]


def multiply(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i + j] += x * y
    return out


def integral(p):
    return sum((x / (i + 1) for i, x in enumerate(p) if i % 2 == 0), F(0))


def direct_quantities(p):
    mass = integral(multiply(p, p))
    dp = [F(i) * p[i] for i in range(1, len(p))] or [F(0)]
    energy = integral(multiply([F(1), F(0), F(-1)], multiply(dp, dp)))
    square = multiply(p, p)
    fourth = integral(multiply(square, square))
    return mass, energy, fourth


def combination(f, g, a, b):
    return [(a * f[i] if i < len(f) else F(0)) + (b * g[i] if i < len(g) else F(0))
            for i in range(max(len(f), len(g)))]


def evaluate(p, x):
    return sum((F(c) * x ** i for i, c in enumerate(p)), F(0))


class InterpolationTrialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.records = {}
        for name in ('legendre_P1_P3', 'centered_caps_5_6'):
            cls.records[name] = json.loads((HERE / 'certificates' / 'interpolation' / (name + '.json')).read_text())
        cls.legendre = cls.records['legendre_P1_P3']

    def test_both_saved_certificates_accept_and_gap_is_exact(self):
        for cert in self.records.values():
            self.assertTrue(adapter.verify(cert, expected_basis=cert['basis']))
            self.assertEqual(F(cert['upper']) - F(cert['lower']), F(1, 10**10))
            self.assertEqual(F(cert['gap']), F(1, 10**10))

    def test_p1_area_normalization_calibration(self):
        mass, energy, fourth = direct_quantities(P1)
        self.assertEqual((mass, energy, fourth), (F(1, 3), F(2, 3), F(1, 5)))
        self.assertEqual(fourth / (4 * mass * energy), F(9, 40))
        # Certificate ratio is pi*K_area = K_probability/4.
        chart = adapter.chart_moments(P1, P3)
        self.assertEqual(evaluate(chart['numerator'], 0) / evaluate(chart['denominator'], 0), F(9, 40))

    def test_independent_legendre_numerator_and_denominator_coefficients(self):
        chart = adapter.chart_moments(P1, P3)
        self.assertEqual(chart['numerator'], [F(1, 5), F(8, 35), F(46, 105), F(48, 385), F(241, 5005)])
        self.assertEqual(chart['mass'], [F(1, 3), F(0), F(1, 7)])
        self.assertEqual(chart['energy'], [F(2, 3), F(0), F(12, 7)])
        self.assertEqual(chart['denominator'], [F(8, 9), F(0), F(8, 3), F(0), F(48, 49)])

    def test_independent_direct_integration_of_both_charts(self):
        for cert in self.records.values():
            basis = [list(map(F, p)) for p in cert['basis']]
            for j, chart in enumerate(cert['charts']):
                for s in [F(-1), F(-3, 7), F(0), F(2, 5), F(1)]:
                    mass, energy, fourth = direct_quantities(combination(basis[j], basis[1-j], F(1), s))
                    self.assertEqual(evaluate(chart['moments']['mass'], s), mass)
                    self.assertEqual(evaluate(chart['moments']['energy'], s), energy)
                    self.assertEqual(evaluate(chart['moments']['numerator'], s), fourth)
                    self.assertEqual(evaluate(chart['moments']['denominator'], s), 4 * mass * energy)
                    self.assertGreater(mass, 0)
                    self.assertGreater(energy, 0)

    def test_two_chart_projective_coverage_including_axes_and_negative_scales(self):
        for cert in self.records.values():
            basis = [list(map(F, p)) for p in cert['basis']]
            for a in range(-3, 4):
                for b in range(-3, 4):
                    if a == b == 0:
                        continue
                    if abs(a) >= abs(b):
                        j, s = 0, F(b, a)
                    else:
                        j, s = 1, F(a, b)
                    self.assertTrue(-1 <= s <= 1)
                    mass, energy, fourth = direct_quantities(combination(*basis, F(a), F(b)))
                    moments = cert['charts'][j]['moments']
                    chart_ratio = evaluate(moments['numerator'], s) / evaluate(moments['denominator'], s)
                    self.assertEqual(chart_ratio, fourth / (4 * mass * energy))
                    self.assertLessEqual(chart_ratio, F(cert['upper']))

    def test_analytical_centering_and_positive_gram_determinants(self):
        for cert in self.records.values():
            f, g = [list(map(F, p)) for p in cert['basis']]
            self.assertEqual(integral(f), 0)
            self.assertEqual(integral(g), 0)
            determinant = integral(multiply(f, f)) * integral(multiply(g, g)) - integral(multiply(f, g)) ** 2
            self.assertGreater(determinant, 0)
        self.assertEqual(integral(adapter.centered_cap(5)), 0)
        self.assertEqual(integral(adapter.centered_cap(6)), 0)

    def test_verifier_never_calls_floating_proposal_or_maximization_search(self):
        with patch.object(adapter, 'propose', side_effect=AssertionError('proposal forbidden in replay')), \
             patch.object(adapter.extrema, 'maximize', side_effect=AssertionError('search forbidden in replay')):
            for cert in self.records.values():
                self.assertTrue(adapter.verify(cert))

    def test_missing_duplicated_or_reordered_chart_is_rejected(self):
        altered = copy.deepcopy(self.legendre)
        altered['charts'].pop()
        self.assertFalse(adapter.verify(altered))
        altered = copy.deepcopy(self.legendre)
        altered['charts'][1] = copy.deepcopy(altered['charts'][0])
        self.assertFalse(adapter.verify(altered))
        altered = copy.deepcopy(self.legendre)
        altered['charts'].reverse()
        self.assertFalse(adapter.verify(altered))

    def test_tampered_upper_gap_or_residual_is_rejected(self):
        for field, value in [('upper', '1/10'), ('lower', '1/10'), ('gap', '0')]:
            altered = copy.deepcopy(self.legendre)
            altered[field] = value
            self.assertFalse(adapter.verify(altered))
        altered = copy.deepcopy(self.legendre)
        altered['upper'] = str(F(altered['upper']) + 1)
        altered['gap'] = str(F(altered['upper']) - F(altered['lower']))
        self.assertFalse(adapter.verify(altered))
        altered = copy.deepcopy(self.legendre)
        altered['charts'][0]['upper_residual_nonpositive']['maximum_upper'] = '-1'
        self.assertFalse(adapter.verify(altered))

    def test_tampered_witness_is_rejected(self):
        for field, value in [('parameter', '2'), ('parameter', '0'), ('chart', True), ('chart', 2)]:
            altered = copy.deepcopy(self.legendre)
            altered['witness'][field] = value
            self.assertFalse(adapter.verify(altered))

    def test_tampered_basis_and_external_basis_binding(self):
        altered = copy.deepcopy(self.legendre)
        altered['basis'][0][1] = '2'
        self.assertFalse(adapter.verify(altered))
        self.assertTrue(adapter.verify(self.legendre, expected_basis=[P1, P3]))
        self.assertFalse(adapter.verify(self.legendre, expected_basis=[P3, P1]))
        other = self.records['centered_caps_5_6']
        self.assertTrue(adapter.verify(other))
        self.assertFalse(adapter.verify(other, expected_basis=[P1, P3]))

    def test_zero_mean_linear_independence_and_exact_inputs_required(self):
        for basis in [[[1], P1], [P1, P1], [[0], P1], [[0, 1.0], P3],
                      [[0, True], P3], [P1], [P1, [0] * 7 + [1]]]:
            with self.assertRaises(ValueError):
                adapter.basis_input(basis)

    def test_nonpositive_denominator_or_wrong_interval_proof_rejected(self):
        for field in ['mass_positive', 'energy_positive']:
            altered = copy.deepcopy(self.legendre)
            altered['charts'][0][field]['maximum_upper'] = '0'
            self.assertFalse(adapter.verify(altered))
        altered = copy.deepcopy(self.legendre)
        altered['charts'][0]['mass_positive']['interval'] = ['0', '1']
        self.assertFalse(adapter.verify(altered))
        altered = copy.deepcopy(self.legendre)
        altered['charts'][0]['moments']['denominator'][0] = '0'
        self.assertFalse(adapter.verify(altered))

    def test_scope_and_area_factor_cannot_be_relabelled(self):
        for field, value in [('scope', 'all_H1_functions_on_S2'),
                             ('ratio', 'K4_for_probability_measure'), ('format', 'unverified')]:
            altered = copy.deepcopy(self.legendre)
            altered[field] = value
            self.assertFalse(adapter.verify(altered))
        altered = copy.deepcopy(self.legendre)
        altered['charts'][0]['moments']['denominator'] = [str(F(x) / 4) for x in altered['charts'][0]['moments']['denominator']]
        self.assertFalse(adapter.verify(altered))

    def test_independent_witness_recalculation(self):
        for cert in self.records.values():
            basis = [list(map(F, p)) for p in cert['basis']]
            j, s = cert['witness']['chart'], F(cert['witness']['parameter'])
            m, e, n = direct_quantities(combination(basis[j], basis[1-j], F(1), s))
            self.assertEqual(F(cert['lower']), n / (4 * m * e))


if __name__ == '__main__':
    unittest.main()
