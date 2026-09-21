"""Independent direct-integration checks of explicit spherical witnesses.

These helpers use Rodrigues polynomials and monomial integrals, not the
implementation's associated-Legendre multiplication matrix.
"""
from copy import deepcopy
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
import sys
import unittest

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import witness as w


def add(a, b):
    c = [F(0)] * max(len(a), len(b))
    for i, v in enumerate(a):
        c[i] += v
    for i, v in enumerate(b):
        c[i] += v
    return c


def scale(a, c):
    return [F(c)*v for v in a]


def mul(a, b):
    c = [F(0)] * (len(a)+len(b)-1)
    for i, v in enumerate(a):
        for j, w in enumerate(b):
            c[i+j] += v*w
    return c


def derivative(a):
    return [i*a[i] for i in range(1, len(a))] or [F(0)]


def integral(a):
    return sum((2*v/F(i+1) for i, v in enumerate(a) if i % 2 == 0), F(0))


def associated_polynomial(l, m):
    """Return (-1)^m D^m P_l, excluding the common (1-t^2)^(m/2)."""
    a = [F(1)]
    for _ in range(l):
        a = mul(a, [-1, 0, 1])
    for _ in range(l+m):
        a = derivative(a)
    return scale(a, F((-1)**m, 2**l*factorial(l)))


def direct_witness(q, coefficients, m, degrees, mean_zero=True):
    q, coefficients = list(map(F, q)), list(map(F, coefficients))
    polynomial, kinetic = [F(0)], [F(0)]
    for l, c in zip(degrees, coefficients):
        b = associated_polynomial(l, m)
        polynomial = add(polynomial, scale(b, c))
        kinetic = add(kinetic, scale(b, c*l*(l+1)))
    weight = [F(1)]
    for _ in range(m):
        weight = mul(weight, [1, 0, -1])
    h_polynomial = add(kinetic, mul(q, polynomial))
    removed_l0 = F(0)
    if m == 0 and mean_zero:
        removed_l0 = integral(h_polynomial)/2
        h_polynomial = add(h_polynomial, [-removed_l0])
    mass = integral(mul(weight, mul(polynomial, polynomial)))
    energy = integral(mul(weight, mul(polynomial, h_polynomial)))
    r = energy/mass
    residual = add(h_polynomial, scale(polynomial, -r))
    residual_mass = integral(mul(weight, mul(residual, residual)))
    orthogonality = integral(mul(weight, mul(residual, polynomial)))
    return {'mass': mass, 'energy': energy, 'rayleigh': r,
            'residual_mass': residual_mass, 'removed_l0': removed_l0,
            'orthogonality': orthogonality}


class IndependentAnalyticCalibration(unittest.TestCase):
    def test_simple_cartesian_counterexample(self):
        result = direct_witness([0, 0, 2], [1, F(-1, 36)], 1, [1, 3])
        self.assertEqual(result['mass'], F(505, 378))
        self.assertEqual(result['energy'], F(1805, 567))
        self.assertEqual(result['rayleigh'], F(722, 303))
        self.assertLess(result['rayleigh'], F(12, 5))
        self.assertEqual(result['orthogonality'], 0)

    def test_projected_residual_analytic_calibration(self):
        projected = direct_witness([0, 1], [1], 0, [1])
        ordinary = direct_witness([0, 1], [1], 0, [1], mean_zero=False)
        self.assertEqual(projected['removed_l0'], F(1, 3))
        self.assertEqual(projected['residual_mass'], F(8, 45))
        self.assertEqual(projected['residual_mass']/projected['mass'], F(4, 15))
        self.assertEqual(ordinary['residual_mass']/ordinary['mass'], F(3, 5))

    def test_full_interval_of_two_harmonic_counterexamples(self):
        for s in [F(-1, 21), F(-1, 36), F(-1, 40), F(1, 40), F(-1, 10)]:
            result = direct_witness([0, 0, 2], [1, s], 1, [1, 3])
            self.assertEqual(result['mass'], F(4, 3)+F(24, 7)*s*s)
            self.assertEqual(result['energy']-F(12, 5)*result['mass'], 16*s*(4+79*s)/35)


class WitnessImplementationReview(unittest.TestCase):
    def simple(self, scale=F(1)):
        return w.rayleigh_certificate([0, 0, 2], [scale, -scale/36], m=1,
                                      degrees=[1, 3], threshold=F(12, 5))

    def test_independent_forms_and_complete_residual(self):
        cases = [([0, 0, 2], [1, F(-1, 36)], 1, [1, 3]),
                 ([1, -2, 3, 0, 1], [2, -1, F(3, 7)], 0, [1, 2, 5]),
                 ([1, 0, -1, 2, 0, 0, 1], [1, F(-2, 3), F(1, 5)], 2, [2, 4, 6])]
        for q, coefficients, m, degrees in cases:
            with self.subTest(m=m):
                expected = direct_witness(q, coefficients, m, degrees)
                cert = w.rayleigh_certificate(q, coefficients, m=m, degrees=degrees)
                residual = w.residual_certificate(cert)
                self.assertTrue(w.verify(cert, expected_q=q, expected_m=m,
                                         expected_mean_zero=True, expected_coefficients=coefficients))
                self.assertTrue(w.verify(residual))
                self.assertEqual(F(cert['norm_squared_divided_by_azimuth_factor']), expected['mass'])
                self.assertEqual(F(cert['energy_divided_by_azimuth_factor']), expected['energy'])
                self.assertEqual(F(cert['rayleigh_quotient']), expected['rayleigh'])
                self.assertEqual(F(residual['residual_norm_squared_divided_by_azimuth_factor']), expected['residual_mass'])
                self.assertEqual(F(residual['removed_l0_coefficient']), expected['removed_l0'])
                self.assertEqual(F(residual['orthogonality_to_trial']), 0)
                self.assertEqual(F(residual['in_support_residual_squared'])+F(residual['outside_support_residual_squared']),
                                 F(residual['normalized_residual_squared']))
                self.assertEqual(max(row['degree'] for row in residual['residual_coefficients']),
                                 max(degrees)+len(q)-1)

    def test_explicit_counterexample_and_angular_factor(self):
        cert = self.simple()
        self.assertTrue(w.verify(cert, expected_q=[0, 0, 2], expected_threshold=F(12, 5)))
        self.assertEqual(cert['azimuth_integral_factor'], 'pi')
        self.assertEqual(cert['status'], 'strict_counterexample')
        self.assertEqual(F(cert['rayleigh_quotient']), F(722, 303))
        self.assertEqual(F(cert['strict_margin']), F(26, 1515))
        m0 = w.rayleigh_certificate([0], [1], m=0)
        self.assertEqual(m0['azimuth_integral_factor'], '2*pi')

    def test_projection_deletes_only_actual_constant_output(self):
        projected = w.rayleigh_certificate([0, 1], [1], degrees=[1])
        ordinary = w.rayleigh_certificate([0, 1], [1], degrees=[1], mean_zero=False)
        r = w.residual_certificate(projected)
        self.assertEqual(r['residual_coefficients'], [{'degree': 2, 'coefficient': '2/3'}])
        self.assertEqual(r['removed_l0_coefficient'], '1/3')
        self.assertEqual(F(r['normalized_residual_squared']), F(4, 15))
        self.assertEqual(F(w.residual_certificate(ordinary)['normalized_residual_squared']), F(3, 5))

    def test_mean_zero_and_nonzero_constraints(self):
        for values, degrees in [([1], [0]), ([1, 1], [0, 1]), ([0], [1])]:
            with self.assertRaises(ValueError):
                w.rayleigh_certificate([0], values, degrees=degrees)
        c = w.rayleigh_certificate([0], [0, 1], degrees=[0, 1])
        self.assertEqual(c['degrees'], [1])
        for m, degrees in [(True, [1]), (1, [0]), (0, [1, 1])]:
            with self.assertRaises(ValueError):
                w.rayleigh_certificate([0], [1]*len(degrees), m=m, degrees=degrees)

    def test_scope_formula_and_input_tampering(self):
        cert = self.simple()
        changes = [('scope', 'all_H1_functions_on_unit_S2'), ('rayleigh_quotient', '0'),
                   ('full_space_spectral_upper', '0'), ('azimuth_integral_factor', '2*pi'),
                   ('degrees', [1, 5]), ('coefficients', ['1', '1/36']),
                   ('basis_convention', 'arbitrary'), ('status', 'proved')]
        for key, value in changes:
            bad = deepcopy(cert)
            bad[key] = value
            self.assertFalse(w.verify(bad), key)
        self.assertFalse(w.verify(cert, expected_q=[0, 0, 1]))
        self.assertFalse(w.verify(cert, expected_m=0))
        self.assertFalse(w.verify(cert, expected_mean_zero=False))
        self.assertFalse(w.verify(cert, expected_coefficients=[1, F(-1, 40)]))
        self.assertFalse(w.verify(cert, expected_threshold=2))
        residual = w.residual_certificate(cert)
        for key, value in [('residual_coefficients', []), ('normalized_residual_squared', '0'),
                           ('scope', 'ground_spectral_lower'), ('removed_l0_coefficient', '1')]:
            bad = deepcopy(residual)
            bad[key] = value
            self.assertFalse(w.verify(bad), key)

    def test_float_proposals_have_only_exact_trial_authority(self):
        proposal = w.mass_normalized_proposal([0, 0, 2], m=1, modes=3)
        proposal['numerical_rayleigh'] = -10**50
        ladder = w.quantize_proposal(proposal, bits=(1, 8, 16), threshold=F(12, 5))
        best = ladder['certificate']
        self.assertTrue(w.verify(best))
        self.assertGreater(F(best['rayleigh_quotient']), 2)
        self.assertLess(F(best['rayleigh_quotient']), F(12, 5))
        for row in ladder['trace']:
            if 'certificate' in row:
                self.assertTrue(w.verify(row['certificate']))
                self.assertLessEqual(F(best['rayleigh_quotient']), F(row['certificate']['rayleigh_quotient']))
        proposal['coefficients_float'] = [0.0]*3
        with self.assertRaises(ValueError):
            w.quantize_proposal(proposal)

    def test_sparse_rejection_and_residual_growth_are_real(self):
        good = self.simple()
        sparse = w.sparsify(good, keep_counts=(1,))
        self.assertEqual(sparse['trace'][0]['status'], 'rejected_no_gain')
        self.assertEqual(sparse['certificate'], good)
        initial = w.rayleigh_certificate([0, 0, 2], [1], m=1, degrees=[1], threshold=F(12, 5))
        enriched = w.residual_enrich(initial, bits=(8, 16))
        self.assertEqual(enriched['added_degrees'], [3])
        self.assertEqual(enriched['status'], 'accepted_improvement')
        self.assertTrue(w.verify(enriched['certificate']))
        self.assertLess(F(enriched['certificate']['rayleigh_quotient']), F(12, 5))
        exact = w.rayleigh_certificate([0], [1], degrees=[1])
        self.assertEqual(w.residual_enrich(exact)['status'], 'support_did_not_grow')

    def test_two_trial_ritz_scope_and_linear_independence(self):
        first = w.rayleigh_certificate([0, 0, 2], [1], m=1, degrees=[1], threshold=F(12, 5))
        second = w.rayleigh_certificate([0, 0, 2], [1], m=1, degrees=[3], threshold=F(12, 5))
        result = w.two_trial_ritz(first, second, bits=20)
        self.assertTrue(w.verify(result))
        self.assertEqual([[F(x) for x in row] for row in result['mass_matrix']], [[F(4, 3), 0], [0, F(24, 7)]])
        self.assertEqual([[F(x) for x in row] for row in result['form_matrix']],
                         [[F(16, 5), F(32, 35)], [F(32, 35), F(1552, 35)]])
        self.assertGreater(F(result['ritz_lower']), F(238, 100))
        self.assertLess(F(result['actual_trial_upper']), F(12, 5))
        self.assertEqual(result['scope'], 'minimum_Rayleigh_quotient_in_the_span_of_the_two_stated_trials_only')
        for key, value in [('scope', 'full_sphere'), ('ritz_lower', '3')]:
            bad = deepcopy(result)
            bad[key] = value
            self.assertFalse(w.verify(bad))
        with self.assertRaises(ValueError):
            w.two_trial_ritz(first, first)

    def test_residual_line_search_is_scale_robust(self):
        quotients = []
        for scale_factor in [F(1), F(10**60), F(1, 10**100)]:
            with self.subTest(scale_factor=scale_factor):
                initial = self.simple(scale_factor)
                line = w.residual_line_search(initial, bits=16)
                self.assertEqual(line['status'], 'accepted_improvement')
                self.assertLess(F(line['certificate']['rayleigh_quotient']), F(initial['rayleigh_quotient']))
                self.assertTrue(w.verify(line['certificate']))
                quotients.append(F(line['certificate']['rayleigh_quotient']))
        self.assertEqual(len(set(quotients)), 1)
        exact = w.rayleigh_certificate([0], [1], degrees=[1])
        self.assertEqual(w.residual_line_search(exact)['status'], 'zero_residual_or_degree_budget')

    def test_enrichment_at_dimension_budget_retains_valid_trial(self):
        initial = w.rayleigh_certificate([0, 1], [1]*64, degrees=list(range(1, 65)))
        enriched = w.residual_enrich(initial)
        self.assertEqual(enriched['status'], 'support_dimension_budget')
        self.assertEqual(enriched['added_degrees'], [])
        self.assertEqual(enriched['certificate'], initial)
        self.assertTrue(w.verify(enriched['certificate']))

    def test_sector_selection_uses_verified_spectrum(self):
        source = json.loads((Path(__file__).parent/'baseline.json').read_text())['old_certificate']['evidence']
        selection = w.select_sector(source, F(12, 5))
        self.assertEqual(selection['selected_m'], 1)
        self.assertTrue(selection['spectral_existence_below_threshold'])
        self.assertTrue(w.verify(selection, expected_q=[0, 0, 2], expected_m=1,
                                 expected_mean_zero=True, expected_threshold=F(12, 5)))
        bad = deepcopy(source)
        bad['sectors'][0]['upper'] = '-100'
        with self.assertRaises(ValueError):
            w.select_sector(bad)

    def test_end_to_end_success_and_honest_not_found(self):
        source = json.loads((Path(__file__).parent/'baseline.json').read_text())['old_certificate']['evidence']
        no = w.find_counterexample([0, 0, 2], F(12, 5), full_certificate=source, modes=2, max_steps=0)
        self.assertEqual(no['status'], 'not_found')
        self.assertTrue(w.verify(no))
        yes = w.find_counterexample([0, 0, 2], F(12, 5), full_certificate=source, modes=2, max_steps=1)
        self.assertEqual(yes['status'], 'explicit_counterexample_verified')
        self.assertTrue(w.verify(yes, expected_q=[0, 0, 2], expected_threshold=F(12, 5)))
        self.assertLess(F(yes['certificate']['rayleigh_quotient']), F(12, 5))
        bad = deepcopy(yes)
        bad['explicit_real_harmonic_expansion'] = []
        self.assertFalse(w.verify(bad))
        with self.assertRaises(ValueError):
            w.find_counterexample([0, 0, 1], F(12, 5), full_certificate=source)


if __name__ == '__main__':
    unittest.main()
