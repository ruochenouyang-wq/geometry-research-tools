"""Analytic calibrations, adversarial scope binding, and full residual replay."""
import json
from copy import deepcopy
from pathlib import Path
import unittest
from unittest.mock import patch

import enriched_trial as e

F = e.F


class EnrichedTrialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parent/'certificates'/'singular_n4.json'
        cls.source = json.loads(path.read_text())
        p = e.propose_trial(e.DEFAULT_EXPONENTS[:10])
        cls.cert = e.trial_certificate(p['powers'], p['coefficients'], cls.source)

    def test_weighted_moment_from_difference_of_plain_moments(self):
        for r in (F(-1, 2), F(-1, 4), F(0), F(7, 4), F(13, 2)):
            self.assertEqual(e.weighted_moment(r), 2/(r+1)-2/(r+3))
        for r in (-1, -2):
            with self.assertRaises(ValueError):
                e.weighted_moment(r)

    def test_domain_does_not_admit_delta_or_non_L2_action(self):
        for s in ('-1', '1/4', '1', '3/2'):
            with self.assertRaises(ValueError):
                e.trial_matrices([0, s])
        self.assertEqual(e.exponents([0, '1500001/1000000']), (F(0), F(1500001, 1000000)))
        self.assertEqual(e.action_terms(0), {F(0): F(2), F(-1, 4): F(-1)})

    def test_exact_input_and_power_budget(self):
        for ss in ([False], [0.0], [0, 2, 2], [2, 0], list(range(2, 19)), [65]):
            with self.assertRaises(ValueError):
                e.trial_matrices(ss)
        with self.assertRaises(ValueError):
            e.trial_statistics([0], [1.0])
        with self.assertRaises(ValueError):
            e.trial_statistics([0], [0])

    def test_constant_radial_trial_has_analytic_rayleigh_and_variance(self):
        stats = e.trial_statistics([0], [1])
        self.assertEqual(F(stats['mass']), F(4, 3))
        self.assertEqual(F(stats['rayleigh']), F(6, 11))
        self.assertEqual(F(stats['residual_squared']), F(172, 605))

    def test_first_fractional_coefficient_cancels_actual_singularity(self):
        terms = {}
        for s, a in ((F(0), F(1)), (F(7, 4), F(-16, 21))):
            for r, c in e.action_terms(s).items():
                terms[r] = terms.get(r, F(0))+a*c
        terms = {r: c for r, c in terms.items() if c}
        self.assertEqual(terms, {F(0): F(2), F(3, 2): F(16, 21), F(7, 4): F(-55, 7)})

    def test_operator_form_matches_independent_weak_gradient_integral(self):
        ss = e.DEFAULT_EXPONENTS[:6]
        M, H, R = e.trial_matrices(ss)
        for i, s in enumerate(ss):
            for j, t in enumerate(ss):
                r = s+t
                leading = (2*s*t*(1/(r-1)-2/(r+1)+1/(r+3)) if s*t else F(0))
                kinetic = leading-r*e.weighted_moment(r)+2/(r+3)+2/(r+1)
                self.assertEqual(H[i][j], kinetic-e.weighted_moment(r-F(1, 4)))
                self.assertEqual(M[i][j], e.weighted_moment(r))
                self.assertEqual(H[i][j], H[j][i])
                self.assertEqual(R[i][j], R[j][i])

    def test_complete_residual_integral_matches_matrix_variance(self):
        ss, v = (F(0), F(7, 4), F(2)), (F(1), F(-16, 21), F(1, 5))
        stats = e.trial_statistics(ss, v)
        mu = F(stats['rayleigh'])
        residual_terms = {}
        for s, a in zip(ss, v):
            for r, c in e.action_terms(s).items():
                residual_terms[r] = residual_terms.get(r, F(0))+a*c
            residual_terms[s] = residual_terms.get(s, F(0))-mu*a
        integral = sum((a*b*(2/(r+t+1)-2/(r+t+3))
                        for r, a in residual_terms.items()
                        for t, b in residual_terms.items()), F(0))
        self.assertEqual(integral/F(stats['mass']), F(stats['residual_squared']))

    def test_trial_scale_does_not_change_rayleigh_or_full_residual(self):
        a = e.trial_statistics([0, '7/4'], [1, '-16/21'])
        b = e.trial_statistics([0, '7/4'], [-21, 16])
        self.assertEqual(a['rayleigh'], b['rayleigh'])
        self.assertEqual(a['residual_squared'], b['residual_squared'])

    def test_second_bound_is_single_real_sector_and_both_radial_parities(self):
        p = e.second_eigenvalue_bound()
        eta = F(p['form_bound']['eta_upper'])
        self.assertGreaterEqual(eta*eta, F(2, 9))
        self.assertLess(eta, 1)
        self.assertEqual(F(p['lower']), (1-eta)*6-F(4, 3)-eta)
        self.assertEqual(p['real_component'], 'cos_phi_only')
        self.assertEqual(p['radial_parities'], 'both_even_and_odd')

    def test_temple_rejects_high_rayleigh_even_if_trial_domain_valid(self):
        with self.assertRaisesRegex(ValueError, 'Temple'):
            e.trial_certificate([2], [1], self.source)

    def test_ten_terms_reaches_goal_with_exact_residual_and_other_sectors(self):
        c = self.cert
        self.assertTrue(e.verify(c, self.source, '1/100000000'))
        self.assertEqual(c['status'], 'target_met')
        self.assertLess(F(c['exact_width']), F(1, 10**9))
        self.assertGreater(F(c['exact_width']), 0)
        mu = F(c['statistics']['rayleigh'])
        sigma2 = F(c['statistics']['residual_squared'])
        beta = F(c['second_eigenvalue_proof']['lower'])
        self.assertEqual(F(c['temple_lower']), mu-sigma2/(beta-mu))
        self.assertEqual(c['lower'], c['temple_lower'])
        self.assertTrue(all(F(row['lower']) > F(c['upper']) for row in c['other_sector_lowers']))
        self.assertGreater(F(c['angular_tail_lower']), F(c['upper']))

    def test_sixteen_terms_delivers_more_than_sixteen_decimal_digits(self):
        p = e.propose_trial(e.DEFAULT_EXPONENTS)
        c = e.trial_certificate(p['powers'], p['coefficients'], self.source, '1/10000000000000000')
        self.assertTrue(e.verify(c))
        self.assertLess(F(c['exact_width']), F(1, 10**17))
        self.assertLessEqual(F(self.cert['lower']), F(c['lower']))
        self.assertGreaterEqual(F(self.cert['upper']), F(c['upper']))

    def test_source_function_mean_constraint_and_m1_presence_are_bound(self):
        for key, value in (('function', dict(e.FUNCTION, amplitude='1')),
                           ('mean_zero', False), ('scope', 'axisymmetric_only')):
            bad = deepcopy(self.source)
            bad[key] = value
            with self.assertRaises(ValueError):
                e.trial_certificate([0], [1], bad)
        wrong = deepcopy(self.source)
        wrong['sectors'].pop()
        with self.assertRaises(ValueError):
            e.trial_certificate([0], [1], wrong)

    def test_full_record_mutations_rejected(self):
        mutations = [
            (('statistics', 'residual_squared'), '0'),
            (('second_eigenvalue_proof', 'lower'), '6'),
            (('second_eigenvalue_proof', 'real_component'), 'cos_and_sin'),
            (('trial', 'radial_parity'), 'all'),
            (('function', 'exponent'), '-1/3'),
            (('geometry',), 'unit_circle'),
            (('source_digest',), 'wrong'),
            (('scope',), 'axisymmetric_only'),
            (('lower',), self.cert['upper']),
            (('full_strong_residual_not_projected_residual',), False)]
        for keys, value in mutations:
            with self.subTest(keys=keys):
                bad = deepcopy(self.cert)
                target = bad
                for key in keys[:-1]:
                    target = target[key]
                target[keys[-1]] = value
                self.assertFalse(e.verify(bad))
        bad = deepcopy(self.cert)
        bad['matrices']['R'][0][0] = '0'
        self.assertFalse(e.verify(bad))

    def test_json_certificate_replay_never_calls_candidate_search(self):
        record = json.loads(json.dumps(self.cert))
        with patch.object(e, 'propose_trial', side_effect=AssertionError('search in verifier')):
            with patch.object(e.direct, 'certify_sector', side_effect=AssertionError('source search')):
                self.assertTrue(e.verify(record, self.source))
        self.assertFalse(e.verify(record, expected_tolerance='1/1000000000'))

    def test_small_budget_preserves_failure_and_uses_source_if_trial_is_worse(self):
        result = e.enrich(self.source, max_terms=1)
        c = result['certificate']
        self.assertEqual(result['status'], 'certified_open')
        self.assertEqual(c['lower'], self.source['lower'])
        self.assertEqual(c['upper'], self.source['upper'])
        self.assertGreater(F(c['statistics']['rayleigh']), F(self.source['upper']))
        self.assertEqual([a['terms'] for a in result['attempts']], [1])

    def test_search_trace_retains_coarse_failures_and_stops_on_target(self):
        result = e.enrich(self.source)
        self.assertEqual([a['terms'] for a in result['attempts']], [1, 3, 6, 10])
        self.assertTrue(all(a['status'] == 'certified_open' for a in result['attempts'][:-1]))
        self.assertEqual(result['status'], 'target_met')
        self.assertTrue(e.verify(result['certificate']))


if __name__ == '__main__':
    unittest.main()
