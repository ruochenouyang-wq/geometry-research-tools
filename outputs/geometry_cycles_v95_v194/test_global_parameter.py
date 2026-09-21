"""Adversarial acceptance checks for cycle 5's global parameter certificates."""
import copy
from fractions import Fraction as F
import json
import unittest
from unittest.mock import patch
import global_parameter as G


class GlobalParameterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.quartic = {'q0': ['0'], 'direction': ['0'],
                       'penalty': ['0', '0', '-1', '0', '1'], 'domain': ['-1', '1']}
        cls.boundary = {'q0': ['0'], 'direction': ['1'], 'penalty': ['0'], 'domain': ['-1', '1']}
        cls.rough = G.branch_and_bound(cls.quartic, max_leaves=2, tolerance=F(1, 1000000))
        cls.cert = G.resume(cls.quartic, cls.rough, max_leaves=8, tolerance=F(1, 1000000))
        cls.edge = G.branch_and_bound(cls.boundary, max_leaves=4)

    def test_three_point_sampling_misses_global_minimum(self):
        # All three values are 2, while the true minimum is 7/4.
        self.assertEqual([2+s**4-s**2 for s in [-1, 0, 1]], [2, 2, 2])
        self.assertLess(F(self.cert['upper']), F(18, 10))
        self.assertLessEqual(F(self.cert['lower']), F(7, 4))
        self.assertGreaterEqual(F(self.cert['upper']), F(7, 4))

    def test_cheap_complete_domain_bound_without_solvers(self):
        obj = {'q0': ['0', '0', '2'], 'direction': ['0', '1'],
               'penalty': ['0', '0', '1/10'], 'domain': ['-2', '2']}
        with patch.object(G.projected, 'full_ground', side_effect=AssertionError('No spectral solve allowed')):
            with patch.object(G.E, 'maximize', side_effect=AssertionError('No extremum solve allowed')):
                cert = G.cheap_objective_enclosure(obj)
                self.assertEqual((cert['lower'], cert['upper']), ('0', '12/5'))
                self.assertEqual(cert['trial_function'], 'x1')
                self.assertTrue(G.verify_cheap_objective(cert, obj))

    def test_cheap_verifier_does_not_call_generator(self):
        cert = G.cheap_objective_enclosure(self.quartic)
        self.assertEqual((cert['lower'], cert['upper']), ('1', '2'))
        with patch.object(G, 'cheap_objective_enclosure', side_effect=AssertionError('Independent verifier')):
            self.assertTrue(G.verify_cheap_objective(cert, self.quartic))

    def test_cheap_constant_case_exact(self):
        obj = {'q0': ['3'], 'direction': ['0'], 'penalty': ['4'], 'domain': ['-1', '1']}
        cert = G.cheap_objective_enclosure(obj)
        self.assertEqual((cert['lower'], cert['upper']), ('9', '9'))
        self.assertTrue(G.verify_cheap_objective(cert))

    def test_cheap_negative_quadratic_selects_x3(self):
        obj = {'q0': ['0', '0', '-2'], 'direction': ['0'], 'penalty': ['0'], 'domain': ['-1', '1']}
        cert = G.cheap_objective_enclosure(obj)
        self.assertEqual(cert['trial_function'], 'x3')
        self.assertEqual(cert['upper'], '4/5')
        self.assertTrue(G.verify_cheap_objective(cert))

    def test_cheap_certificate_mutations_rejected(self):
        cert = G.cheap_objective_enclosure(self.quartic)
        for key, value in [('lower', '3'), ('upper', '0'), ('feasible_parameter', '1'),
                           ('trial_squared_L2_probability_norm', '1'),
                           ('function_space', 'axisymmetric_only')]:
            bad = copy.deepcopy(cert)
            bad[key] = value
            self.assertFalse(G.verify_cheap_objective(bad))
        other = copy.deepcopy(self.quartic)
        other['domain'] = ['-2', '2']
        self.assertFalse(G.verify_cheap_objective(cert, other))

    def test_cheap_nonsymmetric_parameter_interval(self):
        obj = {'q0': ['-1', '2', '-3'], 'direction': ['2', '-1', '1'],
               'penalty': ['1', '-2', '-3', '1'], 'domain': ['1', '2']}
        cert = G.cheap_objective_enclosure(obj)
        self.assertTrue(G.verify_cheap_objective(cert, obj))
        self.assertEqual(cert['feasible_parameter'], '3/2')
        self.assertEqual(cert['penalty_term_intervals'], [['1', '1'], ['-4', '-2'], ['-12', '-3'], ['1', '8']])

    def test_nonconvex_penalty_extrema(self):
        proof = G.penalty_enclosure(self.quartic['penalty'], ['-1', '1'])
        self.assertEqual(F(proof['lower']), F(-1, 4))
        self.assertEqual(F(proof['upper']), 0)

    def test_exact_boundary_minimum(self):
        self.assertEqual(self.edge['lower'], '1')
        self.assertEqual(self.edge['upper'], '1')
        self.assertEqual(self.edge['incumbent_parameter'], '-1')

    def test_symmetric_nonunique_minimizers_retained(self):
        boxes = G.minimizer_enclosure(self.cert)
        self.assertFalse(boxes['uniqueness_claimed'])
        # Rational interval conditions for +/-sqrt(1/2), avoiding float evidence.
        pos = any(F(b) > 0 and max(F(a), 0)**2 <= F(1, 2) <= F(b)**2
                  for a, b in boxes['all_minimizers_in'])
        neg = any(F(a) < 0 and min(F(b), 0)**2 <= F(1, 2) <= F(a)**2
                  for a, b in boxes['all_minimizers_in'])
        self.assertTrue(pos and neg)

    def test_missing_interval_rejected(self):
        bad = copy.deepcopy(self.cert)
        bad['cells'].pop(0)
        self.assertFalse(G.verify(bad))

    def test_overlapping_interval_rejected(self):
        bad = copy.deepcopy(self.cert)
        bad['cells'].append(copy.deepcopy(bad['cells'][0]))
        bad['max_leaves'] = 32
        self.assertFalse(G.verify(bad))

    def test_wrong_seed_objective_rejected(self):
        other = copy.deepcopy(self.quartic)
        other['q0'] = ['1']
        with self.assertRaises(ValueError):
            G.resume(other, self.rough, max_leaves=4)

    def test_wrong_seed_domain_rejected(self):
        other = copy.deepcopy(self.quartic)
        other['domain'] = ['-2', '2']
        with self.assertRaises(ValueError):
            G.resume(other, self.rough, max_leaves=4)

    def test_checkpoint_monotonicity(self):
        self.assertGreaterEqual(F(self.cert['lower']), F(self.rough['lower']))
        self.assertLessEqual(F(self.cert['upper']), F(self.rough['upper']))

    def test_forged_incumbent_rejected(self):
        bad = copy.deepcopy(self.cert)
        bad['upper'] = '0'
        self.assertFalse(G.verify(bad))

    def test_forged_point_parameter_rejected(self):
        bad = copy.deepcopy(self.cert)
        key = next(iter(bad['points']))
        bad['points'][key]['parameter'] = '100'
        self.assertFalse(G.verify(bad))

    def test_wrong_function_space_rejected(self):
        bad = copy.deepcopy(self.edge)
        bad['function_space'] = 'axisymmetric_only'
        self.assertFalse(G.verify(bad))

    def test_boolean_not_exact_number(self):
        bad = copy.deepcopy(self.quartic)
        bad['penalty'] = [True]
        with self.assertRaises(ValueError):
            G.canonical_objective(bad)

    def test_point_outside_domain_rejected(self):
        with self.assertRaises(ValueError):
            G.feasible_point(self.quartic, 5)

    def test_stale_cached_spectrum_rejected(self):
        cache = {'0': G.feasible_point(self.quartic, 0)}
        changed = copy.deepcopy(self.quartic)
        changed['q0'] = ['1']
        with self.assertRaises(ValueError):
            G.cached_point(changed, 0, cache)

    def test_json_roundtrip(self):
        self.assertTrue(G.verify(json.loads(json.dumps(self.cert)), self.quartic))

    def test_all_threshold_outcomes(self):
        lo, hi = F(self.cert['lower']), F(self.cert['upper'])
        expected = ['proved', 'refuted_by_feasible_parameter_spectral_existence', 'undetermined']
        for threshold, status in zip([lo, hi+1, (lo+hi)/2], expected):
            cert = G.threshold_certificate(self.cert, threshold)
            self.assertEqual(cert['status'], status)
            self.assertTrue(G.verify_threshold(cert, self.quartic, threshold))
            self.assertFalse(G.verify_threshold(cert, self.quartic, threshold+1))

    def test_chord_is_stronger_than_lipschitz(self):
        cache = {str(s): G.feasible_point(self.boundary, s) for s in [-1, 0, 1]}
        lip = G.lipschitz_cell(self.boundary, [-1, 1], list(cache.values()))
        chord = G.concavity_cell(self.boundary, [-1, 1], cache['-1'], cache['1'])
        self.assertGreaterEqual(F(chord['lower']), F(lip['lower']))

    def test_actual_spherical_family_keeps_open_budget(self):
        objective = {'q0': ['0', '0', '2'], 'direction': ['0', '1'],
                     'penalty': ['0', '0', '1/10'], 'domain': ['-2', '2']}
        certificate = G.branch_and_bound(objective, max_leaves=1, tolerance=F(1, 10**8))
        self.assertTrue(G.verify(certificate, objective))
        self.assertEqual(certificate['status'], 'certified_bound_open_gap')
        self.assertGreater(F(certificate['gap']), F(certificate['tolerance']))
        self.assertEqual(certificate['minimizer_intervals'], [['-2', '2']])

    def test_penalty_changes_preferred_parameter_region(self):
        objective = {'q0': ['0'], 'direction': ['0', '1'],
                     'penalty': ['0', '0', '1/100'], 'domain': ['-2', '2']}
        edge = G.branch_and_bound(objective, max_leaves=8)
        self.assertTrue(G.verify(edge))
        self.assertEqual(edge['incumbent_parameter'], '-2')
        self.assertLess(F(edge['upper']), 2)
        objective['penalty'][-1] = '1/10'
        middle = G.branch_and_bound(objective, max_leaves=16)
        self.assertTrue(G.verify(middle))
        self.assertEqual(middle['incumbent_parameter'], '0')
        self.assertTrue(all(abs(F(a)) <= F(1, 16) and abs(F(b)) <= F(1, 16)
                            for a, b in middle['minimizer_intervals']))


if __name__ == '__main__':
    unittest.main(verbosity=2)
