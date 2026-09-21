"""Continuum scope, whole-domain coverage, concavity direction, and tamper tests."""
import copy
import unittest
from unittest.mock import patch
from common import F, projected
import parameter_family as p
import witness


class FamilyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.family = p.validate_family([0, 0, 2], [[0, 1]], [[-1, 1]])
        cls.cache = {}
        cls.norms = p.direction_norms(cls.family)
        cls.anchor = p.anchor_cell(cls.family, cache=cls.cache, modes=4, bits=24, norms=cls.norms)
        cls.vertices = p.cell_lower_bound(cls.family, cache=cls.cache, modes=4, bits=24)
        cls.trial = witness.rayleigh_certificate([0, 0, 2], [1], m=1)
        cls.majorant = p.affine_majorant(cls.family, cls.trial)

    def test_v125_cheap_spectral_enclosure(self):
        c = p.cheap_family_enclosure(self.family)
        self.assertTrue(p.verify(c))
        self.assertEqual((c['lower'], c['upper']), ('1', '16/5'))
        for vertex in self.vertices['vertices']:
            self.assertLessEqual(F(c['lower']), F(vertex['spectrum']['lower']))
            self.assertGreaterEqual(F(c['upper']), F(vertex['spectrum']['upper']))
        c['upper'] = '2'
        self.assertFalse(p.verify(c))

    def test_v125_cheap_unconstrained_and_exact_constants(self):
        for mean_zero, expected in [(True, F(5)), (False, F(3))]:
            family = p.validate_family([3], [[0]], [[-1, 1]], mean_zero=mean_zero)
            c = p.cheap_family_enclosure(family)
            self.assertEqual(F(c['lower']), expected)
            self.assertEqual(F(c['upper']), expected)
            self.assertTrue(p.verify(c))
        f = p.validate_family([0], [[1]], [[-2, 3]], mean_zero=False)
        c = p.cheap_family_enclosure(f)
        self.assertEqual((c['lower'], c['upper']), ('-2', '3'))
        self.assertEqual(c['rayleigh_maximizing_parameter'], ['3'])

    def test_v125_reject_float_dimension_scope(self):
        for args in [([0.0], [[0, 1]], [[-1, 1]]), ([0], [[0, 1]], [[1, -1]]),
                     ([0], [[0, 1]]*4, [[-1, 1]]*4), ([0]*8, [[1]], [[-1, 1]])]:
            with self.assertRaises(ValueError):
                p.validate_family(*args)
        c = copy.deepcopy(self.family)
        c['function_space'] = 'all_H1_functions_on_unit_S2'
        self.assertFalse(p.verify(c))

    def test_v125_exact_affine_evaluation(self):
        self.assertEqual(p.potential_at(self.family, [F(1, 3)]), [0, F(1, 3), 2])
        with self.assertRaises(ValueError):
            p.potential_at(self.family, [2])

    def test_v126_norm_and_tamper(self):
        self.assertTrue(p.verify(self.norms))
        self.assertEqual(self.norms['norm_upper'], ['1'])
        c = copy.deepcopy(self.norms)
        c['proofs'][0]['positive']['interval'] = ['0', '1']
        self.assertFalse(p.verify(c))
        c = copy.deepcopy(self.norms)
        c['norm_upper'][0] = '0'
        self.assertFalse(p.verify(c))

    def test_v127_lipschitz_radius(self):
        self.assertTrue(p.verify(self.anchor))
        self.assertEqual(self.anchor['lipschitz_radius'], '1')
        c = copy.deepcopy(self.anchor)
        c['lower'] = c['anchor_proof']['lower']
        self.assertFalse(p.verify(c))

    def test_v127_wrong_cached_family_rejected(self):
        other = p.validate_family([0], [[0, 1]], [[-1, 1]])
        with self.assertRaises(ValueError):
            p.anchor_cell(other, cache=self.cache, modes=4, bits=24)

    def test_v128_coverage_gap_overlap_rejected(self):
        self.assertTrue(p.verify(p.partition_certificate(self.family, [[[-1, 0]], [[0, 1]]])))
        for cells in [[[[-1, 0]]], [[[-1, F(1, 2)]], [[0, 1]]], [[[0, 1]], [[0, 1]]]]:
            with self.assertRaises(ValueError):
                p.partition_certificate(self.family, cells)

    def test_v128_two_dimensional_coverage(self):
        f = p.validate_family([0], [[0, 1], [0, 0, 1]], [[-1, 1], [-1, 1]])
        cells = [[[a, b], [c, d]] for a, b in [(-1, 0), (0, 1)] for c, d in [(-1, 0), (0, 1)]]
        self.assertTrue(p.verify(p.partition_certificate(f, cells)))

    def test_v129_concavity_gain_and_missing_vertex(self):
        self.assertTrue(p.verify(self.vertices))
        self.assertGreater(F(self.vertices['lower']), F(self.anchor['lower']))
        c = copy.deepcopy(self.vertices)
        c['vertices'].pop()
        self.assertFalse(p.verify(c))

    def test_v129_vertex_max_is_not_global_upper(self):
        center_lo = F(self.anchor['anchor_proof']['lower'])
        endpoint_hi = max(F(row['spectrum']['upper']) for row in self.vertices['vertices'])
        self.assertGreater(center_lo, endpoint_hi)
        self.assertNotIn('maximum_upper', self.vertices)

    def test_v129_cache_replay_avoids_search(self):
        with patch.object(projected, 'full_ground', side_effect=AssertionError('cache missed')):
            cached = p.cell_lower_bound(self.family, cache=self.cache, modes=4, bits=24)
        self.assertEqual(cached, self.vertices)

    def test_v130_affine_matches_direct_rayleigh(self):
        self.assertTrue(p.verify(self.majorant))
        for point in [-1, F(1, 7), 1]:
            q = p.potential_at(self.family, [point])
            direct = witness.rayleigh_certificate(q, [1], m=1)
            self.assertEqual(p._affine_at(self.majorant, [point]), F(direct['rayleigh_quotient']))
        c = copy.deepcopy(self.majorant)
        c['slopes'][0] = '1'
        self.assertFalse(p.verify(c))

    def test_v130_no_constant_trial_under_mean_zero(self):
        bad = witness.rayleigh_certificate([0], [1], mean_zero=False)
        with self.assertRaises(ValueError):
            p.affine_majorant(self.family, bad)

    def test_v131_intersection_and_parallel_lines(self):
        plus = p.affine_majorant(self.family, witness.rayleigh_certificate([0], [1, F(1, 8)], m=1))
        minus = p.affine_majorant(self.family, witness.rayleigh_certificate([0], [1, -F(1, 8)], m=1))
        envelope = p.envelope_1d(self.family, [plus, minus, plus])
        self.assertTrue(p.verify(envelope))
        self.assertIn('0', [x['parameter'] for x in envelope['candidates']])
        self.assertEqual(envelope['envelope_maximizer'], '0')
        c = copy.deepcopy(envelope)
        c['candidates'] = [r for r in c['candidates'] if r['parameter'] != '0']
        self.assertFalse(p.verify(c))

    def test_v132_higher_dimension_safe(self):
        f = p.validate_family([0, 0, 2], [[0, 1], [1]], [[-1, 1], [-1, 1]])
        major = p.affine_majorant(f, self.trial)
        envelope = p.envelope_box(f, [major])
        self.assertTrue(p.verify(envelope))
        self.assertEqual(F(envelope['maximum_upper']), F(17, 5))
        with self.assertRaises(ValueError):
            p.envelope_1d(f, [major])

    def test_v133_honest_budget_and_coverage(self):
        result = p.adaptive_family(self.family, tolerance=F(1, 10**10), max_cells=2,
                                  modes=4, bits=24, cache=self.cache)
        self.assertTrue(p.verify(result))
        self.assertEqual(result['status'], 'budget_open')
        self.assertEqual(len(result['cells']), 2)
        c = copy.deepcopy(result)
        c['status'] = 'target_met'
        self.assertFalse(p.verify(c))

    def test_v134_prove_refute_undetermined(self):
        proof = p.universal_inequality(self.family, F(23, 10), lower_proof=self.vertices, modes=4)
        self.assertEqual(proof['status'], 'proved')
        self.assertTrue(p.verify(proof, expected_threshold=F(23, 10)))
        refutation = p.universal_inequality(self.family, F(12, 5), lower_proof=self.vertices, modes=4)
        self.assertEqual(refutation['status'], 'refuted_with_explicit_parameter_and_function')
        self.assertTrue(p.verify(refutation))
        unresolved = p.universal_inequality(self.family, F(self.vertices['lower'])+F(1, 10**20),
                                         lower_proof=self.vertices, modes=4)
        self.assertEqual(unresolved['status'], 'undetermined')
        self.assertTrue(p.verify(unresolved))
        c = copy.deepcopy(refutation)
        c['counterexample']['parameter'] = ['0']
        self.assertFalse(p.verify(c))

    def test_v134_input_threshold_binding(self):
        result = p.universal_inequality(self.family, 2, lower_proof=self.vertices)
        self.assertFalse(p.verify(result, expected_threshold=3))
        other = p.validate_family([0, 0, 3], [[0, 1]], [[-1, 1]])
        self.assertFalse(p.verify(result, expected_family=other))


if __name__ == '__main__':
    unittest.main()
