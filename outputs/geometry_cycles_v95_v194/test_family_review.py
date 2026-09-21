"""Independent cycle-4 review: analytical calibrations and hostile replay.

Run: python3 -B -m unittest -v test_family_review
Expected values below are derived from sphere moments, not from solver output.
"""
import copy
import json
from pathlib import Path
import unittest
from unittest.mock import patch
from fractions import Fraction as F

import parameter_family as P
import witness as W

HERE = Path(__file__).resolve().parent


def family(q0, directions, box, mean_zero=True):
    return P.validate_family(q0, directions, box, mean_zero)


class FamilyIndependentReview(unittest.TestCase):
    def test_constant_shift_exact_spectrum_and_constraint_floor(self):
        # Constants commute with -Delta: constrained lambda=2+c, full lambda=c.
        for zero_mean, floor in [(True, 2), (False, 0)]:
            fam = family([3], [[2]], [[-2, 1]], zero_mean)
            cheap = P.cheap_family_enclosure(fam)
            anchor = P.anchor_cell(fam, modes=2)
            corners = P.cell_lower_bound(fam, modes=2)
            self.assertEqual(F(cheap['lower']), floor-1)
            self.assertEqual(F(cheap['upper']), floor+5)
            self.assertEqual(F(anchor['lower']), floor-1)
            self.assertEqual(F(anchor['upper']), floor+5)
            self.assertEqual(F(corners['lower']), floor-1)
            self.assertEqual(F(corners['minimum_upper']), floor-1)
            self.assertTrue(P.verify(cheap, expected_family=fam))

    def test_two_parameter_constant_shift_and_transfer_of_physical_trial(self):
        fam = family([5], [[2], [-3]], [[-1, 2], [0, 1]])
        # This physical P1 trial was originally certified for a different q.
        trial = W.rayleigh_certificate([7, 0, 1], [1], m=0, degrees=[1])
        line = P.affine_majorant(fam, trial)
        self.assertEqual(F(line['intercept']), 7)
        self.assertEqual(list(map(F, line['slopes'])), [2, -3])
        bound = P.cell_lower_bound(fam, modes=2)
        self.assertEqual(F(bound['lower']), 2)
        self.assertEqual(F(bound['minimum_upper']), 2)
        self.assertEqual(F(P.envelope_box(fam, [line])['maximum_upper']), 11)

    def test_vertex_maximum_fails_as_upper_for_a_concave_family(self):
        fam = family([0], [[0, 1]], [[-1, 1]])
        vertices = P.cell_lower_bound(fam, modes=4)
        # For q=t, u=P1^1-P2^1/20 has exact quotient <2;
        # masses are 4/3,12/5 and cross moment int t P1^1 P2^1=4/5.
        mass = F(4, 3)+F(12, 5)/400
        energy = F(8, 3)+6*F(12, 5)/400-F(2, 20)*F(4, 5)
        self.assertLess(energy/mass, 2)
        self.assertTrue(all(F(row['spectrum']['upper']) < 2 for row in vertices['vertices']))
        center = P.anchor_cell(fam, cell=[[-F(1, 1000), F(1, 1000)]], modes=2)
        self.assertEqual(F(center['anchor_proof']['lower']), 2)
        self.assertNotIn('maximum_upper', vertices)
        bad = copy.deepcopy(vertices)
        bad['maximum_upper'] = max(row['spectrum']['upper'] for row in bad['vertices'])
        self.assertFalse(P.verify(bad))

    def test_one_dimensional_envelope_requires_interior_intersection(self):
        fam = family([0], [[-F(2, 5), 0, 1]], [[-1, 1]])
        lines = [P.affine_majorant(fam, W.rayleigh_certificate([0], [1], m=m, degrees=[1]))
                 for m in (0, 1)]
        # Coordinate squares have E[t²|x3²]=3/5 and E[t²|x1²]=1/5.
        self.assertEqual([list(map(F, row['slopes'])) for row in lines], [[F(1, 5)], [-F(1, 5)]])
        cert = P.envelope_1d(fam, lines+lines[:1])
        self.assertEqual(F(cert['maximum_upper']), 2)
        self.assertEqual(F(cert['envelope_maximizer']), 0)
        endpoints = [r for r in cert['candidates'] if abs(F(r['parameter'])) == 1]
        self.assertTrue(all(F(r['envelope_value']) == F(9, 5) for r in endpoints))
        self.assertTrue(P.verify(cert))

    def test_two_dimensional_min_max_order_counterexample(self):
        fam = family([0], [[0], [-F(2, 5), 0, 1]], [[-1, 1], [-1, 1]])
        lines = [P.affine_majorant(fam, W.rayleigh_certificate([0], [1], m=m, degrees=[1]))
                 for m in (0, 1)]
        # min(2+r/5,2-r/5)=2-|r|/5. Vertex maximum=9/5 is UNSAFE:
        # at r=0 the actual constrained spectral minimum is exactly 2.
        cert = P.envelope_box(fam, lines)
        self.assertEqual(F(cert['maximum_upper']), F(11, 5))
        self.assertGreaterEqual(F(cert['maximum_upper']), 2)
        self.assertIn('Not generally', cert['limitation'])
        bad = copy.deepcopy(cert)
        bad['maximum_upper'] = '9/5'
        self.assertFalse(P.verify(bad))

    def test_three_dimensional_nonuniform_cover_and_exact_tiny_gap(self):
        fam = family([0], [[0], [0], [0]], [[0, 1], [0, 1], [0, 1]])
        cells = [[[0, F(1, 3)], [0, 1], [0, 1]],
                 [[F(1, 3), 1], [0, F(2, 5)], [0, 1]],
                 [[F(1, 3), 1], [F(2, 5), 1], [0, 1]]]
        good = P.partition_certificate(fam, cells)
        self.assertTrue(P.verify(good))
        self.assertEqual(sum(map(F, good['cell_volumes'])), 1)
        gap = copy.deepcopy(cells)
        gap[0][0][1] -= F(1, 10**15)
        with self.assertRaises(ValueError):
            P.partition_certificate(fam, gap)
        with self.assertRaises(ValueError):
            P.partition_certificate(fam, cells+[cells[0]])

    def test_cache_reuse_depends_on_actual_q_and_constraint(self):
        linear = family([0], [[0, 1]], [[-1, 1]])
        quadratic = family([0], [[0, 0, 1]], [[-1, 1]])
        cache = {}
        a = P.anchor_cell(linear, anchor=[0], cache=cache, modes=2)
        # Different families with the SAME actual operator at zero may reuse.
        with patch.object(P.projected, 'full_ground', side_effect=AssertionError('cache miss')):
            b = P.anchor_cell(quadratic, anchor=[0], cache=cache, modes=2)
        self.assertEqual(a['anchor_proof'], b['anchor_proof'])
        P.anchor_cell(linear, anchor=[1], cache=cache, modes=2)
        with self.assertRaises(ValueError):
            P.anchor_cell(quadratic, anchor=[1], cache=cache, modes=2)
        unconstrained = family([0], [[0, 1]], [[-1, 1]], False)
        with self.assertRaises(ValueError):
            P.anchor_cell(unconstrained, anchor=[0], cache=cache, modes=2)

    def test_direction_norm_uses_entire_interval_and_retains_cancellation(self):
        fam = family([0], [[0, 1, 0, -1]], [[-1, 1]])
        cert = P.direction_norms(fam)
        norm = F(cert['norm_upper'][0])
        # sup|t-t³|=2/(3 sqrt(3)); 9/25 < it < 2/5, far below coefficient l1=2.
        self.assertGreater(norm, F(9, 25))
        self.assertLess(norm, F(2, 5))
        bad = copy.deepcopy(cert)
        bad['proofs'][0]['positive']['interval'] = ['0', '1']
        self.assertFalse(P.verify(bad))

    def test_universal_threshold_equality_and_strict_explicit_refutation(self):
        fam = family([0], [[1]], [[-1, 1]])
        lower = P.cell_lower_bound(fam, modes=2)
        yes = P.universal_inequality(fam, 1, lower_proof=lower, modes=2)
        no = P.universal_inequality(fam, F(1001, 1000), lower_proof=lower, modes=2)
        self.assertEqual(yes['status'], 'proved')
        self.assertEqual(no['status'], 'refuted_with_explicit_parameter_and_function')
        c = no['counterexample']
        self.assertEqual(list(map(F, c['parameter'])), [-1])
        self.assertEqual(F(c['trial']['rayleigh_quotient']), 1)
        self.assertTrue(W.verify(c['trial'], expected_q=[-1], expected_mean_zero=True))
        self.assertFalse(P.verify(no, expected_threshold=1))
        bad = copy.deepcopy(no)
        bad['counterexample']['parameter'] = ['1']
        self.assertFalse(P.verify(bad))

    def test_universal_rejects_wrong_lower_evidence_kind_or_subbox(self):
        fam = family([0], [[1]], [[-1, 1]])
        wrong = P.anchor_cell(fam, modes=2)
        sub = P.cell_lower_bound(fam, cell=[[0, 1]], modes=2)
        for proof in (wrong, sub):
            with self.assertRaises(ValueError):
                P.universal_inequality(fam, 10, lower_proof=proof, modes=2)

    def test_closed_parameter_domain_cannot_be_relabelled_open(self):
        fam = family([0], [[1]], [[-1, 1]])
        self.assertEqual(P.potential_at(fam, [-1]), [-1])
        self.assertEqual(P.potential_at(fam, [1]), [1])
        bad = copy.deepcopy(fam)
        bad['parameters_quantifier'] = 'every_real_point_in_open_box'
        self.assertFalse(P.verify(bad))
        with self.assertRaises(ValueError):
            family([0], [[1]], [[1, 1]])

    def test_adaptive_supremum_flat_case_and_honest_open_gap(self):
        flat = family([0], [[0]], [[-1, 1]])
        exact = P.adaptive_family(flat, max_cells=1, tolerance=F(1, 10**8), modes=2)
        self.assertEqual((F(exact['maximum_lower']), F(exact['maximum_upper'])), (2, 2))
        self.assertEqual(exact['status'], 'target_met')
        slope = family([0], [[1]], [[-1, 1]])
        unresolved = P.adaptive_family(slope, max_cells=1, tolerance=F(1, 10**8), modes=2)
        self.assertEqual(F(unresolved['maximum_upper']), 3)
        self.assertEqual(unresolved['status'], 'budget_open')
        bad = copy.deepcopy(unresolved)
        bad['status'] = 'target_met'
        self.assertFalse(P.verify(bad))

    def test_cheap_enclosure_is_search_free_and_input_bound(self):
        fam = family([0, 0, 2], [[0, 1]], [[-1, 1]])
        with (patch.object(P.projected, 'full_ground', side_effect=AssertionError('spectral solve')),
              patch.object(P.exact_extrema, 'maximize', side_effect=AssertionError('extrema search'))):
            cert = P.cheap_family_enclosure(fam)
            self.assertTrue(P.verify(cert, expected_family=fam))
        self.assertEqual((F(cert['lower']), F(cert['upper'])), (1, F(16, 5)))
        self.assertEqual(cert['fixed_trial'], 't')
        bad = copy.deepcopy(cert)
        bad['family']['mean_zero'] = False
        self.assertFalse(P.verify(bad))

    def test_replay_all_saved_certificates_without_spectral_or_extrema_search(self):
        paths = sorted((HERE/'round04'/'certificates').glob('*.json'))
        self.assertGreaterEqual(len(paths), 19)
        with (patch.object(P.projected, 'full_ground', side_effect=AssertionError('spectral search')),
              patch.object(P.exact_extrema, 'maximize', side_effect=AssertionError('extrema search'))):
            for path in paths:
                with self.subTest(path=path.name):
                    cert = json.loads(path.read_text())
                    self.assertTrue(P.verify(cert), path.name)
        stages = json.loads((HERE/'round04'/'STAGES.json').read_text())
        self.assertEqual([s['version'] for s in stages], list(range(125, 135)))
        for stage in stages:
            self.assertTrue(callable(getattr(P, stage['callable'].split('.')[-1])))
            for name in stage['evidence']:
                self.assertTrue((HERE/name).is_file(), name)


if __name__ == '__main__':
    unittest.main()
