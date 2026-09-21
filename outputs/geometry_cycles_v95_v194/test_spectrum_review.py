"""Independent review of ordering, absent modes, strict count and LT scaling."""
import copy
from fractions import Fraction as F
from pathlib import Path
import json
import unittest

import ordered_spectrum as s
from common import projected


def exact_levels(k, shift=0):
    answer = []
    degree = 1
    while len(answer) < k:
        answer += [F(degree*(degree+1)+shift)]*(2*degree+1)
        degree += 1
    return answer[:k]


def make_group(q, m, count):
    return [projected.certify_sector(q, m=m, k=j, modes=max(4, count), bits=40)
            for j in range(1, count+1)]


class OrderedSpectrumIndependentReview(unittest.TestCase):
    def test_empty_evidence_preserves_entire_exact_spectrum(self):
        assembly = s.assemble_sectors([-6], [])
        cert = s.ordered_bounds(s.normalize_request([-6], k=24), assembly)
        levels = exact_levels(24, -6)
        self.assertEqual([F(row['lower']) for row in cert['ordered_intervals']], levels)
        self.assertEqual([F(row['upper']) for row in cert['ordered_intervals']], levels)
        self.assertTrue(cert['uncomputed_modes_may_touch_first_k'])
        self.assertTrue(s.verify(cert))

    def test_missing_radial_modes_keep_correct_indices(self):
        groups = [make_group([0], m, 1) for m in range(3)]
        cert = s.ordered_bounds(s.normalize_request([0], k=15), s.assemble_sectors([0], groups))
        self.assertEqual([F(x['lower']) for x in cert['ordered_intervals']], exact_levels(15))
        placeholders = [row for row in cert['slots'] if row['source'] == 'uncomputed_minmax_placeholder']
        self.assertTrue(any(row['m'] == 0 and row['radial_index'] == 3 for row in placeholders))
        self.assertEqual(cert['floors']['radial'][0]['next_radial_index'], 2)

    def test_missing_angular_modes_keep_correct_indices(self):
        assembly = s.assemble_sectors([0], [make_group([0], 0, 5)])
        cert = s.ordered_bounds(s.normalize_request([0], k=24), assembly)
        self.assertEqual([F(x['upper']) for x in cert['ordered_intervals']], exact_levels(24))
        self.assertEqual(cert['floors']['omitted_azimuth_m_start'], 1)
        self.assertTrue(any(row['m'] == 4 for row in cert['uncomputed_modes_may_touch_first_k']))

    def test_nonconstant_missing_modes_give_valid_wider_enclosures(self):
        q = [-6, 0, 2]
        request = s.normalize_request(q, k=8)
        empty = s.ordered_bounds(request, s.assemble_sectors(q, []))
        groups = [make_group(q, m, 3) for m in range(4)]
        populated = s.ordered_bounds(request, s.assemble_sectors(q, groups))
        for wide, narrow in zip(empty['ordered_intervals'], populated['ordered_intervals']):
            self.assertLessEqual(F(wide['lower']), F(narrow['lower']))
            self.assertGreaterEqual(F(wide['upper']), F(narrow['upper']))
        self.assertEqual(empty['status'], 'certified_bound_open_gap')

    def test_real_angular_multiplicity_and_tampering(self):
        assembly = s.assemble_sectors([0], [make_group([0], 0, 1), make_group([0], 1, 1)])
        self.assertEqual([row['multiplicity'] for row in assembly['entries']], [1, 2])
        self.assertEqual(assembly['represented_eigenvalues_with_multiplicity'], 3)
        cert = s.ordered_bounds(s.normalize_request([0], k=3), assembly)
        cert['assembly']['entries'][1]['multiplicity'] = 1
        self.assertFalse(s.verify(cert))

    def test_zero_and_threshold_ties_are_not_negative(self):
        assembly = s.assemble_sectors([-6], [])
        count = s.count_below(assembly, 0)
        self.assertEqual((count['count_lower'], count['count_upper']), (3, 3))
        self.assertEqual(count['exact_tie_multiplicity'], 5)
        exact = s.count_below(assembly, -4)
        self.assertEqual((exact['count_lower'], exact['count_upper']), (0, 0))
        self.assertEqual(exact['exact_tie_multiplicity'], 3)
        above = s.count_below(assembly, F(-399, 100))
        self.assertEqual((above['count_lower'], above['count_upper']), (3, 3))

    def test_all_negative_trace_can_come_from_omitted_modes(self):
        assembly = s.assemble_sectors([-20], [])
        trace = s.negative_trace(assembly)
        # l=1,2,3 give multiplicities 3,5,7; l=4 is exactly zero.
        expected = 3*18+5*14+7*8
        self.assertEqual(F(trace['lower']), expected)
        self.assertEqual(F(trace['upper']), expected)
        self.assertEqual(F(trace['omitted_negative_trace_lower']), expected)
        self.assertEqual(F(trace['omitted_negative_trace_upper']), expected)
        count = s.count_below(assembly)
        self.assertEqual(count['count_lower'], 15)
        self.assertEqual(count['exact_tie_multiplicity'], 9)

    def test_partial_negative_trace_keeps_omitted_angular_contribution(self):
        assembly = s.assemble_sectors([-20], [make_group([-20], 0, 3)])
        trace = s.negative_trace(assembly)
        self.assertEqual(F(trace['lower']), 180)
        self.assertEqual(F(trace['upper']), 180)
        self.assertEqual(F(trace['omitted_negative_trace_upper']), 140)
        bad = copy.deepcopy(trace)
        bad['omitted_negative_trace_upper'] = '0'
        bad['upper'] = '40'
        self.assertFalse(s.verify(bad))

    def test_unresolved_negative_tail_has_two_sided_contribution(self):
        trace = s.negative_trace(s.assemble_sectors([-20, 0, 2], []))
        self.assertEqual(F(trace['lower']), 3*16+5*12+7*6)
        self.assertEqual(F(trace['upper']), 180)
        self.assertFalse(trace['omitted_negative_tail_resolved'])
        self.assertEqual(trace['status'], 'certified_trace_open_gap')
        self.assertTrue(s.verify(trace))

    def test_lt_area_normalization_factor_four_and_scope(self):
        trace = s.negative_trace(s.assemble_sectors([-20], []))
        quotient = s.lt_quotient(trace, [20])
        self.assertEqual(F(quotient['mean_V_squared']), 400)
        self.assertEqual(F(quotient['rational_divisor']), 1600)
        self.assertEqual(F(quotient['lower']), F(9, 80))
        self.assertTrue(quotient['not_an_upper_bound_for_any_universal_constant'])
        self.assertEqual(quotient['claim'], 'fixed_potential_witness_lower_bound_only')
        self.assertEqual(quotient['orthonormal_constant_relation'], 'k_LT=4*L_trace in dimension two')
        nonconstant = s.lt_quotient(s.negative_trace(s.assemble_sectors([-2, 0, -1], [])), [2, 0, 1])
        self.assertEqual(F(nonconstant['mean_V_squared']), F(83, 15))
        self.assertEqual(F(nonconstant['rational_divisor']), F(332, 15))
        quotient['rational_divisor'] = '400'
        self.assertFalse(s.verify(quotient))

    def test_reject_invalid_lt_potential_and_foreign_operator(self):
        trace = s.negative_trace(s.assemble_sectors([-20], []))
        with self.assertRaises(ValueError):
            s.lt_quotient(trace, [19])
        for potential, q in (([-1], [1]), ([0], [0])):
            with self.assertRaises(ValueError):
                s.lt_quotient(s.negative_trace(s.assemble_sectors(q, [])), potential)
        with self.assertRaises(ValueError):
            s.normalize_request([0], mean_zero=False)

    def test_consecutive_indices_and_bundle_binding(self):
        with self.assertRaises(ValueError):
            s.assemble_sectors([0], [[projected.certify_sector([0], m=0, k=2, modes=3)]])
        with self.assertRaises(ValueError):
            s.assemble_sectors([0], [make_group([0], 1, 1)])
        ordered = s.ordered_bounds(s.normalize_request([-6], k=8), s.assemble_sectors([-6], []))
        cert = s.bundle(ordered, assertion={'index': 4, 'lower_bound': 0})
        self.assertEqual(cert['assertion_result'], 'proved')
        cert['negative_trace'] = s.negative_trace(s.assemble_sectors([-20], []))
        self.assertFalse(s.verify(cert))

    def test_saved_stage_and_budget_evidence_replay(self):
        root = Path(__file__).parent/'round07'
        load = lambda name: json.loads((root/name).read_text())
        self.assertTrue(s.verify_initial(load('v155.json')))
        calibrator = load('v156.json')
        self.assertEqual(calibrator, s.laplace_spectrum(calibrator['k'], calibrator['shift']))
        self.assertEqual(s._assembly(load('v157.json')), load('v157.json'))
        self.assertEqual(s.omitted_floors(load('v157.json')), load('v158.json'))
        for name in ('v159.json', 'v161.json', 'v162.json', 'v163.json', 'v164.json'):
            with self.subTest(name=name):
                self.assertTrue(s.verify(load(name)))
        self.assertTrue(s.verify(load('v160.json')['certificate']))
        budget = load('budget_exhaustion.json')
        self.assertTrue(budget['budget_exhausted_before_target'])
        self.assertEqual(budget['certificate']['status'], 'certified_bound_open_gap')
        for key in ('certificate', 'negative_trace', 'threshold_minus3_count'):
            self.assertTrue(s.verify(budget[key]))

    def test_initial_minmax_covers_every_repeated_level(self):
        cert = s.initial_order_bounds(s.normalize_request([-6, 0, 2], k=8))
        self.assertEqual([F(row['lower']) for row in cert['ordered_intervals']],
                         [F(-4)]*3+[F(0)]*5)
        self.assertEqual([F(row['upper']) for row in cert['ordered_intervals']],
                         [F(-2)]*3+[F(2)]*5)
        self.assertEqual(cert['schur_computations'], 0)
        self.assertTrue(s.verify_initial(cert, [-6, 0, 2], 8, True))
        self.assertFalse(s.verify_initial(cert, expected_q=[-6]))
        self.assertFalse(s.verify_initial(cert, expected_k=7))
        self.assertFalse(s.verify_initial(cert, expected_mean_zero=False))
        cert['ordered_intervals'][1]['lower'] = '0'
        self.assertFalse(s.verify_initial(cert))


if __name__ == '__main__':
    unittest.main()
