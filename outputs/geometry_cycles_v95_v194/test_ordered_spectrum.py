import copy
import json
import unittest
from common import F, projected
import ordered_spectrum as s


class FullOrderedSpectrumTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.free = s.adaptive_spectrum([0], k=8)['certificate']
        cls.shift6 = s.adaptive_spectrum([-6], k=8)['certificate']
        cls.shift20 = s.adaptive_spectrum([-20], k=16)['certificate']
        cls.variable = s.adaptive_spectrum([-6, 0, 2], k=8, modes=6,
            max_modes=8, max_m=3, max_radial=3, bits=28,
            tolerance='1/100000')['certificate']
        cls.low_budget = s.adaptive_spectrum([-6, 0, 2], k=8, modes=4,
            max_modes=4, max_m=0, max_radial=1, bits=20,
            tolerance='1/100000')['certificate']

    def test_full_first_eight_repeated_eigenvalues(self):
        self.assertEqual([F(v['lower']) for v in self.free['ordered_intervals']],
                         [2, 2, 2, 6, 6, 6, 6, 6])
        self.assertTrue(all(v['width'] == '0' for v in self.free['ordered_intervals']))

    def test_exact_calibrator_constant_shift(self):
        self.assertEqual(s.laplace_spectrum(4, -6)['eigenvalues'], ['-4', '-4', '-4', '0'])

    def test_initial_global_bounds_without_schur(self):
        c = s.initial_order_bounds(s.normalize_request([-6, 0, 2], k=8))
        self.assertEqual([(F(x['lower']), F(x['upper'])) for x in c['ordered_intervals']],
                         [(-4, -2)]*3+[(0, 2)]*5)
        self.assertEqual(c['schur_computations'], 0)
        self.assertTrue(s.verify_initial(c, expected_q=[-6, 0, 2], expected_k=8, expected_mean_zero=True))
        self.assertTrue(s.verify(json.loads(json.dumps(c))))

    def test_initial_bounds_bind_q_k_and_projection(self):
        c = s.initial_order_bounds(s.normalize_request([-6, 0, 2], k=8))
        self.assertFalse(s.verify_initial(c, expected_q=[-6, 0, 3]))
        self.assertFalse(s.verify_initial(c, expected_k=7))
        self.assertFalse(s.verify_initial(c, expected_k=True))
        self.assertFalse(s.verify_initial(c, expected_mean_zero=False))
        altered = copy.deepcopy(c)
        altered['ordered_intervals'][0]['lower'] = '-3'
        self.assertFalse(s.verify_initial(altered))
        altered = copy.deepcopy(c)
        altered['request']['k'] = 7
        self.assertFalse(s.verify_initial(altered))

    def test_initial_constant_bounds_are_exact_with_full_multiplicity(self):
        c = s.initial_order_bounds(s.normalize_request([-20], k=16))
        self.assertEqual([x['lower'] for x in c['ordered_intervals']],
                         ['-18']*3+['-14']*5+['-8']*7+['0'])
        self.assertTrue(all(x['width'] == '0' for x in c['ordered_intervals']))

    def test_scope_and_rational_input_validation(self):
        for options in ({'mean_zero': False}, {'mean_zero': 1}, {'k': True},
                        {'max_m': True}, {'tolerance': 0.01}, {'bits': 4}):
            with self.assertRaises(ValueError):
                s.normalize_request([0], **options)
        for q in ([True], [0.1], [1]*8):
            with self.assertRaises(ValueError):
                s.normalize_request(q)

    def test_real_angular_multiplicity(self):
        groups = [[projected.certify_sector([0], m=m, k=1, modes=4, bits=16)] for m in range(3)]
        a = s.assemble_sectors([0], groups)
        self.assertEqual([e['multiplicity'] for e in a['entries']], [1, 2, 2])
        self.assertEqual(a['represented_eigenvalues_with_multiplicity'], 5)

    def test_sector_swap_rejected(self):
        c0 = projected.certify_sector([0], m=0, modes=4, bits=16)
        c1 = projected.certify_sector([0], m=1, modes=4, bits=16)
        with self.assertRaises(ValueError):
            s.assemble_sectors([0], [[c1], [c0]])

    def test_radial_gap_rejected(self):
        c2 = projected.certify_sector([0], m=0, k=2, modes=4, bits=16)
        with self.assertRaises(ValueError):
            s.assemble_sectors([0], [[c2]])

    def test_omitted_radial_and_angular_floors(self):
        c = projected.certify_sector([-6], m=0, modes=4, bits=16)
        floors = s.omitted_floors(s.assemble_sectors([-6], [[c]]))
        self.assertEqual(floors['radial'][0]['all_unenumerated_radial_eigenvalues_at_least'], '0')
        self.assertEqual(floors['all_unenumerated_angular_eigenvalues_at_least'], '-4')

    def test_resource_exhaustion_retains_safe_full_bounds(self):
        self.assertEqual(self.low_budget['status'], 'certified_bound_open_gap')
        self.assertTrue(self.low_budget['uncomputed_modes_may_touch_first_k'])
        for coarse, fine in zip(self.low_budget['ordered_intervals'], self.variable['ordered_intervals']):
            self.assertLessEqual(F(coarse['lower']), F(fine['lower']))
            self.assertGreaterEqual(F(coarse['upper']), F(fine['upper']))

    def test_empty_assembly_still_covers_full_spectrum(self):
        c = s.ordered_bounds(s.normalize_request([0], k=8), s.assemble_sectors([0], []))
        self.assertEqual([x['lower'] for x in c['ordered_intervals']], ['2']*3+['6']*5)
        self.assertTrue(c['first_k_resolved_without_enumerating_all_modes'])

    def test_strict_count_at_repeated_eigenvalue(self):
        c = s.count_below(self.free['assembly'], 2)
        self.assertEqual((c['count_lower'], c['count_upper'], c['exact_tie_multiplicity']), (0, 0, 3))
        c = s.count_below(self.free['assembly'], 6)
        self.assertEqual((c['count_lower'], c['count_upper'], c['exact_tie_multiplicity']), (3, 3, 5))

    def test_zero_cutoff_does_not_count_zero_eigenvalues(self):
        c = s.count_below(self.shift6['assembly'])
        self.assertEqual((c['count_lower'], c['count_upper'], c['exact_tie_multiplicity']), (3, 3, 5))
        trace = s.negative_trace(self.shift6['assembly'])
        self.assertEqual((trace['lower'], trace['upper']), ('12', '12'))

    def test_minus_twenty_complete_negative_trace(self):
        c = s.count_below(self.shift20['assembly'])
        self.assertEqual((c['count_lower'], c['count_upper']), (15, 15))
        trace = s.negative_trace(self.shift20['assembly'])
        self.assertEqual((trace['lower'], trace['upper']), ('180', '180'))
        self.assertGreater(F(trace['omitted_negative_trace_upper']), 0)
        self.assertTrue(trace['omitted_negative_tail_resolved'])

    def test_omitted_nonconstant_negative_tail_is_not_discarded(self):
        trace = s.negative_trace(self.low_budget['assembly'])
        self.assertGreater(F(trace['omitted_negative_trace_upper']), F(trace['omitted_negative_trace_lower']))
        self.assertEqual(trace['status'], 'certified_trace_open_gap')

    def test_unresolved_threshold_cluster_is_explicit(self):
        count = s.count_below(self.low_budget['assembly'], -3)
        self.assertGreater(count['count_upper'], count['count_lower'])
        self.assertTrue(count['unresolved_threshold_clusters'])

    def test_LT_scale_and_fixed_potential_semantics(self):
        q = s.lt_quotient(s.negative_trace(self.shift20['assembly']), [20])
        self.assertEqual((q['lower'], q['upper']), ('9/80', '9/80'))
        self.assertEqual(q['orthonormal_constant_relation'], 'k_LT=4*L_trace in dimension two')
        self.assertTrue(q['not_an_upper_bound_for_any_universal_constant'])

    def test_LT_scope_potential_and_zero_denominator_rejected(self):
        with self.assertRaises(ValueError):
            s.lt_quotient(s.negative_trace(self.free['assembly']), [0])
        with self.assertRaises(ValueError):
            s.lt_quotient(s.negative_trace(self.shift20['assembly']), [19])
        bad = s.assemble_sectors([0, 1], [])
        with self.assertRaises(ValueError):
            s.lt_quotient(s.negative_trace(bad), [0, -1])

    def test_json_replay_and_expected_input(self):
        c = json.loads(json.dumps(self.variable))
        self.assertTrue(s.verify(c, expected_q=[-6, 0, 2]))
        self.assertFalse(s.verify(c, expected_q=[-6, 0, 3]))

    def test_multiplicity_and_missing_sector_tampering_rejected(self):
        c = copy.deepcopy(self.variable)
        c['assembly']['entries'][0]['multiplicity'] = 2
        self.assertFalse(s.verify(c))
        c = copy.deepcopy(self.variable)
        c['assembly']['sector_lists'].pop(0)
        self.assertFalse(s.verify(c))

    def test_count_trace_and_ordering_tampering_rejected(self):
        c = s.count_below(self.shift6['assembly'])
        c['count_upper'] = 8
        self.assertFalse(s.verify(c))
        t = s.negative_trace(self.shift20['assembly'])
        t['upper'] = '54'
        self.assertFalse(s.verify(t))
        c = copy.deepcopy(self.variable)
        c['ordered_intervals'][0]['lower'] = c['ordered_intervals'][0]['upper']
        self.assertFalse(s.verify(c))

    def test_bundle_and_three_way_assertion(self):
        yes = s.bundle(self.free, assertion={'index': 4, 'lower_bound': 6})
        no = s.bundle(self.free, assertion={'index': 4, 'lower_bound': 7})
        maybe = s.bundle(self.low_budget, assertion={'index': 1, 'lower_bound': -3})
        self.assertEqual([x['assertion_result'] for x in (yes, no, maybe)],
                         ['proved', 'refuted', 'undetermined'])
        self.assertTrue(s.verify(yes))
        yes['count'] = s.count_below(self.shift6['assembly'])
        self.assertFalse(s.verify(yes))


if __name__ == '__main__':
    unittest.main()
