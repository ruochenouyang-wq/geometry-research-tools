"""C14 independent review: full operator, complete trace and scope binding."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import unittest

import general_spectrum as g
from support import sectors


class GeneralSpectrumIndependentReview(unittest.TestCase):
    def test_unprojected_l0_and_complete_real_multiplicity(self):
        full = g.analytic_bounds([0], 10, False)
        centered = g.analytic_bounds([0], 9, True)
        self.assertEqual([F(x['lower']) for x in full['ordered']], [0]+[2]*3+[6]*5+[12])
        self.assertEqual([F(x['lower']) for x in centered['ordered']], [2]*3+[6]*5+[12])
        self.assertEqual(full['l0_multiplicity'], 1)
        self.assertEqual(centered['l0_multiplicity'], 0)

    def test_original_c14_both_spaces_and_no_negative_tail(self):
        for mean_zero, count, trace in ((True, 8, 20), (False, 9, 27)):
            model = g.assemble([-7], [], mean_zero)
            counted, traced = g.strict_count(model), g.negative_trace(model)
            self.assertEqual((counted['count_lower'], counted['count_upper']), (count, count))
            self.assertEqual((F(traced['lower']), F(traced['upper'])), (trace, trace))
            self.assertGreater(F(traced['tails']['unlisted_floor']), 0)
            self.assertTrue(g.verify(counted, [-7], mean_zero))
            self.assertTrue(g.verify(traced, [-7], mean_zero))

    def test_missing_radial_and_angular_modes_remain_in_ordering(self):
        group = g.enumerate_sector([-7], 0, 1, False, modes=4)
        model = g.assemble([-7], [group], False)
        cert = g.ordered_merge(model, 10)
        self.assertEqual([F(x['lower']) for x in cert['ordered']], [-7]+[-5]*3+[-1]*5+[5])
        self.assertTrue(any(row['uncomputed'] and row['m'] == 0 and row['j'] == 3 for row in cert['slots']))
        self.assertTrue(any(row['uncomputed'] and row['m'] == 2 for row in cert['slots']))

    def test_nonconstant_potential_couples_constant_mode(self):
        # For u=1-t/2: M=13/12, energy=1/6, potential=-1/3, R=-2/13.
        mass, kinetic, potential = F(13, 12), F(1, 6), F(-1, 3)
        rayleigh = (kinetic+potential)/mass
        self.assertEqual(rayleigh, F(-2, 13))
        a = sectors.assemble([0, 1], 0, 0, 2)
        self.assertEqual(a['A'][0][1], F(2, 3))
        unprojected = g.enumerate_sector([0, 1], 0, 1, False, modes=6, bits=48)
        projected = g.enumerate_sector([0, 1], 0, 1, True, modes=6, bits=48)
        self.assertLess(F(unprojected['eigenvalues'][0]['upper']), rayleigh)
        self.assertGreater(F(projected['eigenvalues'][0]['lower']), 1)
        self.assertEqual(unprojected['degree_start'], 0)
        self.assertEqual(projected['degree_start'], 1)

    def test_zero_is_strictly_excluded_in_both_spaces(self):
        for mean_zero in (False, True):
            model = g.assemble([-2], [], mean_zero)
            count = g.strict_count(model)
            trace = g.negative_trace(model)
            self.assertEqual(count['count_lower'], 0 if mean_zero else 1)
            self.assertEqual(count['count_upper'], count['count_lower'])
            self.assertEqual(count['exact_tie_multiplicity'], 3)
            self.assertEqual(F(trace['upper']), 0 if mean_zero else 2)
        zero = g.strict_count(g.assemble([0], [], False))
        self.assertEqual(zero['count_upper'], 0)
        self.assertEqual(zero['exact_tie_multiplicity'], 1)

    def test_uncomputed_negative_trace_is_not_discarded(self):
        full = g.negative_trace(g.assemble([-20, 0, 1], [], False))
        self.assertEqual(F(full['lower']), 19+3*17+5*13+7*7)
        self.assertEqual(F(full['upper']), 20+3*18+5*14+7*8)
        self.assertEqual(full['lower'], full['omitted_negative_trace_lower'])
        self.assertEqual(full['upper'], full['omitted_negative_trace_upper'])
        self.assertFalse(full['target_met'])
        tampered = copy.deepcopy(full)
        tampered['omitted_negative_trace_upper'] = '0'
        self.assertFalse(g.verify(tampered))

    def test_trace_driven_search_does_not_stop_after_ground_closes(self):
        options = dict(k=1, mean_zero=False, quantities=['ordered', 'count', 'trace'],
                       tolerance=F(1, 10**6), modes=4, max_modes=8,
                       max_m=4, max_radial=5, bits=40)
        result = g.adaptive_spectrum([-20, 0, 1], max_steps=16, **options)
        early = [row for row in result['search_history']
                 if row['decisions']['ordered'] == 'certified_target_met'
                 and row['decisions']['trace'] == 'certified_bound_open']
        self.assertTrue(early)
        self.assertGreater(result['search_history'][-1]['step'], early[0]['step'])
        self.assertTrue(result['certificate']['all_requested_targets_met'])
        self.assertTrue(g.verify(result['certificate']))
        limited = g.adaptive_spectrum([-20, 0, 1], max_steps=1, **options)
        self.assertTrue(limited['certificate']['answers']['ordered']['target_met'])
        self.assertFalse(limited['certificate']['answers']['trace']['target_met'])
        self.assertFalse(limited['certificate']['all_requested_targets_met'])
        self.assertTrue(limited['budget_or_precision_exhausted'])

    def test_zero_budget_has_honest_per_quantity_results(self):
        cert = g.adaptive_spectrum([-20, 0, 1], k=1, max_steps=0)['certificate']
        self.assertFalse(cert['answers']['ordered']['target_met'])
        self.assertFalse(cert['answers']['trace']['target_met'])
        self.assertTrue(cert['answers']['count']['target_met'])
        self.assertTrue(g.verify(cert))

    def test_constant_spectral_shift_preserves_projection(self):
        for mean_zero in (False, True):
            shifted = g.spectral_shift(g.assemble([0], [], mean_zero), 2)
            self.assertEqual(shifted['shifted_q'], ['-2'])
            self.assertIs(shifted['shifted_mean_zero'], mean_zero)
            self.assertEqual(shifted['strict_negative_count_for_shifted_operator']['count_upper'],
                             0 if mean_zero else 1)
            self.assertEqual(F(shifted['shifted_negative_trace_upper']), 0 if mean_zero else 2)
            self.assertTrue(g.verify(shifted))

    def test_small_negative_constant_disproves_unprojected_uniform_lt(self):
        values = []
        for epsilon in (F(1, 2), F(1, 16), F(1, 256)):
            trace = g.negative_trace(g.assemble([-epsilon], [], False))
            quotient = g.fixed_potential_quotient(trace, [epsilon])
            self.assertEqual(F(trace['lower']), epsilon)
            self.assertEqual(F(quotient['lower']), 1/(4*epsilon))
            self.assertFalse(quotient['finite_universal_L2_trace_inequality_compatible_with_scope'])
            self.assertIsNotNone(quotient['unprojected_constant_mode_obstruction'])
            self.assertTrue(quotient['not_a_universal_upper_bound'])
            self.assertTrue(g.verify(quotient))
            values.append(F(quotient['lower']))
        self.assertEqual(values, sorted(values))

    def test_interlacing_order_and_potential_projection_binding(self):
        full = g.ordered_merge(g.assemble([-7], [], False), 10)
        mean_zero = g.ordered_merge(g.assemble([-7], [], True), 9)
        cert = g.projection_interlacing(full, mean_zero, 9)
        for row in cert['rows']:
            self.assertLessEqual(F(row['full_k']['lower']), F(row['projected_k']['lower']))
            self.assertLessEqual(F(row['projected_k']['upper']), F(row['full_k_plus_1']['upper']))
        self.assertTrue(g.verify(cert, expected_q=[-7], expected_k=9))
        with self.assertRaises(ValueError):
            g.projection_interlacing(mean_zero, full)
        with self.assertRaises(ValueError):
            g.projection_interlacing(full, g.ordered_merge(g.assemble([-6], [], True), 9))
        cert['relation'] = 'lambda_k(mean_zero) <= lambda_k(full)'
        self.assertFalse(g.verify(cert))

    def test_nonconstant_interlacing_uses_same_original_q(self):
        q = [0, 1]
        full_groups = [g.enumerate_sector(q, m, 2, False, 5, 40) for m in range(2)]
        projected_groups = [g.enumerate_sector(q, m, 2, True, 5, 40) for m in range(2)]
        full = g.ordered_merge(g.assemble(q, full_groups, False), 5)
        projected = g.ordered_merge(g.assemble(q, projected_groups, True), 4)
        cert = g.projection_interlacing(full, projected, 4)
        self.assertTrue(g.verify(cert, expected_q=q))
        self.assertTrue(cert['nonconstant_q_does_not_have_an_added_constant_eigenvalue'])
        self.assertFalse(g.verify(cert, expected_q=[0]))

    def test_driver_cache_and_scope_are_full_input_bound(self):
        cache = {}
        first = g.research_driver([-7], cache=cache, k=9)
        self.assertTrue(g.verify(first['certificate'], expected_q=[-7]))
        second = g.research_driver([-7], cache=cache, k=9)
        self.assertEqual(second['execution']['cache_hits'], 2)
        bad_key = next(iter(cache))
        cache[bad_key]['request']['q'] = ['-6']
        restored = g.research_driver([-7], cache=cache, k=9)
        self.assertEqual(restored['execution']['cache_rejections'], 1)
        self.assertTrue(g.verify(restored['certificate'], expected_q=[-7]))
        proof = first['certificate']['results'][0]
        self.assertFalse(g.verify(proof, expected_mean_zero=False))
        self.assertFalse(g.verify(proof, expected_k=8))

    def test_missing_sector_wrong_multiplicity_and_scope_rejected(self):
        group = g.enumerate_sector([0], 1, 1, False, 4, 32)
        with self.assertRaises(ValueError):
            g.assemble([0], [group], False)
        model = g.assemble([0], [], False)
        model['real_angular_multiplicity']['m>0'] = 1
        self.assertFalse(g.verify(model))
        cert = g.analytic_bounds([0], 5, False)
        cert['space'] = 'full_mean_zero'
        self.assertFalse(g.verify(cert))

    def test_valid_wrong_format_cache_is_rejected_and_recomputed(self):
        req = g.request([-7], mean_zero=True, k=9)
        for wrong in (g.analytic_bounds([-7], 9, True),
                      g.negative_trace(g.assemble([-7], [], True))):
            cache = {g.digest(req): wrong}
            answer = g.research_driver([-7], spaces=['full_mean_zero'], cache=cache, k=9)
            self.assertEqual(answer['execution']['cache_rejections'], 1)
            self.assertTrue(g.verify(answer['certificate']))
            with self.assertRaises(ValueError):
                g.driver_certificate([req], [wrong])

    def test_every_saved_c14_certificate_replays(self):
        formats = {g.INITIAL, g.SECTOR, g.MODEL, g.ORDERED, g.COUNT, g.TRACE,
                   g.MULTI, g.SHIFT, g.LT, g.INTERLACE, g.DRIVER}
        count = 0
        def walk(value, path):
            nonlocal count
            if isinstance(value, dict):
                if value.get('format') in formats:
                    with self.subTest(path=path):
                        self.assertTrue(g.verify(value))
                    count += 1
                    return
                for key, item in value.items():
                    walk(item, path+'/'+str(key))
            elif isinstance(value, list):
                for i, item in enumerate(value):
                    walk(item, path+'/'+str(i))
        for path in sorted((Path(__file__).parent/'C14').rglob('*.json')):
            walk(json.loads(path.read_text()), str(path.name))
        self.assertGreaterEqual(count, 10)


if __name__ == '__main__':
    unittest.main()
