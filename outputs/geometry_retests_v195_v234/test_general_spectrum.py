import copy
import json
import unittest
from support import F, sectors, digest, old_spectrum
import general_spectrum as g


class GeneralSpectrumTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = g.research_driver([-7])['certificate']
        cls.models = {mean: g.assemble([-7], [], mean) for mean in (False, True)}
        cls.zero = {mean: g.assemble([0], [], mean) for mean in (False, True)}
        cls.linear = {mean: g.adaptive_spectrum([0, 1], k=6 if not mean else 5,
            mean_zero=mean, tolerance='1/1000000', bits=32)['certificate']
            for mean in (False, True)}
        cls.trace_aware = g.adaptive_spectrum([-20, 0, 1], k=1, mean_zero=True,
            quantities=['ordered', 'trace'], tolerance='1/1000000', max_m=4,
            max_radial=6, bits=32)

    def test_original_old_unprojected_request_rejected(self):
        with self.assertRaisesRegex(ValueError, 'mean_zero=True'):
            old_spectrum.adaptive_spectrum([-7], k=9, mean_zero=False)

    def test_original_C14_both_spaces_and_values(self):
        got = {}
        for result in self.original['results']:
            answer = result['answers']
            got[result['request']['space']] = (answer['count']['count_lower'],
                answer['count']['count_upper'], answer['trace']['lower'], answer['trace']['upper'])
        self.assertEqual(got, {'full_mean_zero': (8, 8, '20', '20'),
                               'full_unprojected': (9, 9, '27', '27')})
        self.assertTrue(g.verify(self.original, expected_q=[-7]))

    def test_unprojected_zero_first_nine(self):
        c = g.analytic_bounds([0], 9, False)
        self.assertEqual([F(x['lower']) for x in c['ordered']], [0, 2, 2, 2, 6, 6, 6, 6, 6])
        self.assertEqual(c['l0_multiplicity'], 1)

    def test_mean_zero_zero_first_nine(self):
        c = g.analytic_bounds([0], 9, True)
        self.assertEqual([F(x['lower']) for x in c['ordered']], [2, 2, 2, 6, 6, 6, 6, 6, 12])
        self.assertEqual(c['l0_multiplicity'], 0)

    def test_initial_nonconstant_global_bounds(self):
        c = g.analytic_bounds([0, 1], 4, False)
        self.assertEqual([(F(x['lower']), F(x['upper'])) for x in c['ordered']],
                         [(-1, 1), (1, 3), (1, 3), (1, 3)])
        self.assertEqual(c['schur_computations'], 0)

    def test_scope_and_input_type_validation(self):
        for value in (0, 1, 'false', None):
            with self.assertRaises(ValueError): g.analytic_bounds([0], mean_zero=value)
        for q in ([0.1], [True], [1]*8):
            with self.assertRaises(ValueError): g.analytic_bounds(q)
        with self.assertRaises(ValueError): g.analytic_bounds([0], k=True)

    def test_unprojected_sector_starts_at_l0(self):
        p = g.enumerate_sector([-7], m=0, radial_count=3, mean_zero=False, modes=4, bits=16)
        self.assertEqual(p['degree_start'], 0)
        self.assertEqual([c['lower'] for c in p['eigenvalues']], ['-7', '-5', '-1'])
        self.assertTrue(g.verify(p, expected_mean_zero=False))

    def test_projected_sector_starts_at_l1(self):
        p = g.enumerate_sector([-7], m=0, radial_count=3, mean_zero=True, modes=4, bits=16)
        self.assertEqual(p['degree_start'], 1)
        self.assertEqual([c['lower'] for c in p['eigenvalues']], ['-5', '-1', '5'])

    def test_nonlinear_translation_by_appending_constant_is_invalid(self):
        full = self.linear[False]['answers']['ordered']['ordered'][0]
        mean = self.linear[True]['answers']['ordered']['ordered'][0]
        self.assertLess(F(full['upper']), 0)
        self.assertGreater(F(mean['lower']), 0)
        # q=t couples l=0 to l=1; its average 0 is not an added eigenvalue.
        form = sectors.assemble([0, 1], 0, 0, 2)
        self.assertEqual(form['A'][0][1], F(2, 3))

    def test_full_ground_respects_explicit_mixed_trial(self):
        # u=1-t/2: mean energy=-1/6 and mean mass=13/12.
        upper = F(self.linear[False]['answers']['ordered']['ordered'][0]['upper'])
        self.assertLessEqual(upper, -F(2, 13))

    def test_angular_multiplicities_remain_real(self):
        groups = [g.enumerate_sector([0], m=m, mean_zero=False, modes=4, bits=16) for m in range(3)]
        a = g.assemble([0], groups, False)
        self.assertEqual([e['multiplicity'] for e in a['entries']], [1, 2, 2])

    def test_wrong_projection_sector_rejected(self):
        group = g.enumerate_sector([0], mean_zero=True, modes=4, bits=16)
        with self.assertRaises(ValueError): g.assemble([0], [group], False)

    def test_missing_angular_sector_rejected(self):
        group = g.enumerate_sector([0], m=1, mean_zero=False, modes=4, bits=16)
        with self.assertRaises(ValueError): g.assemble([0], [group], False)

    def test_missing_radial_prefix_rejected(self):
        group = g.enumerate_sector([0], radial_count=2, mean_zero=False, modes=4, bits=16)
        group['eigenvalues'].pop(0)
        with self.assertRaises(ValueError): g.assemble([0], [group], False)

    def test_complete_merge_has_l0_slot(self):
        ordered = g.ordered_merge(self.zero[False], k=9)
        self.assertEqual([x['l'] for x in ordered['slots']].count(0), 1)
        self.assertTrue(ordered['target_met'])

    def test_empty_enumeration_preserves_all_negative_modes(self):
        t = g.negative_trace(self.models[False])
        self.assertEqual((t['lower'], t['upper'], t['omitted_negative_trace_upper']), ('27', '27', '27'))
        self.assertGreater(F(t['tails']['unlisted_floor']), 0)

    def test_exact_zero_ties_excluded_unprojected(self):
        model = g.assemble([-6], [], False)
        c, t = g.strict_count(model), g.negative_trace(model)
        self.assertEqual((c['count_lower'], c['count_upper'], c['exact_tie_multiplicity']), (4, 4, 5))
        self.assertEqual(t['lower'], '18')

    def test_exact_zero_ties_excluded_projected(self):
        model = g.assemble([-6], [], True)
        c, t = g.strict_count(model), g.negative_trace(model)
        self.assertEqual((c['count_lower'], c['count_upper'], c['exact_tie_multiplicity']), (3, 3, 5))
        self.assertEqual(t['lower'], '12')

    def test_qzero_full_zero_is_not_negative(self):
        c = g.strict_count(self.zero[False])
        self.assertEqual((c['count_lower'], c['count_upper'], c['exact_tie_multiplicity']), (0, 0, 1))

    def test_positive_potential_has_no_negative_tail(self):
        t = g.negative_trace(g.assemble([3], [], False))
        self.assertEqual((t['lower'], t['upper']), ('0', '0'))

    def test_trace_aware_search_continues_after_order_closes(self):
        history = self.trace_aware['search_history']
        interim = [h for h in history if h['decisions']['ordered'] == 'certified_target_met'
                   and h['decisions']['trace'] == 'certified_bound_open']
        self.assertTrue(interim)
        self.assertGreater(history[-1]['step'], interim[0]['step'])
        self.assertTrue(self.trace_aware['certificate']['all_requested_targets_met'])

    def test_trace_aware_low_budget_keeps_separate_open_status(self):
        low = g.adaptive_spectrum([-20, 0, 1], k=1, mean_zero=True,
            quantities=['ordered', 'trace'], tolerance='1/1000000', bits=32, max_steps=2)
        self.assertEqual(low['certificate']['decisions'],
                         {'ordered': 'certified_target_met', 'trace': 'certified_bound_open'})
        self.assertTrue(low['budget_or_precision_exhausted'])
        self.assertTrue(g.verify(low['certificate']))

    def test_zero_budget_returns_safe_analytic_intervals(self):
        low = g.adaptive_spectrum([-7, 0, 1], max_steps=0)
        self.assertEqual(low['sector_solves'], 0)
        self.assertTrue(low['budget_or_precision_exhausted'])
        self.assertTrue(g.verify(low['certificate']))

    def test_nonzero_threshold_shift_full(self):
        c = g.spectral_shift(self.zero[False], 6)
        self.assertEqual(c['shifted_q'], ['-6'])
        self.assertFalse(c['shifted_mean_zero'])
        self.assertEqual(c['strict_negative_count_for_shifted_operator']['count_lower'], 4)
        self.assertEqual(c['shifted_negative_trace_lower'], '18')

    def test_nonzero_threshold_shift_projected(self):
        c = g.spectral_shift(self.zero[True], 6)
        self.assertTrue(c['shifted_mean_zero'])
        self.assertEqual(c['strict_negative_count_for_shifted_operator']['count_lower'], 3)
        self.assertEqual(c['shifted_negative_trace_lower'], '12')

    def test_rational_threshold_shift_preserves_spectrum(self):
        c = g.spectral_shift(g.assemble([0], [], False), F(5, 2))
        direct = g.negative_trace(g.assemble(['-5/2'], [], False))
        self.assertEqual(c['shifted_negative_trace_lower'], direct['lower'])

    def test_full_tiny_constant_ratio_diverges(self):
        for eps, expected in ((F(1, 10), F(5, 2)), (F(1, 100), F(25))):
            trace = g.negative_trace(g.assemble([-eps], [], False))
            ratio = g.fixed_potential_quotient(trace, [eps])
            self.assertEqual(F(ratio['lower']), expected)
            self.assertFalse(ratio['finite_universal_L2_trace_inequality_compatible_with_scope'])
            self.assertEqual(ratio['unprojected_constant_mode_obstruction']['limit_as_epsilon_to_zero_positive'], 'positive_infinity')

    def test_projected_tiny_constant_has_zero_trace(self):
        trace = g.negative_trace(g.assemble(['-1/100'], [], True))
        ratio = g.fixed_potential_quotient(trace, ['1/100'])
        self.assertEqual(ratio['lower'], '0')
        self.assertIsNone(ratio['unprojected_constant_mode_obstruction'])

    def test_original_case_LT_ratio_scope_difference(self):
        full = g.fixed_potential_quotient(g.negative_trace(self.models[False]), [7])
        mean = g.fixed_potential_quotient(g.negative_trace(self.models[True]), [7])
        self.assertEqual((full['lower'], mean['lower']), ('27/196', '5/49'))

    def test_LT_wrong_potential_negative_or_zero_rejected(self):
        with self.assertRaises(ValueError): g.fixed_potential_quotient(g.negative_trace(self.models[False]), [6])
        with self.assertRaises(ValueError): g.fixed_potential_quotient(g.negative_trace(self.zero[False]), [0])
        with self.assertRaises(ValueError): g.fixed_potential_quotient(g.negative_trace(g.assemble([0, 1], [], False)), [0, -1])

    def test_nonlinear_projection_interlacing(self):
        full = self.linear[False]['answers']['ordered']
        constrained = self.linear[True]['answers']['ordered']
        proof = g.projection_interlacing(full, constrained)
        self.assertEqual(proof['k'], 5)
        self.assertTrue(g.verify(proof, expected_q=[0, 1]))
        self.assertTrue(proof['nonconstant_q_does_not_have_an_added_constant_eigenvalue'])

    def test_interlacing_wrong_q_and_scope_rejected(self):
        full = self.linear[False]['answers']['ordered']
        constrained = self.linear[True]['answers']['ordered']
        with self.assertRaises(ValueError): g.projection_interlacing(constrained, full)
        wrong = g.ordered_merge(g.assemble([1], [], True), 5)
        with self.assertRaises(ValueError): g.projection_interlacing(full, wrong)

    def test_driver_cache_reuse_and_full_input_binding(self):
        cache = {}
        first = g.research_driver([-7], cache=cache)
        second = g.research_driver([-7], cache=cache)
        self.assertEqual(second['execution']['cache_hits'], 2)
        self.assertTrue(same_json(first['certificate'], second['certificate']))
        changed = g.research_driver([-6], cache=cache)
        self.assertEqual(changed['execution']['cache_hits'], 0)

    def test_valid_cache_for_wrong_scope_is_rejected(self):
        cache = {}
        full = g.research_driver([-7], spaces=['full_unprojected'], cache=cache)
        key = next(iter(cache))
        mean = g.adaptive_spectrum([-7], mean_zero=True)['certificate']
        cache[key] = mean
        result = g.research_driver([-7], spaces=['full_unprojected'], cache=cache)
        self.assertEqual(result['execution']['cache_rejections'], 1)
        self.assertEqual(result['certificate']['results'][0]['answers']['trace']['lower'], '27')

    def test_malformed_cache_is_rejected_and_recomputed(self):
        cache = {}
        g.research_driver([-7], spaces=['full_unprojected'], cache=cache)
        cache[next(iter(cache))]['answers']['trace']['upper'] = '20'
        result = g.research_driver([-7], spaces=['full_unprojected'], cache=cache)
        self.assertEqual(result['execution']['cache_rejections'], 1)
        self.assertTrue(g.verify(result['certificate']))

    def test_valid_other_certificate_kind_in_cache_is_rejected(self):
        cache = {}
        g.research_driver([-7], spaces=['full_unprojected'], cache=cache)
        cache[next(iter(cache))] = g.analytic_bounds([-7], 9, False)
        result = g.research_driver([-7], spaces=['full_unprojected'], cache=cache)
        self.assertEqual(result['execution']['cache_rejections'], 1)
        self.assertTrue(g.verify(result['certificate']))

    def test_driver_quantity_decisions_are_not_collapsed(self):
        result = g.research_driver([-20, 0, 1], spaces=['full_mean_zero'],
            k=1, quantities=['ordered', 'trace'], tolerance='1/1000000', bits=32, max_steps=2)
        self.assertFalse(result['certificate']['all_requested_targets_met'])
        self.assertEqual(result['certificate']['decisions']['full_mean_zero']['trace'], 'certified_bound_open')

    def test_l0_removal_metadata_tampering_rejected(self):
        proof = g.analytic_bounds([0], 9, False)
        proof['l0_multiplicity'] = 0
        self.assertFalse(g.verify(proof))
        proof = g.ordered_merge(self.zero[False], 9)
        proof['slots'] = [s for s in proof['slots'] if s['l'] != 0]
        self.assertFalse(g.verify(proof))

    def test_projection_and_tail_tampering_rejected(self):
        proof = g.negative_trace(self.models[False])
        proof['model']['mean_zero'] = True
        self.assertFalse(g.verify(proof))
        proof = g.negative_trace(self.models[False])
        proof['tails'].pop('unlisted_floor')
        self.assertFalse(g.verify(proof))

    def test_external_q_scope_k_binding_and_json_replay(self):
        proof = g.analytic_bounds([-7], 9, False)
        self.assertTrue(g.verify(json.loads(json.dumps(proof)), expected_q=[-7], expected_mean_zero=False, expected_k=9))
        self.assertFalse(g.verify(proof, expected_q=[-6]))
        self.assertFalse(g.verify(proof, expected_mean_zero=True))
        self.assertFalse(g.verify(proof, expected_k=8))


def same_json(a, b):
    return json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


if __name__ == '__main__':
    unittest.main()
