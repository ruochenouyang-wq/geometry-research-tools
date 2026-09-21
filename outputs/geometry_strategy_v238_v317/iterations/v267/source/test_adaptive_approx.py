from fractions import Fraction as F
import unittest
import adaptive_approx as adaptive

SINGULAR = {'kind':'axis_profile','profile':'abs_power','exponent':'-1/4','amplitude':'-1'}
STEP = {'kind':'axis_profile','profile':'step','amplitude':'3','offset':'-2'}


class ApproximationTests(unittest.TestCase):
    def test_decomposition_matches_full_integrated_certificate(self):
        for ratio in ['1/2', '2/3', '3/4']:
            terms = adaptive.error_terms(3, 5, ratio, '-3/2')
            function = dict(SINGULAR, amplitude='-3/2', offset='9/7', axis=1)
            cert = adaptive.pieces.piecewise_model(function, 5, 3, ratio)
            self.assertEqual(F(terms['error_squared']), F(cert['error_squared']))
            self.assertEqual(F(terms['core_error_squared'])+F(terms['annulus_error_squared']),
                             F(terms['error_squared']))

    def test_step_is_exact_and_bound_to_input(self):
        result = adaptive.solve(STEP)
        self.assertEqual(result['status'], 'met')
        self.assertEqual(result['certificate']['error_squared'], '0')
        self.assertTrue(adaptive.pieces.verify_piecewise(result['certificate'], STEP))

    def test_invalid_budget_and_target_are_rejected(self):
        for kwargs in [{'tolerance':'0'}, {'max_levels':0}, {'max_degree':True}, {'root_ratio':'1'}]:
            with self.assertRaises(ValueError):
                adaptive.solve(SINGULAR, **kwargs)

    def test_level_selection_is_minimal_not_rounded_logarithm(self):
        result = adaptive.solve(SINGULAR, tolerance='3/100', root_ratio='1/2')
        cert = result['certificate']
        self.assertEqual(result['status'], 'met')
        self.assertGreater(cert['levels'], 4)
        previous = adaptive.error_terms(cert['degree'], cert['levels']-1, cert['root_ratio'])
        self.assertGreater(F(previous['error_squared']), F(3,100)**2)
        self.assertLessEqual(result['cost']['level_comparisons'], 7)

    def test_unreachable_fixed_degree_is_proved_before_level_search(self):
        result = adaptive.solve(SINGULAR, max_degree=3, root_ratio='1/2')
        terms = result['diagnostics']['error_decomposition']
        self.assertEqual(result['status'], 'open')
        self.assertEqual(terms['limiting_component'], 'fixed_degree_floor')
        self.assertGreater(F(terms['infinite_levels_error_floor_squared']), F(1,10**8)**2)
        self.assertEqual(result['cost']['level_comparisons'], 0)

    def test_degree_increases_only_until_exact_target_is_attainable(self):
        result = adaptive.solve(SINGULAR, tolerance='1/10000', root_ratio='1/2')
        self.assertEqual(result['status'], 'met')
        self.assertGreater(result['certificate']['degree'], 3)
        self.assertLess(result['certificate']['degree'], 24)
        self.assertEqual(result['cost']['certificate_builds'], 1)
        self.assertTrue(all(a['limiting_component'] == 'fixed_degree_floor'
                            for a in result['attempts'][:-1]))

    def test_ratio_change_escapes_the_single_ratio_precision_floor(self):
        fixed = adaptive.solve(SINGULAR, root_ratio='1/2')
        adaptive_result = adaptive.solve(SINGULAR)
        self.assertEqual(fixed['status'], 'open')
        self.assertEqual(adaptive_result['status'], 'met')
        self.assertNotEqual(adaptive_result['certificate']['root_ratio'], '1/2')
        self.assertEqual(adaptive_result['cost']['certificate_builds'], 1)

    def test_square_root_rounding_is_refined_without_changing_the_model(self):
        result = adaptive.solve(SINGULAR, tolerance='3/100', root_ratio='1/2', sqrt_bits=0)
        cert = result['certificate']
        self.assertEqual(result['status'], 'met')
        self.assertGreater(cert['sqrt_bits'], 0)
        self.assertGreater(adaptive.pieces.sqrt_upper(F(cert['error_squared']), cert['sqrt_bits']-1), F(3,100))
        self.assertEqual(result['cost']['certificate_builds'], 1)

    def test_square_root_budget_is_not_mistaken_for_mathematical_error(self):
        result = adaptive.solve(SINGULAR, tolerance='3/100', root_ratio='1/2',
                                sqrt_bits=0, max_sqrt_bits=0)
        self.assertEqual(result['status'], 'open')
        self.assertEqual(result['diagnostics']['rounding']['reason'], 'sqrt_bit_budget')
        self.assertLess(F(result['certificate']['error_squared']), F(3,100)**2)

    def test_zero_amplitude_has_no_shape_search_and_keeps_source_identity(self):
        source = dict(SINGULAR, amplitude='0', offset='17/11', axis=0)
        result = adaptive.solve(source, tolerance='1/100000000000000000000')
        self.assertEqual(result['status'], 'met')
        self.assertEqual(result['cost']['shape_evaluations'], 0)
        self.assertEqual(result['attempts'], [])
        self.assertTrue(adaptive.pieces.verify_piecewise(result['certificate'], source))

    def test_amplitude_normalizes_target_without_discarding_offset_or_axis(self):
        unit = adaptive.solve(SINGULAR, tolerance='1/10000', root_ratio='1/2')
        source = dict(SINGULAR, amplitude='3/7', offset='29', axis=0)
        scaled = adaptive.solve(source, tolerance='3/70000', root_ratio='1/2')
        for field in ['degree', 'levels', 'root_ratio']:
            self.assertEqual(unit['certificate'][field], scaled['certificate'][field])
        self.assertEqual(F(scaled['certificate']['error_squared']), F(9,49)*F(unit['certificate']['error_squared']))
        self.assertTrue(adaptive.pieces.verify_piecewise(scaled['certificate'], source))

    def test_explicit_shape_cache_reuses_scaled_offset_and_axis_problems(self):
        cache = adaptive.ShapeCache()
        cold = adaptive.solve(SINGULAR, tolerance='1/10000', root_ratio='1/2', shape_cache=cache)
        source = dict(SINGULAR, amplitude='2', offset='17', axis=1)
        warm = adaptive.solve(source, tolerance='1/5000', root_ratio='1/2', shape_cache=cache)
        self.assertGreater(cold['cost']['shape_evaluations'], 0)
        self.assertEqual(warm['cost']['shape_evaluations'], 0)
        self.assertEqual(warm['cost']['shape_cache_hits'], cold['cost']['shape_evaluations'])
        self.assertEqual(warm['cost']['certificate_replays'], 1)
        self.assertTrue(adaptive.pieces.verify_piecewise(warm['certificate'], source))

    def test_shape_cache_is_bounded_and_explicit(self):
        cache = adaptive.ShapeCache(capacity=1)
        cost = {'shape_evaluations':0}
        adaptive._shape_error(0, F(1,2), cache, cost)
        adaptive._shape_error(1, F(1,2), cache, cost)
        adaptive._shape_error(0, F(1,2), cache, cost)
        self.assertEqual(cost['shape_evaluations'], 3)
        self.assertEqual(len(cache._values), 1)
        cache.clear()
        self.assertEqual(len(cache._values), 0)

    def test_budget_exhaustion_preserves_best_examined_candidate(self):
        result = adaptive.solve(SINGULAR, max_shape_evaluations=2)
        self.assertEqual(result['status'], 'open')
        self.assertEqual(result['cost']['shape_evaluations'], 2)
        self.assertEqual(result['diagnostics']['stopping_reason'], 'shape_evaluation_budget')
        self.assertEqual(F(result['certificate']['error_squared']),
                         min(F(a['error_squared']) for a in result['attempts']))
        self.assertTrue(adaptive.pieces.verify_piecewise(result['certificate'], SINGULAR))

    def test_zero_search_budget_never_manufactures_a_certificate(self):
        result = adaptive.solve(SINGULAR, max_shape_evaluations=0)
        self.assertEqual(result['status'], 'open')
        self.assertIsNone(result['certificate'])
        self.assertEqual(result['cost']['certificate_builds'], 0)
        self.assertFalse(result['diagnostics']['target_proved_unattainable'])
        exact = adaptive.solve(STEP, max_shape_evaluations=0)
        self.assertEqual(exact['status'], 'met')

    def test_certificate_slot_budget_is_enforced(self):
        result = adaptive.solve(SINGULAR, max_certificate_slots=13)
        self.assertLessEqual(result['certificate']['expanded_polynomial_coefficient_slots'], 13)
        self.assertTrue(all(a['expanded_coefficient_slots'] <= 13 for a in result['attempts']))

    def test_warm_search_can_run_with_zero_new_shape_budget(self):
        cache = adaptive.ShapeCache()
        first = adaptive.solve(SINGULAR, tolerance='1/10000', root_ratio='1/2', shape_cache=cache)
        warm = adaptive.solve(SINGULAR, tolerance='1/10000', root_ratio='1/2',
                              shape_cache=cache, max_shape_evaluations=0)
        self.assertEqual(first['certificate'], warm['certificate'])
        self.assertEqual(warm['status'], 'met')
        self.assertEqual(warm['cost']['shape_evaluations'], 0)

    def test_minimum_slots_matches_exhaustive_old_certificate_comparison(self):
        ratios = ['1/2', '2/3', '3/4']
        target = F(1,100)
        result = adaptive.solve(SINGULAR, tolerance=str(target), max_degree=6, max_levels=16,
                                root_ratios=ratios, sqrt_bits=40, max_sqrt_bits=40)
        feasible_slots = []
        for ratio in ratios:
            for degree in range(7):
                for levels in range(1,17):
                    cert = adaptive.pieces.piecewise_model(SINGULAR, levels, degree, ratio, 40)
                    if F(cert['error_upper']) <= target:
                        feasible_slots.append(cert['expanded_polynomial_coefficient_slots'])
        self.assertEqual(result['status'], 'met')
        self.assertEqual(result['certificate']['expanded_polynomial_coefficient_slots'], min(feasible_slots))
        self.assertTrue(result['diagnostics']['minimum_slots_within_configured_family'])
        self.assertGreater(result['cost']['shape_candidates_pruned_by_cost_bound'], 0)

    def test_cost_optimality_is_not_claimed_after_budget_exhaustion(self):
        result = adaptive.solve(SINGULAR, tolerance='1/100', max_shape_evaluations=6)
        self.assertFalse(result['diagnostics']['minimum_slots_within_configured_family'])

    def test_rounding_constraint_enters_model_selection(self):
        result = adaptive.solve(SINGULAR, tolerance='3/100', sqrt_bits=6, max_sqrt_bits=6)
        self.assertEqual(result['status'], 'met')
        self.assertEqual(F(result['diagnostics']['reportable_target_at_max_bits']), F(1,64))
        self.assertLessEqual(F(result['certificate']['error_upper']), F(3,100))


if __name__ == '__main__':
    unittest.main()
