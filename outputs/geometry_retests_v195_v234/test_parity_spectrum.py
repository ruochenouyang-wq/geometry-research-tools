import copy
import json
import unittest
from fractions import Fraction as F
import parity_spectrum as p


C18 = {'2,0,0': 1, '0,2,0': 2, '0,0,2': 3}
QUARTIC = {'4,0,0': '1/2', '2,2,0': '-1/3', '0,0,2': 1}


class ParitySpectrumTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.block = p.block_ground(C18, (1, 0, 0), L=5, bits=40)
        cls.result = p.adaptive_ground(C18)

    def test_v225_original_quadratic_range(self):
        c = p.potential_range(C18)
        self.assertEqual((c['lower'], c['upper']), ('1', '3'))
        self.assertTrue(p.verify_range(c, C18))
        self.assertEqual(len(c['Bernstein_coefficients']), 3)

    def test_v225_constant_on_sphere_range(self):
        q = {'2,0,0': 1, '0,2,0': 1, '0,0,2': 1}
        c = p.potential_range(q)
        self.assertEqual((c['lower'], c['upper']), ('1', '1'))

    def test_v225_quartic_bernstein_cross_coefficient(self):
        c = p.potential_range({'2,2,0': 1})
        self.assertEqual((c['lower'], c['upper']), ('0', '1/2'))
        self.assertEqual(c['simplex_degree'], 2)

    def test_v225_corrupt_range_or_identity_rejected(self):
        c = p.potential_range(C18)
        bad = copy.deepcopy(c)
        bad['lower'] = '2'
        self.assertFalse(p.verify_range(bad))
        bad = copy.deepcopy(c)
        bad['sphere_identity']['quotient']['0,0,0'] = '99'
        self.assertFalse(p.verify_range(bad))
        self.assertFalse(p.verify_range(c, {'2,0,0': 1}))

    def test_v226_all_sign_flip_invariance(self):
        for q in (C18, QUARTIC, {}):
            proof = p.reflection_symmetry(q)
            self.assertTrue(proof['all_eight_reflections'])
            self.assertEqual(len(proof['checks']), 8)

    def test_v226_partial_symmetry_is_not_eight_blocks(self):
        proof = p.reflection_symmetry({'1,1,0': 1})
        self.assertFalse(proof['all_eight_reflections'])
        self.assertEqual(sum(c['invariant'] for c in proof['checks']), 4)

    def test_v226_symmetry_respects_exact_sphere_identity(self):
        # x(r²-1)+x² is the same sphere potential as x².
        q = {'3,0,0': 1, '1,2,0': 1, '1,0,2': 1, '1,0,0': -1, '2,0,0': 1}
        self.assertTrue(p.reflection_symmetry(q)['all_eight_reflections'])
        self.assertEqual(p.potential_range(q)['effective_potential'], {'2,0,0': '1'})

    def test_v226_projectors_are_complete_and_idempotent(self):
        polynomial = p.a.polynomial({'0,0,0': 1, '1,0,0': 2, '0,1,1': 3, '1,1,1': 4, '2,2,0': 5})
        pieces = [p.parity_project(polynomial, character) for character in p.PARITIES]
        self.assertEqual(p.a.add(*pieces), polynomial)
        for character, piece in zip(p.PARITIES, pieces):
            self.assertEqual(p.parity_project(piece, character), piece)

    def test_v227_dimensions_cover_each_entire_degree(self):
        for degree in range(9):
            pieces = [p.degree_basis(character, degree) for character in p.PARITIES]
            self.assertEqual(sum(map(len, pieces)), 2*degree+1)
            for character, basis in zip(p.PARITIES, pieces):
                k = sum(character)
                expected = (degree-k)//2+1 if degree >= k and (degree-k) % 2 == 0 else 0
                self.assertEqual(len(basis), expected)

    def test_v227_zero_character_removes_only_constant(self):
        b = p.parity_basis((0, 0, 0), 4)
        self.assertEqual([v['degree'] for v in b], [2, 2, 4, 4, 4])
        self.assertTrue(all(p.a.integrate(v['polynomial']) == 0 for v in b))

    def test_v227_higher_basis_energy_independent_integral(self):
        for b in p.degree_basis((1, 0, 0), 9):
            self.assertEqual(p.a.tangent_energy(b['polynomial'], b['polynomial']), 90*b['mass'])

    def test_v228_original_x_trial_is_18_over_5(self):
        f = p.block_form(C18, (1, 0, 0), 1)
        self.assertEqual(f['A'][0][0]/f['mass'][0], F(18, 5))

    def test_v228_all_cross_blocks_vanish_for_even_potential(self):
        for i, character in enumerate(p.PARITIES):
            for other in p.PARITIES[:i]:
                self.assertTrue(all(v == 0 for row in p.cross_block_form(QUARTIC, character, other) for v in row))

    def test_v228_cross_block_witness_exposes_xy(self):
        matrix = p.cross_block_form({'1,1,0': 1}, (1, 0, 0), (0, 1, 0), 1)
        self.assertEqual(matrix, [[F(1, 15)]])

    def test_v229_actual_next_degree_uses_character(self):
        self.assertEqual(p.next_omitted_degree((1, 0, 0), 4), 5)
        self.assertEqual(p.next_omitted_degree((1, 0, 0), 3), 5)
        self.assertEqual(p.next_omitted_degree((0, 0, 0), 3), 4)
        self.assertEqual(p.next_omitted_degree((1, 1, 1), 3), 5)

    def test_v229_quartic_tail_keeps_every_same_character(self):
        tail = p.tail_couplings(QUARTIC, (1, 0, 0), 3)
        self.assertEqual(tail['first_omitted_degree'], 5)
        self.assertEqual([c['degree'] for c in tail['columns']], [5]*3+[7]*4)

    def test_v229_far_tail_has_zero_exact_coupling(self):
        finite = p.parity_basis((1, 0, 0), 3)
        far = p.degree_basis((1, 0, 0), 9)
        q = p.a.polynomial(QUARTIC)
        self.assertTrue(all(p.a.inner(b['polynomial'], p.a.multiply(q, c['polynomial'])) == 0
                            for b in far for c in finite))

    def test_v230_unperturbed_all_eight_character_floors(self):
        for character in p.PARITIES:
            l = 2 if sum(character) == 0 else sum(character)
            cert = p.block_ground({}, character, L=max(3, l))
            self.assertEqual((cert['lower'], cert['upper']), (str(l*(l+1)),)*2)
            self.assertTrue(p.verify_block(cert, {}, character))

    def test_v230_negative_constant_shift(self):
        cert = p.block_ground({'0,0,0': -5}, (1, 0, 0), L=1)
        self.assertEqual((cert['lower'], cert['upper']), ('-3', '-3'))

    def test_v230_mass_tail_and_endpoint_mutations_fail(self):
        for change in ('mass', 'tail_mass', 'tail_drop', 'tail_start', 'lower'):
            bad = copy.deepcopy(self.block)
            if change == 'mass':
                bad['mass'][0] = '1'
            elif change == 'tail_mass':
                bad['complete_tail']['columns'][0]['mass'] = '1'
            elif change == 'tail_drop':
                bad['complete_tail']['columns'].pop()
            elif change == 'tail_start':
                bad['complete_tail']['first_omitted_degree'] = 8
            else:
                bad['lower'] = '4'
            self.assertFalse(p.verify_block(bad), change)

    def test_v230_input_character_and_scope_bound(self):
        self.assertFalse(p.verify_block(self.block, {'2,0,0': 1}))
        self.assertFalse(p.verify_block(self.block, C18, (0, 1, 0)))
        bad = copy.deepcopy(self.block)
        bad['scope'] = 'all_H1'
        self.assertFalse(p.verify_block(bad))

    def test_v231_all_characters_present_in_final_result(self):
        cert = self.result['certificate']
        self.assertEqual([tuple(c['parity']) for c in cert['characters']], list(p.PARITIES))
        self.assertEqual(sum(c['format'] == p.BLOCK_FORMAT for c in cert['characters']), 3)
        self.assertTrue(p.verify_full(cert, C18))

    def test_v231_missing_or_duplicate_character_fails(self):
        bad = copy.deepcopy(self.result['certificate'])
        bad['characters'].pop()
        self.assertFalse(p.verify_full(bad))
        bad = copy.deepcopy(self.result['certificate'])
        bad['characters'][0] = copy.deepcopy(bad['characters'][1])
        self.assertFalse(p.verify_full(bad))

    def test_v231_analytic_floor_cannot_be_invented(self):
        bad = copy.deepcopy(self.result['certificate'])
        bad['characters'][0]['lower'] = '100'
        self.assertFalse(p.verify_full(bad))

    def test_v232_tail_diagnosis_uses_actual_retests(self):
        diagnosis = p.error_diagnosis(C18, (1, 0, 0), L=3, bits=32, extra_bits=16)
        self.assertEqual(diagnosis['diagnosis'], 'tail_dominated')
        self.assertGreater(F(diagnosis['fixed_space_width_ratio']), F(99, 100))
        self.assertLess(F(diagnosis['larger_space_width_ratio']), F(1, 1000))

    def test_v232_only_the_competitor_was_grown(self):
        actions = self.result['trace'][1:]
        self.assertEqual(len(actions), 1)
        self.assertEqual(actions[0]['parity'], [1, 0, 0])
        self.assertEqual(actions[0]['retained_degree'], 5)
        self.assertEqual(actions[0]['dimension'], 6)

    def test_v233_cache_is_reused_and_bound_to_inputs(self):
        cache = p.FormCache()
        one = p.block_ground(C18, (1, 0, 0), 3, 32, cache)
        two = p.block_ground(C18, (1, 0, 0), 3, 48, cache)
        self.assertEqual(cache.stats['form_builds'], 1)
        self.assertEqual(cache.stats['form_hits'], 1)
        p.block_ground(C18, (0, 1, 0), 3, 24, cache)
        p.block_ground({'2,0,0': 1}, (1, 0, 0), 3, 24, cache)
        self.assertEqual(cache.stats['form_builds'], 3)
        self.assertTrue(p.verify_block(one))
        self.assertTrue(p.verify_block(two))

    def test_v233_returned_cache_copy_cannot_poison_future_form(self):
        cache = p.FormCache()
        finite, _ = cache.form(C18, (1, 0, 0), 3)
        original = finite['A'][0][0]
        finite['A'][0][0] = F(999)
        again, _ = cache.form(C18, (1, 0, 0), 3)
        self.assertEqual(again['A'][0][0], original)

    def test_v234_original_goal_met_without_changing_input_or_width(self):
        r = self.result
        self.assertTrue(p.verify_result(r, C18, F(1, 10**8)))
        self.assertEqual(r['status'], 'target_met')
        self.assertLessEqual(F(r['certificate']['exact_width']), F(1, 10**8))
        self.assertGreaterEqual(F(r['certificate']['lower']), 3)
        self.assertLessEqual(F(r['certificate']['upper']), F(18, 5))
        self.assertTrue(r['within_wall_budget'])

    def test_v234_even_quartic_is_also_accepted(self):
        r = p.adaptive_ground(QUARTIC, target=F(1, 10**7), start_L=3, max_L=7)
        self.assertTrue(p.verify_result(r, QUARTIC, F(1, 10**7)))
        self.assertEqual(r['status'], 'target_met')

    def test_v234_unknown_symmetry_rejected_by_all_spectral_routes(self):
        xy = {'1,1,0': 1}
        for function in (lambda: p.block_ground(xy, (1, 0, 0)),
                         lambda: p.full_ground(xy), lambda: p.adaptive_ground(xy)):
            with self.assertRaises(ValueError):
                function()

    def test_v234_insufficient_degree_budget_remains_open(self):
        r = p.adaptive_ground(C18, start_L=3, max_L=3)
        self.assertTrue(p.verify_result(r, C18))
        self.assertEqual(r['status'], 'certified_open_gap')
        self.assertEqual(r['limit_reason'], 'retained_degree_cap_reached')

    def test_v234_initial_coverage_respects_degree_budget(self):
        r = p.adaptive_ground({'4,0,0': -100}, start_L=1, max_L=1, bits=16)
        self.assertTrue(p.verify_result(r))
        self.assertEqual(r['status'], 'certified_open_gap')
        self.assertTrue(all(e['retained_degree'] <= 1 for e in r['certificate']['characters']
                            if e['format'] == p.BLOCK_FORMAT))

    def test_v234_driver_request_and_status_tampering_rejected(self):
        self.assertFalse(p.verify_result(self.result, C18, F(1, 10**10)))
        bad = copy.deepcopy(self.result)
        bad['requested_width'] = '1/10000000000'
        self.assertFalse(p.verify_result(bad))
        bad = copy.deepcopy(self.result)
        bad['status'] = 'certified_open_gap'
        self.assertFalse(p.verify_result(bad))

    def test_v234_json_roundtrip(self):
        self.assertTrue(p.verify_result(json.loads(json.dumps(self.result)), C18))

    def test_reject_float_coefficients_boolean_parity_and_wrong_sizes(self):
        for function in (lambda: p.potential_range({'2,0,0': 0.5}),
                         lambda: p.block_ground(C18, (True, 0, 0)),
                         lambda: p.adaptive_ground(C18, start_L=3, max_L=2),
                         lambda: p.block_ground(C18, (1, 0, 0), bits=True)):
            with self.assertRaises(ValueError):
                function()


if __name__ == '__main__':
    unittest.main()
