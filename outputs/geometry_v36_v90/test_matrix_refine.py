"""Meaningful exact-proof, conditioning, failure, and continuation tests."""
import copy
import unittest
from fractions import Fraction as F
from unittest.mock import patch
import algebra as a
import dual_v33 as dual
import family_v29 as family
import matrix_refine as refine
import optimize_v34 as old


DIRS = [[0, 1], [0, 0, 1]]
COST = [[2, F(1, 2)], [F(1, 2), 1]]


def fractions(matrix):
    return [list(map(F, row)) for row in matrix]


class CoordinateTests(unittest.TestCase):
    def test_nonsymmetric_transform_and_nondiagonal_cost(self):
        transform = [[2, 1], [F(1, 3), -1]]
        directions, cost, inverse = refine.transform_problem(DIRS, COST, transform)
        _, old_v, _ = dual.setup(DIRS, COST)
        _, new_v, _ = dual.setup(directions, cost)
        for t in [F(-1), F(-2, 7), F(0), F(3, 5), F(1)]:
            self.assertEqual(dual.at(new_v, t), refine.congruence(transform, dual.at(old_v, t)))
        g = [[3, F(1, 5)], [F(1, 5), 2]]
        h = refine.congruence(transform, g)
        self.assertEqual(dual.inner(COST, g), dual.inner(cost, h))
        self.assertEqual(refine.congruence(inverse, h), g)
        wrong_cost = refine.congruence(transform, COST)
        self.assertNotEqual(dual.inner(wrong_cost, h), dual.inner(COST, g))

    def test_dual_coordinate_pullback(self):
        transform = [[2, 1], [0, F(1, 2)]]
        dirs, c, _ = refine.transform_problem(DIRS, COST, transform)
        atom = {'t': '1/3', 'matrix': dual.encode(refine._scale(c, F(1, 2)))}
        transformed = dual.lower_certificate(dirs, c, [atom])
        pulled = {'t': atom['t'], 'matrix': dual.encode(refine.congruence(a.transpose(transform), fractions(atom['matrix'])))}
        original = dual.lower_certificate(DIRS, COST, [pulled])
        self.assertEqual(transformed['lower_bound'], original['lower_bound'])
        self.assertTrue(dual.verify_lower(original))

    def test_singular_and_wrong_size_transform_rejected(self):
        for transform in [[[1, 2], [2, 4]], [[1]], [[1, 0, 2], [0, 1, 3]]]:
            with self.assertRaises(ValueError):
                refine.transform_problem(DIRS, COST, transform)

    def test_exact_decorrelation(self):
        directions = [[0, 1000, 1000], [0, 0, F(1, 1000)]]
        _, v, _ = dual.setup(directions)
        transform = refine.preconditioner(v)
        dirs, _, _ = refine.transform_problem(directions, None, transform)
        _, new_v, _ = dual.setup(dirs)
        width = max(map(len, new_v))
        rows = [p + [F(0)] * (width - len(p)) for p in new_v]
        gram = a.matmul(rows, a.transpose(rows))
        self.assertEqual(gram[0][1], 0)
        self.assertTrue(all(F(1, 4) <= gram[i][i] <= 4 for i in range(2)))


class ConstraintAndDualTests(unittest.TestCase):
    def test_scalar_dominance_is_exact(self):
        _, v, _ = dual.setup([[0, 1]])
        samples = [F(i, 4) for i in range(-4, 5)]
        compressed = refine.effective_samples(v, samples)
        self.assertEqual([t for t, _ in compressed], [F(0)])
        for t in samples:
            self.assertEqual(a.inertia(refine._add(compressed[0][1], refine._scale(dual.at(v, t), -1)))[0], 0)

    def test_merging_identical_a_preserves_objective_and_sum(self):
        atoms = [{'t': '-1/2', 'matrix': [['1/3']]}, {'t': '1/2', 'matrix': [['1/4']]}]
        before = dual.lower_certificate([[0, 1]], [[1]], atoms)
        after = dual.lower_certificate([[0, 1]], [[1]], refine.merge_atoms([[0, 1]], atoms))
        self.assertEqual(len(after['atoms']), 1)
        self.assertEqual(before['lower_bound'], after['lower_bound'])
        self.assertEqual(before['cost_slack'], after['cost_slack'])

    def test_upward_repair_improves_old_capped_scaling(self):
        atoms = [{'t': '0', 'matrix': [['1/100']]}]
        old_lower = dual.repair_dual([[0, 1]], [[1]], atoms)
        new_lower = refine.repair_dual([[0, 1]], [[1]], atoms)
        self.assertEqual(new_lower['lower_bound'], '1/4')
        self.assertGreater(F(new_lower['lower_bound']), 99 * F(old_lower['lower_bound']))
        self.assertTrue(dual.verify_lower(new_lower))

    def test_huge_common_atom_scale_does_not_exhaust_bisection(self):
        atoms = [{'t': '0', 'matrix': [[10**100]]}]
        cert = refine.repair_dual([[0, 1]], [[1]], atoms, F(1, 10**10), fill=False)
        self.assertEqual(cert['lower_bound'], '1/4')
        self.assertTrue(dual.verify_lower(cert))

    def test_nondiagonal_cost_downward_repair_and_slack_fill(self):
        atoms = [{'t': '-1/3', 'matrix': [[4, 2], [2, 3]]},
                 {'t': '2/3', 'matrix': [[2, -1], [-1, 2]]}]
        before_fill = refine.repair_dual(DIRS, COST, atoms, fill=False)
        after_fill = refine.fill_slack(DIRS, COST, before_fill, [F(i, 8) for i in range(-8, 9)])
        self.assertTrue(dual.verify_lower(after_fill))
        self.assertEqual(fractions(after_fill['cost_slack']), [[0, 0], [0, 0]])
        self.assertGreaterEqual(F(after_fill['lower_bound']), F(before_fill['lower_bound']))

    def test_negative_atom_rejected(self):
        with self.assertRaises(ValueError):
            refine.repair_dual([[0, 1]], [[1]], [{'t': '0', 'matrix': [[-1]]}])

    def test_slack_fill_is_bound_to_cost_and_directions(self):
        lower = dual.lower_certificate([[0, 1]], [[1]], [{'t': '0', 'matrix': [[F(1, 2)]]}])
        for dirs, cost in [([[0, 2]], [[1]]), ([[0, 1]], [[2]])]:
            with self.assertRaises(ValueError):
                refine.fill_slack(dirs, cost, lower, [0])

    def test_sparsity_discards_zero_contribution_and_refills(self):
        atoms = [{'t': '0', 'matrix': [['1/2']]}, {'t': '1', 'matrix': [['1/2']]}]
        cert = dual.lower_certificate([[0, 1]], [[1]], atoms)
        sparse, loss = refine.sparsify_dual([[0, 1]], [[1]], cert, 0, [0, 1])
        self.assertEqual(loss, 0)
        self.assertEqual(len(sparse['atoms']), 1)
        self.assertEqual(sparse['lower_bound'], '1/4')

    def test_psd_factor_rounding_handles_rank_deficiency(self):
        matrix = [[F(1, 3), F(1, 5)], [F(1, 5), F(3, 25)]]
        entrywise = [[F(round(x * 100), 100) for x in row] for row in matrix]
        self.assertGreater(a.inertia(entrywise)[0], 0)
        reconstructed = refine.reconstruct_psd(matrix, 1000)
        self.assertEqual(a.inertia(reconstructed)[0], 0)
        self.assertLess(max(abs(matrix[i][j] - reconstructed[i][j]) for i in range(2) for j in range(2)), F(1, 100))

    def test_sample_budget_always_admits_new_separators(self):
        _, v, _ = dual.setup(DIRS)
        lower = dual.lower_certificate(DIRS, COST, [{'t': '0', 'matrix': dual.encode(COST)}])
        old_samples = {F(i, 4) for i in range(-4, 5)}
        new_samples = {F(2, 7), F(-2, 7)}
        kept = refine.update_samples(old_samples, new_samples, v, lower, 9)
        self.assertTrue(new_samples <= kept)
        self.assertEqual(len(kept), 9)
        self.assertTrue(dual.verify_lower(lower))


class BudgetAndProofTests(unittest.TestCase):
    def test_error_budget_accounts_for_objective_scale(self):
        one = refine.error_allocation(F(1, 10**6), F(1, 10), 1)
        large = refine.error_allocation(F(1, 10**6), F(1, 10), 1000)
        self.assertEqual(one['range'], 1000 * large['range'])
        later = refine.error_allocation(F(1, 10**6), F(1, 10), 1, 8)
        self.assertLess(later['search'], one['search'])
        self.assertGreaterEqual(later['search'], F(1, 8 * 10**6))

    def test_zero_round_scalar_exact_proof(self):
        result = refine.refine([[0, 1]], budget=0)
        self.assertEqual(result['certificate']['gap'], '0')
        self.assertTrue(refine.verify(result))
        self.assertTrue(dual.verify(result['certificate']))

    def test_wrong_range_polynomial_is_rejected_and_falls_back(self):
        def wrong(q, epsilon, leaves):
            return family.maximum([F(1)], epsilon, leaves)
        result = refine.refine(DIRS, budget=0, range_evaluator=wrong)
        self.assertTrue(refine.verify(result))
        self.assertEqual(result['trace'][0]['outcome'], 'range_evaluator_fallback')

    def test_injected_legacy_range_evaluator_is_used(self):
        with patch('matrix_refine.family.maximum', wraps=family.maximum) as evaluator:
            result = refine.refine(DIRS, budget=1, range_evaluator=evaluator)
        self.assertGreaterEqual(evaluator.call_count, 2)
        self.assertTrue(refine.verify(result))

    def test_failed_candidate_retains_valid_fallback(self):
        def fail(*args, **kwargs):
            raise ValueError('Injected proposal failure')
        result = refine.refine(DIRS, candidate_solver=fail)
        self.assertTrue(refine.verify(result))
        self.assertEqual(result['trace'][-1]['outcome'], 'proposal_failed')
        self.assertEqual(result['certificate']['status'], 'global_gap_open')

    def test_weak_candidate_does_not_claim_convergence(self):
        def weak(samples, cost, final_mu, **kwargs):
            return [[100.0, 0.0], [0.0, 100.0]], 0.1, {'newton_steps': 0, 'reached_final_barrier': True}
        result = refine.refine(DIRS, budget=1, candidate_solver=weak)
        self.assertTrue(refine.verify(result))
        self.assertEqual(result['certificate']['status'], 'global_gap_open')

    def test_tampered_bound_is_rejected(self):
        result = refine.refine([[0, 1]], budget=0)
        result['certificate']['lower_bound'] = '1'
        self.assertFalse(refine.verify(result))

    def test_search_metadata_has_no_mathematical_authority(self):
        result = refine.refine([[0, 1]], budget=0)
        result['statistics'] = {'claimed_spectral_optimality': True}
        result['conditioning'] = {'transformed_cost': 'bad'}
        self.assertTrue(refine.verify(result))
        self.assertIn('not the sharp spectral inequality', result['certificate']['limitation'])

    def test_invalid_budgets(self):
        for budget in [-1, 25, {'range_leaves': 0}, {'newton_steps': True}, {'max_samples': 8}, {'unknown': 1}, {'precondition': 1}]:
            with self.assertRaises(ValueError):
                refine.refine(DIRS, budget=budget)
        for tolerance in [F(0), F(1, 10**11), F(1, 2)]:
            with self.assertRaises(ValueError):
                refine.refine(DIRS, tolerance=tolerance)


class SearchAndContinuationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.seed = refine.refine(DIRS, tolerance=F(1, 10**5), budget=1)
        cls.refined = refine.refine(DIRS, tolerance=F(1, 10**6), seed_result=cls.seed, budget=8)

    def test_actual_error_improvement(self):
        self.assertTrue(refine.verify(self.seed))
        self.assertTrue(refine.verify(self.refined))
        self.assertLess(F(self.refined['certificate']['gap']), F(self.seed['certificate']['gap']) / 100)
        self.assertLessEqual(F(self.refined['certificate']['gap']), F(1, 10**6))

    def test_seed_bounds_are_independently_monotone(self):
        self.assertGreaterEqual(F(self.refined['certificate']['lower_bound']), F(self.seed['certificate']['lower_bound']))
        self.assertLessEqual(F(self.refined['certificate']['upper_bound']), F(self.seed['certificate']['upper_bound']))
        self.assertTrue(self.refined['statistics']['seed_reused'])
        self.assertTrue(all(row['warm_start'] for row in self.refined['trace'] if row['outcome'] == 'certified'))

    def test_zero_work_reuses_both_exact_bounds(self):
        result = refine.refine(DIRS, tolerance=F(1, 10**7), seed_result=self.seed, budget=0)
        self.assertEqual(result['certificate']['lower_bound'], self.seed['certificate']['lower_bound'])
        self.assertEqual(result['certificate']['upper_bound'], self.seed['certificate']['upper_bound'])

    def test_seed_survives_failed_continuation(self):
        def fail(*args, **kwargs):
            raise ValueError('Continuation failed')
        result = refine.refine(DIRS, seed_result=self.seed, candidate_solver=fail)
        self.assertEqual(result['certificate']['lower_bound'], self.seed['certificate']['lower_bound'])
        self.assertEqual(result['certificate']['upper_bound'], self.seed['certificate']['upper_bound'])
        self.assertTrue(refine.verify(result))

    def test_corrupt_optional_warm_cache_is_ignored(self):
        for cached in [None, 4, {'samples': 9}, {'samples': {'bad': 'cache'}}]:
            seed = copy.deepcopy(self.seed)
            seed['warm_state'] = cached
            result = refine.refine(DIRS, seed_result=seed, budget=0)
            self.assertEqual(result['certificate']['gap'], self.seed['certificate']['gap'])

    def test_internal_preconditioning_does_not_reapply_input_denominator_limit(self):
        cases = [
            [[0, F(1, 97), F(1, 89), F(1, 83)], [0, F(1, 79), F(1, 73), F(1, 71)]],
            [[0, 1, F(1, 999983)], [0, F(1, 999979), 1]],
        ]
        for dirs in cases:
            template = refine._template(dirs, [[1, 0], [0, 1]], F(1, 1000), 4, family.maximum)
            lower = dual.lower_certificate(dirs, [[1, 0], [0, 1]], [{'t': '0', 'matrix': [[0, 0], [0, 0]]}])
            seed = {'certificate': dual.certificate(template, lower, F(1, 1000))}
            result = refine.refine(dirs, seed_result=seed, budget=0)
            self.assertTrue(refine.verify(result))
            self.assertEqual(result['certificate']['upper_bound'], seed['certificate']['upper_bound'])
            fallback = refine.refine(dirs, budget={'rounds': 0, 'range_leaves': 4})
            self.assertTrue(refine.verify(fallback))

    def test_seed_rejects_wrong_problem_and_tampering(self):
        for dirs, cost in [(DIRS, COST), ([[0, 2], [0, 0, 1]], None)]:
            with self.assertRaises(ValueError):
                refine.refine(dirs, cost, seed_result=self.seed)
        tampered = copy.deepcopy(self.seed)
        tampered['certificate']['template']['scale'] = '0'
        with self.assertRaises(ValueError):
            refine.refine(DIRS, seed_result=tampered)

    def test_nondiagonal_cost_closes(self):
        result = refine.refine(DIRS, COST, F(1, 10**6), budget=8)
        self.assertTrue(refine.verify(result))
        self.assertLessEqual(F(result['certificate']['gap']), F(1, 10**6))
        self.assertEqual(result['certificate']['dual']['cost'], dual.encode(COST))

    def test_four_directions_close(self):
        dirs = [[0] * i + [1] for i in range(1, 5)]
        result = refine.refine(dirs, tolerance=F(1, 10**5), budget=8)
        self.assertTrue(refine.verify(result))
        self.assertLessEqual(F(result['certificate']['gap']), F(1, 10**5))
        self.assertEqual(result['statistics']['model_calls'], 0)

    def test_scaled_coordinates_improve_v34_at_same_round_budget(self):
        dirs = [[0, 1000], [0, 0, F(1, 1000)]]
        cost = [[F(1, 10**6), 0], [0, 10**6]]
        baseline = old.solve(dirs, cost, F(1, 10**5), max_rounds=5)
        result = refine.refine(dirs, cost, F(1, 10**5), budget=5)
        self.assertTrue(dual.verify(baseline['certificate']))
        self.assertTrue(refine.verify(result))
        self.assertLess(F(result['certificate']['gap']), F(baseline['certificate']['gap']) / 100)
        self.assertLessEqual(F(result['certificate']['gap']), F(1, 10**5))


if __name__ == '__main__':
    unittest.main()
