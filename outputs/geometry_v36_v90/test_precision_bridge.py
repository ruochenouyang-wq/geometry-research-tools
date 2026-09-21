"""V88 integration, binding, failure honesty and deterministic comparisons."""

import copy
from fractions import Fraction as F
import json
import unittest
from unittest.mock import patch

import dual_v33 as dual
import exact_extrema as extrema
import family_v29 as family
import matrix_refine
import positive_search
import precision_bridge as bridge


DIRS = [[0, 1, F(1, 3), 0, F(1, 5)]]


def legacy_seed(directions=DIRS, weight=((1,),), cost=((1,),), leaves=1):
    dirs, w, _, _, q = family.assemble(directions, weight)
    template = family.build(dirs, w, family.maximum(q, F(1, 10**12), leaves),
                            "poisson_matrix_v30")
    lower = dual.lower_certificate(dirs, cost, [{"t": "0", "matrix": dual.encode(cost)}])
    return dual.certificate(template, lower, F(1, 10**10))


class BackendTests(unittest.TestCase):
    def setUp(self):
        self.backend = bridge.ExactRangeBackend()

    def test_minimum_sign_and_budget_conversion(self):
        q = [0, 1, 0, -1]
        proof = self.backend.minimum(q, F(1, 10**14), 33)
        lower, upper = self.backend.verify_minimum(q, proof)
        self.assertEqual(proof["max_nodes"], 32)
        self.assertEqual(proof["polynomial"], ["0", "-1", "0", "1"])
        self.assertLessEqual(lower, upper)
        self.assertGreaterEqual(lower*lower, F(4, 27))
        self.assertLessEqual(upper*upper, F(4, 27))
        self.assertLessEqual(upper-lower, F(1, 10**14))

    def test_polynomial_domain_and_tampering_binding(self):
        q = [0, 0, 1]
        proofs = [extrema.maximize([0, 0, -2]),
                  extrema.maximize([0, 0, -1], interval=(0, 1)),
                  family.maximum([0, 0, -1])]
        altered = self.backend.minimum(q, F(1, 1000), 1)
        altered["maximum_upper"] = "999"
        proofs.append(altered)
        for proof in proofs:
            with self.assertRaises(ValueError):
                self.backend.verify_minimum(q, proof)

    def test_open_proof_remains_a_valid_minimum_bound(self):
        proof = self.backend.minimum([0, 1, 0, -1], F(1, 10**30), 1)
        self.assertEqual(proof["status"], "open")
        lower, upper = self.backend.verify_minimum([0, 1, 0, -1], proof)
        self.assertLessEqual(lower, upper)

    def test_positive_search_hook_requires_the_matching_backend(self):
        lower = positive_search.lower_certificate([0, 0, 6], [[0, 0, 1]], [F(-4, 5)],
                                                   range_backend=self.backend)
        self.assertTrue(positive_search.verify_lower(lower, self.backend))
        self.assertFalse(positive_search.verify_lower(lower))
        self.assertEqual(lower["range_backend"], "exact_sturm_minimum_v1")
        self.assertEqual(lower["range_proof"]["format"], extrema.FORMAT)
        upper = positive_search.ansatz_upper_certificate([0, 0, 6], [[0, 0, 1]],
                    [{"t": "-1/2", "weight": "1/2"}, {"t": "1/2", "weight": "1/2"}])
        combined = positive_search.optimization_certificate(lower, upper, range_backend=self.backend)
        self.assertTrue(positive_search.verify(combined, self.backend))
        self.assertFalse(positive_search.verify(combined))

    def test_budget_and_input_validation(self):
        for budget in (0, -1, 4098, True, 2.0):
            with self.assertRaises(ValueError):
                self.backend.minimum([0, 1], F(1, 1000), budget)
        with self.assertRaises(ValueError):
            self.backend.minimum([0.1, 1], F(1, 1000), 2)


class FamilyTests(unittest.TestCase):
    def test_same_exact_family_formula_and_new_format(self):
        dirs = [[F(2, 3), 1], [0, 0, 1]]
        weight = [[2, F(1, 2)], [F(1, 2), 1]]
        cert = bridge.precise_family(dirs, weight, F(1, 10**14), 64)
        self.assertTrue(bridge.verify_family(cert))
        self.assertTrue(bridge.verify(json.loads(json.dumps(cert))))
        normalized, w, r, means, q = family.assemble(dirs, weight)
        self.assertEqual(cert["envelope_polynomial"], list(map(str, q)))
        self.assertEqual(cert["poisson_solutions"], [list(map(str, p)) for p in r])
        self.assertEqual(cert["means"], list(map(str, means)))
        self.assertEqual(cert["quadratic_bound"], dual.encode([[F(cert["scale"])*v for v in row] for row in w]))
        self.assertFalse(family.verify(cert))
        self.assertFalse(family.verify_range(cert["range_proof"]))

    def test_constant_direction_and_exact_zero_scale(self):
        cert = bridge.precise_family([[3]], [[1]], max_nodes=0)
        self.assertTrue(bridge.verify_family(cert))
        self.assertEqual(cert["scale"], "0")
        self.assertEqual(cert["range_status"], "exact")

    def test_family_tampering(self):
        original = bridge.precise_family(DIRS, [[1]], max_nodes=16)
        mutations = [
            lambda c: c.__setitem__("scale", "0"),
            lambda c: c["quadratic_bound"][0].__setitem__(0, "0"),
            lambda c: c["poisson_solutions"][0].__setitem__(1, "0"),
            lambda c: c.__setitem__("means", ["99"]),
            lambda c: c.__setitem__("scope", "all geometries"),
            lambda c: c.__setitem__("formal_assistant_checked", 0),
            lambda c: c.__setitem__("max_nodes", True),
            lambda c: c.__setitem__("range_proof", extrema.maximize([0, 0, 1], max_nodes=0)),
        ]
        for mutate in mutations:
            cert = copy.deepcopy(original)
            mutate(cert)
            self.assertFalse(bridge.verify_family(cert))
        cert = copy.deepcopy(original)
        q = list(map(F, cert["envelope_polynomial"]))
        cert["range_proof"] = extrema.maximize(q, interval=(-F(1, 2), F(1, 2)))
        self.assertFalse(bridge.verify_family(cert))

    def test_invalid_weight_and_input_rejected(self):
        for weight in ([[0]], [[-1]], [[1, 0]], [[1.0]]):
            with self.assertRaises((ValueError, TypeError)):
                bridge.precise_family([[0, 1]], weight)
        for kwargs in [dict(tolerance=0), dict(max_nodes=-1), dict(max_nodes=True)]:
            with self.assertRaises(ValueError):
                bridge.precise_family([[0, 1]], [[1]], **kwargs)


class MatrixTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.seed = legacy_seed()
        cls.result = bridge.precise_matrix(cls.seed, F(1, 10**14), 64)

    def test_actual_improvement_and_unchanged_dual(self):
        cert = self.result
        self.assertTrue(bridge.verify_matrix(cert))
        self.assertEqual(cert["dual"], self.seed["dual"])
        self.assertEqual(cert["lower_bound"], self.seed["lower_bound"])
        self.assertEqual(cert["weight"], self.seed["template"]["weight"])
        self.assertLess(F(cert["upper_bound"]), F(self.seed["upper_bound"]))
        self.assertEqual(cert["precision_outcome"], "improved")
        self.assertEqual(cert["range_status"], "enclosed")
        self.assertEqual(cert["status"], "global_gap_open")
        self.assertEqual(F(cert["gap"]), F(cert["upper_bound"])-F(cert["lower_bound"]))
        self.assertIn("Poisson", cert["scope"])
        self.assertIn("not a sharp spectral", cert["limitation"])
        self.assertFalse(dual.verify(cert))

    def test_negative_gain_retains_verified_seed(self):
        seed = legacy_seed(leaves=128)
        cert = bridge.precise_matrix(seed, F(1, 10**14), 0)
        self.assertTrue(bridge.verify_matrix(cert))
        self.assertEqual(cert["range_status"], "open")
        self.assertEqual(cert["precision_outcome"], "negative_gain_seed_retained")
        self.assertLess(F(cert["candidate_upper_improvement"]), 0)
        self.assertEqual(cert["upper_bound"], seed["upper_bound"])
        self.assertEqual(cert["upper_improvement"], "0")
        self.assertEqual(cert["upper_source"], "legacy_seed")
        self.assertEqual(cert["quadratic_bound"], seed["template"]["quadratic_bound"])

    def test_real_matrix_refine_wrapper_and_cost_scaling(self):
        dirs = [[0, 1, F(1, 3)], [0, 0, 1, F(1, 5)]]
        cost = [[2, F(1, 2)], [F(1, 2), 1]]
        seed_result = matrix_refine.refine(dirs, cost, budget=0)
        cert = bridge.precise_matrix(seed_result, F(1, 10**12), 64)
        self.assertTrue(bridge.verify_matrix(cert))
        self.assertEqual(cert["cost"], dual.encode(cost))
        self.assertEqual(cert["lower_bound"], seed_result["certificate"]["lower_bound"])
        self.assertEqual(F(cert["range_tolerance"])*F(cert["weight_cost_multiplier"]), F(1, 10**12))
        self.assertLessEqual(F(cert["range_objective_error_bound"]), F(1, 10**12))
        self.assertGreater(F(cert["upper_improvement"]), 0)

    def test_wrapper_metadata_has_no_authority_and_input_is_not_mutated(self):
        seed = copy.deepcopy(self.seed)
        wrapper = {"certificate": seed, "statistics": "garbage", "warm_state": {"samples": ["bad"]},
                   "cost": [[999]], "directions": [[999]], "status": "globally optimal"}
        before = copy.deepcopy(wrapper)
        cert = bridge.precise_matrix(wrapper, F(1, 10**14), 64)
        self.assertEqual(wrapper, before)
        self.assertEqual(cert, self.result)
        cert["seed_certificate"]["dual"]["lower_bound"] = "999"
        self.assertEqual(wrapper, before)

    def test_scalar_tie_and_exact_gap(self):
        seed = legacy_seed([[0, 1]], leaves=1)
        cert = bridge.precise_matrix(seed, F(1, 10**20), 0)
        self.assertTrue(bridge.verify_matrix(cert))
        self.assertEqual(cert["precision_outcome"], "tied")
        self.assertEqual(cert["gap"], "0")
        self.assertEqual(cert["status"], "global_gap_closed")

    def test_verification_does_not_repeat_synthesis(self):
        with patch.object(bridge, "precise_matrix", side_effect=AssertionError("search")), \
             patch.object(bridge, "precise_family", side_effect=AssertionError("search")), \
             patch.object(extrema, "maximize", side_effect=AssertionError("search")):
            self.assertTrue(bridge.verify_matrix(self.result))
            self.assertTrue(bridge.verify_family(self.result["precise_family"]))

    def test_matrix_binding_and_claim_tampering(self):
        mutations = [
            lambda c: c.__setitem__("directions", [["0", "2"]]),
            lambda c: c.__setitem__("cost", [["2"]]),
            lambda c: c.__setitem__("weight", [["2"]]),
            lambda c: c.__setitem__("lower_bound", "999"),
            lambda c: c.__setitem__("upper_bound", "0"),
            lambda c: c.__setitem__("gap", "0"),
            lambda c: c.__setitem__("candidate_upper_bound", "0"),
            lambda c: c.__setitem__("precision_outcome", "tied"),
            lambda c: c.__setitem__("upper_source", "legacy_seed"),
            lambda c: c.__setitem__("range_tolerance", "1"),
            lambda c: c.__setitem__("range_objective_error_bound", "0"),
            lambda c: c.__setitem__("scope", "sharp spectral optimum"),
            lambda c: c.__setitem__("formal_assistant_checked", 0),
            lambda c: c["dual"].__setitem__("lower_bound", "99"),
        ]
        for mutate in mutations:
            cert = copy.deepcopy(self.result)
            mutate(cert)
            self.assertFalse(bridge.verify_matrix(cert))
        cert = copy.deepcopy(self.result)
        cert["precise_family"] = bridge.precise_family([[0, 2]], [[1]],
                                                       F(cert["range_tolerance"]), cert["max_nodes"])
        self.assertFalse(bridge.verify_matrix(cert))
        cert = copy.deepcopy(self.result)
        cert["precise_family"] = bridge.precise_family(DIRS, [[2]],
                                                       F(cert["range_tolerance"]), cert["max_nodes"])
        self.assertFalse(bridge.verify_matrix(cert))

    def test_invalid_legacy_evidence_rejected(self):
        for seed in (None, {}, {"certificate": {}}, self.result):
            with self.assertRaises(ValueError):
                bridge.precise_matrix(seed)
        seed = copy.deepcopy(self.seed)
        seed["dual"]["cost"] = [["2"]]
        with self.assertRaises(ValueError):
            bridge.precise_matrix(seed)
        self.assertFalse(bridge.verify(None))


def comparison_rows():
    """Actual exact gains and losses; no elapsed times or model calls."""
    rows = []
    for name, leaves, nodes, dirs in [
        ("coarse_legacy_positive_gain", 1, 64, DIRS),
        ("tight_legacy_zero_budget_negative_gain", 128, 0, DIRS),
        ("already_exact_scalar_tie", 1, 0, [[0, 1]]),
    ]:
        seed = legacy_seed(dirs, leaves=leaves)
        cert = bridge.precise_matrix(seed, F(1, 10**14), nodes)
        if not bridge.verify_matrix(cert):
            raise AssertionError("Comparison proof failed")
        rows.append({"case": name, "legacy_leaf_budget": leaves, "exact_node_budget": nodes,
                     "actual_exact_nodes": cert["precise_family"]["range_proof"]["refinement_nodes"],
                     "previous_upper": cert["previous_upper_bound"],
                     "candidate_upper": cert["candidate_upper_bound"], "retained_upper": cert["upper_bound"],
                     "candidate_gain": cert["candidate_upper_improvement"],
                     "retained_gain": cert["upper_improvement"], "outcome": cert["precision_outcome"]})
    return rows


if __name__ == "__main__":
    unittest.main()
