"""Functional, independent replay and adversarial certificate checks."""
from fractions import Fraction as F
import copy
import json
import random
import unittest
from unittest.mock import patch

import uniform_refine as u


class UniformRefineTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.box = u.synthesize()
        cls.global_bound = u.synthesize(all_real=True)

    def test_exact_poisson_correctors_and_cubic_residual(self):
        r, energies = u.parameter_correctors(order=3)
        self.assertEqual(r, [[F(0)], [F(0), F(-1, 2)],
                             [F(0), F(0), F(-1, 24)],
                             [F(0), F(1, 48), F(0), F(-1, 144)]])
        self.assertEqual(energies, [F(0), F(0), F(-1, 6), F(0)])
        exact = {}
        u._sparse_add(exact, 2, [F(1, 30)])
        u._sparse_add(exact, 4, u._mul([1, 0, -1], [3, 0, -4]), F(1, 144))
        u._sparse_add(exact, 5, u._mul([0, 1], [1, 0, -2, 0, 1]), F(1, 288))
        u._sparse_add(exact, 6, [1, 0, -3, 0, 3, 0, -1], F(-1, 2304))
        self.assertEqual(u.local_residual([F(0), F(1)], F(1, 5), r), exact)

    def test_remainder_is_not_discarded_at_high_order(self):
        r, energies = u.parameter_correctors(order=4)
        self.assertEqual(energies[4], F(11, 1080))
        polynomial = u.local_residual([F(0), F(1)], F(1, 6), r)
        self.assertEqual(max(i for i, _ in polynomial), 8)
        factor, reduced = u.extract_even_factor(polynomial)
        self.assertEqual(factor, 4)
        self.assertEqual(polynomial, {(i + factor, j): x for (i, j), x in reduced.items()})
        self.assertEqual(u._evaluate_sparse(polynomial, F(0), F(2, 7)), 0)

    def test_whole_original_box(self):
        self.assertTrue(u.is_proved(self.box))
        self.assertEqual(self.box["leaf_count"], 1)
        self.assertEqual(F(self.box["certified_reduced_lower"]), F(2, 135))
        self.assertEqual(self.box["amplitude_box"], ["-2", "2"])
        self.assertEqual(self.box["proof_amplitude_box"], ["0", "2"])
        self.assertEqual(self.box["extracted_amplitude_power"], 2)
        self.assertIsNone(self.box["tail"])

    def test_full_real_glue_covers_boundaries(self):
        cert = self.global_bound
        self.assertTrue(u.is_proved(cert))
        self.assertEqual(cert["amplitude_box"], ["-5", "5"])
        self.assertEqual(cert["leaf_count"], 5)
        self.assertEqual(F(cert["certified_reduced_lower"]), F(139, 81920))
        radius = F(cert["tail"]["radius"])
        b = F(cert["tail"]["centered_potential_abs_upper"])
        c = F(cert["coefficient"])
        self.assertEqual(c * radius, b)
        for a in [-radius, radius, -radius - 1, radius + 1]:
            self.assertGreaterEqual(c * a * a - b * abs(a), 0)

    def test_failures_are_local_ansatz_failures_with_nonzero_exact_witnesses(self):
        for order in [1, 2]:
            cert = u.synthesize(order=order)
            self.assertTrue(u.verify(cert))
            self.assertFalse(u.is_proved(cert))
            self.assertEqual(cert["status"], "ansatz_obstruction")
            witness = cert["failure_witness"]
            self.assertNotEqual(F(witness["amplitude"]), 0)
            self.assertLess(F(witness["residual_value"]), 0)
            self.assertIn("not a spectral counterexample", witness["meaning"])

    def test_budget_exhaustion_is_honest_unresolved(self):
        cert = u.synthesize(all_real=True, max_leaves=1)
        self.assertTrue(u.verify(cert))
        self.assertFalse(u.is_proved(cert))
        self.assertEqual(cert["status"], "unresolved")
        self.assertIsNone(cert["failure_witness"])

    def test_reflection_is_optional_and_nonodd_inputs_are_not_reflected(self):
        no_reflection = u.synthesize(symmetry=False)
        self.assertTrue(u.is_proved(no_reflection))
        self.assertFalse(no_reflection["reflection_used"])
        self.assertEqual(no_reflection["proof_amplitude_box"], ["-2", "2"])
        quadratic = u.synthesize(direction=[0, 0, 1], coefficient=F(1, 30), order=1)
        self.assertTrue(u.is_proved(quadratic))
        self.assertFalse(quadratic["reflection_used"])

    def test_affine_mean_translation_and_nonsymmetric_box(self):
        cert = u.synthesize(direction=[F(7, 3), 1], amplitude_box=[-1, 2])
        self.assertTrue(u.is_proved(cert))
        self.assertEqual(cert["energy_series"][1], "7/3")
        self.assertFalse(cert["reflection_used"])
        self.assertEqual(cert["proof_amplitude_box"], ["-1", "2"])

    def test_independent_leaf_formula_matches_subdivision(self):
        rng = random.Random(34812)
        for _ in range(12):
            poly = {(i, j): F(rng.randrange(-5, 6), rng.randrange(1, 8))
                    for i in range(4) for j in range(5)}
            box = ((F(-3, 2), F(7, 3)), (F(-1), F(1)))
            coefficients = u.tensor_bernstein(poly, *box)
            self.assertEqual(coefficients, u.replay_tensor_bernstein(poly, *box))
            for axis in [0, 1]:
                children = u.split_tensor(coefficients, axis)
                child_boxes = u._split_box(box, axis)
                for child, child_box in zip(children, child_boxes):
                    self.assertEqual(child, u.replay_tensor_bernstein(poly, *child_box))

    def test_verifier_does_not_use_synthesis_range_or_subdivision(self):
        with patch.object(u, "tensor_bernstein", side_effect=AssertionError("generator disabled")), \
             patch.object(u, "split_tensor", side_effect=AssertionError("generator disabled")):
            self.assertTrue(u.verify(self.global_bound))

    def test_tamper_residual_factor_corrector_and_claim(self):
        mutations = [
            lambda c: c["correctors"][1].__setitem__(1, "-1/3"),
            lambda c: c.__setitem__("extracted_amplitude_power", 4),
            lambda c: c["reduced_residual"][0].__setitem__(2, "1"),
            lambda c: c.__setitem__("coefficient", "1/6"),
            lambda c: c.__setitem__("scope", "every surface"),
            lambda c: c.__setitem__("formal_assistant_checked", True),
            lambda c: c.__setitem__("status", "disproved"),
            lambda c: c.__setitem__("extra_unchecked_claim", "stronger"),
        ]
        for mutate in mutations:
            cert = copy.deepcopy(self.global_bound)
            mutate(cert)
            self.assertFalse(u.verify(cert))

    def test_tamper_partition_and_full_real_glue(self):
        mutations = [
            lambda c: c["tree"].__setitem__("axis", 2),
            lambda c: c["tree"].pop("right"),
            lambda c: c.__setitem__("tree", {"lower": "1", "upper": "2"}),
            lambda c: c.__setitem__("leaf_count", 1),
            lambda c: c["tail"].__setitem__("radius", "4"),
            lambda c: c["tail"].__setitem__("centered_potential_abs_upper", "1/2"),
            lambda c: c.__setitem__("amplitude_box", ["-4", "4"]),
            lambda c: c.__setitem__("proof_amplitude_box", ["0", "4"]),
            lambda c: c.__setitem__("reflection_used", False),
        ]
        for mutate in mutations:
            cert = copy.deepcopy(self.global_bound)
            mutate(cert)
            self.assertFalse(u.verify(cert))
        promoted = copy.deepcopy(self.box)
        promoted["all_real"] = True
        self.assertFalse(u.verify(promoted))

    def test_tamper_failure_and_unresolved_to_success(self):
        failure = u.synthesize(order=2)
        altered = copy.deepcopy(failure)
        altered["failure_witness"]["residual_value"] = "-1"
        self.assertFalse(u.verify(altered))
        altered = copy.deepcopy(failure)
        altered["failure_witness"]["amplitude"] = "100"
        self.assertFalse(u.verify(altered))
        altered = u.synthesize(all_real=True, max_leaves=1)
        altered["status"] = "proved"
        self.assertFalse(u.verify(altered))

    def test_json_round_trip_and_malformed_inputs(self):
        self.assertTrue(u.verify(json.loads(json.dumps(self.global_bound))))
        for malformed in [None, [], {}, {"method": u.METHOD}, {"direction": [0.0, 1.0]}]:
            self.assertFalse(u.verify(malformed))
        for kwargs in [{"coefficient": 0.2}, {"order": True}, {"order": 13},
                       {"amplitude_box": [1, 1]}, {"max_leaves": 0},
                       {"symmetry": "yes"}, {"direction": [0] * 7 + [1]}]:
            with self.assertRaises((ValueError, TypeError)):
                u.synthesize(**kwargs)

    def test_supplementary_random_point_checks_use_exact_fractions(self):
        rng = random.Random(116)
        cert = self.global_bound
        reduced = {(i, j): F(v) for i, j, v in cert["reduced_residual"]}
        for _ in range(64):
            a, t = F(rng.randrange(501), 100), F(rng.randrange(-100, 101), 100)
            self.assertGreaterEqual(u._evaluate_sparse(reduced, a, t), 0)

    def test_log_sobolev_sharp_rule_and_bound_theorem(self):
        cert = u.synthesize_log_sobolev()
        self.assertTrue(u.is_proved(cert))
        self.assertEqual(cert["coefficient"], "1/6")
        self.assertEqual(cert["sharp_uniform_quadratic_coefficient"], "1/6")
        self.assertEqual(cert["moment_series_rule"]["induction_difference_coefficients"], ["0", "4", "4"])
        self.assertFalse(cert["trusted_analytic_theorem"]["proved_by_this_program"])
        self.assertEqual(cert["optimality_proof"]["mass_in_amplitude"], ["1", "0", "1/12"])

    def test_log_sobolev_exact_spectral_counterexample(self):
        cert = u.synthesize_log_sobolev(coefficient=F(1, 7))
        self.assertTrue(u.verify(cert))
        self.assertFalse(u.is_proved(cert))
        self.assertEqual(cert["status"], "disproved")
        self.assertEqual(cert["counterexample"]["amplitude"], "1")
        self.assertEqual(F(cert["counterexample"]["spectral_objective_upper"]), F(-1, 91))
        close = u.synthesize_log_sobolev(coefficient=F(1, 6) - F(1, 10**20))
        self.assertTrue(u.verify(close))
        self.assertLess(F(close["counterexample"]["spectral_objective_upper"]), 0)

    def test_log_sobolev_affine_translation_scaling_and_scope(self):
        cert = u.synthesize_log_sobolev([F(7, 3), -2])
        self.assertTrue(u.is_proved(cert))
        self.assertEqual(cert["mean"], "7/3")
        self.assertEqual(cert["coefficient"], "2/3")
        with self.assertRaises(ValueError):
            u.synthesize_log_sobolev([0, 0, 1])
        self.assertTrue(u.is_proved(u.synthesize_log_sobolev([3])))
        self.assertEqual(u.synthesize_log_sobolev([3], -1)["status"], "disproved")

    def test_log_sobolev_tampering_is_rejected(self):
        cert = u.synthesize_log_sobolev()
        mutations = [
            lambda c: c.__setitem__("coefficient", "1/7"),
            lambda c: c["trusted_analytic_theorem"].__setitem__("entropy_to_dirichlet_constant", "1/2"),
            lambda c: c["trusted_analytic_theorem"].__setitem__("proved_by_this_program", True),
            lambda c: c["moment_series_rule"].__setitem__("induction_difference_coefficients", ["0", "0", "0"]),
            lambda c: c["optimality_proof"].__setitem__("necessary_coefficient", "1/7"),
            lambda c: c.__setitem__("scope", "every Riemannian surface"),
        ]
        for mutate in mutations:
            altered = copy.deepcopy(cert)
            mutate(altered)
            self.assertFalse(u.verify(altered))
        failed = u.synthesize_log_sobolev(coefficient="1/7")
        failed["counterexample"]["spectral_objective_upper"] = "-1"
        self.assertFalse(u.verify(failed))

    def test_evidence_keeps_both_proof_paths_and_failures(self):
        records = u.milestone_records()
        self.assertEqual(len(records), 8)
        self.assertTrue(all(u.verify(row["certificate"]) for row in records))
        self.assertEqual([row["certificate"]["status"] for row in records],
                         ["proved", "ansatz_obstruction", "ansatz_obstruction", "proved",
                          "proved", "ansatz_obstruction", "proved", "disproved"])


if __name__ == "__main__":
    unittest.main()
