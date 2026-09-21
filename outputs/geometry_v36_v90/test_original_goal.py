"""V89 original-goal binding, scope, exact extrema and mutation tests."""
from fractions import Fraction as F
import copy
import json
import unittest

import compiler_v23 as compiler
import original_goal as original
import uniform_refine as uniform


def goal(potential="a*t", penalty="a**2/5", parameters=None, domain=None,
         threshold="0", space="full_sphere"):
    return {"geometry": "unit_sphere", "function_space": space, "eigenvalue_index": 1,
            "parameters": ["a"] if parameters is None else parameters,
            "potential": potential, "penalty": penalty,
            "domain": {"kind": "all_real"} if domain is None else domain,
            "threshold": threshold}


class OriginalGoalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raw = goal()
        cls.cert = original.solve(cls.raw)
        cls.barta = original.solve(cls.raw, rule="barta")

    def test_original_goal_exact_zero_with_both_independent_rules(self):
        for cert in (self.cert, self.barta):
            self.assertTrue(original.verify(cert))
            self.assertEqual(cert["method"], "original_spectral_goal_v89")
            self.assertEqual(cert["algorithm_version"], 89)
            self.assertEqual(cert["status"], "proved")
            self.assertEqual(cert["lower_bound"], "0")
            self.assertEqual(cert["upper_bound"], "0")
            self.assertTrue(cert["infimum_exact"])
            self.assertEqual(cert["candidate_parameters"], ["0"])
            self.assertEqual(cert["compiled_goal"]["source"], self.raw)
        self.assertEqual(self.cert["derived_quadratic"]["hessian"], [["1/15"]])
        self.assertEqual(self.barta["derived_quadratic"]["hessian"], [["0"]])

    def test_original_box_exact_zero_and_normalized_domain(self):
        raw = goal(domain={"kind": "box", "box": [["-2", "2"]]})
        cert = original.solve(raw)
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert["lower_bound"], "0")
        self.assertEqual(cert["upper_bound"], "0")
        self.assertEqual(cert["compiled_goal"]["problem"]["domain"], raw["domain"])

    def test_axisymmetric_and_full_sphere_scopes_are_bound(self):
        axisymmetric = original.solve(goal(space="axisymmetric"))
        self.assertTrue(original.verify(axisymmetric))
        self.assertEqual(axisymmetric["operator_binding"]["function_space"], "axisymmetric")
        self.assertEqual(self.cert["operator_binding"]["function_space"], "full_sphere")
        self.assertIn("probability", self.cert["operator_binding"]["measure"])

    def test_fixed_polynomial_remainder_and_affine_parameter_mixing(self):
        raw = goal("3+2*t+t**2+a*(1+2*t)+b*(-2+t)",
                   "(2+2*a+b)**2/6-a+2*b+a**2+b**2", parameters=["a", "b"], threshold="3")
        cert = original.solve(raw)
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert["amplitude_offset"], "2")
        self.assertEqual(cert["amplitude_slopes"], ["2", "1"])
        self.assertEqual(cert["parameter_constants"], ["1", "-2"])
        self.assertEqual(cert["remainder_polynomial"], ["3", "0", "1"])
        self.assertEqual(cert["remainder_lower"], "3")
        self.assertEqual(cert["derived_quadratic"],
                         {"constant": "3", "linear": ["0", "0"], "hessian": [["2", "0"], ["0", "2"]]})
        self.assertEqual(cert["lower_bound"], "3")
        self.assertEqual(cert["upper_bound"], "122/33")
        self.assertFalse(cert["infimum_exact"])

    def test_nonlinear_fixed_potential_range_has_complete_evidence(self):
        cert = original.solve(goal("t**4-t**2+a*t", "a**2/6+1/4"))
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert["remainder_lower"], "-1/4")
        self.assertEqual(cert["remainder_range"]["coordinate"], "t_squared")
        self.assertEqual(cert["status"], "proved")
        self.assertEqual(cert["lower_bound"], "0")

    def test_non_even_fixed_polynomial_remainder_uses_full_range(self):
        cert = original.solve(goal("t**3+a*t", "a**2/6+1"))
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert["remainder_range"]["coordinate"], "t")
        self.assertEqual(cert["remainder_range"]["intervals"][0][0], "-1")
        self.assertEqual(cert["remainder_range"]["intervals"][-1][1], "1")
        self.assertEqual(cert["status"], "proved")

    def test_rank_deficient_affine_directions_and_passive_parameters(self):
        cert = original.solve(goal("(a+2*b)*t+3*c", "(a+2*b)**2/6-3*c+c**2",
                                   parameters=["a", "b", "c"]))
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert["status"], "proved")
        self.assertEqual(cert["lower_bound"], "0")
        self.assertEqual(cert["upper_bound"], "0")
        self.assertEqual(cert["amplitude_slopes"], ["1", "2", "0"])

    def test_zero_parameter_problem_and_constant_trial_normalization(self):
        cert = original.solve(goal("2", "3", parameters=[], threshold="5"))
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert["lower_bound"], "5")
        self.assertEqual(cert["upper_bound"], "5")
        self.assertEqual(cert["candidate_parameters"], [])
        self.assertEqual(cert["constant_trial"]["probability_mass"], "1")
        self.assertEqual(cert["constant_trial"]["dirichlet_energy"], "0")

    def test_indefinite_box_quadratic_is_minimized_over_entire_box(self):
        domain = {"kind": "box", "box": [["1", "2"]]}
        cert = original.solve(goal(penalty="-a**2", domain=domain, threshold="-5"))
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert["lower_bound"], "-14/3")
        self.assertEqual(cert["relaxation_minimizer"], ["2"])
        self.assertEqual(cert["upper_bound"], "-4")
        self.assertEqual(cert["status"], "proved")

    def test_spectral_refutation_uses_original_constant_trial(self):
        cert = original.solve(goal(penalty="-a**2"))
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert["status"], "refuted")
        self.assertIsNone(cert["lower_bound"])
        self.assertEqual(cert["upper_bound"], "-1")
        self.assertEqual(cert["constant_trial"]["trial_power_coefficients"], ["1"])
        self.assertEqual(cert["constant_trial"]["spectral_rayleigh_upper"], "0")
        self.assertEqual(cert["constant_trial"]["original_penalty"], "-1")

    def test_unbounded_lower_relaxation_is_not_called_spectral_refutation(self):
        cert = original.solve(goal(penalty="a**2/7"))
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert["status"], "unresolved")
        self.assertIsNone(cert["lower_bound"])
        self.assertEqual(cert["upper_bound"], "0")
        self.assertEqual(cert["lower_diagnostic"]["reason"], "negative_hessian_direction")

    def test_nullspace_linear_relaxation_diagnostic(self):
        cert = original.solve(goal(penalty="a**2/6+a", threshold="-10"))
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert["lower_diagnostic"]["reason"], "linear_term_in_hessian_nullspace")
        self.assertEqual(cert["status"], "unresolved")
        self.assertIsNone(cert["lower_bound"])

    def test_old_threshold_semantics_and_strict_refutation(self):
        positive_threshold = original.solve(goal(threshold="1"))
        self.assertEqual(positive_threshold["status"], "refuted")
        self.assertEqual(positive_threshold["upper_bound"], "0")
        self.assertEqual(original.solve(goal(threshold="-1"))["status"], "proved")
        unresolved = original.solve(goal("t", "0", parameters=[]))
        self.assertEqual(unresolved["status"], "unresolved")
        self.assertEqual(unresolved["upper_bound"], "0")
        self.assertEqual(unresolved["threshold"], "0")

    def test_unsupported_inputs_fail_explicitly(self):
        invalid = [goal("a*t**2"), goal("a**2*t"), goal("t**7+a*t"),
                   goal(space="arbitrary_functions"), goal(domain={"kind": "half_line"})]
        for raw in invalid:
            with self.assertRaises((ValueError, TypeError)):
                original.solve(raw)
        for field, value in [("geometry", "arbitrary_surface"), ("eigenvalue_index", 2)]:
            raw = goal()
            raw[field] = value
            with self.assertRaises(ValueError):
                original.solve(raw)
        with self.assertRaises(ValueError):
            original.solve(goal(), rule="unverified_local_energy")
        with self.assertRaises(ValueError):
            original.solve(goal(), coefficient=0.2)

    def test_reject_unproved_requested_uniform_coefficients(self):
        with self.assertRaises(ValueError):
            original.solve(goal(), coefficient="1/7")
        with self.assertRaises(ValueError):
            original.solve(goal(), coefficient="1/6", rule="barta")

    def test_classical_theorem_trust_is_explicit_and_separate(self):
        theorem = self.cert["trust"]["classical_theorem"]
        self.assertEqual(theorem["entropy_to_dirichlet_constant"], "1")
        self.assertFalse(theorem["proved_by_this_program"])
        self.assertIsNone(self.barta["trust"]["classical_theorem"])
        self.assertFalse(self.cert["formal_assistant_checked"])

    def test_constant_trial_is_reintegrated_in_the_original_coordinates(self):
        raw = goal("1+t**2+a*(2-t)", "(a-3)**2", domain={"kind": "box", "box": [["1", "2"]]})
        cert = original.solve(raw)
        self.assertTrue(original.verify(cert))
        compiled = compiler.compile_goal(raw)
        for witness in cert["constant_trial_candidates"]:
            a = F(witness["parameters"][0])
            self.assertGreaterEqual(a, 1)
            self.assertLessEqual(a, 2)
            expected = F(4, 3) + 2 * a + (a - 3) ** 2
            self.assertEqual(F(witness["objective_upper"]), expected)
        self.assertEqual(cert["compiled_goal"], compiled)

    def test_tamper_original_goal_operator_domain_and_threshold(self):
        mutations = [
            lambda c: c["compiled_goal"]["source"].__setitem__("potential", "2*a*t"),
            lambda c: c["compiled_goal"]["source"].__setitem__("penalty", "a**2/7"),
            lambda c: c["compiled_goal"]["problem"].__setitem__("domain", {"kind": "box", "box": [["1", "2"]]}),
            lambda c: c["operator_binding"].__setitem__("geometry", "sphere_radius_2"),
            lambda c: c["operator_binding"].__setitem__("function_space", "axisymmetric"),
            lambda c: c["operator_binding"].__setitem__("eigenvalue_index", 2),
            lambda c: c.__setitem__("threshold", "1"),
            lambda c: c.__setitem__("relation", ">"),
        ]
        for mutate in mutations:
            cert = copy.deepcopy(self.cert)
            mutate(cert)
            self.assertFalse(original.verify(cert))

    def test_tamper_uniform_lemma_and_refuse_finite_box_promotion(self):
        cert = copy.deepcopy(self.barta)
        cert["uniform_lemma"] = uniform.synthesize()
        self.assertTrue(uniform.is_proved(cert["uniform_lemma"]))
        self.assertFalse(original.verify(cert))
        cert = copy.deepcopy(self.cert)
        cert["uniform_lemma"] = uniform.synthesize_log_sobolev(coefficient="1/7")
        self.assertFalse(original.verify(cert))
        cert = copy.deepcopy(self.cert)
        cert["uniform_lemma"]["trusted_analytic_theorem"]["entropy_to_dirichlet_constant"] = "1/2"
        self.assertFalse(original.verify(cert))

    def test_tamper_remainder_and_derived_quadratic(self):
        mutations = [
            lambda c: c.__setitem__("remainder_lower", "1"),
            lambda c: c["remainder_range"].__setitem__("maximum_upper", "-1"),
            lambda c: c["remainder_range"].__setitem__("polynomial", ["1"]),
            lambda c: c.__setitem__("amplitude_offset", "1"),
            lambda c: c.__setitem__("amplitude_slopes", ["0"]),
            lambda c: c.__setitem__("parameter_constants", ["1"]),
            lambda c: c["derived_quadratic"].__setitem__("constant", "1"),
            lambda c: c["derived_quadratic"].__setitem__("hessian", [["0"]]),
            lambda c: c.__setitem__("lower_bound", "1"),
        ]
        for mutate in mutations:
            cert = copy.deepcopy(self.cert)
            mutate(cert)
            self.assertFalse(original.verify(cert))

    def test_tamper_upper_witness_cannot_smuggle_local_energy_as_spectral_upper(self):
        mutations = [
            lambda c: c.__setitem__("upper_bound", "-1/6"),
            lambda c: c["constant_trial"].__setitem__("objective_upper", "-1/6"),
            lambda c: c["constant_trial"].__setitem__("trial_power_coefficients", ["1", "-1/2"]),
            lambda c: c["constant_trial"].__setitem__("probability_mass", "2"),
            lambda c: c.__setitem__("candidate_parameters", ["10"]),
            lambda c: c.__setitem__("infimum_exact", False),
            lambda c: c.__setitem__("status", "refuted"),
            lambda c: c["trust"]["classical_theorem"].__setitem__("proved_by_this_program", True),
        ]
        for mutate in mutations:
            cert = copy.deepcopy(self.cert)
            mutate(cert)
            self.assertFalse(original.verify(cert))

    def test_json_roundtrip_and_invalid_records(self):
        self.assertTrue(original.verify(json.loads(json.dumps(self.cert))))
        for malformed in [None, [], {}, {"method": original.METHOD}]:
            self.assertFalse(original.verify(malformed))
        extra = copy.deepcopy(self.cert)
        extra["unchecked_claim"] = "all geometries"
        self.assertFalse(original.verify(extra))


if __name__ == "__main__":
    unittest.main()
