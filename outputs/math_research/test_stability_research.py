import copy
from fractions import Fraction as F
import importlib.util
from pathlib import Path
import unittest

import stability_research as s
from math_research import load_json

ROOT = Path(__file__).resolve().parent
NUMPY = importlib.util.find_spec("numpy") is not None


def load(name):
    return load_json(ROOT / "examples" / name)


def log_certificate():
    return {"plan": {"family": "log_quadratic", "degree": 2, "log_indices": [0], "explanation": "测试证书"},
            "coefficients": ["1", "1"], "grams": [[["1", "0", "0"], ["0", "1", "-1"], ["0", "-1", "2"]]]}


def quartic_certificate():
    return {"plan": {"family": "polynomial", "degree": 4, "log_indices": [], "explanation": "测试证书"},
            "coefficients": ["2", "0", "0", "0", "2"],
            "grams": [[["1", "0", "-1"], ["0", "0", "0"], ["-1", "0", "1"]],
                      [["7", "0", "-7", "0"], ["0", "11", "0", "-7"], ["-7", "0", "11", "0"], ["0", "-7", "0", "7"]]]}


class ExactStabilityTests(unittest.TestCase):
    def test_log_certificate(self):
        self.assertTrue(s.verify_certificate(load("stability_no_polynomial.json"), log_certificate()))

    def test_quartic_certificate_with_singular_gram(self):
        self.assertTrue(s.verify_certificate(load("stability_quartic.json"), quartic_certificate()))

    def test_change_coefficient_rejected(self):
        cert = log_certificate()
        cert["coefficients"][0] = "100"
        self.assertFalse(s.verify_certificate(load("stability_no_polynomial.json"), cert))

    def test_missing_variable_coercivity_rejected(self):
        cert = log_certificate()
        cert["coefficients"][1] = "0"
        self.assertFalse(s.verify_certificate(load("stability_no_polynomial.json"), cert))

    def test_change_vector_field_rejected(self):
        p = load("stability_no_polynomial.json")
        p["vector_field"][1][0]["coefficient"] = "1"
        self.assertFalse(s.verify_certificate(p, log_certificate()))

    def test_change_log_variable_rejected(self):
        cert = log_certificate()
        cert["plan"]["log_indices"] = [1]
        self.assertFalse(s.verify_certificate(load("stability_no_polynomial.json"), cert))

    def test_gram_tampering_rejected(self):
        cert = quartic_certificate()
        cert["grams"][1][0][0] = "8"
        self.assertFalse(s.verify_certificate(load("stability_quartic.json"), cert))

    def test_positive_diagonal_is_not_psd(self):
        self.assertFalse(s.exact_psd([[1, 2], [2, 1]]))

    def test_tiny_negative_eigenvalue_is_rejected_exactly(self):
        self.assertFalse(s.exact_psd([[1, 0], [0, -F(1, 10**90)]]))

    def test_zero_pivot_and_symmetry(self):
        self.assertTrue(s.exact_psd([[0, 0], [0, 1]]))
        self.assertFalse(s.exact_psd([[0, 1], [1, 1]]))
        self.assertFalse(s.exact_psd([[1, 1], [0, 1]]))
        self.assertFalse(s.exact_psd([[1, 2]]))

    def test_psd_outer_products(self):
        vectors = [[F(1, 3), -2, 5], [4, F(3, 7), -1]]
        matrix = [[sum(v[i] * v[j] for v in vectors) for j in range(3)] for i in range(3)]
        self.assertTrue(s.exact_psd(matrix))

    def test_rref_rational_reconstruction(self):
        rows = [[F(1), F(2), F(3)], [F(2), F(4), F(6)]]
        reduced = s.rref(rows, [F(7), F(14)])
        self.assertEqual(reduced[1], [0])
        self.assertIsNone(s.rref(rows, [F(7), F(15)]))

    def test_reject_non_equilibrium_and_changed_domain(self):
        p = load("stability_no_polynomial.json")
        p["vector_field"][0].append({"powers": [0, 0], "coefficient": "1"})
        with self.assertRaises(ValueError):
            s.parse_problem(p)
        p = load("stability_no_polynomial.json")
        p["domain"] = "small ball"
        with self.assertRaises(ValueError):
            s.parse_problem(p)

    def test_reject_malformed_plans(self):
        for indices in [[0, 0], [2], [True]]:
            plan = log_certificate()["plan"]
            plan["log_indices"] = indices
            with self.assertRaises(ValueError):
                s.normalize_plan(plan, 2)

    def test_multivariate_lie_derivative_by_hand(self):
        n, f = s.parse_problem(load("stability_quartic.json"))
        v = {(4, 0): F(1), (0, 4): F(1)}
        self.assertEqual(s.lie(v, f), {(6, 0): F(-4), (0, 6): F(-4)})


@unittest.skipUnless(NUMPY, "数值搜索测试需要 NumPy；上面的精确验证测试无此依赖。")
class SearchTests(unittest.TestCase):
    def test_representation_switch_succeeds(self):
        p = load("stability_no_polynomial.json")
        result = s.run_research(p, iterations=2500)
        self.assertEqual([r["status"] for r in result["history"]], ["search_failed", "search_failed", "certified"])
        self.assertEqual(result["certificate"]["plan"]["family"], "log_quadratic")
        self.assertTrue(s.verify_certificate(p, result["certificate"]))

    def test_degree_expansion_succeeds(self):
        p = load("stability_quartic.json")
        result = s.run_research(p, iterations=2500)
        self.assertEqual(result["certificate"]["plan"]["degree"], 4)
        self.assertTrue(s.verify_certificate(p, result["certificate"]))

    def test_unstable_system_not_certified(self):
        result = s.run_research(load("stability_unstable.json"), rounds=5, iterations=500)
        self.assertEqual(result["status"], "not_solved")
        self.assertNotIn("certificate", result)

    def test_manual_gpt_plan_uses_same_verifier(self):
        p = load("stability_no_polynomial.json")
        result = s.run_research(p, plans=[log_certificate()["plan"]])
        self.assertEqual(len(result["history"]), 1)
        self.assertTrue(s.verify_certificate(p, result["certificate"]))

    def test_three_dimensional_stable_system(self):
        p = {"title": "三维线性稳定基准", "dimension": 3,
             "vector_field": [[{"powers": list(s.unit(3, i)), "coefficient": str(-i - 1)}] for i in range(3)]}
        result = s.run_research(p, rounds=1)
        self.assertTrue(s.verify_certificate(p, result["certificate"]))


if __name__ == "__main__":
    unittest.main()
