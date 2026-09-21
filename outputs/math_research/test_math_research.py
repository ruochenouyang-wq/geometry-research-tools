import copy
from fractions import Fraction as F
import json
from pathlib import Path
import random
import tempfile
import unittest
from unittest.mock import patch
import urllib.error

import math_research as m


def problem(kind="sum", coefficients=None):
    return m.Problem.parse({"kind": kind, "coefficients": coefficients or [0, 0, 0, 1]})


def formula(coefficients):
    return {"coefficients": m.encoded(coefficients), "explanation": "测试候选"}


def divisor(value):
    return {"divisor": str(value), "explanation": "测试候选"}


class ProofTests(unittest.TestCase):
    def test_sum_cubes_known_formula(self):
        self.assertEqual(m.certify(problem(), formula([0, 0, F(1, 4), F(1, 2), F(1, 4)]))["status"], "certified")

    def test_wrong_base_constant_cannot_pass_induction(self):
        result = m.certify(problem(), formula([1, 0, F(1, 4), F(1, 2), F(1, 4)]))
        self.assertEqual(result["status"], "refuted")
        self.assertEqual(result["counterexample"]["n"], 0)

    def test_finite_tests_do_not_prove_identity(self):
        p = problem(coefficients=[1])
        perturbation = [F(1)]
        for n in range(9):
            perturbation = m.multiply(perturbation, [-F(n), F(1)])
        candidate = formula(m.add([0, 1], perturbation))
        self.assertEqual(m.sample_checks(p, candidate, range(9))["status"], "tested_only")
        proof = m.certify(p, candidate)
        self.assertEqual(proof["status"], "refuted")
        self.assertEqual(proof["counterexample"]["n"], 9)

    def test_zero_and_constant_sums(self):
        for p, q in [([0], [0]), ([7], [0, 7]), ([F(1, 3)], [0, F(1, 3)])]:
            spec = problem(coefficients=m.encoded(p))
            self.assertEqual(m.certify(spec, formula(q))["status"], "certified")

    def test_random_rational_sums_against_direct_arithmetic(self):
        rng = random.Random(31)
        for _ in range(24):
            coeff = [F(rng.randint(-7, 7), rng.randint(1, 5)) for _ in range(rng.randint(1, 8))]
            p = problem(coefficients=m.encoded(coeff))
            run = m.research(p, m.baseline_propose, rounds=10)
            self.assertTrue(m.verify_run(run))
            q = m.polynomial(run["history"][-1]["candidate"]["coefficients"])
            for n in [0, 1, 17, 41]:
                self.assertEqual(m.evaluate(q, n), sum((m.evaluate(coeff, k) for k in range(1, n + 1)), F(0)))

    def test_sum_degree_limit(self):
        p = problem(coefficients=[0] * m.MAX_DEGREE + [1])
        candidate = m.baseline_propose(p, [], m.MAX_DEGREE + 1)
        self.assertEqual(m.certify(p, candidate)["status"], "certified")

    def test_fixed_divisor_is_30(self):
        p = problem("fixed_divisor", [0, -1, 0, 0, 0, 1])
        result = m.certify(p, divisor(30))
        self.assertEqual(result["status"], "certified")
        self.assertEqual(result["certificate"]["forward_differences"], ["0", "0", "30", "150", "240", "120"])

    def test_valid_divisor_is_not_automatically_maximal(self):
        p = problem("fixed_divisor", [0, -1, 0, 0, 0, 1])
        self.assertEqual(m.certify(p, divisor(6))["status"], "feasible_not_optimal")

    def test_false_divisor_returns_checkable_witness(self):
        p = problem("fixed_divisor", [0, -1, 0, 1])
        result = m.certify(p, divisor(12))
        w = result["counterexample"]
        self.assertNotEqual(int(m.evaluate(p.coefficients, w["n"])) % 12, 0)

    def test_negative_inputs_and_negative_constant(self):
        for coeff, expected in [([-12], 12), ([0, -1, 1], 2), ([0, -1, 0, 1], 6)]:
            p = problem("fixed_divisor", coeff)
            self.assertEqual(m.certify(p, divisor(expected))["status"], "certified")
            self.assertEqual(m.sample_checks(p, divisor(expected), range(-100, 101))["status"], "tested_only")

    def test_random_fixed_divisors_against_enumeration(self):
        rng = random.Random(91)
        for _ in range(50):
            coeff = [rng.randint(-20, 20) for _ in range(rng.randint(1, 10))]
            if not any(coeff):
                continue
            p = problem("fixed_divisor", coeff)
            expected = m.reduce(m.gcd, (abs(int(m.evaluate(coeff, n))) for n in range(-50, 51)), 0)
            self.assertEqual(m.certify(p, divisor(expected))["status"], "certified")


class WorkflowTests(unittest.TestCase):
    def test_failure_feedback_reaches_next_proposal(self):
        calls = []
        def proposer(p, history, turn):
            calls.append(copy.deepcopy(history))
            return formula([0, 1]) if turn == 1 else formula([0, 0, F(1, 4), F(1, 2), F(1, 4)])
        run = m.research(problem(), proposer)
        self.assertEqual(run["status"], "certified")
        self.assertEqual(calls[1][0]["status"], "refuted")
        self.assertIn("counterexample", calls[1][0])

    def test_duplicate_normalization_and_budget(self):
        run = m.research(problem(), lambda *args: formula([0, 1]), rounds=3)
        self.assertEqual(run["status"], "budget_exhausted")
        self.assertEqual([x["status"] for x in run["history"]], ["refuted", "duplicate", "duplicate"])

    def test_invalid_candidate_recorded(self):
        run = m.research(problem(), lambda *args: {"coefficients": [0.25], "explanation": "x"}, rounds=1)
        self.assertEqual(run["history"][0]["status"], "invalid_candidate")

    def test_forged_success_rejected(self):
        run = m.research(problem(), m.baseline_propose)
        self.assertTrue(m.verify_run(run))
        run["history"][-1]["candidate"]["coefficients"][0] = "1"
        self.assertFalse(m.verify_run(run))

    def test_tampered_certificate_rejected(self):
        run = m.research(problem(), m.baseline_propose)
        run["history"][-1]["certificate"]["base_value"] = "99"
        self.assertFalse(m.verify_run(run))

    def test_changed_problem_rejected(self):
        run = m.research(problem(), m.baseline_propose)
        run["problem"]["coefficients"] = ["2"]
        self.assertFalse(m.verify_run(run))

    def test_problem_rejects_unsupported_or_ambiguous_input(self):
        invalid = [{"kind": "anything", "coefficients": [1]},
                   {"kind": "fixed_divisor", "coefficients": [0]},
                   {"kind": "fixed_divisor", "coefficients": ["1/2"]},
                   {"kind": "sum", "coefficients": [1], "domain": "n>10"}]
        for obj in invalid:
            with self.assertRaises(ValueError):
                m.Problem.parse(obj)

    def test_numbers_reject_code_floats_and_bad_fractions(self):
        for value in [True, 0.1, "1/0", "1e10", "__import__('os')", "1" * 101]:
            with self.assertRaises(ValueError):
                m.exact(value)

    def test_serialization_and_reverification(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "run.json"
            m.save_json(path, m.research(problem(), m.baseline_propose))
            self.assertTrue(m.verify_run(m.load_json(path)))


class APITests(unittest.TestCase):
    def response(self):
        return {"id": "mock_response", "status": "completed", "usage": {"output_tokens": 42},
                "output": [{"type": "reasoning", "summary": []}, {"type": "message", "content": [
                    {"type": "output_text", "text": json.dumps(formula([0, 0, F(1, 4), F(1, 2), F(1, 4)]))}]}]}

    def test_extract_completed_response(self):
        self.assertEqual(m.extract_response(self.response())["coefficients"][-1], "1/4")

    def test_refusal_incomplete_and_missing_text(self):
        invalid = [{"status": "incomplete"}, {"status": "completed", "output": []},
                   {"status": "completed", "output": [{"type": "message", "content": [{"type": "refusal"}]}]}]
        for obj in invalid:
            with self.assertRaises(m.ProviderError):
                m.extract_response(obj)

    @patch.dict("os.environ", {"OPENAI_API_KEY": "test-key"})
    @patch("urllib.request.urlopen")
    def test_mock_api_round_trip_and_request(self, mocked):
        mocked.return_value.__enter__.return_value.read.return_value = json.dumps(self.response()).encode()
        proposer = m.GPTProposer("test-model")
        result = m.research(problem(), proposer, provider_name="gpt", model="test-model")
        self.assertTrue(m.verify_run(result))
        request = mocked.call_args.args[0]
        body = json.loads(request.data)
        self.assertEqual(request.full_url, "https://api.openai.com/v1/responses")
        self.assertEqual(body["model"], "test-model")
        self.assertFalse(body["store"])
        self.assertTrue(body["text"]["format"]["strict"])
        self.assertNotIn("test-key", json.dumps(result))
        self.assertEqual(result["api_usage"][0]["usage"]["output_tokens"], 42)

    @patch.dict("os.environ", {"OPENAI_API_KEY": "test-key"})
    @patch("urllib.request.urlopen")
    def test_http_error_stops_without_retry(self, mocked):
        mocked.side_effect = urllib.error.HTTPError("https://api.openai.com", 429, "quota", None, None)
        run = m.research(problem(), m.GPTProposer("test-model"))
        self.assertEqual(run["status"], "provider_error")
        self.assertEqual(mocked.call_count, 1)
        self.assertEqual(len(run["history"]), 1)

    @patch.dict("os.environ", {}, clear=True)
    def test_missing_key_is_explicit(self):
        with self.assertRaisesRegex(ValueError, "OPENAI_API_KEY"):
            m.GPTProposer("test-model")


if __name__ == "__main__":
    unittest.main()
