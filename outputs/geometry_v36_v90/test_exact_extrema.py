"""Focused exact arithmetic, adversarial replay and precision comparisons."""

import copy
from fractions import Fraction as F
import json
from math import comb
import random
import unittest
from unittest.mock import patch

import exact_extrema as ex


def multiply(a, b):
    result = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return ex.trim(result)


def from_roots(roots):
    result = [F(1)]
    for root in roots:
        result = multiply(result, [-F(root), F(1)])
    return result


def integrate(p):
    return [F(0)] + [F(c, i+1) for i, c in enumerate(p)]


class ExactExtremaTests(unittest.TestCase):
    def checked(self, q, **options):
        cert = ex.maximize(q, **options)
        self.assertTrue(ex.verify(cert))
        self.assertTrue(ex.verify(json.loads(json.dumps(cert))))
        self.assertLessEqual(F(cert["maximum_lower"]), F(cert["maximum_upper"]))
        return cert

    def test_constants_affine_and_singleton(self):
        for q, domain, answer in [([7], (-5, 9), 7), ([2, -3], (-2, 4), 8),
                                  ([0], (-1, 1), 0), ([1, 2, 3], (F(2, 3), F(2, 3)), F(11, 3))]:
            with self.subTest(q=q, domain=domain):
                c = self.checked(q, interval=domain, max_nodes=0)
                self.assertEqual(c["status"], "exact")
                self.assertEqual(F(c["maximum_upper"]), answer)

    def test_one_sided_sturm_counts_at_every_boundary(self):
        sf, _ = ex.square_free(from_roots([-1, -1, 0, 0, 0, 1, 2]))
        seq = ex.sturm_sequence(sf)
        self.assertEqual(len(sf)-1, 4)
        for a, b, answer in [(-1, 1, 1), (-2, 2, 3), (-2, 3, 4),
                              (0, 1, 0), (-1, 0, 0), (F(-1, 2), F(3, 2), 2)]:
            self.assertEqual(ex.count_open(seq, F(a), F(b)), answer)

    def test_sturm_negative_leading_remainders(self):
        for roots in [(-2, -1, 0, 1, 3), (F(-4, 5), F(1, 7), F(5, 6))]:
            p = from_roots(roots)
            for sign in (1, -1):
                sequence = ex.sturm_sequence([sign*c for c in p])
                for i in range(-20, 20):
                    a, b = F(i, 10), F(i+3, 10)
                    self.assertEqual(ex.count_open(sequence, a, b), sum(a < r < b for r in roots))

    def test_repeated_rational_flat_peak(self):
        for degree in (4, 8, 24, 32):
            q = [-F(comb(degree, k))*(-F(1, 3))**(degree-k) for k in range(degree+1)]
            c = self.checked(q, max_nodes=0)
            self.assertEqual(c["status"], "exact")
            self.assertEqual(c["maximum_upper"], "0")
            self.assertEqual(c["rational_critical_points"], [
                {"at": "1/3", "value": "0", "derivative_multiplicity": degree-1}])

    def test_stationary_inflection_and_endpoint_roots(self):
        c = self.checked([0, 0, 0, 1], max_nodes=0)
        self.assertEqual(c["maximum_upper"], "1")
        self.assertEqual(c["rational_critical_points"][0]["derivative_multiplicity"], 2)
        q = integrate(from_roots([-1, -1, 0, 1, 1, 1]))
        c = self.checked(q, suggestions=[0], tolerance=F(1, 10**10), max_nodes=12)
        self.assertEqual(c["interior_critical_count"], 1)
        self.assertEqual([r["at"] for r in c["rational_critical_points"]], ["-1", "0", "1"])
        self.assertEqual(c["status"], "exact")

    def test_even_reduction_on_all_interval_positions(self):
        q = [0, 0, 1, 0, -1]
        for domain, image, maximum in [((-1, 1), ["0", "1"], F(1, 4)),
                                       ((-2, -1), ["1", "4"], F(0)),
                                       ((0, F(1, 2)), ["0", "1/4"], F(3, 16)),
                                       ((-1, F(1, 4)), ["0", "1"], F(1, 4))]:
            c = self.checked(q, interval=domain, max_nodes=0)
            self.assertEqual(c["coordinate"], "x_squared")
            self.assertEqual(c["working_interval"], image)
            self.assertEqual(F(c["maximum_upper"]), maximum)
            self.assertEqual(c["status"], "exact")

    def test_irrational_critical_points(self):
        # max(x-x^3)=2/(3 sqrt(3)). Squaring checks the algebraic truth exactly.
        c = self.checked([0, 1, 0, -1], tolerance=F(1, 10**18), max_nodes=128)
        lo, hi = F(c["maximum_lower"]), F(c["maximum_upper"])
        self.assertLessEqual(lo*lo, F(4, 27))
        self.assertGreaterEqual(hi*hi, F(4, 27))
        self.assertLessEqual(hi-lo, F(1, 10**18))
        self.assertEqual(c["interior_critical_count"], 2)
        self.assertTrue(c["isolation_complete"])
        self.assertEqual(c["status"], "enclosed")

    def test_multiple_peaks_and_complete_critical_coverage(self):
        roots = [F(-4, 5), F(-1, 3), F(1, 11), F(1, 11), F(1, 2), F(4, 5)]
        q = integrate([-v for v in from_roots(roots)])
        true_max = max(ex.evaluate(q, x) for x in [-1, 1]+roots)
        c = self.checked(q, suggestions=roots, tolerance=F(1, 10**20), max_nodes=64)
        self.assertEqual(c["interior_critical_count"], len(set(roots)))
        self.assertEqual(c["status"], "exact")
        self.assertEqual(F(c["maximum_lower"]), true_max)
        self.assertEqual(F(c["maximum_upper"]), true_max)

    def test_close_roots_require_honest_isolation(self):
        roots = [F(1, 3), F(1, 3)+F(1, 2**35), F(2, 3)]
        q = integrate(from_roots(roots))
        c = self.checked(q, tolerance=F(10**6), max_nodes=0)
        self.assertTrue(c["target_met"])
        self.assertFalse(c["isolation_complete"])
        self.assertEqual(c["status"], "open")
        self.assertTrue(c["budget_exhausted"])
        c2 = self.checked(q, suggestions=roots, tolerance=F(1, 10**20), max_nodes=16)
        self.assertEqual(c2["status"], "exact")

    def test_zero_budget_and_successive_refinement(self):
        previous = None
        for budget in (0, 1, 2, 4, 8, 16, 32):
            c = self.checked([0, 1, 0, -1], tolerance=F(1, 10**20), max_nodes=budget)
            self.assertLessEqual(c["refinement_nodes"], budget)
            lo, hi = map(F, (c["maximum_lower"], c["maximum_upper"]))
            self.assertLessEqual(lo*lo, F(4, 27))
            self.assertGreaterEqual(hi*hi, F(4, 27))
            if previous is not None:
                # Global enclosure refinement is checked for this example;
                # monotonicity is not a general promise of intersected local bounds.
                self.assertLessEqual(hi-lo, previous)
            previous = hi-lo
            self.assertEqual(c["status"], "enclosed" if hi-lo <= F(1, 10**20)
                             and c["isolation_complete"] else "open")

    def test_interval_newton_contraction_encloses_irrational_root(self):
        data = ex._prepare([0, 1, 0, -1], (-1, 1))
        lower, upper = ex._newton_contraction(data, F(1, 2), F(3, 4))
        self.assertLessEqual(lower*lower, F(1, 3))
        self.assertGreaterEqual(upper*upper, F(1, 3))
        self.assertLessEqual(upper-lower, F(1, 8))
        self.assertIsNone(ex._newton_contraction(data, F(-1), F(1)))

    def test_float_suggestions_are_only_untrusted_proposals(self):
        q = integrate(from_roots([F(1, 3), F(-2, 5), F(3, 4)]))
        hints = [float("nan"), float("inf"), -10.0, 0.333333333333, 0.750000001, -0.4]
        c = self.checked(q, suggestions=hints, tolerance=F(1, 10**12), max_nodes=64)
        answer = max(ex.evaluate(q, x) for x in [-1, 1, F(1, 3), F(-2, 5), F(3, 4)])
        self.assertLessEqual(F(c["maximum_lower"]), answer)
        self.assertGreaterEqual(F(c["maximum_upper"]), answer)
        for point in c["rational_critical_points"]:
            self.assertEqual(ex.evaluate(ex.derivative(q), F(point["at"])), 0)

    def test_degree32_monotone(self):
        q = [F(1), F(1)] + [F(0)]*30 + [F(1, 100)]
        c = self.checked(q, max_nodes=0)
        self.assertEqual(c["interior_critical_count"], 0)
        self.assertEqual(c["maximum_upper"], "201/100")
        self.assertEqual(c["status"], "exact")

    def test_known_rational_extrema_randomized(self):
        rng = random.Random(9421)
        for _ in range(20):
            roots = [F(rng.randint(-6, 6), 7) for _ in range(rng.randint(2, 7))]
            q = integrate(from_roots(roots))
            c = self.checked(q, suggestions=roots, tolerance=F(1, 10**15), max_nodes=64)
            answer = max(ex.evaluate(q, x) for x in [-1, 1]+roots)
            self.assertEqual(F(c["maximum_lower"]), answer)
            self.assertEqual(F(c["maximum_upper"]), answer)
            self.assertEqual(c["status"], "exact")

    def test_all_bound_components_contain_known_stationary_values(self):
        roots = [F(-2, 3), F(1, 5), F(3, 5)]
        q = integrate(from_roots(roots))
        c = self.checked(q, max_nodes=0)
        for method, bounds in c["partition"][0]["bounds"].items():
            lo, hi = map(F, bounds)
            for root in roots:
                self.assertLessEqual(lo, ex.evaluate(q, root), method)
                self.assertGreaterEqual(hi, ex.evaluate(q, root), method)

    def test_replay_does_not_repeat_search(self):
        c = self.checked([0, 1, 0, -1], max_nodes=8)
        with patch.object(ex, "maximize", side_effect=AssertionError("Search invoked")):
            self.assertTrue(ex.verify(c))

    def test_tampering_and_gaps_are_rejected(self):
        original = self.checked([0, 1, 0, -1], max_nodes=8)
        mutations = []
        for key, value in [("maximum_upper", "0"), ("maximum_lower", "999"),
                           ("error_bound", "0"), ("status", "exact"),
                           ("target_met", not original["target_met"]), ("interior_critical_count", 1),
                           ("refinement_nodes", 0), ("coordinate", "x_squared"),
                           ("formal_assistant_checked", True), ("partition_cells", True),
                           ("max_nodes", 8.0)]:
            c = copy.deepcopy(original)
            c[key] = value
            mutations.append(c)
        c = copy.deepcopy(original)
        c["partition"].pop(1)
        mutations.append(c)
        c = copy.deepcopy(original)
        c["partition"][0]["interval"][1] = "-999/1000"
        mutations.append(c)
        c = copy.deepcopy(original)
        c["partition"][0]["root_count"] += 1
        mutations.append(c)
        c = copy.deepcopy(original)
        c["sturm_sequence"][0][0] = "999"
        mutations.append(c)
        c = copy.deepcopy(original)
        c["derivative_gcd"] = ["2"]
        mutations.append(c)
        c = copy.deepcopy(original)
        c["critical_remainder"] = ["0"]
        mutations.append(c)
        c = copy.deepcopy(original)
        index = next(i for i, cell in enumerate(c["partition"]) if cell["root_count"])
        c["partition"][index]["exact_root"] = "1/2"
        mutations.append(c)
        c = copy.deepcopy(original)
        c["partition"][index]["bounds"]["stationary_taylor"][1] = "0"
        mutations.append(c)
        c = copy.deepcopy(original)
        c["sample_witness"]["value"] = "1"
        mutations.append(c)
        c = copy.deepcopy(original)
        c["unexpected"] = "discard me"
        mutations.append(c)
        for i, c in enumerate(mutations):
            self.assertFalse(ex.verify(c), f"Mutation {i} accepted")

    def test_boundary_root_cannot_be_omitted(self):
        q = integrate(from_roots([-F(1, 2), F(0), F(1, 2)]))
        c = self.checked(q, suggestions=[0, -F(1, 2), F(1, 2)], max_nodes=8)
        altered = copy.deepcopy(c)
        altered["rational_critical_points"] = []
        self.assertFalse(ex.verify(altered))

    def test_malformed_inputs_and_certificates(self):
        for q in [[], [0]*34, [0.5, 1], [True, 1], ["NaN"], None]:
            with self.assertRaises((TypeError, ValueError)):
                ex.maximize(q)
        for options in [dict(tolerance=0), dict(tolerance=-1), dict(max_nodes=-1),
                        dict(max_nodes=4097), dict(max_nodes=True), dict(interval=(2, -1)),
                        dict(interval=(0.0, 1)), dict(suggestions=[0]*257)]:
            with self.assertRaises((TypeError, ValueError)):
                ex.maximize([1, 1], **options)
        for cert in [None, [], {}, {"format": ex.FORMAT}, True]:
            self.assertFalse(ex.verify(cert))


def comparison_rows():
    """Small deterministic width/node comparison; deliberately no elapsed times."""
    import family_v29 as legacy
    cases = [
        ("cubic_irrational", [0, 1, 0, -1]),
        ("flat_peak_one_third", [-F(comb(8, k))*(-F(1, 3))**(8-k) for k in range(9)]),
        ("two_unequal_peaks", integrate([-v for v in from_roots([F(-4, 5), F(-1, 4), F(2, 3)])])),
        ("even_octic_flat_peak", [F(-1, 81), 0, F(4, 27), 0, F(-2, 3), 0, F(4, 3), 0, -1]),
    ]
    rows = []
    for name, q in cases:
        old = legacy.maximum(q, F(1, 10**12), max_leaves=33)
        new = ex.maximize(q, tolerance=F(1, 10**12), max_nodes=32)
        if not legacy.verify_range(old) or not ex.verify(new):
            raise AssertionError("Comparison certificate failed replay")
        old_nodes = 0 if old["intervals"] is None else len(old["intervals"])-1
        rows.append(dict(case=name, tolerance="1/1000000000000", allowed_refinements=32,
                         legacy_nodes=old_nodes, exact_nodes=new["refinement_nodes"],
                         legacy_width=str(F(old["maximum_upper"])-F(old["maximum_lower"])),
                         exact_width=new["error_bound"], legacy_status=old["status"],
                         exact_status=new["status"]))
    return rows


class ComparisonTests(unittest.TestCase):
    def test_recorded_comparison(self):
        rows = comparison_rows()
        self.assertEqual(len(rows), 4)
        self.assertTrue(all(F(r["exact_width"]) <= F(r["legacy_width"]) for r in rows))
        self.assertTrue(any(F(r["exact_width"]) < F(r["legacy_width"]) for r in rows))


if __name__ == "__main__":
    unittest.main()
