"""V88: bind exact extrema to positive functions and Poisson matrix envelopes.

This is a new certificate format. Legacy verifiers are used only for legacy
evidence, never to accept an exact-Sturm range as a legacy Bernstein proof.
"""

from copy import deepcopy
from fractions import Fraction as F

import exact_extrema as extrema
import family_v29 as family
import dual_v33 as dual


FAMILY_METHOD = "precise_poisson_matrix_v88"
MATRIX_METHOD = "precise_global_matrix_v88"
FAMILY_SCOPE = "full unit sphere; all real amplitude vectors; fixed polynomial Poisson exponential ansatz"
MATRIX_SCOPE = "minimum trace(C G) over Poisson matrix envelopes on the full interval [-1,1]"
LIMITATION = "Only Poisson-ansatz optimality is enclosed; this is not a sharp spectral-inequality optimum."


def _rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, F)):
        raise ValueError("Exact integers, Fraction or rational strings required")
    return F(value)


def _budget(tolerance, max_nodes):
    tolerance = _rational(tolerance)
    if tolerance <= 0:
        raise ValueError("Strictly positive rational tolerance required")
    if type(max_nodes) is not int or not 0 <= max_nodes <= extrema.MAX_NODES:
        raise ValueError("Refinement budget must be an integer in 0..4096")
    return tolerance


def _same(actual, expected):
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, dict):
        return actual.keys() == expected.keys() and all(_same(actual[k], expected[k]) for k in expected)
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(_same(a, b) for a, b in zip(actual, expected))
    return actual == expected


class ExactRangeBackend:
    """Trusted positive_search hook with exact polynomial AND domain binding."""

    name = "exact_sturm_minimum_v1"

    def minimum(self, p, tolerance, max_leaves):
        if type(max_leaves) is not int or not 1 <= max_leaves <= extrema.MAX_NODES+1:
            raise ValueError("Leaf budget must be an integer in 1..4097")
        negative = [-v for v in extrema.polynomial(p)]
        return extrema.maximize(negative, interval=(-1, 1), tolerance=tolerance,
                                max_nodes=max_leaves-1)

    def verify_minimum(self, p, proof):
        negative = [-v for v in extrema.polynomial(p)]
        if (not extrema.verify(proof)
                or proof["polynomial"] != list(map(str, negative))
                or proof["interval"] != ["-1", "1"]):
            raise ValueError("Exact range is not bound to this polynomial on [-1,1]")
        return -F(proof["maximum_upper"]), -F(proof["maximum_lower"])


def _assemble(directions, weight):
    if not isinstance(directions, (list, tuple)) or not 1 <= len(directions) <= 4:
        raise ValueError("One to four polynomial directions required")
    directions = [extrema.polynomial(p) for p in directions]
    if not isinstance(weight, (list, tuple)):
        raise ValueError("Rational weight matrix required")
    weight = [[_rational(v) for v in row] for row in weight]
    # Reuse exactly the established Poisson identity and SPD weight check.
    return family.assemble(directions, weight)


def _family_record(directions, weight, proof):
    dirs, w, r, means, polynomial = _assemble(directions, weight)
    if (not extrema.verify(proof)
            or proof["polynomial"] != list(map(str, polynomial))
            or proof["interval"] != ["-1", "1"]):
        raise ValueError("Extremum certificate does not match the Poisson envelope")
    lower, upper = F(proof["maximum_lower"]), F(proof["maximum_upper"])
    if upper < 0:
        raise ValueError("A rank-one envelope has nonnegative maximum")
    return {
        "method": FAMILY_METHOD,
        "directions": [list(map(str, p)) for p in dirs], "weight": dual.encode(w),
        "poisson_solutions": [list(map(str, p)) for p in r], "means": list(map(str, means)),
        "envelope_polynomial": list(map(str, polynomial)), "interval": ["-1", "1"],
        "scale": str(upper), "scale_lower": str(lower), "scale_upper": str(upper),
        "scale_error_bound": str(upper-lower),
        "quadratic_bound": dual.encode([[upper*v for v in row] for row in w]),
        "range_proof": proof, "range_status": proof["status"],
        "range_tolerance": proof["tolerance"], "max_nodes": proof["max_nodes"],
        "scope": FAMILY_SCOPE, "limitation": LIMITATION, "formal_assistant_checked": False,
    }


def precise_family(directions, weight, tolerance=F(1, 10**12), max_nodes=256):
    """Certify the same Poisson family with a separately verified Sturm range.

    An open range remains a valid conservative quadratic bound. Its width and
    incomplete isolation are retained explicitly in ``range_proof``.
    """
    tolerance = _budget(tolerance, max_nodes)
    dirs, w, _, _, polynomial = _assemble(directions, weight)
    proof = extrema.maximize(polynomial, tolerance=tolerance, max_nodes=max_nodes)
    return _family_record(dirs, w, proof)


def verify_family(certificate):
    """Recompute Poisson identities, SPD weight, bound polynomial and range."""
    try:
        if not isinstance(certificate, dict) or certificate.get("method") != FAMILY_METHOD:
            return False
        expected = _family_record(certificate["directions"], certificate["weight"],
                                  certificate["range_proof"])
        return _same(certificate, expected)
    except (ArithmeticError, AttributeError, IndexError, KeyError, TypeError, ValueError):
        return False


def _legacy_seed(seed_result):
    if not isinstance(seed_result, dict):
        raise ValueError("A verified legacy matrix certificate or result wrapper is required")
    seed = seed_result.get("certificate", seed_result)
    if not isinstance(seed, dict) or not dual.verify(seed):
        raise ValueError("Invalid legacy primal/dual certificate")
    # Result statistics, seconds, warm samples and solver flags have no authority.
    return seed


def _matrix_record(seed, candidate, tolerance, max_nodes):
    tolerance = _budget(tolerance, max_nodes)
    seed = _legacy_seed(seed)
    if not verify_family(candidate):
        raise ValueError("Invalid precise Poisson candidate")
    dirs, _, cost = dual.setup(seed["dual"]["directions"], seed["dual"]["cost"])
    canonical_dirs = [list(map(str, p)) for p in dirs]
    if (candidate["directions"] != canonical_dirs
            or seed["template"]["directions"] != canonical_dirs
            or candidate["weight"] != seed["template"]["weight"]):
        raise ValueError("Candidate must use the seed's original directions and retained best weight")
    weight = [list(map(F, row)) for row in candidate["weight"]]
    multiplier = dual.inner(cost, weight)
    if multiplier <= 0:
        raise ValueError("Positive trace(C W) required")
    range_tolerance = tolerance/multiplier
    if candidate["range_tolerance"] != str(range_tolerance) or candidate["max_nodes"] != max_nodes:
        raise ValueError("Candidate range budget is not bound to objective tolerance")
    lower = F(seed["dual"]["lower_bound"])
    previous = F(seed["upper_bound"])
    new_upper = dual.inner(cost, [list(map(F, row)) for row in candidate["quadratic_bound"]])
    if lower > min(previous, new_upper):
        raise ValueError("Contradictory verified primal and dual bounds")
    use_new = new_upper <= previous
    upper = min(new_upper, previous)
    retained = candidate["quadratic_bound"] if use_new else seed["template"]["quadratic_bound"]
    outcome = "improved" if new_upper < previous else "tied" if new_upper == previous else "negative_gain_seed_retained"
    return {
        "method": MATRIX_METHOD, "directions": canonical_dirs, "cost": dual.encode(cost),
        "weight": candidate["weight"], "seed_certificate": seed, "dual": seed["dual"],
        "precise_family": candidate, "quadratic_bound": retained,
        "lower_bound": str(lower), "upper_bound": str(upper), "gap": str(upper-lower),
        "previous_upper_bound": str(previous), "candidate_upper_bound": str(new_upper),
        "upper_improvement": str(previous-upper), "candidate_upper_improvement": str(previous-new_upper),
        "precision_outcome": outcome, "upper_source": "precise_family" if use_new else "legacy_seed",
        "weight_cost_multiplier": str(multiplier), "range_tolerance": str(range_tolerance),
        "range_objective_error_bound": str(multiplier*F(candidate["scale_error_bound"])),
        "range_status": candidate["range_status"], "tolerance": str(tolerance), "max_nodes": max_nodes,
        "status": "global_gap_closed" if upper-lower <= tolerance else "global_gap_open",
        "objective": "trace(C G); G dominates (1-t^2) v(t)v(t)^T for all t in [-1,1]",
        "scope": MATRIX_SCOPE, "limitation": LIMITATION, "formal_assistant_checked": False,
    }


def precise_matrix(seed_result, tolerance=F(1, 10**12), max_nodes=256):
    """Re-enclose the seed's best weight, retaining its verified dual lower.

    Accept a matrix_refine result wrapper or a legacy dual_v33 certificate.
    The tighter of the old and new uppers is retained. Negative candidate gain
    is reported explicitly. No weight optimization or new dual search occurs.
    """
    tolerance = _budget(tolerance, max_nodes)
    seed = deepcopy(_legacy_seed(seed_result))
    weight = [list(map(F, row)) for row in seed["template"]["weight"]]
    cost = [list(map(F, row)) for row in seed["dual"]["cost"]]
    multiplier = dual.inner(cost, weight)
    candidate = precise_family(seed["template"]["directions"], weight,
                               tolerance=tolerance/multiplier, max_nodes=max_nodes)
    return _matrix_record(seed, candidate, tolerance, max_nodes)


def verify_matrix(certificate):
    """Replay old dual/primal and new range; check directions, cost and gap.

    No search or ``precise_family``/``precise_matrix`` synthesis is repeated.
    A closed gap concerns the Poisson envelope optimization only.
    """
    try:
        if not isinstance(certificate, dict) or certificate.get("method") != MATRIX_METHOD:
            return False
        expected = _matrix_record(certificate["seed_certificate"], certificate["precise_family"],
                                  certificate["tolerance"], certificate["max_nodes"])
        return _same(certificate, expected)
    except (ArithmeticError, AttributeError, IndexError, KeyError, TypeError, ValueError):
        return False


def verify(certificate):
    if not isinstance(certificate, dict):
        return False
    if certificate.get("method") == FAMILY_METHOD:
        return verify_family(certificate)
    if certificate.get("method") == MATRIX_METHOD:
        return verify_matrix(certificate)
    return False
