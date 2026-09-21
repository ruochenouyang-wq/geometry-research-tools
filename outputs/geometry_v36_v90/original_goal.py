"""V89: apply verified scalar uniform lemmas to the original structured goal.

The parameter directions must be affine in the unit-sphere coordinate t.
The fixed base potential may have degree up to six. No pointwise local-energy
upper value is ever used as a spectral upper bound: the upper witness is u=1.
"""
from fractions import Fraction as F
from itertools import product
import copy

import algebra
import compiler_v23 as compiler
import family_v29 as family
import uniform_refine as uniform


METHOD = "original_spectral_goal_v89"


def _exact(value):
    if isinstance(value, bool) or not isinstance(value, (str, int, F)):
        raise ValueError("Use an exact rational coefficient")
    return F(value)


def _trim(values):
    values = list(values)
    while len(values) > 1 and not values[-1]:
        values.pop()
    return values or [F(0)]


def _encode(values):
    return list(map(str, values))


def _data(problem):
    problem = compiler.normalized(problem)
    q0 = list(map(F, problem["q0"]))
    directions = [list(map(F, row)) for row in problem["directions"]]
    if any(any(row[2:]) for row in directions):
        raise ValueError("V89 only supports parameter directions affine in t")
    penalty = problem["penalty"]
    constant = F(penalty["constant"])
    linear = list(map(F, penalty["linear"]))
    hessian = [list(map(F, row)) for row in penalty["hessian"]]
    return q0, directions, constant, linear, hessian


def _floor(constant, linear, hessian, box):
    """A failure of this relaxation is never called a spectral refutation."""
    try:
        low, where = algebra.quadratic_min(constant, linear, hessian, box)
        return low, where, {"kind": "finite_exact_quadratic_minimum"}
    except ValueError:
        if box is not None:
            raise
        inertia = algebra.inertia(hessian)
        if inertia[0]:
            reason = "negative_hessian_direction"
        else:
            try:
                algebra.rectangular(hessian, [-x for x in linear], len(linear))
            except ValueError:
                reason = "linear_term_in_hessian_nullspace"
            else:
                raise
        return None, None, {"kind": "quadratic_relaxation_unbounded_below",
                            "reason": reason, "hessian_inertia": inertia,
                            "meaning": "This lower relaxation gives no finite floor; the spectral goal is not thereby refuted."}


def _candidate_points(problem, relaxation_minimizer, trial_quadratic):
    dimensions = len(problem["directions"])
    box = ([list(map(F, row)) for row in problem["domain"]["box"]]
           if problem["domain"]["kind"] == "box" else None)
    candidates = []

    def add(point):
        point = list(map(F, point))
        if len(point) != dimensions:
            raise ValueError("Candidate dimension mismatch")
        if box is not None and any(not lo <= x <= hi for x, (lo, hi) in zip(point, box)):
            raise ValueError("Candidate outside the original parameter domain")
        if point not in candidates:
            candidates.append(point)

    if relaxation_minimizer is not None:
        add(relaxation_minimizer)
    if box is None:
        add([F(0)] * dimensions)
        for i in range(dimensions):
            for sign in (-1, 1):
                point = [F(0)] * dimensions
                point[i] = F(sign)
                add(point)
    else:
        add([min(max(F(0), lo), hi) for lo, hi in box])
        add([(lo + hi) / 2 for lo, hi in box])
        for point in product(*box):
            add(point)
    try:
        _, trial_minimizer = algebra.quadratic_min(*trial_quadratic, box)
        add(trial_minimizer)
    except ValueError:
        # Finitely many candidate points suffice for valid upper witnesses.
        # Failure to minimize their quadratic is not an original lower claim.
        pass
    return candidates


def _constant_trial(problem, parameters):
    """Reassemble the original potential and integrate the original u=1 trial."""
    q0, directions, constant, linear, hessian = _data(problem)
    parameters = list(map(F, parameters))
    if len(parameters) != len(directions):
        raise ValueError("Trial parameter dimension mismatch")
    if problem["domain"]["kind"] == "box":
        for value, (lo, hi) in zip(parameters, problem["domain"]["box"]):
            if not F(lo) <= value <= F(hi):
                raise ValueError("Trial lies outside the original parameter domain")
    degree = max([len(q0)] + [len(row) for row in directions])
    potential = q0 + [F(0)] * (degree - len(q0))
    for amplitude, direction in zip(parameters, directions):
        for j, value in enumerate(direction):
            potential[j] += amplitude * value
    potential = _trim(potential)
    mean = sum((value / (j + 1) for j, value in enumerate(potential) if j % 2 == 0), F(0))
    penalty = algebra.value(constant, linear, hessian, parameters)
    return {"parameters": _encode(parameters), "trial_power_coefficients": ["1"],
            "original_potential": _encode(potential),
            "probability_mass": "1", "dirichlet_energy": "0",
            "potential_mean": str(mean), "original_penalty": str(penalty),
            "spectral_rayleigh_upper": str(mean), "objective_upper": str(mean + penalty),
            "meaning": "The original constant trial gives an upper bound on the original objective infimum."}


def _build(compiled, coefficient, rule, lemma, remainder_range):
    if not compiler.verify(compiled):
        raise ValueError("Invalid original structured goal compilation")
    coefficient = _exact(coefficient)
    if rule not in ("log_sobolev", "barta"):
        raise ValueError("Choose log_sobolev or barta")
    if not uniform.is_proved(lemma):
        raise ValueError("A verified successful uniform spectral lemma is required")
    if lemma["direction"] != ["0", "1"] or lemma["coefficient"] != str(coefficient):
        raise ValueError("Uniform lemma is not bound to the normalized linear coordinate and coefficient")
    if rule == "log_sobolev":
        if lemma["method"] != uniform.LOG_SOBOLEV_METHOD:
            raise ValueError("The requested analytic rule does not match the lemma")
        trust = {"classical_theorem": copy.deepcopy(lemma["trusted_analytic_theorem"]),
                 "other_rules": ["spectral monotonicity under pointwise potential domination",
                                 "exact parameter quadratic minimization", "Rayleigh variational upper bound"],
                 "formal_assistant_checked": False}
    else:
        if lemma["method"] != uniform.METHOD or lemma["all_real"] is not True:
            raise ValueError("Barta transport requires an all-real lemma, not a finite amplitude box")
        trust = {"classical_theorem": None,
                 "other_rules": ["positive-function quadratic-form identity",
                                 "full tensor Bernstein cover and all-real tail gluing",
                                 "spectral monotonicity under pointwise potential domination",
                                 "exact parameter quadratic minimization", "Rayleigh variational upper bound"],
                 "formal_assistant_checked": False}
    problem = compiled["problem"]
    q0, directions, constant, linear, hessian = _data(problem)
    dimensions = len(directions)
    offset = q0[1] if len(q0) > 1 else F(0)
    slopes = [row[1] if len(row) > 1 else F(0) for row in directions]
    parameter_constants = [row[0] for row in directions]
    remainder = q0[:]
    if len(remainder) > 1:
        remainder[1] = F(0)
    remainder = _trim(remainder)
    if not family.verify_range(remainder_range):
        raise ValueError("Invalid full-interval remainder range")
    if remainder_range["polynomial"] != _encode([-x for x in remainder]):
        raise ValueError("Remainder range is not bound to minus the actual remainder")
    remainder_lower = -F(remainder_range["maximum_upper"])
    new_constant = constant + remainder_lower - coefficient * offset ** 2
    new_linear = [linear[i] + parameter_constants[i] - 2 * coefficient * offset * slopes[i]
                  for i in range(dimensions)]
    new_hessian = [[hessian[i][j] - 2 * coefficient * slopes[i] * slopes[j]
                    for j in range(dimensions)] for i in range(dimensions)]
    box = ([list(map(F, row)) for row in problem["domain"]["box"]]
           if problem["domain"]["kind"] == "box" else None)
    low, where, diagnostic = _floor(new_constant, new_linear, new_hessian, box)
    mean_q0 = sum((value / (j + 1) for j, value in enumerate(q0) if j % 2 == 0), F(0))
    trial_quadratic = (constant + mean_q0,
                       [linear[i] + parameter_constants[i] for i in range(dimensions)], hessian)
    points = _candidate_points(problem, where, trial_quadratic)
    trials = [_constant_trial(problem, point) for point in points]
    best = min(range(len(trials)), key=lambda i: F(trials[i]["objective_upper"]))
    trial = trials[best]
    upper = F(trial["objective_upper"])
    threshold = F(compiled["threshold"])
    if low is not None and low > upper:
        raise ArithmeticError("Certified lower bound exceeds the original Rayleigh upper witness")
    status = ("proved" if low is not None and low >= threshold else
              "refuted" if upper < threshold else "unresolved")
    return {"method": METHOD, "algorithm_version": 89, "compiled_goal": compiled,
            "rule": rule, "coefficient": str(coefficient), "uniform_lemma": lemma,
            "operator_binding": {"geometry": "unit_sphere", "laplacian": "nonnegative -Delta",
                                 "function_space": problem["function_space"], "eigenvalue_index": 1,
                                 "measure": "uniform sphere probability measure; Rayleigh normalization invariant",
                                 "scope_transfer": "The full-sphere lower lemma also bounds its axisymmetric restriction; u=1 belongs to both spaces."},
            "amplitude_offset": str(offset), "amplitude_slopes": _encode(slopes),
            "parameter_constants": _encode(parameter_constants),
            "remainder_polynomial": _encode(remainder), "remainder_range": remainder_range,
            "remainder_lower": str(remainder_lower),
            "derived_quadratic": {"constant": str(new_constant), "linear": _encode(new_linear),
                                  "hessian": [_encode(row) for row in new_hessian]},
            "lower_bound": None if low is None else str(low),
            "relaxation_minimizer": None if where is None else _encode(where),
            "lower_diagnostic": diagnostic,
            "upper_bound": str(upper), "candidate_parameters": trial["parameters"],
            "constant_trial": trial, "constant_trial_candidates": trials,
            "candidate_index": best,
            "infimum_exact": low is not None and low == upper,
            "infimum_gap": None if low is None else str(upper - low),
            "threshold": str(threshold), "relation": ">=", "status": status,
            "trust": trust, "formal_assistant_checked": False}


def solve(raw, coefficient=None, rule="log_sobolev"):
    """Return the original-goal certificate, with honest unresolved diagnostics.

    Unsupported geometry, index, language or nonaffine parameter directions
    raise ValueError. A valid but insufficient lemma gives unresolved, not a
    claimed proof. The optional coefficient must itself have an accepted
    all-real uniform lemma; failed local-energy candidates are rejected.
    """
    compiled = compiler.compile_goal(raw)
    q0, _, _, _, _ = _data(compiled["problem"])
    if rule == "log_sobolev":
        coefficient = F(1, 6) if coefficient is None else _exact(coefficient)
        lemma = uniform.synthesize_log_sobolev(coefficient=coefficient)
    elif rule == "barta":
        coefficient = F(1, 5) if coefficient is None else _exact(coefficient)
        lemma = uniform.synthesize(coefficient=coefficient, order=3, all_real=True)
    else:
        raise ValueError("Choose log_sobolev or barta")
    if not uniform.is_proved(lemma):
        raise ValueError("The requested coefficient has no accepted all-real lemma in this rule")
    remainder = q0[:]
    if len(remainder) > 1:
        remainder[1] = F(0)
    remainder = _trim(remainder)
    range_proof = family.maximum([-x for x in remainder], epsilon=F(1, 10**6), max_leaves=64)
    cert = _build(compiled, coefficient, rule, lemma, range_proof)
    if not verify(cert):
        raise ArithmeticError("Generated original-goal certificate failed verification")
    return cert


def verify(certificate):
    """Recompile the original source and replay every binding and accepted lemma.

    True means the evidence and status are valid, including unresolved/refuted.
    To accept the original universal inequality, also require status=='proved'.
    """
    try:
        return certificate == _build(certificate["compiled_goal"], certificate["coefficient"],
                                     certificate["rule"], certificate["uniform_lemma"],
                                     certificate["remainder_range"])
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError, ArithmeticError,
            RecursionError, OverflowError):
        return False


if __name__ == "__main__":
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("goal", help="Path to the original structured goal JSON")
    parser.add_argument("--rule", choices=["log_sobolev", "barta"], default="log_sobolev")
    parser.add_argument("--coefficient")
    args = parser.parse_args()
    with open(args.goal, encoding="utf-8") as handle:
        raw_goal = json.load(handle)
    print(json.dumps(solve(raw_goal, coefficient=args.coefficient, rule=args.rule),
                     ensure_ascii=False, sort_keys=True, indent=2))
