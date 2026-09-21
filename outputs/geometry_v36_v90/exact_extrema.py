"""Exact univariate maximum certificates, independent of the legacy format.

Coefficients are in ascending power order. Only Python's standard library is
used. ``maximize`` searches; ``verify`` replays the mathematical evidence, not
the search. See PROOF_EXACT_EXTREMA.md for scope, guarantees and resource limits.
"""

from fractions import Fraction as F
from math import comb, isfinite


FORMAT = "exact_sturm_extrema_v1"
MAX_DEGREE = 32
MAX_NODES = 4096


def _rational(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, F)):
        raise ValueError("Use integers, Fraction, or rational strings")
    return F(value)


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p or [F(0)]


def polynomial(raw):
    if not isinstance(raw, (list, tuple)) or not 1 <= len(raw) <= MAX_DEGREE + 1:
        raise ValueError("One to 33 rational coefficients required")
    return trim([_rational(x) for x in raw])


def evaluate(p, x):
    result = F(0)
    for c in reversed(p):
        result = result*x + c
    return result


def derivative(p):
    return trim([i*p[i] for i in range(1, len(p))])


def divide(a, b):
    """Exact Euclidean polynomial division in Q[x]."""
    a, b = trim(map(F, a)), trim(map(F, b))
    if b == [0]:
        raise ZeroDivisionError("Zero divisor polynomial")
    quotient = [F(0)] * max(1, len(a)-len(b)+1)
    while a != [0] and len(a) >= len(b):
        shift, ratio = len(a)-len(b), a[-1]/b[-1]
        quotient[shift] += ratio
        for i, c in enumerate(b):
            a[i+shift] -= ratio*c
        a = trim(a)
    return trim(quotient), a


def monic(p):
    p = trim(map(F, p))
    return [c/p[-1] for c in p] if p != [0] else p


def gcd(a, b):
    a, b = trim(a), trim(b)
    while b != [0]:
        a, b = b, divide(a, b)[1]
    return monic(a)


def square_free(p):
    p = trim(p)
    if len(p) == 1:
        return [F(1)], [F(0) if p == [0] else F(1)]
    common = gcd(p, derivative(p))
    quotient, remainder = divide(p, common)
    if remainder != [0]:
        raise ArithmeticError("Inexact square-free division")
    return monic(quotient), common


def _positive_normalize(p):
    # Positive rescaling preserves every Sturm sign, including negative leads.
    return [F(c)/abs(p[-1]) for c in p]


def sturm_sequence(p):
    p = trim(p)
    if len(p) <= 1:
        return []
    result = [_positive_normalize(p), _positive_normalize(derivative(p))]
    while len(result[-1]) > 1:
        remainder = divide(result[-2], result[-1])[1]
        if remainder == [0]:
            raise ValueError("Sturm input must be square free")
        result.append(_positive_normalize([-v for v in remainder]))
    return result


def _side_sign(p, x, side):
    """Sign immediately to one side, even if p(x)=0."""
    order = 0
    while p != [0]:
        value = evaluate(p, x)
        if value:
            return (1 if value > 0 else -1) * (side if order % 2 else 1)
        p = derivative(p)
        order += 1
    return 0


def variations(sequence, x, side):
    if side not in (-1, 1):
        raise ValueError("Side must be -1 or +1")
    signs = [_side_sign(p, x, side) for p in sequence]
    signs = [s for s in signs if s]
    return sum(a != b for a, b in zip(signs, signs[1:]))


def count_open(sequence, left, right):
    """Number of distinct roots in the OPEN interval (left, right)."""
    if not left < right:
        raise ValueError("Strict interval required")
    count = variations(sequence, left, 1)-variations(sequence, right, -1)
    if count < 0:
        raise ArithmeticError("Negative Sturm count")
    return count


def bernstein_bounds(p, left, right):
    """Convex-hull enclosure on a closed rational interval."""
    n, width = len(p)-1, right-left
    power = [sum((p[j]*comb(j, k)*left**(j-k)*width**k
                  for j in range(k, n+1)), F(0)) for k in range(n+1)]
    values = [sum((power[j]*F(comb(i, j), comb(n, j))
                   for j in range(i+1)), F(0)) for i in range(n+1)]
    return min(values), max(values)


def stationary_bounds(p, left, right):
    """Enclose p(r) for every stationary r in [left,right].

    If p(m+z)=sum a_k z^k and p'(m+z)=0, then
    p(m+z)=a_0-sum_{k>=2}(k-1)*a_k*z^k. This removes the
    linear term exactly and gives a quadratic-order error enclosure.
    """
    middle, radius = (left+right)/2, (right-left)/2
    lower = upper = evaluate(p, middle)
    for k in range(2, len(p)):
        coefficient = -(k-1)*sum((p[j]*comb(j, k)*middle**(j-k)
                                  for j in range(k, len(p))), F(0))
        amount = coefficient*radius**k
        if k % 2:
            lower -= abs(amount)
            upper += abs(amount)
        else:
            lower += min(F(0), amount)
            upper += max(F(0), amount)
    return lower, upper


def _strings(values):
    return [str(v) for v in values]


def _prepare(raw, interval):
    original = polynomial(raw)
    if not isinstance(interval, (list, tuple)) or len(interval) != 2:
        raise ValueError("A pair of rational interval endpoints is required")
    a, b = map(_rational, interval)
    if a > b:
        raise ValueError("Interval must be ordered")
    if len(original) > 1 and not any(original[1::2]):
        q = trim(original[::2])
        left = F(0) if a <= 0 <= b else min(a*a, b*b)
        right = max(a*a, b*b)
        coordinate = "x_squared"
    else:
        q, left, right, coordinate = original, a, b, "x"
    d = derivative(q)
    sf, common = square_free(d)
    seq = sturm_sequence(sf)
    rem = divide(q, sf)[1] if len(sf) > 1 else q
    return dict(original=original, interval=(a, b), q=q, left=left,
                right=right, coordinate=coordinate, d=d, sf=sf,
                common=common, sequence=seq, remainder=rem)


def _exact_root(data, left, right, count, suggestions=()):
    if count != 1:
        return None
    sf = data["sf"]
    candidates = [(left+right)/2]
    if len(sf) == 2:
        candidates.append(-sf[0]/sf[1])
    candidates.extend(suggestions)
    for x in candidates:
        if left < x < right and evaluate(sf, x) == 0:
            return x
    return None


def _cell(data, left, right, exact_root=None):
    count = count_open(data["sequence"], left, right)
    result = {"interval": _strings((left, right)), "root_count": count,
              "kind": "empty" if not count else "isolated" if count == 1 else "cluster",
              "exact_root": None, "bounds": None, "critical_value_bounds": None}
    if exact_root is not None:
        if count != 1 or not left < exact_root < right or evaluate(data["sf"], exact_root):
            raise ValueError("Invalid rational critical point")
        result["exact_root"] = str(exact_root)
    if count:
        bounds = {"bernstein": bernstein_bounds(data["q"], left, right),
                  "stationary_taylor": stationary_bounds(data["q"], left, right),
                  "critical_remainder": bernstein_bounds(data["remainder"], left, right)}
        if exact_root is not None:
            value = evaluate(data["q"], exact_root)
            bounds["rational_root"] = (value, value)
        lower = max(v[0] for v in bounds.values())
        upper = min(v[1] for v in bounds.values())
        if lower > upper:
            raise ArithmeticError("Contradictory critical-value bounds")
        result["bounds"] = {k: _strings(v) for k, v in bounds.items()}
        result["critical_value_bounds"] = _strings((lower, upper))
    return result


def _multiplicity(p, x):
    if p == [0]:
        return None  # The identically zero derivative is handled separately.
    result = 0
    while p != [0] and evaluate(p, x) == 0:
        result += 1
        p = derivative(p)
    return result


def _summary(data, cells):
    left, right, q = data["left"], data["right"], data["q"]
    samples = {left, right}
    boundary = {left, right}
    exact_roots = set()
    for cell in cells:
        a, b = map(F, cell["interval"])
        samples.update((a, (a+b)/2, b))
        boundary.update((a, b))
        if cell["exact_root"] is not None:
            exact_roots.add(F(cell["exact_root"]))
    if data["d"] != [0]:
        exact_roots.update(x for x in boundary if evaluate(data["d"], x) == 0)
    samples.update(exact_roots)
    # Ties choose the smallest coordinate for deterministic serialization.
    at = max(sorted(samples), key=lambda x: evaluate(q, x))
    witness = {"coordinate": data["coordinate"], "at": str(at), "value": str(evaluate(q, at))}
    lower = F(witness["value"])
    upper = max(evaluate(q, left), evaluate(q, right))
    for x in exact_roots:
        upper = max(upper, evaluate(q, x))
    lower_source = {"kind": "rational_coordinate", **witness}
    for i, cell in enumerate(cells):
        if cell["root_count"]:
            lo, hi = map(F, cell["critical_value_bounds"])
            upper = max(upper, hi)
            if lo > lower:
                lower = lo
                lower_source = {"kind": "critical_interval", "cell_index": i, "value": str(lo)}
    if lower > upper:
        raise ArithmeticError("Contradictory global bounds")
    roots = [{"at": str(x), "value": str(evaluate(q, x)),
              "derivative_multiplicity": _multiplicity(data["d"], x)} for x in sorted(exact_roots)]
    # Interior roots on split boundaries must not be lost or double counted.
    boundary_count = sum(left < x < right for x in boundary
                         if data["d"] != [0] and evaluate(data["d"], x) == 0)
    covered_count = sum(c["root_count"] for c in cells) + boundary_count
    total_count = count_open(data["sequence"], left, right) if left < right else 0
    if covered_count != total_count:
        raise ArithmeticError("Incomplete critical-point coverage")
    return lower, upper, witness, lower_source, roots, total_count


def _assemble(data, cells, tolerance, budget):
    lower, upper, witness, source, roots, total_count = _summary(data, cells)
    isolated = all(c["root_count"] <= 1 for c in cells)
    target_met = upper-lower <= tolerance
    met = isolated and target_met
    status = "exact" if met and lower == upper else "enclosed" if met else "open"
    nodes = max(0, len(cells)-1)
    return {
        "format": FORMAT, "coefficient_order": "ascending",
        "polynomial": _strings(data["original"]), "interval": _strings(data["interval"]),
        "coordinate": data["coordinate"], "working_polynomial": _strings(data["q"]),
        "working_interval": _strings((data["left"], data["right"])),
        "derivative": _strings(data["d"]), "derivative_gcd": _strings(data["common"]),
        "critical_square_free": _strings(data["sf"]),
        "critical_remainder": _strings(data["remainder"]),
        "sturm_sequence": [_strings(p) for p in data["sequence"]],
        "partition": cells, "rational_critical_points": roots,
        "interior_critical_count": total_count, "constant_polynomial": len(data["q"]) == 1,
        "maximum_lower": str(lower), "maximum_upper": str(upper),
        "error_bound": str(upper-lower), "sample_witness": witness, "lower_source": source,
        "tolerance": str(tolerance), "target_met": target_met,
        "isolation_complete": isolated, "status": status,
        "max_nodes": budget, "refinement_nodes": nodes, "partition_cells": len(cells),
        "budget_exhausted": not met and nodes >= budget,
        "formal_assistant_checked": False,
    }


def _validate_budget(tolerance, budget):
    tolerance = _rational(tolerance)
    if tolerance <= 0:
        raise ValueError("Strictly positive rational tolerance required")
    if type(budget) is not int or not 0 <= budget <= MAX_NODES:
        raise ValueError("max_nodes must be an integer in 0..4096")
    return tolerance


def _suggestions(raw, data):
    """Rationalize optional ORIGINAL-coordinate hints; all claims are rechecked."""
    if raw is None:
        return []
    if not isinstance(raw, (list, tuple)) or len(raw) > 256:
        raise ValueError("At most 256 optional suggestions allowed")
    result = set()
    for value in raw:
        if isinstance(value, float):
            if not isfinite(value):
                continue
            value = F(value).limit_denominator(2**30)
        else:
            value = _rational(value)
        if data["interval"][0] <= value <= data["interval"][1]:
            if data["coordinate"] == "x_squared":
                value *= value
            if data["left"] < value < data["right"]:
                result.add(value)
    return sorted(result)


def _newton_contraction(data, left, right):
    """Exact interval Newton proposal; replay still uses Sturm root counts.

    Rounding out to a dyadic grid controls denominator growth. A contraction
    is optional and need not be trusted by the eventual verifier.
    """
    dlo, dhi = bernstein_bounds(derivative(data["sf"]), left, right)
    if dlo <= 0 <= dhi:
        return None
    middle, width = (left+right)/2, right-left
    value = evaluate(data["sf"], middle)
    candidates = (middle-value/dlo, middle-value/dhi)
    lower, upper = max(left, min(candidates)), min(right, max(candidates))
    if lower > upper:
        raise ArithmeticError("Newton proposal lost the certified root")
    if lower == upper and evaluate(data["sf"], lower) == 0:
        return lower, upper
    bits = min(512, 2*max(4, width.denominator.bit_length()-width.numerator.bit_length()+1))
    scale = 2**bits
    lower = max(left, F((lower*scale).__floor__(), scale))
    upper = min(right, F((upper*scale).__ceil__(), scale))
    if lower <= upper and upper-lower <= width/2:
        return lower, upper
    return None


def maximize(polynomial, interval=(-1, 1), tolerance=F(1, 10**12),
             max_nodes=256, suggestions=None):
    """Return a replayable enclosure of max p(x) on a closed rational interval.

    ``max_nodes`` counts interval refinements, not arithmetic operations or time.
    A valid ``open`` certificate is returned if this budget is insufficient.
    Finite floating suggestions select split points only; they are never proof.
    """
    tolerance = _validate_budget(tolerance, max_nodes)
    data = _prepare(polynomial, interval)
    hints = _suggestions(suggestions, data)
    left, right = data["left"], data["right"]
    if left == right:
        return _assemble(data, [], tolerance, max_nodes)
    count = count_open(data["sequence"], left, right)
    exact = _exact_root(data, left, right, count, hints)
    cells = [_cell(data, left, right, exact)]
    while True:
        cert = _assemble(data, cells, tolerance, max_nodes)
        if cert["status"] != "open" or len(cells)-1 >= max_nodes:
            return cert
        clusters = [i for i, c in enumerate(cells) if c["root_count"] > 1]
        if clusters:
            index = max(clusters, key=lambda i: (cells[i]["root_count"],
                        F(cells[i]["interval"][1])-F(cells[i]["interval"][0])))
        else:
            candidates = [i for i, c in enumerate(cells)
                          if c["root_count"] and c["exact_root"] is None]
            if not candidates:
                raise ArithmeticError("Open bound with no refinable critical cell")
            index = max(candidates, key=lambda i: F(cells[i]["critical_value_bounds"][1]))
        a, b = map(F, cells[index]["interval"])
        # A certified isolated root can often be bracketed quadratically faster
        # by interval Newton. Two added boundaries cost two refinement nodes.
        contraction = (_newton_contraction(data, a, b)
                       if cells[index]["root_count"] == 1 else None)
        if contraction is not None:
            lo, hi = contraction
            if lo == hi and a < lo < b and evaluate(data["sf"], lo) == 0:
                cells[index] = _cell(data, a, b, lo)
                continue
            breaks = sorted({a, lo, hi, b})
            extra_nodes = len(breaks)-2
            if 0 < extra_nodes <= max_nodes-(len(cells)-1):
                children = []
                for lo, hi in zip(breaks, breaks[1:]):
                    count = count_open(data["sequence"], lo, hi)
                    exact = _exact_root(data, lo, hi, count, hints)
                    children.append(_cell(data, lo, hi, exact))
                cells[index:index+1] = children
                continue
        inside_hints = [h for h in hints if a < h < b]
        middle = min(inside_hints, key=lambda h: abs(h-(a+b)/2)) if inside_hints else (a+b)/2
        if middle in hints:
            hints.remove(middle)
        children = []
        for lo, hi in ((a, middle), (middle, b)):
            count = count_open(data["sequence"], lo, hi)
            exact = _exact_root(data, lo, hi, count, hints)
            children.append(_cell(data, lo, hi, exact))
        cells[index:index+1] = children


def verify(cert):
    """Reject malformed, incomplete, altered, or falsely labeled certificates.

    The verifier recomputes algebra, one-sided Sturm counts, critical ranges,
    endpoint/witness values, global bounds and success status. It never calls
    ``maximize`` and needs no floating-point information or external state.
    """
    try:
        if not isinstance(cert, dict) or cert.get("format") != FORMAT:
            return False
        tolerance = _validate_budget(cert["tolerance"], cert["max_nodes"])
        data = _prepare(cert["polynomial"], cert["interval"])
        partition = cert["partition"]
        if not isinstance(partition, list) or len(partition) > cert["max_nodes"]+1:
            return False
        if data["left"] == data["right"]:
            if partition:
                return False
        elif not partition:
            return False
        cells, last = [], data["left"]
        for raw in partition:
            if not isinstance(raw, dict) or not isinstance(raw["interval"], list) or len(raw["interval"]) != 2:
                return False
            left, right = map(_rational, raw["interval"])
            if left != last or not left < right <= data["right"]:
                return False
            exact = None if raw["exact_root"] is None else _rational(raw["exact_root"])
            cells.append(_cell(data, left, right, exact))
            last = right
        if last != data["right"]:
            return False
        expected = _assemble(data, cells, tolerance, cert["max_nodes"])
        # Strict canonical types avoid True==1 and 1.0==1 accepting mutations.
        return _same_canonical(cert, expected)
    except (ArithmeticError, IndexError, KeyError, TypeError, ValueError, OverflowError):
        return False


def _same_canonical(actual, expected):
    if type(actual) is not type(expected):
        return False
    if isinstance(expected, dict):
        return actual.keys() == expected.keys() and all(
            _same_canonical(actual[k], expected[k]) for k in expected)
    if isinstance(expected, list):
        return len(actual) == len(expected) and all(
            _same_canonical(a, b) for a, b in zip(actual, expected))
    return actual == expected
