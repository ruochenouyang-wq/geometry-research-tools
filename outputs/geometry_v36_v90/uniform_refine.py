"""Exact, reusable scalar spectral bounds from parameter-dependent positive functions.

Only the standard library is used.  A successful certificate concerns every
amplitude in its stated domain and the full unit sphere.  ``verify`` checks the
record including an honest unresolved outcome; use ``is_proved`` for acceptance.
"""
from fractions import Fraction as F
from math import comb


METHOD = "parameter_poisson_tensor_bernstein_v1"
LOG_SOBOLEV_METHOD = "unit_sphere_affine_log_sobolev_v1"
MAX_ORDER = 12
MAX_LEAVES = 2048


def _fraction(value):
    if isinstance(value, bool) or not isinstance(value, (str, int, F)):
        raise ValueError("Use exact rational strings, integers or Fraction")
    result = F(value)
    if result.numerator.bit_length() > 4096 or result.denominator.bit_length() > 4096:
        raise ValueError("Rational input budget exceeded")
    return result


def _trim(poly):
    result = list(poly)
    while len(result) > 1 and not result[-1]:
        result.pop()
    return result or [F(0)]


def _add(left, right):
    result = [F(0)] * max(len(left), len(right))
    for i, value in enumerate(left):
        result[i] += value
    for i, value in enumerate(right):
        result[i] += value
    return _trim(result)


def _scale(poly, scalar):
    return _trim([value * scalar for value in poly])


def _mul(left, right):
    result = [F(0)] * (len(left) + len(right) - 1)
    for i, value in enumerate(left):
        for j, other in enumerate(right):
            result[i + j] += value * other
    return _trim(result)


def _derivative(poly):
    return _trim([i * poly[i] for i in range(1, len(poly))])


def _mean(poly):
    return sum((value / (i + 1) for i, value in enumerate(poly) if i % 2 == 0), F(0))


def _eval(poly, point):
    result = F(0)
    for value in reversed(poly):
        result = result * point + value
    return result


def laplacian(poly):
    """H0 = -( (1-t^2) d/dt ) d/dt, in power coefficients."""
    return _add(_mul([0, 2], _derivative(poly)),
                _scale(_mul([1, 0, -1], _derivative(_derivative(poly))), -1))


def poisson_zero_mean(target):
    """Solve H0 r = target exactly; choose r(0)=0 to fix the constant."""
    target = _trim(target)
    if _mean(target):
        raise ValueError("Poisson target must have zero spherical mean")
    remainder = target[:]
    result = [F(0)] * len(target)
    for j in range(len(target) - 1, 0, -1):
        value = remainder[j] / (j * (j + 1))
        result[j] = value
        remainder[j] = F(0)
        if j >= 2:
            remainder[j - 2] += j * (j - 1) * value
    if any(remainder) or laplacian(result) != target:
        raise ArithmeticError("Poisson identity failed")
    return _trim(result)


def parameter_correctors(direction=(0, 1), order=3):
    """Return r[0..order] and formal local-energy constants e[0..order].

    For s=sum(a**j*r[j](t),j=1..order), all local-energy coefficients
    through the requested order are constants in t.  Higher terms remain
    explicit and must still be certified, not discarded as a Taylor error.
    """
    if type(order) is not int or not 1 <= order <= MAX_ORDER:
        raise ValueError("Order must be an integer from 1 to 12")
    direction = _trim([_fraction(value) for value in direction])
    if len(direction) > 7:
        raise ValueError("Direction degree at most six")
    mean = _mean(direction)
    result = [[F(0)], poisson_zero_mean(_add(_scale(direction, -1), [mean]))]
    energies = [F(0), mean]
    for n in range(2, order + 1):
        square = [F(0)]
        for i in range(1, n):
            square = _add(square, _mul(_derivative(result[i]),
                                      _derivative(result[n - i])))
        forcing = _mul([1, 0, -1], square)
        average = _mean(forcing)
        result.append(poisson_zero_mean(_add(forcing, [-average])))
        energies.append(-average)
    return result, energies


def _sparse_add(target, a_degree, polynomial, scale=F(1)):
    for t_degree, value in enumerate(polynomial):
        key = (a_degree, t_degree)
        target[key] = target.get(key, F(0)) + scale * value
        if not target[key]:
            del target[key]


def local_residual(direction, coefficient, correctors):
    """Exact polynomial L = a(p-mean p)+c*a^2+H0*s-(1-t^2)*s_t^2."""
    result = {}
    _sparse_add(result, 1, _add(direction, [-_mean(direction)]))
    _sparse_add(result, 2, [coefficient])
    for n in range(1, len(correctors)):
        _sparse_add(result, n, laplacian(correctors[n]))
        for m in range(1, len(correctors)):
            product = _mul([1, 0, -1], _mul(_derivative(correctors[n]),
                                          _derivative(correctors[m])))
            _sparse_add(result, n + m, product, F(-1))
    return result


def extract_even_factor(polynomial):
    """Extract the largest known even amplitude power, never divide numerically."""
    if not polynomial:
        return 0, {}
    factor = min(key[0] for key in polynomial)
    factor -= factor % 2
    return factor, {(i - factor, j): value for (i, j), value in polynomial.items()}


def _encode_sparse(polynomial):
    return [[i, j, str(value)] for (i, j), value in sorted(polynomial.items())]


def _evaluate_sparse(polynomial, amplitude, t):
    return sum((value * amplitude ** i * t ** j
                for (i, j), value in polynomial.items()), F(0))


def tensor_bernstein(polynomial, amplitude_box, t_box=(F(-1), F(1))):
    """Exact tensor Bernstein coefficients after affine mapping to [0,1]^2."""
    n = max((i for i, _ in polynomial), default=0)
    m = max((j for _, j in polynomial), default=0)
    alo, ahi = map(F, amplitude_box)
    tlo, thi = map(F, t_box)
    if not alo < ahi or not tlo < thi:
        raise ValueError("Nonempty intervals required")
    affine = [[F(0)] * (m + 1) for _ in range(n + 1)]
    for (i, j), value in polynomial.items():
        for k in range(i + 1):
            av = value * comb(i, k) * alo ** (i - k) * (ahi - alo) ** k
            for ell in range(j + 1):
                affine[k][ell] += av * comb(j, ell) * tlo ** (j - ell) * (thi - tlo) ** ell
    first = [[sum((affine[i][j] * F(comb(k, i), comb(n, i))
                   for i in range(k + 1)), F(0))
              for j in range(m + 1)] for k in range(n + 1)]
    return [[sum((first[k][j] * F(comb(ell, j), comb(m, j))
                  for j in range(ell + 1)), F(0))
             for ell in range(m + 1)] for k in range(n + 1)]


def replay_tensor_bernstein(polynomial, amplitude_box, t_box=(F(-1), F(1))):
    """Independent direct formula, used by the verifier on every leaf box.

    This does not call the synthesizer's affine transform, power-to-Bernstein
    conversion or de Casteljau subdivision. It starts with the exact residual
    polynomial and the independently reconstructed rational leaf intervals.
    """
    degrees = (max((i for i, _ in polynomial), default=0),
               max((j for _, j in polynomial), default=0))
    tables = []
    for degree, endpoints in zip(degrees, (amplitude_box, t_box)):
        lo, hi = endpoints
        tables.append({(power, index): sum((F(comb(power, j) * comb(index, j), comb(degree, j))
                                            * lo ** (power - j) * (hi - lo) ** j
                                            for j in range(min(power, index) + 1)), F(0))
                       for power in range(degree + 1) for index in range(degree + 1)})
    return [[sum((value * tables[0][(i, k)] * tables[1][(j, ell)]
                  for (i, j), value in polynomial.items()), F(0))
             for ell in range(degrees[1] + 1)] for k in range(degrees[0] + 1)]


def _split_line(coefficients):
    level = list(coefficients)
    left, right = [level[0]], [level[-1]]
    while len(level) > 1:
        level = [(level[i] + level[i + 1]) / 2 for i in range(len(level) - 1)]
        left.append(level[0])
        right.append(level[-1])
    return left, list(reversed(right))


def split_tensor(coefficients, axis):
    """de Casteljau bisection; axis 0=amplitude, 1=t."""
    if axis not in (0, 1):
        raise ValueError("Unknown tensor axis")
    if axis == 1:
        pairs = [_split_line(row) for row in coefficients]
        return [pair[0] for pair in pairs], [pair[1] for pair in pairs]
    transposed = [list(column) for column in zip(*coefficients)]
    left, right = split_tensor(transposed, 1)
    return [list(column) for column in zip(*left)], [list(column) for column in zip(*right)]


def _bounds(coefficients):
    flat = [value for row in coefficients for value in row]
    return min(flat), max(flat)


def _leaf_record(coefficients):
    lo, hi = _bounds(coefficients)
    return {"lower": str(lo), "upper": str(hi)}


def _split_box(box, axis):
    midpoint = (box[axis][0] + box[axis][1]) / 2
    left, right = list(box), list(box)
    left[axis] = (box[axis][0], midpoint)
    right[axis] = (midpoint, box[axis][1])
    return tuple(left), tuple(right)


def _witness(polynomial, box, factor=0):
    amplitudes = (box[0][0], sum(box[0]) / 2, box[0][1])
    points = (box[1][0], sum(box[1]) / 2, box[1][1])
    value, _, amplitude, t = min((_evaluate_sparse(polynomial, a, t), -abs(a), a, t)
                                 for a in amplitudes for t in points)
    if value >= 0:
        return None
    if factor and not amplitude:
        endpoint = next(a for a in amplitudes if a)
        amplitude = endpoint
        # A strict negative value at zero has a nonzero negative neighbor.
        # Continuity motivates the search; an exact negative value accepts it.
        for _ in range(128):
            value = _evaluate_sparse(polynomial, amplitude, t)
            if value < 0:
                break
            amplitude /= 2
        else:
            return None
    return {"amplitude": str(amplitude), "t": str(t), "reduced_value": str(value),
            "residual_value": str(value * amplitude ** factor),
            "meaning": "This local-energy ansatz fails; this is not a spectral counterexample."}


def _setup(direction, coefficient, amplitude_box, order, symmetry, all_real):
    direction = _trim([_fraction(value) for value in direction])
    coefficient = _fraction(coefficient)
    if coefficient <= 0:
        raise ValueError("Positive quadratic coefficient required")
    if type(symmetry) is not bool or type(all_real) is not bool:
        raise ValueError("Boolean flags required")
    correctors, energies = parameter_correctors(direction, order)
    centered = _add(direction, [-_mean(direction)])
    envelope = sum(map(abs, centered), F(0))
    if all_real:
        # A constant direction is handled by any nonempty inner interval.
        radius = max(envelope / coefficient, F(1))
        domain = (-radius, radius)
    else:
        if not isinstance(amplitude_box, (tuple, list)) or len(amplitude_box) != 2:
            raise ValueError("Two exact amplitude endpoints required")
        domain = tuple(map(_fraction, amplitude_box))
        if not domain[0] < domain[1]:
            raise ValueError("Amplitude interval must be nonempty")
    residual = local_residual(direction, coefficient, correctors)
    factor, reduced = extract_even_factor(residual)
    # Centered odd directions satisfy H(a) unitarily equivalent to H(-a)
    # after subtracting a*mean. Check both the input and the actual residual.
    reflected = (symmetry and not any(centered[::2]) and domain[0] == -domain[1]
                 and all((i + j) % 2 == 0 for i, j in reduced))
    proof_domain = (F(0), domain[1]) if reflected else domain
    return {"direction": direction, "coefficient": coefficient, "correctors": correctors,
            "energies": energies, "domain": domain, "proof_domain": proof_domain,
            "reflected": reflected, "factor": factor, "reduced": reduced,
            "residual": residual, "envelope": envelope, "all_real": all_real}


def _base_record(data, order, symmetry):
    tail = None
    if data["all_real"]:
        tail = {"centered_potential_abs_upper": str(data["envelope"]),
                "radius": str(data["domain"][1]),
                "rule": "lambda_1-a*mean(p) >= -B*abs(a); c*abs(a)>=B outside radius"}
    return {"method": METHOD, "direction": list(map(str, data["direction"])),
            "coefficient": str(data["coefficient"]), "order": order,
            "requested_symmetry": symmetry, "all_real": data["all_real"],
            "amplitude_box": list(map(str, data["domain"])),
            "proof_amplitude_box": list(map(str, data["proof_domain"])),
            "reflection_used": data["reflected"],
            "correctors": [list(map(str, p)) for p in data["correctors"]],
            "energy_series": list(map(str, data["energies"])),
            "extracted_amplitude_power": data["factor"],
            "residual": _encode_sparse(data["residual"]),
            "reduced_residual": _encode_sparse(data["reduced"]),
            "tail": tail,
            "claim": "lambda_1(-Delta+a*p(t)) >= a*mean(p)-coefficient*a^2",
            "scope": "full unit sphere; all real amplitudes" if data["all_real"] else
                     "full unit sphere; every amplitude in the stated closed box",
            "formal_assistant_checked": False}


def synthesize(direction=(0, 1), coefficient=F(1, 5), amplitude_box=(-2, 2),
               order=3, max_leaves=256, symmetry=True, all_real=False):
    """Propose and exactly certify a lower bound, or retain a verified failure.

    ``all_real=True`` selects its inner radius from the cheap potential tail
    estimate; a successful finite box alone never becomes an all-real claim.
    The parameter order controls a finite exact polynomial, not an asymptotic
    approximation. Adaptive selection affects speed, never acceptance.
    """
    if type(max_leaves) is not int or not 1 <= max_leaves <= MAX_LEAVES:
        raise ValueError("Leaf budget must be between 1 and 2048")
    data = _setup(direction, coefficient, amplitude_box, order, symmetry, all_real)
    root_coefficients = tensor_bernstein(data["reduced"], data["proof_domain"])
    tree = _leaf_record(root_coefficients)
    leaves = [(tree, root_coefficients, (data["proof_domain"], (F(-1), F(1))), 0)]
    failure = None
    while True:
        failing = [i for i, (_, coeff, _, _) in enumerate(leaves) if _bounds(coeff)[0] < 0]
        if not failing:
            status = "proved"
            break
        worst = min(failing, key=lambda i: _bounds(leaves[i][1])[0])
        node, coefficients, box, depth = leaves[worst]
        failure = _witness(data["reduced"], box, data["factor"])
        if failure is not None:
            status = "ansatz_obstruction"
            break
        if len(leaves) >= max_leaves or depth >= 128:
            status = "unresolved"
            break
        choices = []
        for axis in (0, 1):
            left, right = split_tensor(coefficients, axis)
            choices.append((min(_bounds(left)[0], _bounds(right)[0]), axis, left, right))
        # Prefer the split giving the better certified minimum. Tie-breaking
        # alternates with geometric widths to avoid a permanently flat axis.
        _, axis, left, right = max(choices, key=lambda item:
                                  (item[0], box[item[1]][1] - box[item[1]][0]))
        left_box, right_box = _split_box(box, axis)
        left_node, right_node = _leaf_record(left), _leaf_record(right)
        node.clear()
        node.update({"axis": axis, "left": left_node, "right": right_node})
        leaves.pop(worst)
        leaves.extend([(left_node, left, left_box, depth + 1),
                       (right_node, right, right_box, depth + 1)])
    result = _base_record(data, order, symmetry)
    result.update({"tree": tree, "leaf_count": len(leaves), "max_leaves": max_leaves,
                   "status": status, "failure_witness": failure,
                   "certified_reduced_lower": str(min(_bounds(c)[0] for _, c, _, _ in leaves)),
                   "certified_reduced_upper": str(max(_bounds(c)[1] for _, c, _, _ in leaves))})
    if not verify(result):
        raise ArithmeticError("Generated record failed exact verification")
    return result


def verify(certificate):
    """Reconstruct all identities and the entire cover, including failure status."""
    try:
        if not isinstance(certificate, dict):
            return False
        if certificate.get("method") == LOG_SOBOLEV_METHOD:
            return verify_log_sobolev(certificate)
        order = certificate["order"]
        symmetry = certificate["requested_symmetry"]
        max_leaves = certificate["max_leaves"]
        if type(max_leaves) is not int or not 1 <= max_leaves <= MAX_LEAVES:
            return False
        data = _setup(certificate["direction"], certificate["coefficient"],
                      certificate["amplitude_box"], order, symmetry, certificate["all_real"])
        expected = _base_record(data, order, symmetry)
        stack = [(certificate["tree"], (data["proof_domain"], (F(-1), F(1))), 0)]
        count, bounds, visited = 0, [], 0
        while stack:
            node, box, depth = stack.pop()
            visited += 1
            if visited > 2 * MAX_LEAVES - 1 or depth > 128 or not isinstance(node, dict):
                return False
            if "axis" in node:
                if set(node) != {"axis", "left", "right"} or type(node["axis"]) is not int:
                    return False
                axis = node["axis"]
                if axis not in (0, 1):
                    return False
                left_box, right_box = _split_box(box, axis)
                stack.extend([(node["left"], left_box, depth + 1),
                              (node["right"], right_box, depth + 1)])
            else:
                coeff = replay_tensor_bernstein(data["reduced"], box[0], box[1])
                if node != _leaf_record(coeff):
                    return False
                count += 1
                bounds.append(_bounds(coeff))
        if not 1 <= count <= max_leaves:
            return False
        lower = min(pair[0] for pair in bounds)
        upper = max(pair[1] for pair in bounds)
        failure = certificate["failure_witness"]
        if lower >= 0:
            status = "proved"
            if failure is not None:
                return False
        elif failure is not None:
            amplitude, t = _fraction(failure["amplitude"]), _fraction(failure["t"])
            if not data["proof_domain"][0] <= amplitude <= data["proof_domain"][1] or not -1 <= t <= 1:
                return False
            value = _evaluate_sparse(data["reduced"], amplitude, t)
            residual_value = value * amplitude ** data["factor"]
            if value >= 0 or residual_value >= 0:
                return False
            exact = {"amplitude": str(amplitude), "t": str(t), "reduced_value": str(value),
                     "residual_value": str(residual_value),
                     "meaning": "This local-energy ansatz fails; this is not a spectral counterexample."}
            if failure != exact:
                return False
            status = "ansatz_obstruction"
        else:
            status = "unresolved"
        expected.update({"tree": certificate["tree"], "leaf_count": count,
                         "max_leaves": max_leaves, "status": status,
                         "failure_witness": failure,
                         "certified_reduced_lower": str(lower),
                         "certified_reduced_upper": str(upper)})
        return certificate == expected
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError, OverflowError,
            ArithmeticError, RecursionError):
        return False


def is_proved(certificate):
    """The sole acceptance predicate for the spectral inequality."""
    return verify(certificate) and certificate["status"] == "proved"


def _log_sobolev_record(direction, coefficient):
    direction = _trim([_fraction(value) for value in direction])
    if len(direction) > 2:
        raise ValueError("The log-Sobolev moment rule is restricted to affine p(t)=b+d*t")
    constant = direction[0]
    slope = direction[1] if len(direction) == 2 else F(0)
    sharp = slope * slope / 6
    coefficient = sharp if coefficient is None else _fraction(coefficient)
    # Universal induction, not a finite numerical check of the moment series:
    # (2k+3)(2k+2)-6(k+1) = 4k+4k^2 >= 0 for every integer k>=0.
    induction = _add(_mul([F(3), F(2)], [F(2), F(2)]), [F(-6), F(-6)])
    if any(value < 0 for value in induction):
        raise ArithmeticError("The universal moment induction failed")
    witness = None
    if coefficient < sharp:
        amplitude = F(1)
        for _ in range(8192):
            mass = 1 + amplitude ** 2 * slope ** 2 / 12
            centered_numerator = -amplitude ** 2 * slope ** 2 / 6
            objective_upper = centered_numerator / mass + coefficient * amplitude ** 2
            if objective_upper < 0:
                break
            amplitude /= 2
        else:
            raise ArithmeticError("Unable to materialize a strict rational counterexample")
        witness = {"amplitude": str(amplitude),
                   "trial_power_coefficients": ["1", str(-amplitude * slope / 2)],
                   "mass": str(mass), "centered_energy_numerator": str(centered_numerator),
                   "spectral_objective_upper": str(objective_upper),
                   "meaning": "An exact Rayleigh trial disproves the requested all-real spectral bound."}
    return {"method": LOG_SOBOLEV_METHOD, "direction": list(map(str, direction)),
            "coefficient": str(coefficient), "mean": str(constant),
            "sharp_uniform_quadratic_coefficient": str(sharp),
            "status": "proved" if coefficient >= sharp else "disproved",
            "scope": "full unit sphere; all real amplitudes; affine coordinate potentials only",
            "claim": "lambda_1(-Delta+a*p(t)) >= a*mean(p)-coefficient*a^2",
            "trusted_analytic_theorem": {
                "name": "Sharp logarithmic Sobolev inequality on the unit sphere",
                "source": "https://arxiv.org/pdf/1210.1853v1",
                "location": "Corollary 2, page 3; uniform probability measure; dimension d=2",
                "entropy_to_dirichlet_constant": "1",
                "statement": "Ent_mu(|u|^2) <= integral_S2 |grad u|^2 dmu",
                "proved_by_this_program": False},
            "moment_series_rule": {
                "identity": "integral exp(z*t) dmu = sinh(z)/z, continued as 1 at z=0",
                "domination": "(2*k+1)! >= 6^k*k! for every integer k>=0",
                "induction_difference_coefficients": list(map(str, induction)),
                "consequence": "sinh(z)/z <= exp(z^2/6) for every real z"},
            "variational_rule": "Dirichlet >= entropy; Gibbs inequality gives energy >= -log moment",
            "optimality_proof": {
                "trial": "u(a,t)=1-a*d*t/2, where d is the affine slope",
                "mass_in_amplitude": ["1", "0", str(slope ** 2 / 12)],
                "centered_energy_numerator_in_amplitude": ["0", "0", str(-slope ** 2 / 6)],
                "limit_rule": "Divide any valid quadratic bound by a^2 and let nonzero a tend to zero",
                "necessary_coefficient": str(sharp)},
            "counterexample": witness, "formal_assistant_checked": False}


def synthesize_log_sobolev(direction=(0, 1), coefficient=None):
    """An explicitly separate analytic rule, dependent on a classical theorem.

    The default coefficient is sharp: d^2/6 for p(t)=b+d*t. This rule neither
    calls nor claims to follow from the tensor Bernstein verifier. A smaller
    requested coefficient produces an exact spectral counterexample.
    """
    result = _log_sobolev_record(direction, coefficient)
    if not verify_log_sobolev(result):
        raise ArithmeticError("Log-Sobolev rule failed verification")
    return result


def verify_log_sobolev(certificate):
    """Check theorem binding, universal induction and exact Rayleigh arithmetic.

    Trust in the cited log-Sobolev theorem is explicit. This is not a proof
    assistant checking that theorem or its analytic consequences from axioms.
    """
    try:
        expected = _log_sobolev_record(certificate["direction"], certificate["coefficient"])
        if expected != certificate:
            return False
        witness = certificate["counterexample"]
        if witness is not None:
            amplitude = F(witness["amplitude"])
            trial = list(map(F, witness["trial_power_coefficients"]))
            p = list(map(F, certificate["direction"]))
            centered = _add(p, [-_mean(p)])
            trial_squared = _mul(trial, trial)
            mass = _mean(trial_squared)
            kinetic = _mean(_mul([1, 0, -1], _mul(_derivative(trial), _derivative(trial))))
            potential = amplitude * _mean(_mul(centered, trial_squared))
            upper = (kinetic + potential) / mass + F(certificate["coefficient"]) * amplitude ** 2
            if (str(mass) != witness["mass"] or str(kinetic + potential) != witness["centered_energy_numerator"]
                    or str(upper) != witness["spectral_objective_upper"] or upper >= 0):
                return False
        return True
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError, ArithmeticError):
        return False


def milestone_records():
    """Deterministic small evidence set, retaining both failures and advances."""
    requests = [
        ("linear_1_over_4", dict(coefficient=F(1, 4), order=1)),
        ("linear_1_over_5_obstruction", dict(coefficient=F(1, 5), order=1)),
        ("quadratic_1_over_5_obstruction", dict(coefficient=F(1, 5), order=2)),
        ("cubic_box_1_over_5", dict(coefficient=F(1, 5), order=3)),
        ("cubic_all_real_1_over_5", dict(coefficient=F(1, 5), order=3,
                                        all_real=True, max_leaves=256)),
        ("quartic_all_real_1_over_6_obstruction", dict(coefficient=F(1, 6), order=4,
                                                     all_real=True)),
    ]
    result = [{"name": name, "certificate": synthesize(**options)} for name, options in requests]
    result.extend([
        {"name": "log_sobolev_all_real_sharp_1_over_6", "certificate": synthesize_log_sobolev()},
        {"name": "below_sharp_1_over_7_spectral_counterexample",
         "certificate": synthesize_log_sobolev(coefficient=F(1, 7))}])
    return result


if __name__ == "__main__":
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coefficient", default="1/5")
    parser.add_argument("--radius", default="2")
    parser.add_argument("--order", type=int, default=3)
    parser.add_argument("--max-leaves", type=int, default=256)
    parser.add_argument("--all-real", action="store_true")
    parser.add_argument("--log-sobolev", action="store_true")
    parser.add_argument("--milestones", action="store_true")
    args = parser.parse_args()
    result = milestone_records() if args.milestones else (
        synthesize_log_sobolev(coefficient=args.coefficient) if args.log_sobolev else synthesize(
            coefficient=args.coefficient, amplitude_box=(-F(args.radius), F(args.radius)),
            order=args.order, max_leaves=args.max_leaves, all_real=args.all_real))
    print(json.dumps(result, ensure_ascii=False, sort_keys=True, indent=2))
