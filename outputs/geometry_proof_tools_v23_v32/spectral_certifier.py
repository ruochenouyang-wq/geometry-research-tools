#!/usr/bin/env python3
"""Exact two-sided spectral bounds, including the infinite omitted tail.

Model: -Delta_{S^2} + q(cos(theta)), RESTRICTED to axisymmetric functions.
q is a rational polynomial. k is one-based and includes the constant mode.
Uses exact Legendre assembly, rational inertia, and a Schur tail bound.
Python 3.9+, standard library only. No numerical eigensolver is trusted.
"""
import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import re
import sys

MODEL = "unit_sphere_axisymmetric_polynomial_potential_v1"
SCOPE = "axisymmetric functions on the unit two-sphere; eigenvalues indexed from 1 including the constant mode at q=0"
MAX_MODES = 32


def rational(raw):
    if type(raw) is int:
        raw = str(raw)
    if not isinstance(raw, str) or len(raw) > 2000 or not re.fullmatch(r"-?\d+(?:/[1-9]\d*)?", raw):
        raise ValueError("需要整数或分数字符串；不接受浮点数。")
    return F(raw)


def potential(raw):
    if not isinstance(raw, list) or not 1 <= len(raw) <= 7:
        raise ValueError("势函数 q 的系数按升幂排列，次数最多为 6。")
    q = [rational(c) for c in raw]
    if any(abs(c) > 1000 or c.denominator > 10**9 for c in q):
        raise ValueError("势函数系数绝对值最多 1000，分母最多 10^9。")
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q


def check_sizes(k, modes, bits=None):
    if type(k) is not int or type(modes) is not int or not 1 <= k <= modes <= MAX_MODES:
        raise ValueError("必须满足 1 <= k <= modes <= 32；k 从 1 开始，包括零模态。")
    if bits is not None and (type(bits) is not int or not 8 <= bits <= 64):
        raise ValueError("bits 必须在 8..64 之间。")


def x_times(series):
    """x P_n = (n+1)/(2n+1) P_{n+1} + n/(2n+1) P_{n-1}."""
    out = {}
    for n, c in series.items():
        out[n + 1] = out.get(n + 1, F(0)) + c * F(n + 1, 2 * n + 1)
        if n:
            out[n - 1] = out.get(n - 1, F(0)) + c * F(n, 2 * n + 1)
    return {n: c for n, c in out.items() if c}


def q_times_p(q, n):
    result, term = {}, {n: F(1)}
    for c in q:
        for i, v in term.items():
            result[i] = result.get(i, F(0)) + c * v
        term = x_times(term)
    return {i: c for i, c in result.items() if c}


def assemble(q, modes):
    """A: finite form, M: diagonal mass, G=C M_tail^-1 C^T (exact)."""
    degree = len(q) - 1
    mass = [F(2, 2 * i + 1) for i in range(modes + degree)]
    a = [[F(0)] * modes for _ in range(modes)]
    c = [[F(0)] * degree for _ in range(modes)]
    for j in range(modes + degree):
        for i, coefficient in q_times_p(q, j).items():
            if i < modes:
                if j < modes:
                    a[i][j] += mass[i] * coefficient
                else:
                    c[i][j - modes] = mass[i] * coefficient
    for i in range(modes):
        a[i][i] += i * (i + 1) * mass[i]
    if any(a[i][j] != a[j][i] for i in range(modes) for j in range(modes)):
        raise ArithmeticError("内部错误：算子矩阵不对称。")
    g = [[sum((c[i][j] * c[l][j] / mass[modes + j] for j in range(degree)), F(0))
          for l in range(modes)] for i in range(modes)]
    radius = sum((abs(x) for x in q[1:]), F(0))
    qmin, qmax = q[0] - radius, q[0] + radius
    return {"A": a, "M": mass[:modes], "G": g, "qmin": qmin, "qmax": qmax,
            "beta": modes * (modes + 1) + qmin}


def inertia(matrix):
    """Exact symmetric 1x1/2x2 elimination; returns negative, zero, positive."""
    a = [[F(x) for x in row] for row in matrix]
    size = len(a)
    if any(len(row) != size for row in a) or any(a[i][j] != a[j][i] for i in range(size) for j in range(size)):
        raise ValueError("惯性检查要求精确对称方阵。")
    negative, zero, positive = 0, 0, 0
    while a:
        n = len(a)
        pivot = next((i for i in range(n) if a[i][i]), None)
        if pivot is not None:
            order = [pivot] + [i for i in range(n) if i != pivot]
            a = [[a[i][j] for j in order] for i in order]
            d = a[0][0]
            negative += d < 0
            positive += d > 0
            a = [[a[i][j] - a[i][0] * a[0][j] / d for j in range(1, n)] for i in range(1, n)]
            continue
        pair = next(((i, j) for i in range(n) for j in range(i + 1, n) if a[i][j]), None)
        if pair is None:
            zero += n
            break
        # All remaining diagonal entries are zero. A nonzero off-diagonal
        # creates the invertible block [[0,b],[b,0]], of inertia (1,0,1).
        order = list(pair) + [i for i in range(n) if i not in pair]
        a = [[a[i][j] for j in order] for i in order]
        b = a[0][1]
        negative += 1
        positive += 1
        a = [[a[i][j] - (a[i][0] * a[1][j] + a[i][1] * a[0][j]) / b
              for j in range(2, n)] for i in range(2, n)]
    return [negative, zero, positive]


def shifted(data, value, tail=False):
    n = len(data["A"])
    if tail and value >= data["beta"]:
        raise ValueError("下界候选必须严格小于尾部下界 beta。")
    return [[data["A"][i][j] - (value * data["M"][i] if i == j else 0)
             - (data["G"][i][j] / (data["beta"] - value) if tail else 0)
             for j in range(n)] for i in range(n)]


def upper_holds(data, k, endpoint):
    count = inertia(shifted(data, endpoint))
    return count[0] + count[1] >= k


def lower_holds(data, k, endpoint):
    return endpoint < data["beta"] and inertia(shifted(data, endpoint, tail=True))[0] < k


def outward_decimal(number, places=12, upper=False):
    factor = 10**places
    scaled = number * factor
    integer = -((-scaled.numerator) // scaled.denominator) if upper else scaled.numerator // scaled.denominator
    sign = "-" if integer < 0 else ""
    whole, fraction = divmod(abs(integer), factor)
    return f"{sign}{whole}.{fraction:0{places}d}"


def certify(q, k=1, modes=8, bits=36):
    check_sizes(k, modes, bits)
    q = potential([str(c) for c in q])
    data = assemble(q, modes)
    # The first k unperturbed modes give an analytic Ritz upper bound.
    low, high = data["qmin"] - 1, data["qmax"] + (k - 1) * k + 1
    for _ in range(bits):
        mid = (low + high) / 2
        if upper_holds(data, k, mid):
            high = mid
        else:
            low = mid
    upper = high
    low, high = data["qmin"] - 1, upper
    # qmin is a true operator lower bound, but the more conservative Schur
    # sufficient condition need not hold at qmin-1 for strong coupling.
    # Establish the invariant explicitly before bisecting. It eventually
    # holds since -low*M dominates and G/(beta-low) tends to zero.
    distance = F(1)
    while not lower_holds(data, k, low):
        distance *= 2
        low = data["qmin"] - distance
    for _ in range(bits):
        mid = (low + high) / 2
        if lower_holds(data, k, mid):
            low = mid
        else:
            high = mid
    lower = low
    evidence = {"analytic_potential_lower": str(data["qmin"]), "tail_lower": str(data["beta"]),
                "lower_schur_inertia": inertia(shifted(data, lower, tail=True)),
                "upper_finite_inertia": inertia(shifted(data, upper))}
    result = {"model": MODEL, "q_coefficients": [str(c) for c in q], "eigenvalue_index": k,
              "modes": modes, "bisection_steps": bits, "lower": str(lower), "upper": str(upper),
              "evidence": evidence,
              "decimal_enclosure": [outward_decimal(lower), outward_decimal(upper, upper=True)],
              "exact_width": str(upper - lower),
              "scope": SCOPE,
              "verification": "exact_rational_with_analytic_tail", "formal_assistant_checked": False}
    if not verify(result):
        raise ArithmeticError("内部错误：生成的谱证书没有通过复核。")
    return result


def verify(certificate):
    if not isinstance(certificate, dict) or certificate.get("model") != MODEL:
        return False
    if certificate.get("scope") != SCOPE or certificate.get("verification") != "exact_rational_with_analytic_tail" or certificate.get("formal_assistant_checked") is not False:
        return False
    try:
        q = potential(certificate["q_coefficients"])
        k, modes = certificate["eigenvalue_index"], certificate["modes"]
        check_sizes(k, modes)
        lo, hi = rational(certificate["lower"]), rational(certificate["upper"])
        if lo > hi:
            return False
        data = assemble(q, modes)
        if not lower_holds(data, k, lo) or not upper_holds(data, k, hi):
            return False
        expected = {"analytic_potential_lower": str(data["qmin"]), "tail_lower": str(data["beta"]),
                    "lower_schur_inertia": inertia(shifted(data, lo, tail=True)),
                    "upper_finite_inertia": inertia(shifted(data, hi))}
        if certificate.get("evidence") != expected:
            return False
        if certificate.get("decimal_enclosure") != [outward_decimal(lo), outward_decimal(hi, upper=True)]:
            return False
        if certificate.get("exact_width") != str(hi - lo):
            return False
        return True
    except (ValueError, KeyError, TypeError, ZeroDivisionError):
        return False


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description="几何谱工具：有限维矩阵 + 无限维尾部的精确上下界认证。")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("bound")
    run.add_argument("--q", default="0,1", help="q(t) 的有理系数，按常数、t、t²...排列")
    run.add_argument("--k", type=int, default=1)
    run.add_argument("--modes", type=int, default=8)
    run.add_argument("--bits", type=int, default=36)
    run.add_argument("--out", required=True)
    check = sub.add_parser("verify")
    check.add_argument("certificate")
    demo = sub.add_parser("demo")
    demo.add_argument("--out-dir", default=str(Path(__file__).resolve().parent / "results"))
    args = parser.parse_args(argv)
    try:
        if args.command == "verify":
            path = Path(args.certificate)
            if path.stat().st_size > 1_000_000:
                raise ValueError("证书文件超过大小限制。")
            passed = verify(json.loads(path.read_text(encoding="utf-8")))
            print("谱上下界证书复核通过。" if passed else "谱上下界证书复核未通过。")
            return 0 if passed else 1
        if args.command == "demo":
            cases = [("sphere_first_positive", ["0"], 2, 6),
                     ("linear_potential_coarse", ["0", "1"], 1, 2),
                     ("linear_potential_fine", ["0", "1"], 1, 8),
                     ("quadratic_potential", ["1", "-2", "3"], 1, 10),
                     ("second_eigenvalue", ["0", "1"], 2, 8)]
            results = []
            for name, q, k, modes in cases:
                cert = certify(q, k, modes)
                write_json(Path(args.out_dir) / (name + ".json"), cert)
                results.append({"case": name, "k": k, "modes": modes, "interval": cert["decimal_enclosure"], "verified": verify(cert)})
                print(name + ": [" + ", ".join(cert["decimal_enclosure"]) + "]")
            write_json(Path(args.out_dir) / "summary.json", results)
            return 0
        result = certify(potential(args.q.split(",")), args.k, args.modes, args.bits)
        write_json(args.out, result)
        print("已认证区间：[" + ", ".join(result["decimal_enclosure"]) + "]；证书：" + args.out)
        return 0
    except (ValueError, OSError, TypeError) as exc:
        print("错误：" + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
