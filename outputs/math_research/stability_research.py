#!/usr/bin/env python3
"""Representation search + numerical SOS synthesis + exact certificate checking.

Only synthesis requires NumPy. Verification uses Python's exact fractions.
No generated code is executed. A failed search never means instability.
"""
import argparse
from fractions import Fraction as F
from itertools import product, combinations
import json
import os
import sys
import time

from math_research import exact, GPTProposer, ProviderError, load_json, save_json


def clean(p):
    return {e: F(c) for e, c in p.items() if c}


def plus(*polys):
    out = {}
    for p in polys:
        for e, c in p.items():
            out[e] = out.get(e, F(0)) + c
    return clean(out)


def scale(p, c):
    return clean({e: v * c for e, v in p.items()})


def times(p, q):
    out = {}
    for e, c in p.items():
        for f, d in q.items():
            g = tuple(a + b for a, b in zip(e, f))
            out[g] = out.get(g, F(0)) + c * d
    return clean(out)


def deriv(p, i):
    out = {}
    for e, c in p.items():
        if e[i]:
            f = list(e)
            f[i] -= 1
            out[tuple(f)] = c * e[i]
    return clean(out)


def lie(p, field):
    return plus(*(times(deriv(p, i), f) for i, f in enumerate(field)))


def value(p, point):
    return sum((c * product_value([x ** e for x, e in zip(point, powers)]) for powers, c in p.items()), F(0))


def product_value(values):
    out = F(1)
    for x in values:
        out *= x
    return out


def monomials(n, degree, minimum=None):
    minimum = degree if minimum is None else minimum
    return sorted((e for e in product(range(degree + 1), repeat=n) if minimum <= sum(e) <= degree), key=lambda e: (sum(e), e))


def unit(n, i, power=1):
    return tuple(power if j == i else 0 for j in range(n))


def radial(n, power):
    out = {(0,) * n: F(1)}
    base = {unit(n, i, 2): F(1) for i in range(n)}
    for _ in range(power):
        out = times(out, base)
    return out


def poly_json(p):
    return [{"powers": list(e), "coefficient": str(c)} for e, c in sorted(p.items()) if c]


def parse_poly(terms, n):
    if not isinstance(terms, list) or len(terms) > 150:
        raise ValueError("多项式必须是最多含 150 项的列表。")
    out = {}
    for term in terms:
        if not isinstance(term, dict) or set(term) != {"powers", "coefficient"}:
            raise ValueError("每一项必须包含 powers 和 coefficient。")
        e = term["powers"]
        if not isinstance(e, list) or len(e) != n or any(type(x) is not int or not 0 <= x <= 8 for x in e) or sum(e) > 8:
            raise ValueError("幂次必须为非负整数，且总次数不超过 8。")
        e = tuple(e)
        out[e] = out.get(e, F(0)) + exact(term["coefficient"])
    return clean(out)


def parse_problem(obj):
    if not isinstance(obj, dict) or set(obj) != {"title", "dimension", "vector_field"}:
        raise ValueError("稳定性问题必须且只能包含 title、dimension、vector_field；定义域固定为整个实空间。")
    n = obj["dimension"]
    if type(n) is not int or not 2 <= n <= 3:
        raise ValueError("首版支持 2 或 3 个变量。")
    if not isinstance(obj["title"], str) or len(obj["title"]) > 300:
        raise ValueError("title 必须为不超过 300 字符的字符串。")
    raw = obj["vector_field"]
    if not isinstance(raw, list) or len(raw) != n:
        raise ValueError("向量场分量数不正确。")
    field = [parse_poly(p, n) for p in raw]
    if any(p.get((0,) * n, 0) for p in field):
        raise ValueError("原点必须是平衡点。")
    return n, field


PLAN_SCHEMA = {"type": "object", "properties": {
    "family": {"type": "string", "enum": ["polynomial", "log_quadratic"]},
    "degree": {"type": "integer"}, "log_indices": {"type": "array", "items": {"type": "integer"}},
    "explanation": {"type": "string"}},
    "required": ["family", "degree", "log_indices", "explanation"], "additionalProperties": False}


def normalize_plan(plan, n):
    if not isinstance(plan, dict) or set(plan) != set(PLAN_SCHEMA["required"]):
        raise ValueError("研究计划字段不符合约定。")
    if plan["family"] not in ("polynomial", "log_quadratic"):
        raise ValueError("不支持的函数表示。")
    if type(plan["degree"]) is not int or plan["degree"] not in (2, 4, 6):
        raise ValueError("可搜索次数为 2、4、6。")
    indices = plan["log_indices"]
    if not isinstance(indices, list) or any(type(i) is not int or not 0 <= i < n for i in indices) or len(indices) != len(set(indices)):
        raise ValueError("log_indices 必须是不重复的变量下标。")
    if plan["family"] == "polynomial" and indices:
        raise ValueError("多项式计划不能含对数项。")
    if plan["family"] == "log_quadratic" and (not indices or plan["degree"] != 2):
        raise ValueError("对数模板的 degree 固定为 2，且至少指定一个对数变量。")
    if not isinstance(plan["explanation"], str) or len(plan["explanation"]) > 3000:
        raise ValueError("计划说明过长或类型错误。")
    return dict(plan, log_indices=sorted(indices))


def plan_prompt(problem, history):
    brief = [{k: r[k] for k in ("plan", "status", "feedback", "counterexample") if k in r} for r in history[-8:]]
    return ("研究下述多项式微分方程的原点是否全局渐近稳定。你负责选择 Lyapunov 函数表示，"
            "数值引擎求系数，精确验证器判断证书。搜索失败不说明不稳定，也不证明该表示不存在。"
            "polynomial: 搜索 2/4/6 次多项式；齐次奇次向量场采用齐次 V，其余采用含低阶项的 V。"
            "log_quadratic: V=sum_i a_i*log(1+x_i^2)（i在log_indices中）+sum_i a_i*x_i^2（其余i），a_i>0。"
            "如果多项式搜索连续失败，考虑乘法耦合所提示的对数表示。不要重复相同计划或改变原问题。"
            "返回计划 JSON；explanation 给出简短数学动机，不能自报已证明。\n问题："
            + json.dumps(problem, ensure_ascii=False) + "\n历史：" + json.dumps(brief, ensure_ascii=False)
            + "\nSchema：" + json.dumps(PLAN_SCHEMA, ensure_ascii=False))


def baseline_plans(n):
    for degree in (2, 4):
        yield {"family": "polynomial", "degree": degree, "log_indices": [], "explanation": "扩大多项式搜索空间。"}
    for count in range(1, n + 1):
        for indices in combinations(range(n), count):
            yield {"family": "log_quadratic", "degree": 2, "log_indices": list(indices), "explanation": "转换函数表示，检验对数项能否消除乘法耦合。"}
    yield {"family": "polynomial", "degree": 6, "log_indices": [], "explanation": "扩大到六次多项式。"}


def template(n, field, plan):
    """Return affine polynomials to certify; every target is sum c_i*T_i + T0."""
    zero = (0,) * n
    if plan["family"] == "polynomial":
        degree = plan["degree"]
        degrees = {sum(e) for f in field for e in f}
        homogeneous = len(degrees) == 1 and next(iter(degrees)) % 2 == 1
        vbasis = monomials(n, degree, degree if homogeneous else 2)
        pieces = [{e: F(1)} for e in vbasis]
        vdot = [scale(lie(p, field), -1) for p in pieces]
        max_d = max((sum(e) for p in vdot for e in p), default=2)
        lower_v = radial(n, degree // 2 if homogeneous else 1)
        lower_d = radial(n, max_d // 2 if homogeneous else 1)
        target_v = (pieces, scale(lower_v, -1), monomials(n, degree // 2, degree // 2 if homogeneous else 1))
        target_d = (vdot, scale(lower_d, -1), monomials(n, (max_d + 1) // 2, (max_d + 1) // 2 if homogeneous else 1))
        return {"coefficient_count": len(pieces), "v_basis": vbasis, "homogeneous": homogeneous,
                "targets": [target_v, target_d], "lower_v": lower_v, "lower_d": lower_d,
                "denominator": {zero: F(1)}}
    logs = set(plan["log_indices"])
    factors = {i: {zero: F(1), unit(n, i, 2): F(1)} for i in logs}
    denominator = {zero: F(1)}
    for factor in factors.values():
        denominator = times(denominator, factor)
    pieces = []
    for i in range(n):
        factor = {zero: F(1)}
        for j, p in factors.items():
            if j != i:
                factor = times(factor, p)
        if i not in logs:
            factor = denominator
        pieces.append(scale(times(times({unit(n, i): F(1)}, field[i]), factor), -2))
    lower = radial(n, 1)
    support = set(lower) | {e for p in pieces for e in p}
    # A sparse sufficient basis; omission only causes missed certificates.
    basis = sorted({tuple(x // 2 for x in e) for e in support if all(x % 2 == 0 for x in e) and sum(e) > 0})
    return {"coefficient_count": n, "v_basis": [], "homogeneous": False,
            "targets": [(pieces, scale(lower, -1), basis)], "lower_d": lower,
            "denominator": denominator}


def build_affine(spec):
    count = spec["coefficient_count"]
    blocks, width = [], count
    for _, _, basis in spec["targets"]:
        entries = [(i, j) for i in range(len(basis)) for j in range(i, len(basis))]
        blocks.append({"offset": width, "basis": basis, "entries": entries})
        width += len(entries)
    if width > 600:
        raise ValueError("当前表示超过 600 个搜索变量的资源限制。")
    rows, rhs = [], []
    for (pieces, constant, _), block in zip(spec["targets"], blocks):
        basis, entries = block["basis"], block["entries"]
        support = set(constant) | {e for p in pieces for e in p}
        support |= {tuple(a + b for a, b in zip(basis[i], basis[j])) for i, j in entries}
        for e in sorted(support):
            row = [F(0)] * width
            for i, p in enumerate(pieces):
                row[i] = p.get(e, F(0))
            for k, (i, j) in enumerate(entries):
                if tuple(a + b for a, b in zip(basis[i], basis[j])) == e:
                    row[block["offset"] + k] = F(-1 if i == j else -2)
            rows.append(row)
            rhs.append(-constant.get(e, F(0)))
    return rows, rhs, blocks, width


def exact_psd(matrix):
    """Exact symmetric Schur elimination, including singular PSD matrices."""
    a = [[F(x) for x in row] for row in matrix]
    n = len(a)
    if not n or any(len(row) != n for row in a) or any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        return False
    for k in range(n):
        pivot = a[k][k]
        if pivot < 0:
            return False
        if pivot == 0:
            if any(a[k][j] for j in range(k + 1, n)):
                return False
            continue
        for i in range(k + 1, n):
            for j in range(i, n):
                a[i][j] -= a[i][k] * a[k][j] / pivot
                a[j][i] = a[i][j]
    return True


def matrices_from(vector, blocks):
    matrices = []
    for block in blocks:
        size = len(block["basis"])
        matrix = [[F(0)] * size for _ in range(size)]
        for k, (i, j) in enumerate(block["entries"]):
            matrix[i][j] = matrix[j][i] = vector[block["offset"] + k]
        matrices.append(matrix)
    return matrices


def rref(rows, rhs):
    a = [list(row) + [b] for row, b in zip(rows, rhs)]
    width = len(rows[0])
    pivots, r = [], 0
    for col in range(width):
        pivot = next((i for i in range(r, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        factor = a[r][col]
        a[r] = [x / factor for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][col]:
                factor = a[i][col]
                a[i] = [x - factor * y for x, y in zip(a[i], a[r])]
        pivots.append(col)
        r += 1
        if r == len(a):
            break
    if any(not any(row[:-1]) and row[-1] for row in a):
        return None
    return a[:r], pivots


def recover_exact(numeric, reduced, blocks, count, positive_coefficients):
    if reduced is None:
        return None
    rows, pivots = reduced
    free = [j for j in range(len(numeric)) if j not in pivots]
    for denominator in (10, 100, 1000, 10000, 1000000):
        vector = [F(0)] * len(numeric)
        for j in free:
            vector[j] = F(float(numeric[j])).limit_denominator(denominator)
        for row, pivot in zip(rows, pivots):
            vector[pivot] = row[-1] - sum((row[j] * vector[j] for j in free), F(0))
        if positive_coefficients and any(c < 1 for c in vector[:count]):
            continue
        matrices = matrices_from(vector, blocks)
        if all(exact_psd(matrix) for matrix in matrices):
            return vector, matrices
    return None


def synthesize(n, field, plan, iterations=2500, seconds=12):
    import numpy as np
    start = time.monotonic()
    spec = template(n, field, plan)
    rows, rhs, blocks, width = build_affine(spec)
    reduced = rref(rows, rhs)
    if reduced is None:
        return {"status": "search_failed", "feedback": "当前证书的系数恒等式无解；这不是系统不稳定的证明。"}
    a, b = np.array(rows, dtype=float), np.array(rhs, dtype=float)
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        return {"status": "search_failed", "feedback": "数值范围过大，无法完成此次证书搜索。"}
    weights = np.ones(width)
    for block in blocks:
        for k, (i, j) in enumerate(block["entries"]):
            weights[block["offset"] + k] = 1 if i == j else np.sqrt(2)
    aw = a / weights
    try:
        projection = np.linalg.pinv(aw, rcond=1e-11)
    except np.linalg.LinAlgError:
        return {"status": "search_failed", "feedback": "数值分解失败，没有产生可验证证书。"}
    def affine(x):
        return x + projection @ (b - aw @ x)
    positive_coefficients = plan["family"] == "log_quadratic"
    def cone(x):
        out = x.copy()
        if positive_coefficients:
            out[:spec["coefficient_count"]] = np.maximum(out[:spec["coefficient_count"]], 1.0)
        for block in blocks:
            matrix = np.zeros((len(block["basis"]), len(block["basis"])))
            for k, (i, j) in enumerate(block["entries"]):
                pos = block["offset"] + k
                matrix[i, j] = matrix[j, i] = out[pos] / weights[pos]
            eigenvalues, eigenvectors = np.linalg.eigh(matrix)
            matrix = (eigenvectors * np.maximum(eigenvalues, 0)) @ eigenvectors.T
            for k, (i, j) in enumerate(block["entries"]):
                pos = block["offset"] + k
                out[pos] = matrix[i, j] * weights[pos]
        return out
    x = np.zeros(width)
    best, best_residual = None, float("inf")
    for step in range(1, iterations + 1):
        z = cone(x)
        x += affine(2 * z - x) - z
        if not np.all(np.isfinite(x)):
            return {"status": "search_failed", "feedback": "数值迭代溢出，没有产生可验证证书。"}
        if step % 25 == 0:
            candidate = affine(z) / weights
            residual = float(np.max(np.abs(aw @ z - b)))
            if residual < best_residual:
                best, best_residual = candidate, residual
            recovered = recover_exact(candidate, reduced, blocks, spec["coefficient_count"], positive_coefficients) if residual < 1e-5 else None
            if recovered:
                vector, matrices = recovered
                certificate = {"plan": plan, "coefficients": [str(c) for c in vector[:spec["coefficient_count"]]],
                    "grams": [[[str(x) for x in row] for row in matrix] for matrix in matrices]}
                return {"status": "certified", "certificate": certificate, "iterations": step,
                        "elapsed_seconds": round(time.monotonic() - start, 3), "search_variables": width}
            if time.monotonic() - start > seconds:
                break
    result = {"status": "search_failed", "iterations": step, "elapsed_seconds": round(time.monotonic() - start, 3),
              "search_variables": width, "best_numerical_residual": best_residual,
              "feedback": "预算内没有找到能精确复核的证书。可换表示或提高预算；不能据此断言不稳定。"}
    if best is not None:
        candidate_coeff = [F(float(c)).limit_denominator(1000) for c in best[:spec["coefficient_count"]]]
        for point in product(range(-3, 4), repeat=n):
            if not any(point):
                continue
            for pieces, constant, _ in spec["targets"]:
                target = plus(constant, *(scale(p, c) for p, c in zip(pieces, candidate_coeff)))
                measured = value(target, point)
                if measured < 0:
                    result["counterexample"] = {"point": list(point), "candidate_coefficients": [str(c) for c in candidate_coeff],
                        "violated_margin_value": str(measured), "scope": "仅反驳该近似候选的所需裕量，不反驳系统稳定性"}
                    return result
    return result


def gram_polynomial(matrix, basis):
    out = {}
    for i, row in enumerate(matrix):
        for j, coefficient in enumerate(row):
            e = tuple(a + b for a, b in zip(basis[i], basis[j]))
            out[e] = out.get(e, F(0)) + coefficient
    return clean(out)


def verify_certificate(problem, certificate):
    """Recompute V, its derivative, positivity, and all SOS identities exactly."""
    n, field = parse_problem(problem)
    if not isinstance(certificate, dict) or set(certificate) != {"plan", "coefficients", "grams"}:
        return False
    plan = normalize_plan(certificate["plan"], n)
    spec = template(n, field, plan)
    coeff_raw = certificate["coefficients"]
    if not isinstance(coeff_raw, list) or len(coeff_raw) != spec["coefficient_count"]:
        return False
    coefficients = [exact(c) for c in coeff_raw]
    if plan["family"] == "log_quadratic" and any(c < 1 for c in coefficients):
        return False
    grams = certificate["grams"]
    if not isinstance(grams, list) or len(grams) != len(spec["targets"]):
        return False
    for raw, (pieces, constant, basis) in zip(grams, spec["targets"]):
        if not isinstance(raw, list) or len(raw) != len(basis) or any(not isinstance(row, list) or len(row) != len(basis) for row in raw):
            return False
        matrix = [[exact(c) for c in row] for row in raw]
        if not exact_psd(matrix):
            return False
        expected = plus(constant, *(scale(p, c) for p, c in zip(pieces, coefficients)))
        if expected != gram_polynomial(matrix, basis):
            return False
    return True


def describe_certificate(problem, certificate):
    n, field = parse_problem(problem)
    spec = template(n, field, certificate["plan"])
    coefficients = [exact(c) for c in certificate["coefficients"]]
    if certificate["plan"]["family"] == "polynomial":
        return {"V_polynomial": poly_json(dict(zip(spec["v_basis"], coefficients))),
                "V_lower_bound": poly_json(spec["lower_v"]), "negative_derivative_lower_bound": poly_json(spec["lower_d"])}
    return {"V_terms": [{"variable": i, "coefficient": str(c),
                         "function": "log(1+x_i^2)" if i in certificate["plan"]["log_indices"] else "x_i^2"} for i, c in enumerate(coefficients)],
            "positive_denominator": poly_json(spec["denominator"]),
            "negative_derivative_numerator_lower_bound": poly_json(spec["lower_d"])}


def run_research(problem, provider="baseline", model=None, rounds=6, iterations=2500, seconds=12, plans=None):
    n, field = parse_problem(problem)
    if not 1 <= rounds <= 20 or not 25 <= iterations <= 100000 or not 1 <= seconds <= 60:
        raise ValueError("rounds 范围为 1..20，iterations 为 25..100000，seconds 为 1..60。")
    proposer = GPTProposer(model) if provider == "gpt" else None
    schedule = iter(plans if plans is not None else baseline_plans(n))
    history, seen = [], set()
    result = {"version": "0.2.0", "problem": problem, "provider": "manual" if plans is not None else provider, "model": model if provider == "gpt" else None,
              "status": "not_solved", "novelty": "not_assessed", "formal_assistant_checked": False,
              "budget": {"rounds": rounds, "iterations_per_round": iterations, "seconds_per_round_soft": seconds}, "history": history}
    for turn in range(1, rounds + 1):
        try:
            plan = proposer.request_structured(plan_prompt(problem, history), PLAN_SCHEMA) if proposer else next(schedule)
            plan = normalize_plan(plan, n)
            key = (plan["family"], plan["degree"], tuple(plan["log_indices"]))
            if key in seen:
                history.append({"round": turn, "plan": plan, "status": "duplicate", "feedback": "该表示已经搜索过，请换表示。"})
                continue
            seen.add(key)
            found = synthesize(n, field, plan, iterations, seconds)
            record = {"round": turn, "plan": plan, **found}
            history.append(record)
            if found["status"] == "certified":
                if not verify_certificate(problem, found["certificate"]):
                    raise ArithmeticError("独立复核失败，拒绝输出稳定性结论。")
                result["status"] = "certified_global_asymptotic_stability"
                result["certificate"] = found["certificate"]
                result["readable_result"] = describe_certificate(problem, found["certificate"])
                break
        except StopIteration:
            break
        except ProviderError as exc:
            history.append({"round": turn, "status": "provider_error", "feedback": str(exc)})
            result["status"] = "provider_error"
            break
        except ValueError as exc:
            history.append({"round": turn, "status": "invalid_plan", "feedback": str(exc)})
    result["api_usage"] = proposer.usage if proposer else []
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description="非线性稳定性研究：换表示、搜索 SOS 证书、精确复核。")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run")
    run.add_argument("problem")
    run.add_argument("--provider", choices=["baseline", "gpt"], default="baseline")
    run.add_argument("--model", default=os.environ.get("OPENAI_MODEL"))
    run.add_argument("--rounds", type=int, default=6)
    run.add_argument("--iterations", type=int, default=2500)
    run.add_argument("--seconds", type=float, default=12)
    run.add_argument("--plans", help="手动从 GPT 获得的计划列表 JSON")
    run.add_argument("--out", required=True)
    prompt = sub.add_parser("prompt")
    prompt.add_argument("problem")
    prompt.add_argument("--history")
    verify = sub.add_parser("verify")
    verify.add_argument("result")
    args = parser.parse_args(argv)
    try:
        if args.command == "verify":
            result = load_json(args.result)
            ok = verify_certificate(result["problem"], result["certificate"])
            print("全局渐近稳定性证书：精确复核通过。" if ok else "证书复核未通过。")
            return 0 if ok else 1
        problem = load_json(args.problem)
        parse_problem(problem)
        if args.command == "prompt":
            previous = load_json(args.history) if args.history else {"problem": problem, "history": []}
            if previous["problem"] != problem:
                raise ValueError("历史记录属于另一个问题。")
            print(plan_prompt(problem, previous["history"]))
            return 0
        if args.plans and args.provider == "gpt":
            raise ValueError("--plans 与 --provider gpt 不能同时使用。")
        plans = load_json(args.plans) if args.plans else None
        if plans is not None and not isinstance(plans, list):
            raise ValueError("--plans 文件必须包含计划列表。")
        result = run_research(problem, args.provider, args.model, args.rounds, args.iterations, args.seconds, plans)
        save_json(args.out, result)
        print(f"状态：{result['status']}；完成 {len(result['history'])} 轮；结果：{args.out}")
        return 0 if result["status"] == "certified_global_asymptotic_stability" else 1
    except ImportError:
        print("搜索需要 NumPy；精确复核不需要 NumPy。请按 README 使用运行环境。", file=sys.stderr)
        return 2
    except (ValueError, TypeError, KeyError, OSError) as exc:
        print("错误：" + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
