#!/usr/bin/env python3
"""GPT-assisted conjecture search with exact, domain-specific verification.

Python 3.9+, standard library only. Coefficients use ascending powers.
The baseline is deterministic and does not call or simulate GPT.
"""
import argparse
from dataclasses import dataclass
from fractions import Fraction as F
from functools import reduce
import json
from math import comb, gcd
import os
from pathlib import Path
import random
import re
import sys
import urllib.error
import urllib.request

VERSION = "0.1.0"
MAX_DEGREE = 20
MAX_CANDIDATE_DEGREE = 32


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p or [F(0)]


def exact(value):
    if type(value) is int:
        value = str(value)
    if not isinstance(value, str) or not re.fullmatch(r"-?\d{1,100}(?:/[1-9]\d{0,99})?", value):
        raise ValueError("系数必须是整数或分数字符串，例如 3、'1/2'；不接受浮点数。")
    return F(value)


def polynomial(values, limit=MAX_CANDIDATE_DEGREE):
    if not isinstance(values, list) or not 1 <= len(values) <= limit + 1:
        raise ValueError("多项式系数列表为空或超出次数限制。")
    return trim([exact(x) for x in values])


def evaluate(p, n):
    out = F(0)
    for c in reversed(p):
        out = out * n + c
    return out


def add(a, b):
    out = [F(0)] * max(len(a), len(b))
    for p in (a, b):
        for i, c in enumerate(p):
            out[i] += c
    return trim(out)


def multiply(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def shift(p, amount):
    out = [F(0)] * len(p)
    for i, c in enumerate(p):
        for j in range(i + 1):
            out[j] += c * comb(i, j) * amount ** (i - j)
    return trim(out)


def encoded(p):
    return [str(x) for x in trim(p)]


def display(p, variable="n"):
    terms = []
    for i in reversed(range(len(p))):
        c = p[i]
        if not c:
            continue
        power = "" if i == 0 else variable if i == 1 else f"{variable}^{i}"
        term = str(abs(c)) if not power else ("" if abs(c) == 1 else f"({abs(c)})*") + power
        terms.append(("-" if c < 0 else "+", term))
    if not terms:
        return "0"
    sign, term = terms[0]
    return ("-" if sign == "-" else "") + term + "".join(f" {s} {t}" for s, t in terms[1:])


@dataclass(frozen=True)
class Problem:
    kind: str
    coefficients: tuple
    title: str = "数学研究实验"

    @classmethod
    def parse(cls, obj):
        if not isinstance(obj, dict) or set(obj) - {"kind", "coefficients", "title"}:
            raise ValueError("问题只能包含 kind、coefficients、title。")
        if obj.get("kind") not in ("sum", "fixed_divisor"):
            raise ValueError("目前支持 sum 和 fixed_divisor 两类问题。")
        p = polynomial(obj.get("coefficients"), MAX_DEGREE)
        if obj["kind"] == "fixed_divisor":
            if any(c.denominator != 1 for c in p):
                raise ValueError("固定除数问题要求整数系数。")
            if p == [0]:
                raise ValueError("零多项式的所有正整数都是公因子，没有最大的正公因子。")
        title = obj.get("title", "数学研究实验")
        if not isinstance(title, str) or len(title) > 300:
            raise ValueError("title 必须是长度不超过 300 的字符串。")
        return cls(obj["kind"], tuple(p), title)

    def json(self):
        return {"kind": self.kind, "coefficients": encoded(self.coefficients), "title": self.title}

    def statement(self):
        if self.kind == "sum":
            return f"求 Q(n)，使所有整数 n≥0 都有 sum(k=1..n, {display(self.coefficients, 'k')}) = Q(n)。"
        return f"求最大的正整数 m，使所有整数 n 都有 m | ({display(self.coefficients)})。"


def normalize_candidate(problem, obj):
    field = "coefficients" if problem.kind == "sum" else "divisor"
    if not isinstance(obj, dict) or set(obj) != {field, "explanation"}:
        raise ValueError(f"候选必须且只能包含 {field} 和 explanation。")
    if not isinstance(obj["explanation"], str) or len(obj["explanation"]) > 4000:
        raise ValueError("explanation 必须是长度不超过 4000 的字符串。")
    if problem.kind == "sum":
        value = encoded(polynomial(obj[field]))
    else:
        v = exact(obj[field])
        if v.denominator != 1 or v <= 0:
            raise ValueError("divisor 必须为正整数。")
        value = str(v)
    return {field: value, "explanation": obj["explanation"]}


def sample_checks(problem, candidate, points):
    """Finite checks only. Success is deliberately not a proof status."""
    points = sorted(set(points))
    if any(type(n) is not int or (problem.kind == "sum" and n < 0) for n in points):
        raise ValueError("检验点不属于问题的定义域。")
    p = problem.coefficients
    sums = {0: F(0)}
    if problem.kind == "sum":
        q = polynomial(candidate["coefficients"])
        total = F(0)
        for n in range(1, max(points, default=0) + 1):
            total += evaluate(p, n)
            if n in points:
                sums[n] = total
    for index, n in enumerate(points):
        if problem.kind == "sum":
            actual, predicted = sums[n], evaluate(q, n)
            witness = {"n": n, "actual_sum": str(actual), "candidate_value": str(predicted)}
            failed = actual != predicted
        else:
            actual, m = int(evaluate(p, n)), int(candidate["divisor"])
            witness = {"n": n, "polynomial_value": str(actual), "divisor": str(m), "remainder": str(actual % m)}
            failed = actual % m != 0
        if failed:
            return {"status": "refuted", "checked": index + 1, "counterexample": witness}
    return {"status": "tested_only", "checked": len(points), "points": points}


def newton_coefficients(p):
    row = [evaluate(p, n) for n in range(len(p))]
    out = []
    while row:
        out.append(row[0])
        row = [b - a for a, b in zip(row, row[1:])]
    return out


def certify(problem, candidate):
    """Recompute proof obligations from the trusted problem and candidate data."""
    p = problem.coefficients
    if problem.kind == "sum":
        q = polynomial(candidate["coefficients"])
        residual = add(add(q, [-c for c in shift(q, -1)]), [-c for c in p])
        cert = {"method": "polynomial_induction", "domain": "integers n >= 0",
                "base_value": str(evaluate(q, 0)), "residual_coefficients": encoded(residual)}
        if evaluate(q, 0) == 0 and residual == [0]:
            return {"status": "certified", "certificate": cert}
        # The error is a polynomial of degree <= max(deg(Q), deg(P)+1).
        # Hence a wrong identity must fail at one of these D+1 points.
        degree_bound = max(len(q) - 1, len(p))
        result = sample_checks(problem, candidate, range(degree_bound + 1))
        if result["status"] != "refuted":
            raise ArithmeticError("内部错误：证明检查与反例搜索不一致。")
        result["failed_obligations"] = cert
        return result
    m = int(candidate["divisor"])
    differences = newton_coefficients(p)
    if any(c.denominator != 1 for c in differences):
        raise ArithmeticError("内部错误：整数系数多项式出现非整数差分。")
    # P(n) = sum_j Δ^j P(0) * binomial(n,j); binomial(n,j) is an
    # integer for every integer n. The finite-difference transform is
    # invertible over integers, so the gcd is also the maximal divisor.
    maximal = reduce(gcd, (abs(int(c)) for c in differences), 0)
    cert = {"method": "integer_newton_basis", "domain": "all integers n",
            "forward_differences": encoded(differences), "gcd": str(maximal)}
    if any(int(c) % m for c in differences):
        result = sample_checks(problem, candidate, range(len(p)))
        if result["status"] != "refuted":
            raise ArithmeticError("内部错误：差分检查与反例搜索不一致。")
        return result
    if m != maximal:
        return {"status": "feasible_not_optimal", "certificate": cert,
                "feedback": "这个数整除所有取值，但不是最大公因子。"}
    return {"status": "certified", "certificate": cert}


def interpolate(values):
    result = [F(0)]
    for i, value in enumerate(values):
        basis, denominator = [F(1)], 1
        for j in range(len(values)):
            if i != j:
                basis = multiply(basis, [-F(j), F(1)])
                denominator *= i - j
        result = add(result, [c * value / denominator for c in basis])
    return result


def baseline_propose(problem, history, round_number):
    if problem.kind == "sum":
        # Increase the interpolation degree after failures. This is a
        # transparent classical baseline, not an LLM result.
        degree = min(round_number, len(problem.coefficients))
        values, total = [F(0)], F(0)
        for n in range(1, degree + 1):
            total += evaluate(problem.coefficients, n)
            values.append(total)
        return {"coefficients": encoded(interpolate(values)),
                "explanation": f"确定性基线：用 n=0..{degree} 的求和值作精确插值。"}
    values = [abs(int(evaluate(problem.coefficients, n))) for n in range(round_number + 2)]
    m = reduce(gcd, values, 0) or 1
    return {"divisor": str(m), "explanation": "确定性基线：取逐步扩大的有限样本的最大公因子。"}


def candidate_schema(problem):
    field = "coefficients" if problem.kind == "sum" else "divisor"
    field_schema = {"type": "array", "items": {"type": "string"}} if field == "coefficients" else {"type": "string"}
    return {"type": "object", "properties": {field: field_schema, "explanation": {"type": "string"}},
            "required": [field, "explanation"], "additionalProperties": False}


def prompt_for(problem, history):
    return ("你是数学猜想生成器。根据精确问题和失败记录提出一个候选。"
            "explanation 只写简短的数学依据。验证器独立判断真假；不要自报已证明。"
            "系数按常数项、一次项、二次项递增排列，只用整数或分数字符串，禁止浮点数。"
            "不要重复被反例推翻的候选。你可以改变公式，不能改变原问题的定义域或前提。\n"
            + problem.statement() + "\n问题：" + json.dumps(problem.json(), ensure_ascii=False)
            + "\n最近记录：" + json.dumps(history[-6:], ensure_ascii=False)
            + "\n返回符合以下 JSON Schema 的 JSON 对象：" + json.dumps(candidate_schema(problem), ensure_ascii=False))


class ProviderError(RuntimeError):
    pass


def extract_response(body):
    if not isinstance(body, dict) or body.get("status") != "completed":
        raise ProviderError("API 响应未完成；本次没有得到可检验候选。")
    texts = []
    for item in body.get("output", []):
        if item.get("type") != "message":
            continue
        for content in item.get("content", []):
            if content.get("type") == "refusal":
                raise ProviderError("模型未提供候选。")
            if content.get("type") == "output_text":
                texts.append(content.get("text", ""))
    if not texts:
        raise ProviderError("API 响应没有候选文本。")
    try:
        return json.loads("".join(texts))
    except (ValueError, TypeError) as exc:
        raise ProviderError("API 返回的候选不是有效 JSON。") from exc


class GPTProposer:
    def __init__(self, model, max_output_tokens=3000):
        self.key = os.environ.get("OPENAI_API_KEY")
        if not self.key:
            raise ValueError("未设置 OPENAI_API_KEY。可先用 baseline 模式，或用 prompt/check 手动连接 GPT。")
        if not model:
            raise ValueError("GPT 模式需要 --model 指定你账户可用且支持结构化输出的模型。")
        self.model, self.max_output_tokens = model, max_output_tokens
        self.usage = []

    def __call__(self, problem, history, round_number):
        return self.request_structured(prompt_for(problem, history), candidate_schema(problem))

    def request_structured(self, prompt, schema):
        payload = {"model": self.model, "input": [{"role": "user", "content": prompt}],
                   "store": False, "max_output_tokens": self.max_output_tokens,
                   "text": {"format": {"type": "json_schema", "name": "math_candidate",
                                       "strict": True, "schema": schema}}}
        request = urllib.request.Request("https://api.openai.com/v1/responses",
            data=json.dumps(payload).encode(),
            headers={"Content-Type": "application/json", "Authorization": "Bearer " + self.key})
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                raw = response.read(2_000_001)
                if len(raw) > 2_000_000:
                    raise ProviderError("API 响应超出大小限制。")
                body = json.loads(raw)
        except urllib.error.HTTPError as exc:
            raise ProviderError(f"API 请求失败，HTTP {exc.code}。检查模型权限、额度和请求参数。") from exc
        except (urllib.error.URLError, TimeoutError, OSError, ValueError) as exc:
            raise ProviderError("API 连接或响应解析失败；为避免重复计费，没有自动重试。") from exc
        if not isinstance(body, dict):
            raise ProviderError("API 响应不是对象。")
        self.usage.append({"response_id": body.get("id"), "usage": body.get("usage")})
        return extract_response(body)


def research(problem, proposer, rounds=8, seed=0, provider_name="custom", model=None):
    if type(rounds) is not int or not 1 <= rounds <= 100:
        raise ValueError("rounds 必须在 1..100 之间。")
    rng = random.Random(seed)
    points = set(range(9)) | {rng.randint(9, 250) for _ in range(8)}
    if problem.kind == "fixed_divisor":
        points |= {-x for x in points}
    history, seen, status = [], set(), "budget_exhausted"
    for turn in range(1, rounds + 1):
        record = {"round": turn}
        try:
            raw = proposer(problem, history, turn)
            candidate = normalize_candidate(problem, raw)
        except ProviderError as exc:
            history.append({"round": turn, "status": "provider_error", "error": str(exc)})
            status = "provider_error"
            break
        except ValueError as exc:
            history.append({"round": turn, "status": "invalid_candidate", "error": str(exc)})
            continue
        record["candidate"] = candidate
        identity = json.dumps({k: v for k, v in candidate.items() if k != "explanation"}, sort_keys=True)
        if identity in seen:
            record.update(status="duplicate", feedback="请提出与已有候选不同的公式或除数。")
        else:
            seen.add(identity)
            experiment = sample_checks(problem, candidate, points)
            record["experiment"] = experiment
            result = experiment if experiment["status"] == "refuted" else certify(problem, candidate)
            record.update(result)
        history.append(record)
        if record["status"] == "certified":
            status = "certified"
            break
    return {"version": VERSION, "problem": problem.json(), "provider": provider_name,
            "model": model, "seed": seed, "round_limit": rounds,
            "status": status, "novelty": "not_assessed", "formal_assistant_checked": False,
            "history": history, "api_usage": getattr(proposer, "usage", [])}


def verify_run(run):
    """Ignore self-reported truth; independently reconstruct the obligations."""
    problem = Problem.parse(run["problem"])
    if run.get("status") != "certified" or not run.get("history"):
        return False
    record = run["history"][-1]
    candidate = normalize_candidate(problem, record["candidate"])
    checked = certify(problem, candidate)
    return (record.get("status") == "certified" and checked["status"] == "certified"
            and record.get("certificate") == checked["certificate"])


def save_json(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def load_json(path):
    if Path(path).stat().st_size > 5_000_000:
        raise ValueError("输入文件超过 5MB 限制。")
    return json.loads(Path(path).read_text(encoding="utf-8"))


def main(argv=None):
    parser = argparse.ArgumentParser(description="GPT 数学研究原型：猜想、反例、精确证明检查。")
    sub = parser.add_subparsers(dest="command", required=True)
    run = sub.add_parser("run", help="自动研究一个问题")
    run.add_argument("problem")
    run.add_argument("--provider", choices=["baseline", "gpt"], default="baseline")
    run.add_argument("--model", default=os.environ.get("OPENAI_MODEL"))
    run.add_argument("--rounds", type=int, default=8)
    run.add_argument("--seed", type=int, default=0)
    run.add_argument("--max-output-tokens", type=int, default=3000)
    run.add_argument("--out", required=True)
    prompt = sub.add_parser("prompt", help="生成可粘贴到 GPT 的研究提示")
    prompt.add_argument("problem")
    prompt.add_argument("--history", help="先前的 run 或 check 结果")
    check = sub.add_parser("check", help="检查从 GPT 得到的候选 JSON")
    check.add_argument("problem")
    check.add_argument("candidate")
    check.add_argument("--out", required=True)
    verify = sub.add_parser("verify", help="重新检查保存的证明记录")
    verify.add_argument("result")
    args = parser.parse_args(argv)
    try:
        if args.command == "verify":
            passed = verify_run(load_json(args.result))
            print("证明记录复核通过。" if passed else "证明记录复核未通过。")
            return 0 if passed else 1
        problem = Problem.parse(load_json(args.problem))
        if args.command == "prompt":
            history = []
            if args.history:
                previous = load_json(args.history)
                if previous.get("problem") != problem.json():
                    raise ValueError("历史记录属于另一个问题。")
                history = previous["history"]
            print(prompt_for(problem, history))
            return 0
        if args.command == "check":
            candidate = load_json(args.candidate)
            result = research(problem, lambda *_: candidate, rounds=1, provider_name="manual")
        else:
            if not 256 <= args.max_output_tokens <= 32000:
                raise ValueError("max-output-tokens 必须在 256..32000 之间。")
            proposer = GPTProposer(args.model, args.max_output_tokens) if args.provider == "gpt" else baseline_propose
            result = research(problem, proposer, args.rounds, args.seed, args.provider,
                              args.model if args.provider == "gpt" else None)
        save_json(args.out, result)
        labels = {"certified": "精确证明检查通过", "budget_exhausted": "轮数耗尽，尚未完成证明", "provider_error": "GPT 接入失败"}
        print(f"{labels[result['status']]}；已记录 {len(result['history'])} 轮。结果：{args.out}")
        return 0 if result["status"] == "certified" else 1
    except (ValueError, KeyError, TypeError, OSError) as exc:
        print("错误：" + str(exc), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
