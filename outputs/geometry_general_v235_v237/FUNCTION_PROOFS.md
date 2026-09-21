# Exact source functions and remainder evidence

`function_models.py` supplies the source-function front end for the three spectral bridges V235–V237. It does not add a fourth version and does not itself transfer any spectral inequality. Its certificates assert only a polynomial approximation and a proved function-norm error. All integrations use the probability measure `dμ=dσ/(4π)` on the unit sphere, with certificate scope `unit_S2_probability_measure`.

## Public interface and source binding

The public entry points are `normalize_function(function)`, `analytic_model(function, order=4, sqrt_bits=40)`, `l2_model(function, degree=8, sqrt_bits=40)`, `sqrt_upper(value_squared, bits=40)`, and `verify_model(certificate, expected_function=None, expected_polynomial=None)`.

Every accepted numerical input is an integer, `Fraction`, or exact rational string. Floats and booleans are rejected. A polynomial is an exponent-to-coefficient dictionary such as `{"0,0,1":"1/4"}`; the three exponents refer to `(x,y,z)`. Polynomial degree is at most 24. Keys must be nonnegative integer triples, and duplicate triples produced by parsing distinct keys are rejected. Unknown source fields are rejected.

An analytic source has the form

```python
{"kind": "analytic_sum", "polynomial": {}, "terms": [
    {"function": "exp", "argument": {"0,0,1": "1/4"}, "coefficient": "1"},
    {"function": "sin", "argument": {"1,1,0": "1/4"}, "coefficient": "1"},
    {"function": "log1p", "argument": {"0,0,1": "1/4"}, "coefficient": "1"}
]}
```

The defaults are `polynomial={}`, `terms=[]`, and, in each term, `argument={}` and `coefficient=1`. The kind and nonempty term's function are required. At most 12 input terms are accepted. Each original `log1p` term must first satisfy the argument-radius condition `R<1`, including a zero-coefficient term. Only after that domain check are identical function/argument pairs combined with their signed coefficients and zero sums removed. Thus cancellation cannot hide an undefined logarithm. These canonicalization rules are also applied to an expected source during verification.

An axis source has the form

```python
{"kind": "axis_profile", "axis": 2, "profile": "abs_power",
 "exponent": "-1/4", "amplitude": "-1", "offset": "0"}
```

It represents `offset + amplitude*g(x_axis)`, where `g(t)=1_{t>0}` for `step`, or `g(t)=|t|^alpha` for `abs_power`. Defaults are axis 2, amplitude 1, offset 0, and, for `abs_power`, exponent 1. A step source does not accept an exponent field. The power requires `alpha > -1/2`, including when its amplitude is zero; this is an explicit input-family restriction. Amplitudes and offsets may have either sign.

## Uniform analytic remainder

For an argument polynomial `a`, set `R=sum_e |a_e|`. Each coordinate has absolute value at most 1 on the unit sphere, so `|a(x)| <= R`. The approximation is the Taylor polynomial through order `n`, composed with `a`, with all coefficients and polynomial products computed exactly.

For `exp`, the first omitted absolute term is `R^(n+1)/(n+1)!`. If `rho=R/(n+2)<1`, every subsequent term ratio is at most `rho`, giving the complete-tail bound

`B_exp = R^(n+1) / ((n+1)! * (1-rho))`.

Otherwise Taylor's real-variable remainder gives `exp(R)*R^(n+1)/(n+1)!`. The implementation bounds `exp(R)` by `3^ceil(R)`. The elementary series proves `e<3`: for `k>=2`, `k!>=2^(k-1)`, with strict inequality for `k>=3`, hence the tail after `1+1` is strictly less than 1. The envelope is therefore a rigorous rational upper bound. An explicit envelope budget rejects `ceil(R)>4096` in this fallback branch.

For `sin`, all real derivatives have absolute value at most 1, so Taylor's remainder is bounded by

`B_sin = R^(n+1)/(n+1)!`.

For `log1p`, the implementation requires `R<1`. Absolute convergence of the entire remaining series gives

`B_log = sum_(k>=n+1) R^k/k <= R^(n+1)/((n+1)*(1-R))`.

The source's exact polynomial part contributes no error. For signed coefficients `c_j`, the final bound is `sum_j |c_j| B_j`, by the triangle inequality. Each certificate includes every source term, its argument bound, Taylor coefficients, composed polynomial, and complete-tail parameters. Deleted tail evidence is rejected by replay.

Both `order` and expanded approximation degree are at most 24. The implementation does not generate powers beyond the highest nonzero retained Taylor coefficient: for example, order 4 for `sin(x^8)` retains degrees 8 and 24 and does not attempt the unused degree-32 power. Polynomial multiplication has an explicit 262144-pair budget and input dictionaries have a 4096-monomial budget. Bounds remain valid when conservative; no claim of best approximation is made. `sqrt_bits` is validated and recorded but is unused for this `Linf` certificate.

## Exact probability-L2 projection

Rotation invariance implies that any coordinate `t=x_axis` has marginal measure `dt/2` on `[-1,1]`. Consequently the same one-variable formulas apply to all three axes.

For nonnegative integer `k`, the profile moments are

`E[t^k * 1_(t>0)] = 1/(2*(k+1))`,

and

`E[t^k * |t|^alpha] = 1/(alpha+k+1)` for even `k`, and 0 for odd `k`.

For `q=b+A*g`, each source moment is `b*E[t^k]+A*E[t^k*g]`. Its exact squared norm is

`||q||_2^2 = b^2 + 2*b*A*E[g] + A^2*E[g^2]`,

where `E[g^2]=1/2` for the step and `1/(2*alpha+1)` for an absolute power. The signed offset/amplitude cross term is retained. The strict integrability threshold ensures the denominator is positive and the source belongs to L2.

The exact Legendre recurrence constructs `P_l` with probability norm `||P_l||_2^2=1/(2l+1)`. Set `b_l=<q,P_l>` and `a_l=(2l+1)*b_l`. The returned polynomial is

`p_d=sum_(l=0)^d a_l P_l`,

and orthogonality gives the complete, exact error identity

`||q-p_d||_2^2 = ||q||_2^2 - sum_(l=0)^d (2l+1)*b_l^2`.

This difference accounts for every omitted mode, without truncating or estimating an infinite tail by sampled values. The certificate records all monomial moments, all coefficients of every retained Legendre polynomial, all inner products and norm contributions, and both exact squared norms. The result is embedded along the requested xyz axis. It is the full orthogonal projection onto polynomials of that coordinate of degree at most `d`.

The value at coordinate zero is an a.e. convention. The step is zero there. A negative absolute power is assigned the finite value zero there as a representative of its L2 class. The set in question is the equator and has sphere-area measure zero. At exponent zero, the profile is instead 1 everywhere, including the equator. The certificate records the resulting value of `q`, including its offset. No uniform bound is asserted for a singular absolute-power source.

## Outward square root

For an exact nonnegative squared error `s` and `D=2^bits`, compute `k=floor(sqrt(floor(s*D^2)))`. Increment `k` if `k^2/D^2<s`. Then `k/D` is the smallest nonnegative multiple of `1/D` whose square is at least `s`. Integer square root and exact rational comparison suffice; no floating-point rounding enters the bound. `bits` may be any integer from 0 through 256.

Thus `error_squared` is the exact probability-L2 error squared, while `error_upper` is an outward rational bound on its square root. The area-measure squared L2 error is `4π*error_squared`, and the area-measure norm is `sqrt(4π)` times the probability norm. An `Linf` bound is unchanged by the choice of measure.

## Independent exact reference case

For `q=-|z|^(-1/4)`, the source mean is `-4/3` and its probability squared norm is 2. The exact projection errors are:

| Degree | Exact squared probability-L2 error |
| --- | --- |
| 4 | `5202/43681` |
| 8 | `14450/159201` |
| 12 | `4118450/54066609` |
| 24 | `1335263295935522/24084192784943929` |

For each of these projections, `||p_d+4/3||_2^2 = 2/9-error_squared`. This is a function-norm identity usable by the separate fixed-reference bridge; it is not a spectral assertion. The degree-24 polynomial has a largest coefficient of approximately `1.674e7`, beyond the frozen wide-potential backend's `1e6` input-coefficient budget. The front end intentionally supports this valid function certificate. A downstream backend must report its own input-budget refusal rather than changing the source or claiming a spectral result.

## Replay, tests, and limits

`verify_model` rebuilds the entire certificate from its normalized original function and requested order/degree/bits, and compares the full canonical payload. Optional expected-function and expected-polynomial arguments bind it to a caller's requested problem. Wrong axes, functions, exponents, amplitudes, offsets, norms, polynomials, error values, deleted proof rows, added fields, and inexact inputs are rejected. These are executable rational derivations using the stated analytic facts, not formal proof-assistant objects; the certificates explicitly state `formal_assistant_checked=False` and `spectral_transfer_claimed=False`.

The self-tests in `test_function_models.py` cover exact analytic coefficients and full tails, signed combinations, admissibility boundaries, axis binding, independent Rodrigues coefficients, direct polynomial norm integration, all four negative-quarter reference errors, degree-24 backend separation, square-root enclosure, and payload tampering. A release review found that the first draft combined logarithms before checking their domain; the implementation now checks every raw logarithm first, and a regression covers both cancellation and zero multiplication for `log1p(2z)` and `log1p(-1)`, as well as accepted cancellation for `log1p(z/2)`. The corrected suite contains 29 tests. The separate reviewer maintains its own test results and saved-evidence replay; root maintains the full integrated suite.

The corrected self-test run `python3 -B -m unittest test_function_models -v` passed all 29 tests. The final implementation SHA-256 is `7ed1540035eec59da95bca7a7dd59172568e124446fc93ca0ba6bf825ec4b7f3`.

The module is a bounded exact source model, with explicit degree and arithmetic-work budgets. It does not claim to accept all analytic or L2 functions, to compute optimal error bounds, or to establish any global sharp Gagliardo–Nirenberg constant.
