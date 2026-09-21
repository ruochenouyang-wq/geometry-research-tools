# Structure-adapted function approximation

This is one scheme: exact local probability-L2 polynomial projection on a geometric partition, with exact reuse of the shape of each annulus. The two supported original sources are `step(t)=1_(t>0)` and `|t|^(-1/4)`, including rational signed amplitude and constant offset, along any of the three unit-sphere coordinates. The target singular example is `q=-|t|^(-1/4)`.

The step is represented exactly by two constants separated at zero. For the singular source, the same 25 nominal coefficient slots that a degree-24 global polynomial uses give a certified function-error bound of approximately `0.02402252489`, compared with the actually replayed frozen global projection's `0.23546014554`. A predeclared larger-budget trial reaches the rigorous bound `1833/1099511627776 < 1e-8`. These are function-norm results. This module supplies neither a spectral enclosure nor a global polynomial accepted by the old spectral backend.

## Interface and exact source scope

```python
source = {"kind": "axis_profile", "axis": 2, "profile": "abs_power",
          "exponent": "-1/4", "amplitude": "-1", "offset": "0"}
certificate = piecewise_model(source, levels=6, degree=3, root_ratio="1/2", sqrt_bits=40)
verify_piecewise(certificate, expected_function=source)
evaluate(certificate, "1/16")
```

`normalize_function` shares the established source defaults: axis 2, amplitude 1, offset 0, and absolute-power exponent 1. Since this implementation supports only exponent `-1/4`, omitting the exponent is rejected as unsupported. It is never silently interpreted as the singular source. The step does not accept an exponent field. Other profiles, extra fields, floating-point numbers, booleans used as numbers, and mismatched expected sources are rejected.

The arithmetic limits are degree 0–32, 1–64 annuli (`levels`), square-root precision 0–256 bits, and rational fourth-root grid ratio `0<r<1` with denominator at most 4096. These are front-end work budgets, independent of the old degree-24 or coefficient-size restrictions of a spectral engine. A valid certificate makes no assertion that a requested error target has been met unless its exact `error_upper` does so.

The certificate scope is `unit_S2_probability_measure`, and the norm is `L2_probability`. For `dμ=dσ/(4π)`, a coordinate has law `dt/2`; rotation invariance handles every axis. An even residual therefore has squared sphere-probability norm equal to its integral on `[0,1]`. Area-L2 squared error is `4π` times the displayed squared error. No area/probability factor is hidden in the comparison.

At axis zero the original step is defined to be zero, and the singular profile is assigned zero as an a.e. representative. This changes only the equator, a set of sphere measure zero. The piecewise approximant at zero need not equal that source representative. The evaluator handles the step's zero value separately; nonzero seams use the cell beginning at that absolute coordinate, with the final endpoint included.

## Exact cell proof

Write `q(t)=b+A*t^(-1/4)` for `t>0`. Every cell has endpoints `a=L^4` and `c=R^4` for rational nonnegative `L<R`. The substitution `t=x^4` yields, for every nonnegative integer `k`,

`integral_a^c t^k*t^(-1/4) dt = 4*(R^(4k+3)-L^(4k+3))/(4k+3)`.

Also `integral_a^c t^(-1/2) dt = 2*(R^2-L^2)`. Thus all source moments, including the norm of the singular part, are rational. The complete squared source norm on a cell is

`b^2*(c-a) + 2*b*A*(4/3)*(R^3-L^3) + 2*A^2*(R^2-L^2)`.

The signed amplitude/offset cross term is retained. Let `s=(t-a)/(c-a)` and let `L_l(s)=P_l(2s-1)` be the shifted Legendre polynomial. Its cell norm squared is `(c-a)/(2l+1)`. Exact source moments give `h_l=integral q*L_l dt` and projection coefficients `v_l=(2l+1)*h_l/(c-a)`. The returned local polynomial is `sum v_l*L_l(s)`.

For the base annulus, the implementation separately expands the actual returned polynomial into powers of `t`, integrates its square and its cross product with the source, and checks both against the Legendre norm sum. It records the full error

`E_cell = integral q^2 dt - 2*integral q*p dt + integral p^2 dt`.

No sample points, numerical quadrature, or truncated residual tail are used. The verifier reconstructs the source, cover, basis, coefficients, every moment and norm, and the outward square root. Whole-payload comparison rejects missing core cells, gaps, altered fourth roots, altered polynomials or source parameters, changed norms, missing projection rows, and incorrect error claims.

## Geometric cover, core, and shape reuse

The positive half has a core `[0,r^(4J)]` and `J` annuli `[r^(4(j+1)),r^(4j)]` for `j=0,...,J-1`, reflected evenly to the negative half. The core uses its exact mean as a constant approximation; each annulus uses degree `d`. The expanded coefficient-slot count is `1+J*(d+1)`. These are stored local polynomial coefficients, not independent optimization variables: every annulus follows from the same shape containing `d+1` coefficients, plus the single core constant fixed by its exact mean. The certificate distinguishes `expanded_polynomial_coefficient_slots` from `shape_coefficient_count`. The step has two constants. Knots are fixed by `r,J` and are not free optimization parameters.

For the unit-amplitude singular power the core squared error is exactly `(2/9)*r^(2J)`, so it is never discarded as an insignificant interval. Its constant approximant is `(4/3)*r^(-J)` before applying the source amplitude and offset.

Only one degree-`d` shape on `[r^4,1]` is integrated from scratch. For another annulus set `t=R^4*x`, where `R=r^j`. Then `t^(-1/4)=R^(-1)*x^(-1/4)` and `dt=R^4 dx`. Consequently the nonconstant projection coefficients scale by `A/R`, the constant coefficient additionally receives the offset, and squared residual integrals scale by `A^2*R^2`. The implementation checks the resulting source norm minus projection norm against this scaled exact residual in every cell. Independent tests also compare transported cells with separate direct integrations.

If `e_d(r)` is the base annulus's unit-amplitude squared residual, the complete error is

`A^2 * [(2/9)*r^(2J) + e_d(r)*(1-r^(2J))/(1-r^2)]`.

This identity is checked against the sum of all returned cells. It explains both refinement and saturation: at fixed `r,d`, merely increasing `J` approaches the positive squared-error floor `A^2*e_d(r)/(1-r^2)`. Higher local degree or a narrower annulus can reduce that floor. This is why the report retains large-budget trials that do not reach a target instead of treating core refinement as arbitrary-precision convergence.

The outward square root uses integer square root on a dyadic grid and an exact comparison. `error_upper` is the smallest multiple of `2^(-sqrt_bits)` whose square encloses the exact squared error. It is an upper bound on the whole-domain error, including the core.

## Fixed-budget and precision experiments

The frozen `geometry_general_v235_v237/function_models.py` was actually run and verified at degree 24 on the identical normalized originals. For the singular source its exact squared error is `1335263295935522/24084192784943929`, and its norm upper bound is `258891167901/1099511627776`. It has 25 nominal global coefficient slots, but only 13 nonzero coefficients because the source is even. Accordingly both 25-slot and 13-slot piecewise comparisons are shown. They compare different approximation spaces, not a universal cost-equivalence theorem.

For the step, that same frozen degree-24 implementation returns exact squared error `457028729521/70368744177664` and upper bound `676039/8388608`, approximately `0.08059012890`. The two-piece model has exactly zero error with two constant values. The original source and the a.e. convention are unchanged.

All rows below use the original singular source, and the displayed decimals are only readable summaries of rational certificates. The fixed grid is `r=1/2` unless stated otherwise.

| Annuli J | Local degree d | Expanded coefficient slots | Certified L2 upper bound, approximately |
| --- | --- | --- | --- |
| 1 | 1 | 3 | 0.25111893713 |
| 2 | 2 | 7 | 0.12604986678 |
| 4 | 2 | 13 | 0.05470571802 |
| 4 | 3 | 17 | 0.03726886990 |
| 5 | 3 | 21 | 0.02719299716 |
| 6 | 3 | 25 | 0.02402252489 |
| 8 | 4 | 41 | 0.01192917017 |
| 12 | 6 | 85 | 0.00336161010 |
| 16 | 8 | 145 | 0.00101305584 |
| 24 | 8 | 217 | 0.00101303031 |

The smallest row is worse than the frozen reference and is retained. At the same 25 slots (`J=6,d=3`), the three preselected grid ratios `1/2,3/4,7/8` give respectively `0.02402252489`, `0.08390591652`, and `0.21156414703`. At 217 slots (`J=24,d=8`), they give respectively `0.00101303031`, `0.00047300422`, and `0.01912436464`. None of these larger rows reaches `1e-6`.

Two further bounded experiments were declared before execution, retaining the unsuccessful one:

| Grid ratio | J | d | Slots | Exact error upper | Target `1e-8` |
| --- | --- | --- | --- | --- | --- |
| `1/2` | 32 | 24 | 801 | `75073/549755813888` ≈ `1.36557e-7` | Not reached |
| `2/3` | 48 | 20 | 1009 | `1833/1099511627776` ≈ `1.66710e-9` | Reached |

The first certificate serialized to 754819 bytes and took about 0.050 s to build and 0.054 s to replay in the bounded local run. The second serialized to 1669967 bytes and took about 0.097 s to build and 0.105 s to replay. These are measured arithmetic wall times for one run and actual serialized bytes, not language-model token counts or a general speed claim. Timings precede the final metadata naming clarification; byte counts bind the final format. The complete certificate repeats cell-specific evidence, so its storage is larger than the reused shape alone. Root saves the full JSON certificates separately.

## Reproduction and limitations

Run `python3 -B -m unittest test_piecewise_models -v` from this directory. All 32 tests passed on the final implementation. The tests cover source defaults and binding, all axes, exact step values, core inclusion, probability normalization, independent binomial/Rodrigues computations, direct residual integration, exact shape transport, signed scaling and offset cross terms, the true frozen baseline, error refinement, the unsuccessful precision attempt, the successful bounded attempt, and altered evidence.

`python3 -B piecewise_models.py` prints a compact reproducible budget table. The saved recipes below bind selected source/parameter/error records to the SHA-256 of the complete canonical certificate. A test reconstructs and replays each recipe and checks the exact error, serialized bytes, and digest. These recipes are not substitutes for the separately saved full certificates; they provide a compact check that regeneration has not changed the evidence.

This model is piecewise and may jump at knots. It is not a globally smooth trial eigenfunction or a global polynomial potential. Its small L2 function error does not imply a comparably narrow spectral interval, and it is not passed to the frozen polynomial spectral backend as if it were one polynomial. Any downstream spectral transfer needs its own correctly scoped theorem and verifier. No sharp GN constant or open problem is claimed solved.

<!-- SAVED_RECIPES_BEGIN -->
```json
[
  {
    "certificate_bytes": 1101,
    "certificate_sha256": "bda52d2349c826108e927d41e91e8cfe25093d718d1b746e2b0d690d6b5352f8",
    "error_squared": "0",
    "error_upper": "0",
    "function": {
      "amplitude": "1",
      "axis": 2,
      "kind": "axis_profile",
      "offset": "0",
      "profile": "step"
    },
    "name": "step_exact",
    "parameters": {
      "degree": 3,
      "levels": 4,
      "root_ratio": "1/2",
      "sqrt_bits": 40
    }
  },
  {
    "certificate_bytes": 6483,
    "certificate_sha256": "5d408618d597f2ad6ddab8b435a9d42bf1f77f8c6eaac5c2ced4ab78a7c328a3",
    "error_squared": "2155873/720373500",
    "error_upper": "3759348317/68719476736",
    "function": {
      "amplitude": "-1",
      "axis": 2,
      "exponent": "-1/4",
      "kind": "axis_profile",
      "offset": "0",
      "profile": "abs_power"
    },
    "name": "matched_13_slots",
    "parameters": {
      "degree": 2,
      "levels": 4,
      "root_ratio": "1/2",
      "sqrt_bits": 40
    }
  },
  {
    "certificate_bytes": 10536,
    "certificate_sha256": "48008a0b06537cfafc2734b207af727d46c0c1cd0697d7d96da2361829b36c36",
    "error_squared": "1166159779/2020788000000",
    "error_upper": "13206522725/549755813888",
    "function": {
      "amplitude": "-1",
      "axis": 2,
      "exponent": "-1/4",
      "kind": "axis_profile",
      "offset": "0",
      "profile": "abs_power"
    },
    "name": "matched_25_slots",
    "parameters": {
      "degree": 3,
      "levels": 6,
      "root_ratio": "1/2",
      "sqrt_bits": 40
    }
  },
  {
    "certificate_bytes": 754819,
    "certificate_sha256": "ad6e3afe4e85f1ae3879c696e7dae553c87b370cd6318aaeaaf34e7efdddf98d",
    "error_squared": "1217613909675454313249284736205448680624044440486251645229081386356493462726473698995645489539/65295901168031228882524689510604937321254438927626933660672000000000000000000000000000000000000000000000000",
    "error_upper": "75073/549755813888",
    "function": {
      "amplitude": "-1",
      "axis": 2,
      "exponent": "-1/4",
      "kind": "axis_profile",
      "offset": "0",
      "profile": "abs_power"
    },
    "name": "precision_first_unmet",
    "parameters": {
      "degree": 24,
      "levels": 32,
      "root_ratio": "1/2",
      "sqrt_bits": 40
    }
  },
  {
    "certificate_bytes": 1669967,
    "certificate_sha256": "b33186d18d64676b5057840d5f13a0cc5298aa1fa3b156d2067d00a6ce32e525",
    "error_squared": "205443203169454653909126829059848725630476362211850548107196592229534351666047516985787620461037388707808751953047836552353098266/73983902566065940562876685923175268448173428914887783658928351084162472052248469377912555780192764799786451274086601233648252673447132110595703125",
    "error_upper": "1833/1099511627776",
    "function": {
      "amplitude": "-1",
      "axis": 2,
      "exponent": "-1/4",
      "kind": "axis_profile",
      "offset": "0",
      "profile": "abs_power"
    },
    "name": "precision_fallback_met",
    "parameters": {
      "degree": 20,
      "levels": 48,
      "root_ratio": "2/3",
      "sqrt_bits": 40
    }
  }
]
```
<!-- SAVED_RECIPES_END -->
