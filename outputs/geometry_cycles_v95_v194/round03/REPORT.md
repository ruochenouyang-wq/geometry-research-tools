# Round 03: V115–V124 broad-band and exponential potentials

The frozen baseline rejected `q=t^12` before solving it. This round provides a separate degree-24 polynomial path and an original-problem exponential-potential certificate. The new verifier rebuilds the polynomial, range proof, associated masses, finite form, complete radial coupling and angular tail. It does not pass high-degree certificates through the V94 degree-six validator.

## Actual callable increments

[STAGES.json](STAGES.json) lists the ten callable increments and their concrete evidence files. `wide_potential.py` exposes cheap full-sphere ground bounds (V115) and broad-band multiplication/assembly (V116), Sturm and bounded-cost Bernstein range evidence (V117–118), input-bound evidence replay (V119), complete Schur kernel evidence (V120), sector and full-sphere spectral certification (V121–122), rational exponential approximation (V123), and uniform minmax perturbation to the original exponential problem (V124).

Polynomial inputs have ascending rational coefficients, degree at most 24, coefficient absolute value at most 10^6 and denominator at most 2^2048. Integers, Fractions and canonical rational strings are supported; floats and booleans are rejected. The denominator allowance includes ordinary rational Taylor parameters through degree 24, such as a=1/10. Exceptional rational parameters whose Taylor coefficients exceed this explicit input budget are rejected, rather than rounded. No string is executed.

Spectral inputs allow integer m in 0..64, 1<=k<=modes<=64 and 8<=bits<=160. The experiments below use only 8–12 modes. Full-sphere `max_m` and `max_modes` are finite search budgets; exhausting either can produce a valid open enclosure. A target status is determined from the final exact rational width, not from the requested bisection count.

## Proof structure and scope

V115 now computes a standalone full-sphere ground enclosure, not merely wider input acceptance. `cheap_ground_enclosure` uses `(2 if mean_zero else 0)+qlo` for the lower bound. For the upper bound it integrates the original polynomial exactly against the real spherical harmonics z and x, and also 1 when mean-zero is not required, taking the smallest Rayleigh quotient. The radial squared factors are t^2, 1-t^2 and 1, with respective radial masses 2/3, 4/3 and 2; each function's azimuthal factor cancels within its quotient. The upper bound therefore comes from genuine full-sphere admissible functions. It requires no finite matrix, inertia search or numerical quadrature. Its separate verifier reconstructs the original moments, trial admissibility, range evidence and both endpoints without invoking the producer or a spectral solver. For q=t^24, the certified mean-zero bound is [2,451/225]; without the constraint it is [0,1/25]. A constant potential gives an exact enclosure. Bound and scope mutations are rejected.

The complete associated-Legendre expansion applies every power of t before projecting out a basis degree. In particular, removal of the constant mode does not remove intermediate multiplication paths through that mode. The nonorthonormal basis has exact mass `2(l+m)!/((2l+1)(l-m)!)`. All d possible omitted radial degrees are included for a degree-d potential, including their masses and complete columns into the retained window.

For a verified whole-domain range qlo<=q<=qhi, the omitted radial operator lies between the diagonal Laplacian plus qlo and plus qhi. Below its lower spectral bound, inversion reverses this form ordering. Exact Schur-complement inertia therefore provides the lower and upper endpoint tests. Search uses rational banded inertia; certificate replay uses the separate dense rational inertia implementation. The new certificate contains the entire kernel evidence, which replay reconstructs and compares. This is executable checking under the classical form-ordering, Schur-complement and minmax theorems, not a formal-assistant proof.

A sector certificate describes one azimuthal sector only. `full_ground(..., mean_zero=True)` combines consecutive sectors with the rigorous remaining angular bound `(next_m)*(next_m+1)+qlo`. Its scope is all zero-mean H1 functions on the unit two-sphere. The mean-zero m=0 operator is a compression to a fixed function space; multiplication by a nonconstant potential is not assumed to preserve zero mean.

The Sturm range encloses maxima of both q and -q and verifies the full closed domain [-1,1]. Bernstein cells use exact convex-hull bounds on a complete dyadic partition. `potential_range(strategy='auto')` chooses degree>8 Bernstein, otherwise Sturm: a deterministic degree-based cost heuristic, with no claim that its cheap bound is always equally tight. Explicit strategy selection is available. Reusing a range proof performs evidence replay and exact coefficient arithmetic, but no new root isolation or extrema search; a test disables both range generation and `maximize` during replay. Full-sphere construction searches for one range proof and reuses it across every sector.

For exp(a*t), |a|<=3, the degree-n Taylor polynomial has rational remainder bound
`delta = 3^ceil(|a|) * |a|^(n+1)/(n+1)!` on [-1,1], n<=24. The certificate explicitly declares Taylor's theorem, exponential monotonicity/multiplication, and e<3; the latter follows from the factorial-series tail dominated strictly by the geometric series with sum 1. Replay recomputes delta from a and n and does not trust a client-supplied error. The full polynomial spectral interval [L,U] becomes [L-delta,U+delta] for the original exponential potential by minmax on the same fixed mean-zero space. The final certificate embeds and binds both proofs and the original a.

## Bounded experiment results

All displayed endpoints are outward decimal enclosures; exact rationals are preserved in the linked evidence directory.

| Potential | Full mean-zero eigenvalue enclosure | Target 1e-8 |
|---|---|---|
| t^12 | [2.015146253687, 2.015146374703] | certified open |
| (1+t)^12 / 4096 | [2.027631817601, 2.027631878853] | certified open |
| exp(t) | [3.048191488429, 3.048191589548] | certified open |
| exp(-t) | [3.048191488429, 3.048191589548] | certified open |

The exponential reflection intervals are identical as exact rationals, not merely matching decimal displays. A degree-two q=2t^2 result intersects the frozen V94 full-sphere certificate. A two-mode, max_m=0 run remains a valid open bound with an unresolved angular contribution. The degree-24 kernel evidence explicitly contains all 24 omitted coupling columns. Tampering the Taylor delta to zero is rejected.

Recorded generation time for each displayed problem was approximately 0.36, 1.03, 1.61 and 1.61 seconds in this local run. These are small bounded examples, not general runtime guarantees. `results/certs/evaluation.json` contains the raw enclosures, timings and evidence paths. No model API, network call, token claim or expensive historical benchmark is involved.

## Reproduction and limitations

From `outputs/geometry_cycles_v95_v194`:

```sh
python3 -B -m unittest test_wide_potential -q
```

The suite has 56 tests: input budgets, cheap full-sphere Rayleigh and range bounds for both projections, independent low-degree and associated-mass moments, high-degree full multiplication, range coverage/binding/reuse, all-tail evidence, scope separation, altered masses/columns/bounds, old degree-two overlap, finite-budget openness, Taylor and original-potential binding, and reflection agreement. Existing evidence can be replayed with `verify_cheap_ground`, `verify_range`, `verify_kernel`, `verify_sector`, `verify_full`, `verify_approximation` or `verify_exponential`, according to its format. Input/assembly records and the evaluation/tamper diagnostic are records rather than spectral certificates.

This round expands the supported potential class. It does not prove universal Gagliardo–Nirenberg sharpness, provide explicit optimizer functions, guarantee a requested tolerance within finite budgets, support arbitrary nonpolynomial functions, or remove classical analytic assumptions. High degree, large coefficients and wide ranges can make the rational Schur work expensive or the enclosure loose. Improving that precision is separate work.
