Independent mathematical and certificate review completed on 2026-09-11.

Result: **15 independent tests passed; all saved global certificates, threshold certificates and ten stage records replayed successfully. No unresolved mathematical blocker found.** The implementation was not edited by the reviewer. Combined implementation and review tests for rounds 4 and 5 passed: 74 tests in 2.571 seconds. Timings are from one local correctness run, without claims about performance or model capability.

Reviewed implementation SHA-256:

`global_parameter.py: 384cbc2748c2e9628d1af3db0c8cf6a3eb1221c9269155cdddb9acfcb8ceef02`

Independent test SHA-256:

`test_global_review.py: 367dc3d77b2267d4ae1e329e41f24d27264ed7985e627ffcd6a2a324e97ef41f`

Run from `outputs/geometry_cycles_v95_v194`:

```text
python3 -B -m unittest -v test_global_review
Ran 15 tests in 1.260s
OK

python3 -B -m unittest -v test_parameter_family test_family_review test_global_parameter test_global_review
Ran 74 tests in 2.571s
OK
```

The objective is the minimum over a **closed** parameter interval of `C(q0+s*direction)+penalty(s)`, where C uses the full real mean-zero sphere function space. The potential degree is at most 6 and the scalar penalty degree at most 4. Every cell lower bound is valid throughout that cell. The global lower is the minimum over a complete interval cover, optionally strengthened by a valid older global lower. The global upper is a certified upper at an actual feasible parameter. Sampling alone supplies the latter, never the former.

The analytical baseline was reproduced independently: for zero potential and penalty `s^4-s²`, the objective is

`2+s^4-s² = 7/4+(s²-1/2)²`.

Its exact minimum is `7/4` at both `s=+sqrt(1/2)` and `s=-sqrt(1/2)`, although the three samples `-1,0,1` all give 2. At one-cell budget the implementation correctly returns `[7/4,2]` with an open gap. A mutation asserting lower 2 and gap zero is rejected. Exact rational square comparisons confirm that the saved minimizer enclosure retains both irrational minimizers; it makes no uniqueness claim.

The concavity direction was checked separately from that quartic. Taking the constant-shift operator `C(s)=2+s` and penalty `(s-1/4)²` gives the exact objective `2+(s+1/4)²`. The whole-cell concavity construction recovers lower 2, and bounded subdivision finds the actual rational minimizer `s=-1/4` with upper 2. It correctly minimizes the sum of a lower spectral chord and the penalty, rather than a chord of the possibly nonconcave full objective.

| Version | Callable capability checked | Independent mathematical evidence and result |
|---|---|---|
| V135 | `cheap_objective_enclosure` / `verify_cheap_objective` | `q=2t²+s*t`, penalty `s²/10`, `s in [-2,2]` gives global-minimum enclosure `[0,12/5]`. The upper is the actual midpoint coordinate trial `x1`; it is independently matched to the associated-Legendre Rayleigh integral. Spectral and extremum searches are prohibited during the test. |
| V136 | `penalty_enclosure` | The quartic has exact minimum `-1/4`. Equal numerical extrema on symmetric but different intervals do not permit reusing a wrong-interval certificate; changed-polynomial evidence is also rejected. |
| V137 | `lipschitz_cell` | For `C(s)=2+s`, penalty `s²`, center anchor over `[-1,1]` gives honest coarse lower 1, below the true minimum `7/4`. Wrong anchor coordinates are rejected. |
| V138 | `feasible_point` | The closed boundary `s=-1` in the constant-shift objective realizes exact global minimum 1. Saved sphere midpoint evidence is replayed against its original operator and penalty. |
| V139 | `branch_and_bound` | The quartic sampling trap remains an open gap at insufficient budget. Missing, duplicate and width `10^-15` gap mutations in a covering partition are rejected. |
| V140 | `concavity_cell` | The exact quadratic completion gives lower 2; swapped endpoint labels fail input binding. Saved chord-polynomial proof is replayed against its original interval. |
| V141 | `cached_point` | Reuse across two objectives at an identical actual point value is accepted. A changed penalty at the same nonzero point, or an unconstrained spectral proof, is rejected. |
| V142 | `minimizer_enclosure` | Both irrational quartic minimizers remain enclosed. A flat objective retains the entire domain when lower equals upper; equality is not used for strict exclusion. False uniqueness or missing-minimizer claims fail verification. |
| V143 | `resume` | Resume preserves every old point record, does not lower the prior global lower and does not worsen the incumbent upper. Wrong objective and digest mutations are rejected. |
| V144 | `threshold_certificate` | Threshold equal to a certified lower is proved; threshold equal to a strictly larger upper remains undetermined; threshold strictly above that upper is a spectral-existence refutation. The explicit-function flag remains false and cannot be relabelled true. |

The new V135 is a real computational capability. Its coefficient and monomial interval bounds cover every real parameter, while the upper uses a single admissible function and parameter. In particular, its upper bounds the **parameter minimum**, and need not bound the value at every parameter. This differs from the uniform per-parameter upper returned by round 4's `cheap_family_enclosure`. A nonsymmetric negative parameter interval and negative even polynomial coefficients were independently calibrated: at `s=-2`, `q=5+2t-5t²` with penalty 12, the coordinate `x3` trial has objective 16; conservative interval lower is `-17`. The negative result is retained as a coarse but valid enclosure.

The cover is recomputed exactly from sorted, contained, nondegenerate cells with matching endpoints and the original domain endpoints. No open-domain option is implemented. A minimum attained at the closed left boundary is deliberately retained, and attempted open-domain or changed-domain records are rejected. Cell exclusion uses strict `cell_lower > feasible_upper`; this is necessary to avoid deleting a tied minimizer. The returned set encloses every minimizer, not necessarily only minimizers, and a wide enclosure is not evidence of uniqueness.

Saved replay covered the nine top-level full global records:

- `global_boundary.json`, `global_quartic.json`, `global_sphere.json`;
- `open_budget.json`, `resumed.json`;
- `stage_139.json`, `stage_143.json`;
- `transition_alpha_1_10.json`, `transition_alpha_1_100.json`.

It also covered the three threshold certificates in `stage_144.json`, the full global certificate embedded in `stage_142.json`, and every individual V135–V144 stage's evidence with its original objective, polynomial and interval bindings. `projected.full_ground` and `exact_extrema.maximize` were disabled during replay, and `branch_and_bound` was disabled for full-global/threshold verification. V141's persisted cache count equality and reused point were checked; the original historical call count itself is not reconstructed from a trace.

No false mathematical claim or unresolved implementation defect was found in this review. Finite budgets, coarse interval arithmetic and nonunique minimizers remain explicit limitations. The review depends on the sphere min–max principle, bounded-multiplier estimate and concavity of the infimum of affine Rayleigh forms; it is an independent test and rational-certificate audit, not a proof-assistant formalization or a new theorem claim.
