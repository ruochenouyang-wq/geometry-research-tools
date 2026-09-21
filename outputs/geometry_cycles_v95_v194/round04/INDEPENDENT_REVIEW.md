Independent mathematical and certificate review completed on 2026-09-11.

Result: **14 independent tests passed; all 19 saved family certificates replayed successfully. No unresolved mathematical blocker found.** The implementation was not edited by the reviewer. The complete implementation plus review suites for rounds 4 and 5 also passed: 74 tests in 2.571 seconds. These timings describe one local validation run, not a benchmark or a model comparison.

Reviewed implementation SHA-256:

`parameter_family.py: caa137fc6f9587d6ad4a80c1d8b5a8868e37e5b52189b87b98a2c564d4d03b17`

Independent test SHA-256:

`test_family_review.py: a145fe591017a5e94f998fe666f7a3691a9270c36aecd6e357db8bb0551e4ad7`

Run from `outputs/geometry_cycles_v95_v194`:

```text
python3 -B -m unittest -v test_family_review
Ran 14 tests in 0.104s
OK

python3 -B -m unittest -v test_parameter_family test_family_review test_global_parameter test_global_review
Ran 74 tests in 2.571s
OK
```

The object under review is the lowest Rayleigh value on the stated **full** sphere function space, for a fixed affine polynomial potential family. When `mean_zero=True` the unperturbed floor is 2; otherwise it is 0. All parameter boxes are closed, have rational endpoints and dimension 1–3, and cover every real parameter in the box. Nothing here certifies an open parameter domain, uniqueness of an optimizing parameter, or an unrestricted sharp functional-inequality constant.

For every admissible fixed function, its Rayleigh value is affine in the parameters. Taking an infimum over those functions makes the lowest spectral value concave. Consequently, the **minimum of all vertex lower bounds** is a uniform whole-box lower bound; the **minimum of vertex upper bounds** is an upper bound on the parameter minimum. The maximum of vertex values is only a lower bound on the parameter maximum. The review explicitly tests this distinction with `q(s,t)=s*t`: the center has exact constrained lowest value 2, while both endpoints have values strictly below 2.

The upper-envelope test is stronger than a schema check. With potential `q(s,r,t)=r*(t²-2/5)` and two actual coordinate trials, the affine Rayleigh functions are `2+r/5` and `2-r/5`. Their minimum has maximum 2 at `r=0`; its values at every box vertex are only `9/5`. V132 returns the safe, deliberately looser `min(max(line_i))=11/5`. A forged upper `9/5` is rejected. V131 computes the exact one-dimensional envelope maximum, including the interior intersection and parallel or duplicate lines.

| Version | Callable capability checked | Independent mathematical evidence and result |
|---|---|---|
| V125 | `cheap_family_enclosure` | Constant potential shifts give exact constrained/full-space floors. For `2t²+s*t`, `s in [-1,1]`, the whole-family enclosure is `[1,16/5]`, with the upper obtained from the actual trial `u=t`. Generation and verification run with spectral and extremum searches disabled. This is a substantive bound, not just a representation version. |
| V126 | `direction_norms` | `t-t³` has exact norm `2/(3 sqrt(3))`, certified between `9/25` and `2/5`, well below coefficient absolute sum 2; wrong interval evidence is rejected. |
| V127 | `anchor_cell` | Constant shift families produce exact endpoint range via the farthest-distance radius; same-point cache reuse is accepted only when the actual potential and function constraint agree. |
| V128 | `partition_certificate` | Three-dimensional nonuniform T-shaped partition with shared faces verifies. Removing a width `10^-15` slab or duplicating a cell is rejected exactly. |
| V129 | `cell_lower_bound` | Exact one- and two-parameter constant shifts calibrate minima. The `s*t` family rejects an added false vertex-maximum upper claim. |
| V130 | `affine_majorant` | Transfer of a physical `P1` trial originally certified for a different potential yields exact intercept 7 and slopes `[2,-3]`. This reuse is valid because the actual function, rather than its old spectral value, is transferred. |
| V131 | `envelope_1d` | Opposite slopes intersect inside the interval; exact maximum 2 exceeds both endpoint envelope values `9/5`. Duplicate lines remain valid. |
| V132 | `envelope_box` | Explicit two-parameter min/max-order counterexample above verifies the safe upper `11/5`, retaining the documented looseness. |
| V133 | `adaptive_family` | Flat family closes at `[2,2]`; a nonconstant shift family at one-cell budget honestly reports `budget_open`. Altering that status is rejected. |
| V134 | `universal_inequality` | For `q=s`, `s in [-1,1]`, threshold 1 is proved inclusively; `1001/1000` is refuted at `s=-1` by an actual Rayleigh trial with quotient 1. Wrong threshold, parameter, evidence type and subbox evidence are rejected. |

The coverage argument was inspected independently: a finite collection of closed subboxes contained in the parent, with pairwise disjoint interiors and exact total parent volume, covers the parent. Otherwise its complement relative to the parent contains a nonempty relatively open neighborhood of positive volume. Thus the exact volume and interior-overlap checks cover boundaries as well as interiors; they are not a sampling argument.

The tests also cover closed-domain endpoint acceptance and attempted open-domain relabelling, wrong projected/unconstrained cache reuse, fixed-function affine transfer, constant-spectrum equality cases, explicit counterexample binding, and a falsely closed optimization budget.

All 19 files in `round04/certificates/*.json` were independently loaded and verified while `projected.full_ground` and `exact_extrema.maximize` were patched to fail if invoked. Every V125–V134 stage has a callable API and an existing evidence path. The saved outcomes include coarse bounds and unresolved budget outcomes; the review does not convert those into claimed improvements.

One robustness defect was reported directly to the owner: V134 accepted any generically valid lower proof at its entry, then could access a missing `vertices` field on a nonvertex proof. It did not produce a false certificate because final assembly already enforced the correct kind. The owner changed the entry to require a whole-box vertex-concavity certificate and raise `ValueError`; the independent regression `test_universal_rejects_wrong_lower_evidence_kind_or_subbox` now passes. No implementation change was made by the reviewer.

These are independently checked rational computations and tests of certificate scope, with the usual sphere min–max, multiplier and concavity facts as mathematical dependencies. They are not a proof-assistant formalization, a proof of the entire search history, or a new theorem claim.
