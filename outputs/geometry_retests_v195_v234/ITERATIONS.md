# V195–V234：40项实际能力增量

每项链接到可调用数学能力及已保存证据；实现和独立审查已完成。

| 版本 | 原题 | 能力 | 调用 | 证据 |
|---|---|---|---|---|
| V195 | C12 | Validate all global moments of a real trial and a justified Poincare lower bound | `global_constraints.trial_enclosure` | [记录](C12/stages/V195.json) |
| V196 | C12 | Assemble every real harmonic including constant and sine/cosine modes | `global_constraints.full_harmonics` | [记录](C12/stages/V196.json) |
| V197 | C12 | Mass-weighted nullspace of complete cross-m and cross-component moments | `global_constraints.constraint_nullspace` | [记录](C12/stages/V197.json) |
| V198 | C12 | Congruence of the full coupled energy and non-diagonal mass matrix | `global_constraints.reduced_forms` | [记录](C12/stages/V198.json) |
| V199 | C12 | Transform every finite-to-tail harmonic coupling through the global nullspace | `global_constraints.transformed_tail` | [记录](C12/stages/V199.json) |
| V200 | C12 | Full-sphere constrained Schur spectral enclosure with entire infinite tail | `global_constraints.certify` | [记录](C12/stages/V200.json) |
| V201 | C12 | Snap to an exact Laplace level only when both exact inertia inequalities hold | `global_constraints.exact_level` | [记录](C12/stages/V201.json) |
| V202 | C12 | Transfer verified spectra between identical global constraint row spaces | `global_constraints.equivalent_constraints` | [记录](C12/stages/V202.json) |
| V203 | C12 | Increase the complete harmonic head while preserving all constraint support | `global_constraints.adaptive` | [记录](C12/stages/V203.json) |
| V204 | C12 | Prove a global weighted inequality or provide an explicit rational violating function | `global_constraints.inequality` | [记录](C12/stages/V204.json) |
| V205 | C14 | Global analytic ordered bounds include l=0 with multiplicity one in the unprojected space. | `general_spectrum.analytic_bounds` | [记录](C14/v205.json) |
| V206 | C14 | Distinct radial enumeration for actual full and mean-zero operators. | `general_spectrum.enumerate_sector` | [记录](C14/v206.json)、[记录](C14/v206_projected.json) |
| V207 | C14 | Full kth eigenvalue merge with uncomputed angular and radial min-max intervals. | `general_spectrum.ordered_merge` | [记录](C14/v207.json) |
| V208 | C14 | Strict threshold counting with the constant mode and exact threshold ties. | `general_spectrum.strict_count` | [记录](C14/v208.json) |
| V209 | C14 | Complete negative trace including every potentially negative omitted mode. | `general_spectrum.negative_trace` | [记录](C14/v209.json) |
| V210 | C14 | Trace-aware adaptive scheduling continues after requested ordered eigenvalues close. | `general_spectrum.adaptive_spectrum` | [记录](C14/v210.json)、[记录](C14/trace_schedule.json) |
| V211 | C14 | Projection-preserving threshold shift computes both count and sum(T-lambda)+. | `general_spectrum.spectral_shift` | [记录](C14/v211.json) |
| V212 | C14 | Scope-bound fixed-V trace quotient with standard-area factor 4 and constant-mode obstruction. | `general_spectrum.fixed_potential_quotient` | [记录](C14/v212.json) |
| V213 | C14 | Verified codimension-one interlacing for the same nonconstant potential. | `general_spectrum.projection_interlacing` | [记录](C14/v213.json) |
| V214 | C14 | One driver gives separate requested-quantity decisions and verified input-bound cache reuse. | `general_spectrum.research_driver` | [记录](C14/v214.json)、[记录](C14/cache_execution.json) |
| V215 | C16 | Sparse powers of y and exact probability mean kernel | `compact_gn.mean_certificate` | [记录](C16/kernel_evidence.json) |
| V216 | C16 | Centered sparse basis and explicit nonzero scale representative | `compact_gn.centered_function` | [记录](C16/kernel_evidence.json) |
| V217 | C16 | Closed pair mass Gram kernel for centered powers | `compact_gn.mass_form` | [记录](C16/kernel_evidence.json) |
| V218 | C16 | Closed exact sphere Dirichlet pair kernel including constant handling | `compact_gn.energy_form` | [记录](C16/kernel_evidence.json) |
| V219 | C16 | Sparse exponent-sum quartic convolution and exact moments | `compact_gn.quartic_moment` | [记录](C16/kernel_evidence.json)、[记录](C16/benchmark.json) |
| V220 | C16 | Probability GN ratio and original 2^n moment scaling | `compact_gn.ratio` | [记录](C16/original_cases.json) |
| V221 | C16 | Complete sparse Euler-Lagrange residual and exact improvement derivative | `compact_gn.residual` | [记录](C16/residuals.json) |
| V222 | C16 | Global two-chart projective-plane maximum with compressed moment kernels | `compact_gn.plane` | [记录](C16/planes.json) |
| V223 | C16 | Symbolic rational family, strict integer monotonicity and exact asymptotic limit | `compact_gn.family_law` | [记录](C16/family.json) |
| V224 | C16 | Unified evaluator and verified original monomial bridge for n4/16/64 | `compact_gn.bridge_original` | [记录](C16/original_cases.json)、[记录](C16/baseline.json) |
| V225 | C18 | Exact simplex Bernstein potential range on the sphere | `parity_spectrum.potential_range / verify_range` | [记录](C18/v225_evidence.json) |
| V226 | C18 | Exact sign-flip symmetry and eight character projections | `parity_spectrum.reflection_symmetry / parity_project` | [记录](C18/v226_evidence.json) |
| V227 | C18 | Complete character harmonic bases and exact positive masses | `parity_spectrum.degree_basis / parity_basis` | [记录](C18/v227_evidence.json) |
| V228 | C18 | Coupled forms within a character and exact cross-block checks | `parity_spectrum.block_form / cross_block_form` | [记录](C18/v228_evidence.json) |
| V229 | C18 | Complete same-character tail columns and correct first omitted degree | `parity_spectrum.tail_couplings / next_omitted_degree` | [记录](C18/v229_evidence.json) |
| V230 | C18 | Infinite-tail Schur enclosure for a reflection character | `parity_spectrum.block_ground / verify_block / analytic_floor` | [记录](C18/v230_evidence.json) |
| V231 | C18 | All-eight-character global ground merge with valid pruning | `parity_spectrum.full_ground / verify_full` | [记录](C18/v231_evidence.json) |
| V232 | C18 | Measured precision-versus-tail diagnosis and targeted same-parity refinement | `parity_spectrum.error_diagnosis` | [记录](C18/v232_evidence.json) |
| V233 | C18 | Input-bound exact finite-form and potential-range reuse | `parity_spectrum.FormCache / block_ground(cache=...)` | [记录](C18/v233_evidence.json) |
| V234 | C18 | Requested-width full-space driver with competitor-only refinement and honest budgets | `parity_spectrum.adaptive_ground / verify_result` | [记录](C18/v234_evidence.json) |
