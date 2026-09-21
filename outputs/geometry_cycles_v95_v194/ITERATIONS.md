# V95–V194：一百项可执行增量

每行对应一次有保存证据的能力改进；实现按十个模块组织。原题和弱点见各轮 task.json，最终验收见 COMPLETED.json。

| 版本 | 能力 | 可调用入口 | 证据 |
|---|---|---|---|
| V95 | 精确任意有理球谐Rayleigh证书 | `witness.rayleigh_certificate` | [运行证据](round01/results/V95.json) |
| V96 | 从完整谱证书选择候选方位扇区 | `witness.select_sector` | [运行证据](round01/results/V96.json) |
| V97 | 质量归一化最低本征向量提议 | `witness.mass_normalized_proposal` | [运行证据](round01/results/V97.json) |
| V98 | 逐级有理量化且保留最佳商 | `witness.quantize_proposal` | [运行证据](round01/results/V98.json) |
| V99 | 质量贡献排序的稀疏截断与精确接受 | `witness.sparsify` | [运行证据](round01/results/V99.json) |
| V100 | 完整投影残差与所有势乘法高阶 | `witness.residual_certificate` | [运行证据](round01/results/V100.json) |
| V101 | 残差驱动新增试探支持 | `witness.residual_enrich` | [运行证据](round01/results/V101.json) |
| V102 | 二维Ritz精确下界及实际试探上界 | `witness.two_trial_ritz` | [运行证据](round01/results/V102.json) |
| V103 | 负残差方向的驻点及有理线搜索 | `witness.residual_line_search` | [运行证据](round01/results/V103.json) |
| V104 | 自动构造显式球谐违反者 | `witness.find_counterexample` | [运行证据](round01/results/V104.json) |
| V105 | Choose bit, radial or angular work from a verified bound | `precision.diagnose_precision` | [运行证据](round02/results/v105.json) |
| V106 | Scale bisection bits to a rational absolute target | `precision.choose_bits` | [运行证据](round02/results/v106.json) |
| V107 | Use full projected variance and bound the same m second eigenvalue | `precision.temple_certificate` | [运行证据](round02/results/v107.json) |
| V108 | Refine a weak separator and keep Rayleigh only if the gap still fails | `precision.recover_gap` | [运行证据](round02/results/v108.json) |
| V109 | Partition by radial window, retaining every projected tail coefficient | `precision.residual_schedule` | [运行证据](round02/results/v109.json) |
| V110 | Increase dyadic coefficient precision and apply exact inverse steps | `precision.refine_coefficients` | [运行证据](round02/results/v110.json) |
| V111 | Intersect verified Schur, Temple and previously accepted sector bounds | `precision.intersect_bounds` | [运行证据](round02/results/v111.json) |
| V112 | Revalidate q/projection binding and inherit exact bounds monotonically | `precision.checkpoint` | [运行证据](round02/results/v112.json) |
| V113 | Prune sectors by exact lower bounds against the best certified upper | `precision.winning_sectors` | [运行证据](round02/results/v113.json) |
| V114 | Target absolute full-S2 projected ground width with honest finite-budget exits | `precision.full_ground` | [运行证据](round02/results/v114.json) |
| V115 | Cheap certified full-sphere ground enclosure for rational potentials through degree 24 | `wide_potential.cheap_ground_enclosure; wide_potential.verify_cheap_ground` | [运行证据](round03/results/certs/degree24_cheap_mean_zero.json) |
| V116 | Complete broad-band associated-Legendre multiplication and assembly | `wide_potential.q_times_basis; wide_potential.assemble` | [运行证据](round03/results/certs/associated_degree12_assembly.json) |
| V117 | Whole-domain Sturm polynomial range evidence | `wide_potential.sturm_range; wide_potential.verify_range` | [运行证据](round03/results/certs/degree12_sturm_range.json) |
| V118 | Bounded-cost Bernstein enclosure and explicit range strategy dispatch | `wide_potential.bernstein_range; wide_potential.potential_range` | [运行证据](round03/results/certs/asymmetric_degree12_bernstein_range.json) |
| V119 | Input-bound reusable range evidence without a new extrema search | `wide_potential.reuse_range; wide_potential.range_bounds` | [运行证据](round03/results/certs/reused_sturm_sector.json) |
| V120 | Wide-band Schur kernel with all tail columns and associated masses | `wide_potential.kernel_certificate; wide_potential.verify_kernel` | [运行证据](round03/results/certs/degree24_schur_kernel.json) |
| V121 | Certified high-degree azimuth-sector spectrum | `wide_potential.certify_sector; wide_potential.verify_sector` | [运行证据](round03/results/certs/degree12_sector.json) |
| V122 | Full mean-zero sphere combination including all angular tail sectors | `wide_potential.full_ground; wide_potential.verify_full` | [运行证据](round03/results/certs/degree12_full.json) |
| V123 | Uniform rational Taylor enclosure of exp(a*t), |a|<=3, degree<=24 | `wide_potential.exponential_approximation; wide_potential.verify_approximation` | [运行证据](round03/results/certs/exp_t_taylor.json) |
| V124 | Original exponential-potential spectrum via uniform minmax perturbation | `wide_potential.exponential_ground; wide_potential.perturb_certificate; wide_potential.verify_exponential` | [运行证据](round03/results/certs/exp_plus.json) |
| V125 | Compute a search-free whole-family spectral enclosure from coefficient boxes and a fixed admissible harmonic. | `parameter_family.cheap_family_enclosure` | [运行证据](round04/certificates/v125_family.json) |
| V126 | Certify multiplier L-infinity bounds by rational polynomial extrema. | `parameter_family.direction_norms` | [运行证据](round04/certificates/v126_norms.json) |
| V127 | Propagate a verified point spectrum to every point of a whole cell. | `parameter_family.anchor_cell` | [运行证据](round04/certificates/v127_anchor.json) |
| V128 | Prove exact rational box coverage with disjoint interiors. | `parameter_family.partition_certificate` | [运行证据](round04/certificates/v128_partition.json) |
| V129 | Use concavity to certify a uniform lower from every vertex. | `parameter_family.cell_lower_bound` | [运行证据](round04/certificates/v129_concavity.json) |
| V130 | Produce exact affine Rayleigh majorants from fixed explicit functions. | `parameter_family.affine_majorant` | [运行证据](round04/certificates/v130_center_majorant.json) |
| V131 | Certify the maximum of a one-dimensional minimum affine envelope. | `parameter_family.envelope_1d` | [运行证据](round04/certificates/v131_intersection_envelope.json) |
| V132 | Provide a safe multi-parameter cell upper from affine majorants. | `parameter_family.envelope_box` | [运行证据](round04/certificates/v132_two_parameter_upper.json) |
| V133 | Adaptively subdivide the global-supremum interval using reusable anchor proofs. | `parameter_family.adaptive_family` | [运行证据](round04/certificates/v133_budget_open.json) |
| V134 | Resolve a universally quantified weighted inequality or exhibit parameter and function. | `parameter_family.universal_inequality` | [运行证据](round04/certificates/v134_proved.json) |
| V135 | Cheap whole-family interval lower bound and feasible spherical-harmonic upper bound | `global_parameter.cheap_objective_enclosure` | [运行证据](round05/stage_135.json) |
| V136 | Nonconvex polynomial penalty interval extrema | `global_parameter.penalty_enclosure` | [运行证据](round05/stage_136.json) |
| V137 | Certified spectral Lipschitz lower bound on a complete interval | `global_parameter.lipschitz_cell` | [运行证据](round05/stage_137.json) |
| V138 | Feasible rational parameter with full-sphere spectral incumbent | `global_parameter.feasible_point` | [运行证据](round05/stage_138.json) |
| V139 | Best-lower-first complete-domain branch and bound | `global_parameter.branch_and_bound` | [运行证据](round05/stage_139.json) |
| V140 | Concavity-chord spectral bound plus globally minimized penalty | `global_parameter.concavity_cell` | [运行证据](round05/stage_140.json) |
| V141 | Validated cache and endpoint/midpoint reuse | `global_parameter.cached_point` | [运行证据](round05/stage_141.json) |
| V142 | Strict cell exclusion and enclosure of every minimizer | `global_parameter.minimizer_enclosure` | [运行证据](round05/stage_142.json) |
| V143 | Input-bound continuation with monotone global bounds | `global_parameter.resume` | [运行证据](round05/stage_143.json) |
| V144 | Uniform threshold proof/refutation/undetermined and full replay | `global_parameter.threshold_certificate` | [运行证据](round05/stage_144.json) |
| V145 | Analytic full-sphere harmonic-exclusion lower bound and actual harmonic Ritz upper | `constraints.analytic_exclusion_bound` | [运行证据](round06/evidence/v145.json) |
| V146 | Assemble constrained sectors with degree start max(m,cutoff+1) | `constraints.sector_assembly` | [运行证据](round06/evidence/v146.json) |
| V147 | Constraint-aware lower bounds for radial and omitted angular tails | `constraints.tail_floors` | [运行证据](round06/evidence/v147.json) |
| V148 | Exact Schur eigenvalue certificate for a low-harmonic constrained sector | `constraints.certify_exclusion_sector` | [运行证据](round06/evidence/v148.json) |
| V149 | Full sphere aggregation with complete harmonic exclusion | `constraints.full_exclusion_ground` | [运行证据](round06/evidence/v149.json) |
| V150 | Infinite-radial-space analytic variational bracket for arbitrary finite moments | `constraints.analytic_moment_bound` | [运行证据](round06/evidence/v150.json) |
| V151 | Mass-weighted exact constraint nullspace with redundancy removal | `constraints.mass_nullspace` | [运行证据](round06/evidence/v151.json) |
| V152 | Generalized reduced energy, non-diagonal mass, and transformed tail couplings | `constraints.constrained_forms` | [运行证据](round06/evidence/v152.json) |
| V153 | Infinite-dimensional arbitrary finite-moment single-sector spectrum | `constraints.certify_moment_sector` | [运行证据](round06/evidence/v153.json) |
| V154 | Three-way constrained weighted inequality with bound quantifier | `constraints.constrained_inequality` | [运行证据](round06/evidence/v154.json) |
| V155 | Global initial ordered min-max intervals with full multiplicity and no Schur computation. | `ordered_spectrum.initial_order_bounds` | [运行证据](round07/v155.json) |
| V156 | Exact unperturbed and constant-shift ordered-spectrum calibration. | `ordered_spectrum.laplace_spectrum` | [运行证据](round07/v156.json) |
| V157 | Verified sector eigenvalue assembly with 1/2 angular multiplicity. | `ordered_spectrum.assemble_sectors` | [运行证据](round07/v157.json) |
| V158 | Independent radial and azimuthal floors for unenumerated modes. | `ordered_spectrum.omitted_floors` | [运行证据](round07/v158.json) |
| V159 | Full kth order statistics with analytical missing-mode placeholders. | `ordered_spectrum.ordered_bounds` | [运行证据](round07/v159.json) |
| V160 | Adaptive sector/radial enumeration with honest resource exhaustion. | `ordered_spectrum.adaptive_spectrum` | [运行证据](round07/v160.json) |
| V161 | Strict full spectral counts and explicit threshold clusters. | `ordered_spectrum.count_below` | [运行证据](round07/v161.json) |
| V162 | Negative trace intervals include all potentially negative omitted modes. | `ordered_spectrum.negative_trace` | [运行证据](round07/v162.json) |
| V163 | Nonnegative-potential validation and rational pi-scaled fixed-potential LT quotient. | `ordered_spectrum.lt_quotient` | [运行证据](round07/v163.json) |
| V164 | Bound certificate bundle, independent replay and tri-state eigenvalue assertion. | `ordered_spectrum.bundle` | [运行证据](round07/v164.json) |
| V165 | Exact mass, energy and quartic moments through trial degree 36 | `gn_variation.moments` | [运行证据](round08/results.json) |
| V166 | Exact zero-mean centering and nonzero rational scale normalization | `gn_variation.normalize` | [运行证据](round08/results.json) |
| V167 | Bound trial to exact area-normalized pi*K GN quotient | `gn_variation.ratio` | [运行证据](round08/results.json) |
| V168 | Complete projected Euler-Lagrange residual including cubic high degrees | `gn_variation.residual` | [运行证据](round08/results.json) |
| V169 | Exact directional derivative and certified positive local direction | `gn_variation.direction` | [运行证据](round08/results.json) |
| V170 | Global 2D projective-plane GN maximum via two complete charts | `gn_variation.plane` | [运行证据](round08/high_degree_plane.json) |
| V171 | Finite-budget monotone accepted plane refinements | `gn_variation.refine` | [运行证据](round08/refinement.json) |
| V172 | Legendre degree-budget projection with full inside/outside residual evidence | `gn_variation.project_residual` | [运行证据](round08/results.json) |
| V173 | Centered-cap exact trial generation and adaptive target/budget screening | `gn_variation.screen_caps` | [运行证据](round08/caps.json) |
| V174 | Compile GN trial square to degree24 potential and replay full-spectrum diagnostics | `gn_variation.spectral_diagnostic` | [运行证据](round08/spectral_diagnostics.json) |
| V175 | Immediate full-space enclosure with exact linear trial upper bound | `anisotropic.cheap_ground_enclosure / verify_cheap` | [运行证据](round09/v175_evidence.json) |
| V176 | Exact normalized spherical moments | `anisotropic.sphere_moment / integrate` | [运行证据](round09/v176_evidence.json) |
| V177 | Rational real solid harmonic construction | `anisotropic.solid_harmonic` | [运行证据](round09/v177_evidence.json) |
| V178 | Exact Gram and tangential gradient forms | `anisotropic.gram_energy` | [运行证据](round09/v178_evidence.json) |
| V179 | Coupled finite form includes cross-azimuth entries | `anisotropic.coupled_form` | [运行证据](round09/v179_evidence.json) |
| V180 | Complete finite-to-tail coupling columns | `anisotropic.tail_couplings` | [运行证据](round09/v180_evidence.json) |
| V181 | Degree-graded lower and upper Schur comparisons | `anisotropic.CoupledKernel.matrix` | [运行证据](round09/v181_evidence.json) |
| V182 | Whole mean-zero sphere ground certification and exact replay | `anisotropic.certify / verify` | [运行证据](round09/v182_evidence.json) |
| V183 | Sphere-identity potential reduction improves global range | `anisotropic.sphere_reduce / potential_range` | [运行证据](round09/v183_evidence.json) |
| V184 | Exact proper-rotation certificate transfer | `anisotropic.transfer_rotation / verify_rotation` | [运行证据](round09/v184_evidence.json) |
| V185 | Exact sphere-dilation transfer of full mean-zero ground bounds | `sphere_scale.scale_ground` | [运行证据](round10/certificates/V185_R2_laplace.json) |
| V186 | Ambient polynomial coordinate pullback and operator scaling | `sphere_scale.physical_pullback` | [运行证据](round10/certificates/V186_ambient_pullback.json) |
| V187 | Natural-area explicit trial mass, Dirichlet and wide-potential energy scaling | `sphere_scale.scale_trial` | [运行证据](round10/certificates/V187_normalized_t2.json) |
| V188 | Preserve ordered eigenvalue multiplicities and scale the complete negative trace | `sphere_scale.scale_ordered` | [运行证据](round10/certificates/V188_scaled_minus6_unit_spectrum.json) |
| V189 | Uniform radius perturbation enclosure with sign-aware interval products | `sphere_scale.radius_cell` | [运行证据](round10/certificates/V189_laplace_radius_cell.json) |
| V190 | Complete-coverage adaptive global minimization over positive radii | `sphere_scale.adaptive_minimum` | [运行证据](round10/certificates/V190_initial_nonconstant_global.json) |
| V191 | Verified cache continuation retaining old incumbents and monotone parent bounds | `sphere_scale.resume_minimum` | [运行证据](round10/certificates/V191_checkpoint.json) |
| V192 | Explicit scaled function counterexamples with threshold binding | `sphere_scale.transfer_counterexample` | [运行证据](round10/certificates/V192_explicit_R2_counterexample.json) |
| V193 | Choose and bind suitable exact backends by degree and required scaled accuracy | `sphere_scale.backend_ground` | [运行证据](round10/certificates/V193_precision_backend.json) |
| V194 | Unified radius research jobs with content-bound proof replay and honest tri-state decisions | `sphere_scale.research_job` | [运行证据](round10/certificates/V194_proved_interval_inequality.json) |
