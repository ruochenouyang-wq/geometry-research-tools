# 完整对照结果

每项顺序运行三次，取中位数。所有组的重复运行得到相同的精确数学输出。以下是结构化任务或部件实验，模型调用次数为 0；不是模型从自然语言原题出发的端到端成绩。

环境：Python 3.9.6，macOS-26.6.2-arm64-arm-64bit。计时包括适用的表达式编译、约化、求解与证明重建，不含文件保存。

| 算例 | 版本 | 结果 | 参数点记录 | 中位秒数 |
| --- | --- | --- | ---: | ---: |
| v23_expression_compilation | V23 | verified_component | — | 0.000113 |
| v24_full_sphere_ground | V24 | proved | 1 | 0.007568 |
| v25_without_elimination | V25 | proved | 195 | 3.031632 |
| v25_with_elimination | V25 | proved | 3 | 0.017502 |
| v26_full_box | V26 | proved | 13 | 0.057311 |
| v26_half_box | V26 | proved | 7 | 0.037411 |
| v27_four_parameters | V27 | proved | 27 | 0.606039 |
| v27_one_effective_parameter | V27 | proved | 3 | 0.027497 |
| v28_optimize_then_decide | V28 | refuted | 37 | 0.236356 |
| v28_targeted_decision | V28 | refuted | 5 | 0.023724 |
| v29_scalar_linear_lemma | V29 | scale_enclosed | — | 0.000193 |
| v29_scalar_quadratic_lemma | V29 | scale_enclosed | — | 0.000330 |
| v29_scalar_degree_six_lemma | V29 | scale_enclosed | — | 0.003987 |
| v30_matrix_lemma | V30 | scale_enclosed | — | 0.000512 |
| v30_nonidentity_weight | V30 | scale_enclosed | — | 0.002934 |
| v30_scale_budget_open | V30 | scale_gap_open | — | 0.001487 |
| v31_affine_dominance_transfer | V31 | proved | 0 | 0.004776 |
| v32_uniform_linear | V32 | proved | 0 | 0.001046 |
| v32_no_lemmas_uniform_linear | V32 | unresolved | 65 | 0.563285 |
| v32_uniform_two_directions | V32 | proved | 0 | 0.002757 |
| v32_no_lemmas_uniform_two_directions | V32 | unresolved | 94 | 1.759342 |
| v32_counterexample | V32 | refuted | 3 | 0.018350 |
| v32_too_weak_lemma | V32 | unresolved | 43 | 0.169084 |
| v32_four_parameter_counterexample | V32 | refuted | 5 | 0.037616 |
| v32_generic_asymmetric | V32 | proved | 3 | 0.035093 |
| v32_no_lemmas_generic_asymmetric | V32 | proved | 3 | 0.033509 |

参数点记录对应 `points`，包含复用结果，不等于真实谱求解次数；实际求解调用数为原始统计中的 `spectral_searches`，其中也包含自适应加阶后的再次求解。`verified_component` 是通过编译核验的部件输出，不是整道研究题已解决。`scale_gap_open` 是引理常数仍有优化余量，不能与整题的 unresolved 混为一谈。原始重复时间、界与统计见 `results/benchmark.json`。

迁移批次另含 16 个由仿射改参、平移和非负六次势构造的目标；复用一个线性族引理，全部 proved、全部零谱点。通用引理生成和批次运行分别计时，原始记录见 `results/transfer_batch.json`。这些探针用于检查已声明的迁移规律，不宣称是未见题族上的能力评测。
