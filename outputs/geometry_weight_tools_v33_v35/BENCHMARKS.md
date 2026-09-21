# 完整实测：V33–V35

Sequential development probes. Time includes construction/search AND explicit replay of the final certificate and diagnostics. Result file I/O is excluded. Not model or held-out evaluation.

每项顺序运行三次，取中位数。计时包含构造/搜索以及最终证书和诊断的显式重放，不含保存结果文件。环境：3.9.6，macOS-26.6.2-arm64-arm-64bit。全部为开发算例，模型调用次数为 0。

## 矩阵系数目标对照

以下上、下界属于固定 Poisson 形式内的 trace(C G) 全局优化，数值为近似显示；精确有理数见原始 JSON。

| 方向/代价 | 单位权重目标值 | 优化后下界 | 优化后上界 | 全局差距 | 单位权重秒数 | 优化秒数 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| linear | 0.250000000 | 0.249999444 | 0.250000000 | 5.56e-07 | 0.000301 | 0.005501 |
| odd_even | 0.500000000 | 0.337363467 | 0.337367823 | 4.36e-06 | 0.000784 | 0.065470 |
| nonorthogonal | 2.551275651 | 1.445401420 | 1.445409460 | 8.04e-06 | 0.004619 | 0.211351 |
| quartic | 0.500000000 | 0.297405954 | 0.297407425 | 1.47e-06 | 0.001596 | 0.066262 |
| degree_six | 0.500000000 | 0.281519855 | 0.281528107 | 8.25e-06 | 0.002612 | 0.032009 |
| three_directions | 0.937500000 | 0.453161187 | 0.453171113 | 9.93e-06 | 0.001770 | 0.223917 |
| four_directions | 1.250000000 | 0.516027556 | 0.516033842 | 6.29e-06 | 0.003118 | 1.152336 |
| weighted_cost | 0.750000000 | 0.595479672 | 0.595486736 | 7.06e-06 | 0.000741 | 0.217211 |

三个及四个方向使用 t、t²、t³（及 t⁴）。默认 C=I，weighted_cost 使用 C=[[2,1/2],[1/2,1]]。非正交算例的方向为 t+t² 与 2t−t²。系数目标下降不表示新矩阵在所有参数方向上都更小。

## 预算失败与原题判断

| 算例 | 结果 | 中位秒数 | 诊断 |
| --- | --- | ---: | --- |
| v34_one_round | global_gap_open | 0.014747 | — |
| v34_range_budget | global_gap_open | 0.021087 | — |
| v32_anisotropic_goal | unresolved | 2.170818 | — |
| v35_anisotropic_goal | proved | 0.004307 | poisson_envelope_fits |
| v32_mapped_goal | unresolved | 2.789523 | — |
| v35_mapped_goal | proved | 0.005069 | poisson_envelope_fits |
| v35_weak_poisson_ansatz | unresolved | 0.151711 | poisson_envelope_incompatible |
| v35_counterexample | refuted | 0.022671 | poisson_envelope_incompatible |
| v32_generic | proved | 0.040336 | — |
| v35_generic | proved | 0.041802 | poisson_envelope_incompatible |

`one_round` 的矩阵全局误差约为 0.00417302；`range_budget` 的误差约为 0.0388088。两者仍是有效界，但未达到目标精度。

原题对照最多使用 32 个区域叶节点；弱构造和反例案例为 16 个。两类未解决→已证明的对照不能据此计算精确加速倍数。一般算例启用新步骤后稍慢，保留在表中。

原始三轮耗时、界和统计：`results/benchmark.json`。每项完整证明对象：`results/proofs/`。
