# 全部顺序对照实验

每项运行 3 次，时间取中位数。计时包括搜索与最终独立核验，不含 JSON 保存。环境：Python 3.9.6；macOS-26.6.2-arm64-arm-64bit。同一案例的三次精确上下界与状态一致。

点数是不同参数点数量，谱搜索次数允许大于点数（同一点自适应多个 N），也允许小于点数（等价轨道复用）。不同目标或不同容限之间的秒数不构成速度比较。

| 算例 | 状态 | 全局差距（近似） | 谱点 | 谱搜索 | 中位秒数 |
| --- | --- | ---: | ---: | ---: | ---: |
| v11_fixed12_tilt | epsilon_global | 0.00053622131 | 23 | — | 0.219638 |
| v13_tilt | epsilon_global | 0.00056201734 | 23 | 31 | 0.075398 |
| v14_tilt | epsilon_global | 0.00056201734 | 23 | 23 | 0.072884 |
| v15_tilt | epsilon_global | 0.00056201734 | 17 | 18 | 0.064270 |
| v22_tilt | epsilon_global | 0.00057568921 | 17 | 16 | 0.055405 |
| v13_fixed4_floor | global_gap_open | 5.121027e-05 | 33 | 33 | 0.078390 |
| v16_independent_moments | epsilon_global | 5.411588e-05 | 35 | 100 | 0.842040 |
| v16_cuts | epsilon_global | 5.411588e-05 | 34 | 97 | 0.827865 |
| v17_3d | epsilon_global | 0.083913115 | 119 | 64 | 0.983337 |
| v18_3d | epsilon_global | 0.092249632 | 41 | 22 | 0.400242 |
| v22_3d | epsilon_global | 0.093616819 | 41 | 22 | 0.435995 |
| v17_4d | global_gap_open | 0.11910369 | 389 | 204 | 6.658438 |
| v18_4d | epsilon_global | 0.095132346 | 229 | 121 | 5.625378 |
| v22_4d | epsilon_global | 0.096499533 | 229 | 121 | 6.325610 |
| v19_three_to_two | epsilon_global | 0.00079256225 | 9 | 6 | 0.166949 |
| v20_symmetric | epsilon_global | 1.9644659e-06 | 2 | 2 | 0.007477 |
| v20_asymmetric | epsilon_global | 8.3304991e-05 | 6 | 6 | 0.019193 |
| v20_budget_open | global_gap_open | 0.11577343 | 2 | 3 | 0.010722 |
| v20_two_design_two_uncertainty | epsilon_global | 2.8597147e-05 | 8 | 8 | 0.103958 |
| v21_target_1_2 | epsilon_global | 0.00058456125 | 7 | 8 | 0.058387 |
| v21_target_9_10 | epsilon_global | 4.6141816e-06 | 12 | 33 | 0.168642 |
| v21_target_2 | certified_infeasible | — | — | — | 0.000055 |
| v21_feasibility_not_found | feasibility_not_found | — | 1 | 1 | 0.001355 |
| v15_extra_weak | epsilon_global | 3.3403192e-06 | 3 | 2 | 0.005695 |
| v22_extra_weak | epsilon_global | 1.7012194e-05 | 3 | 2 | 0.005591 |
| v15_extra_asymmetric_quadratic | epsilon_global | 0.00016182447 | 11 | 12 | 0.060122 |
| v22_extra_asymmetric_quadratic | epsilon_global | 0.00017549634 | 11 | 10 | 0.049967 |
| v15_extra_even_quartic_2d | epsilon_global | 0.00093745169 | 77 | 80 | 1.971688 |
| v22_extra_even_quartic_2d | epsilon_global | 0.00095112357 | 77 | 80 | 1.960996 |
| v15_extra_spectral_action | epsilon_global | 0.00069031468 | 19 | 13 | 0.052854 |
| v22_extra_spectral_action | epsilon_global | 0.00070398655 | 19 | 12 | 0.051967 |

原始数据：`results/benchmark.json`。逐份证书：`results/certificates/`。
