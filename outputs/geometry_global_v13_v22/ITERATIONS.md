# 已完成的十个新增版本：V13–V22

每版均有实质算法改动、数学论证、针对性测试和运行记录。公共底层合并复用，V19–V21 是新增问题类型的独立入口。

| 版本 | 实质改动 | 实现 | 本轮证据 | 状态 |
| --- | --- | --- | --- | --- |
| V13 | 按点自适应谱空间 | point_oracle.py | 强势点 20t：固定 N=4 点误差 0.2901；自适应到 N=8 后 5.44×10⁻⁷；见 `results/certificates/v13_tilt.json` | 完成 |
| V14 | 反射和常数平移复用 | point_oracle.py | 33 个平移等价势的谱搜索由 33 次降为 1 次，全部仍核验原势证书；见 `results/certificates/v14_tilt.json` | 完成 |
| V15 | 中点单纯形区域下界 | simplex_v15.py / global_core.py | 双势阱同容限谱点由 23 降为 17；见 `results/certificates/v15_tilt.json` | 完成 |
| V16 | 联合矩约束排除 | global_core.py | t、t² 方差关系排除 1 格，谱点 35→34；收益有限；见 `results/certificates/v16_cuts.json` | 完成 |
| V17 | 三、四参数张量 Bernstein | bernstein_v17.py / global_core.py | 3 参数全局证书；4 参数原始方案在 64 叶预算下仍未闭合；见 `results/certificates/v17_3d.json` | 完成 |
| V18 | 不增加谱点的多项式细化 | bernstein_v17.py | 3 参数所需谱点 119→41；4 参数在同样 64 叶上限内闭合；见 `results/certificates/v18_4d.json` | 完成 |
| V19 | 精确消去线性等式 | equality_v19.py | 3 参数化为 2 个自由参数并把候选提升回原约束；见 `results/certificates/v19_three_to_two.json` | 完成 |
| V20 | 连续不确定性的鲁棒对偶界 | robust_v20.py | 连续盒不确定性；含 2 设计 + 2 不确定参数的完整证明；见 `results/certificates/v20_two_design_two_uncertainty.json` | 完成 |
| V21 | 期望值约束的谱 Lagrange 界 | constrained_v21.py | m_t≥1/2 的可行函数与全局下界；另测不可行和未找到可行解；见 `results/certificates/v21_target_1_2.json` | 完成 |
| V22 | 按误差来源调度三种工作 | global_core.py | 组合谱细化、代数细化、参数划分；高维仍比固定 V18 稍慢；见 `results/certificates/v22_tilt.json` | 完成 |

每版测试对应 `test_new_versions.py` 中的 `V13...` 至 `V22...` 测试类。数学条件逐版写在 `PROOF_V13_V22.md`。

“完成”指研发和验证任务完成，并不表示每个算例都达到容限。未闭合区间、不可行目标、未找到可行函数与提速不足均保留在结果中。V11–V12 的原始成果目录未改写。
