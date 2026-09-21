# 十个新增版本 V23–V32

十版均有代码、逐版数学论证、针对性测试和重复实验。完成指研发任务完成，不意味着任意输入都可解决。

| 版本 | 实质改动 | 验证证据 | 状态 |
| --- | --- | --- | --- |
| V23 | 表达式精确编译 | 提取仿射势和二次代价，拒绝非法模型；`results/proofs/v23_expression_compilation.json` | 完成 |
| V24 | 完整球面基态约化 | 用方位角 L² 包络论证完整函数空间的精确约化；`results/proofs/v24_full_sphere_ground.json` | 完成 |
| V25 | 无谱作用变量消元 | 同目标参数点记录 195→3，含核验耗时约 173 倍改善；`results/proofs/v25_with_elimination.json` | 完成 |
| V26 | 参数反射商域 | 同目标参数点记录 13→7，约 1.53 倍改善；`results/proofs/v26_half_box.json` | 完成 |
| V27 | 势函数秩压缩 | 4 参数压到 1 参数；参数点记录 27→3，约 22 倍改善；`results/proofs/v27_one_effective_parameter.json` | 完成 |
| V28 | 目标阈值驱动搜索 | 同一否证任务参数点记录 37→5，约 9.96 倍改善；`results/proofs/v28_targeted_decision.json` | 完成 |
| V29 | Poisson–正函数引理生成 | 自动生成覆盖全部实参数的谱下界；`results/proofs/v29_scalar_quadratic_lemma.json` | 完成 |
| V30 | 秩一矩阵条件消元 | 全参数矩阵约束化为一维最值；保留系数预算失败；`results/proofs/v30_matrix_lemma.json` | 完成 |
| V31 | 引理仿射迁移与势单调性 | 16 个构造性迁移目标均零谱点完成；`results/proofs/v31_affine_dominance_transfer.json` | 完成 |
| V32 | 有条件的证明路线搜索 | 自动选择约化、引理或数值路线，核验回原题；`results/proofs/v32_uniform_two_directions.json` | 完成 |

对应测试类为 `test_proof_versions.py` 中的 V23–V32。所有新证明规则见 `PROOF_V23_V32.md`。模型对照状态单独写在 `MODEL_EVALUATION_STATUS.json`，本轮没有将算法实验记为模型能力成绩。
