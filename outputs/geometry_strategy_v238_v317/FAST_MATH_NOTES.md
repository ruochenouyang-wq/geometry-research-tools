# 独立精确性能组件：v248–v257

本组件仅替换重复计算，不改变数学问题、证书格式、候选序列或最后的独立旧验证。`fast_math.py` 不修改或 monkeypatch 任何历史模块。所有数值运算使用 `Fraction`；缓存内容是有界、进程内的不可变元组、分数或规范 JSON 文本。返回给调用者的模型和 source 均为独立副本。

## 可用接口

- `piecewise_model(function, levels=4, degree=3, root_ratio='1/2', sqrt_bits=40)`：旧模型同参数、同逐字段结果；`verify_piecewise` 仍是旧验证器。
- `evaluate_many(certificate, coordinates)`：先拥有证书副本并用旧验证器完整重建一次，然后计算全部精确坐标。`evaluate(certificate, t)` 为单点包装。端点、接缝和轴零处保持旧行为。
- `matrix_assembly(function, m=0, mean_zero=True, modes=8, near_tail=0)` 与 `Kernel(function, m=0, mean_zero=True, modes=8, sqrt_bits=40, near_tail=0)`：旧矩阵数据和公开计算方法兼容。Kernel 内部固定矩阵和 Schur 项是构造时的不可变快照；公开 `data` 供读取，不用于修改内核。`evidence()` 返回独立副本。
- `certify_sector(...)`、`full_ground(...)`：参数及返回证书与旧版本一致。`full_ground` 调用旧 `direct.full_certificate`，该聚合器完整重建每个 sector；没有跳过最终旧验证。
- `SourceContext(source)` 或 `validate_source(source)`：一次完整旧 `direct.verify_full`，绑定固定 singular mean-zero 原问题且显式包含 m=0、m=1。持有冻结规范文本；`.source()` 返回副本。不接受修改已有 context 的操作。
- `trial_matrices(powers)`、`trial_statistics(powers, coefficients)`、`propose_trial(powers, precision_bits=112, iterations=14)`：旧参数兼容；数学字段逐项等值。`elapsed_seconds` 是本次真实耗时，当然不要求相等。
- `trial_certificate(powers, coefficients, source, tolerance=...)`：source 可为普通证书或已验证 `SourceContext`。普通证书每次重新验证；context 可跨同源尝试复用。
- `enrich(source, max_terms=16, tolerance=..., precision_bits=112, iterations=14)`：与旧候选序列及停止规则一致。即使接受 context，结束时仍用旧 `enriched.verify` 完整重放最佳证书。
- `reset_counters(clear_caches=False)` 和 `counters()`：公开测量计数。`clear_caches=True` 清除本组件缓存，不清除历史模块的 Legendre/radial 缓存。

## 每步真实记录

每一步实际修改源文件后运行当时的完整测试，随后立即调用 `campaign.record` 保存当时源码；不存在完工后用同一份源码填充多个版本。v248 原记录保持不动；恢复工作从 v249 开始。每步详情见 `iterations/vNNN/record.json`，汇总见 `metrics/fast_math.json`。

| 步骤 | 实际改动 | 固定输入的精确对照与计数 |
|---|---|---|
| v248 | 不可变单位形状缓存 | 同形状多个层数复用1次形状积分；原步骤已存在 |
| v249 | 局部坐标积分 | 6阶单元全字段相等；全局多项式组合8→0；平方系数峰值位长234→103 |
| v250 | 一次验证的批量评估 | 四点评估旧证书重建4→1；值、接缝、正负、0完全一致；拒绝损坏证书和浮点输入 |
| v251 | 规范源与原始矩表 | m1、N6、near6：旧moment540、normalize541；新唯一矩32、normalize1；完整矩阵相同 |
| v252 | 对称装配和正确奇偶零项 | 同一 singular 输入唯一矩再降到17，矩请求540→222；积分108→78；非偶 step 的奇耦合保留 |
| v253 | 固定 Schur 项与 Kernel | 16次移位外积3456→构造36、查询0；全部矩阵和惯性精确一致，sector/full证书旧重放通过 |
| v254 | 不可变 source context | max_terms6：完整source验证5→2（初始验证与最终重放）；证书逐字段相等 |
| v255 | 幂和矩及前缀 Gram | 10项Gram旧weighted_moment调用1231→41个唯一幂和；嵌套1/3/6/10不重算已有前缀；H双方向精确检查 |
| v256 | LDL前缀扩展与因子复用 | 嵌套1/3/6/10：枢轴20→10，下三角新项63→45；重复提案只求解；LDL重构、解和候选系数精确一致 |
| v257 | 融合对称M/H/R统计 | 10个非零系数：3次全矩阵二次型共300项→55个系数积和165次Gram收缩；完整enrich证书及每次尝试宽度相同 |

## 数学保持性

局部积分采用 `s=(t-a)/(b-a)` 的精确二项式换元，残差仍由真实返回多项式的平方积分和交叉积分独立核算。径向积的奇偶性是 `l+k`；abs_power 为偶函数。step 一般不为偶函数，因此只删除其“同奇偶且不同阶”的常数偶部分正交项；奇耦合不会误删。step 平方的奇部分只在 `a*(a+2*b)=0` 时消失。

Schur 重排仅把固定的 `column_i*column_j/mass` 提前计算；每次移位仍使用精确分母和严格正尾条件。Gram 使用精确幂和矩，并对 H 的两个方向独立计算后比较。LDL 新行使用标准块扩展，所有枢轴要求严格为正。融合二次型仅利用已检查的对称性，保留完整强算子 Gram R；未以投影残差替换完整残差。

## 测量限制与保留的负收益

记录内所有秒数均为当前机器上的三次串行探针，不是全局性能承诺。共享旧 Legendre/radial 缓存、候选原始矩冷缓存、Gram冷/热、LDL冷/热分别注明；构造成本与移位查询分开报告。

v251 固定输入的候选暖矩表中位数约1.118ms，稍慢于候选冷测量1.023ms，保留此测量波动。v256 候选完整冷构造（Gram和LDL）约2.406ms，比旧版已热Gram但重算LDL的约1.112ms慢；候选全热约0.793ms。二者缓存条件不同，不能将该冷暖差异解释为算法普遍退步或普遍获益。

v257 暖统计约0.265ms→0.153ms；相同六项、112位、14次迭代、每次fresh context且结束旧重放的enrich约19.76ms→9.73ms。只有这些明确固定输入受到支持。

最终测试12项覆盖逐字段模型、批评估、矩表、多个m/mean_zero/near组合、step反例、Schur和惯性、旧sector/full/enriched重放、source变异隔离、Gram不可变前缀、LDL重构、精确解、不同符号/稀疏统计与非法输入。所有旧目录只读；Python均以`-B`运行。
