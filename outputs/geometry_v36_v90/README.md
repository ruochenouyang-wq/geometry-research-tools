# 几何研究工具 V36–V90

这是供 GPT 调用的数学研究工具原型：把反复展开公式、搜索候选、估计误差和检查证明交给确定性程序，模型负责选问题、选工具和解释结果。当前数学范围是**单位二维球面上的基态谱不等式、相关一元多项式极值和 Poisson 矩阵包络**，不是通用几何定理证明器。

本轮从 V35 新增 55 个有代码、测试和证据对应的功能里程碑：V36–V85 为用户要求的 50 次迭代，再完成 V86–V90。V36–V55 优化可见文本 token；达到事先声明的工程关口后，V56–V90 转向误差界、全局优化与集成。台账见 [VERSIONS.md](VERSIONS.md)，完整结果见 [REPORT.md](REPORT.md)。它们是同一个工具包的功能演进，并非 55 份独立历史快照。

## 能做什么

| 入口 | 用途 | 结论范围 |
|---|---|---|
| `extrema` | 有理多项式的全区间最大值，处理重根、平坦峰和严格误差 | 指定有理闭区间、次数不超过 32 |
| `positive` | 搜索更好的正函数，收紧固定势的 Barta 谱下界 | 对数多项式基最多 6 维；另算真实 Rayleigh 谱上界 |
| `refine` | 改善矩阵条件数、全局对偶下界，并精确重估上界 | Poisson 包络函数族内的全局误差 |
| `uniform` | 证明全振幅不等式，或给出可核验的障碍/真实反例 | 全球面；默认 log-Sobolev 路径限仿射坐标势 |
| `research` | 把统一不等式用于用户提交的原始参数化目标 | 势的每个参数方向为 `b+d*t`；基础势可到六次 |
| `submit`、`solve` | 保留 V35 的多项式参数方向搜索和已验证引理复用 | 原 V35 支持范围 |
| `inspect`、`checkpoint`、`restore` | 只读取需要的证据，保留完整任务绑定并恢复会话 | 引用依赖原来的本地存储目录 |

最明显的数学改进之一是

\[
\lambda_1(-\Delta_{S^2}+a t)\ge-\frac{a^2}{6},\qquad a\in\mathbb R,
\]

且该族中的系数 `1/6` 最优。此路径明确采用经典球面 log-Sobolev 定理，不声称程序从公理证明了该定理。另有正函数和有理张量 Bernstein 路径证明全参数系数 `1/5`。[证明与来源](PROOF_UNIFORM_REFINE.md)

## 运行

数学核心和证书核验只需 Python 标准库。本次基准使用 Python 3.9.6。带 token 预算的接口另外需要固定版本的真实文本编码器：

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-token.txt
.venv/bin/python run.py call examples/v90/original_goal_sharp.json --store .geometry-store
```

返回包含经过核验的状态、向外取整的界和证据引用。`wire_tokens` 包含完整已知返回包装；预算容纳不下必需信息时明确返回失败。两套编码器的词表已随包保存；安装依赖以后，例子不需调用模型或访问网络。

持续多轮交互使用：

```sh
.venv/bin/python run.py serve --store .geometry-store
```

每行提交一个 JSON 请求；例如：

```json
{"op":"uniform","spec":{"coefficient":"1/7"}}
```

这会得到 `disproved`，证据含实际 Rayleigh 反例 `-1/91`。完整调用说明及可复制给模型的工具规则见 [MODEL_TOOL_GUIDE.md](MODEL_TOOL_GUIDE.md)。

## 重放与复现

```sh
.venv/bin/python -m unittest discover -s . -p 'test_*.py'
.venv/bin/python run.py verify results/service_certificates
.venv/bin/python benchmark_precision.py --verify results/math_benchmark.json
.venv/bin/python benchmark_tokens.py
.venv/bin/python version_token_experiments.py
.venv/bin/python benchmark_service.py
```

重新测数学耗时使用 `benchmark_precision.py`，它会覆盖本包对应的基准结果文件。历史 V35 的对照输入及证书已放在 `fixtures/v35`，解压后可单独复现。`results/final_validation.json` 记录发布前检查；`MANIFEST.json` 记录文件指纹。

`verified` 表示软件对有理证据和规则绑定作了重放；也可以核验一个诚实的 `unresolved` 结果。只有对应目标的成功状态才表示证明成功。函数族上界不能用作真实谱上界，且本包没有通过 Lean 等形式化内核。

**尚未做 GPT-5.6/Astra 的实际任务对照。** 本轮可见文本计数、毫秒级算子基准和开发 agent 的使用，均不能证明“5.6 半小时胜过 Astra 八小时”。参见 [MODEL_EVALUATION_STATUS.md](MODEL_EVALUATION_STATUS.md)。
