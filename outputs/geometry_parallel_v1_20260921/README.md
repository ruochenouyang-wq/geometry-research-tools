# 三个并列的几何研究工具版本

本轮比较三条不同路线。它们从同一份只读基线出发，分别扩展数学能力、复用参数证明、减少不必要计算。结果和选择理由见 [比较报告](REPORT.md)，评测规则见 [计划](PLAN.md)。

**谱值**可理解为数学算子的能量层级。**证书**是保存原问题、数学界限和核对所需数据的记录。**重放**是按原数学规则重新核对证书，不依赖“成功”状态文字。

| 版本 | 目的 | 代码 |
|---|---|---|
| A | 认证负谱个数、重数和负谱和区间；重数是同一谱值的独立方向数量 | [a_multispectrum/variant.py](a_multispectrum/variant.py) |
| B | 证明固定指数下整段幅度都满足精度，再快速回答其中的查询 | [b_parameter/variant.py](b_parameter/variant.py) |
| C | 根据目标误差检查足够强的谱界，减少固定精度的重复二分 | [c_adaptive/variant.py](c_adaptive/variant.py) |

三者不是前后继承的版本。A、B 的普通单点任务明确回退原基线；新增能力通过专项接口调用。C 对当前未优化的任务也保留原路线。各版都仍受原函数与证明前提约束。

## 本轮实测

冻结后的四系统对照已完成：共同 8 题各重复两次全部达标，C 总壁钟成本相对基线下降 33.9%；B 在 16 次参数查询中含准备约快 5.16 倍，加入逐题独立冷重放后约快 1.87 倍；A 的负谱计数 6/6 成功，谱和精度 4/6 达标。39 次证书攻击全部拒绝。这是公开开发对照，没有未见题或真实模型 token 结论。

详细精度代价、失败和路线选择见 [比较报告](REPORT.md)；原始数据见 [最终结果](evaluation/final_20260922T013006751900Z/report.json)。完整比较的复现命令见 [评估说明](evaluation/README.md)。

## 运行普通任务

使用 Python 3.12，在整个项目根目录运行（必须保留相邻历史模块）：

```sh
python3 -B outputs/geometry_parallel_v1_20260921/run.py --variant c < outputs/geometry_parallel_v1_20260921/examples/spectrum.jsonl > response.jsonl
```

`--variant` 可选 `baseline`、`a`、`b`、`c`。输入为 JSONL，即每行一个 JSON 任务。标准输出逐行返回结果，数学参数使用精确分数字符串。公开示例已经用于开发，不可当作新的未见题。

## 使用 A 的负谱接口

```sh
python3 -B outputs/geometry_parallel_v1_20260921/run.py --variant a < outputs/geometry_parallel_v1_20260921/examples/negative_spectrum.jsonl > negative_result.jsonl
```

严格负谱只计小于零的谱值；零模不计入。请分别查看计数是否确定、负谱和是否达到宽度目标。计数成功不能自动解释为和的精度成功。有符号和是负数；其绝对值和由同一上下界反向取负得到，数学含义在证书中注明。

## 使用 B 的参数接口

**参数证书库**是一组覆盖整个幅度区间的证明，不是只保存若干计算点。先准备一份库：

```sh
python3 -B outputs/geometry_parallel_v1_20260921/run.py --variant b --operation prepare-family < outputs/geometry_parallel_v1_20260921/examples/family.jsonl > family_result.jsonl
```

结果中的 `bank` 为库；若为空，准备未完成。对 `--operation query-family`，每行输入 `{"bank":完整库,"task":原单点任务}`。这一路会重放完整证明，便于独立核验。

程序内的重复查询可使用 `prepare_context(bank)`：先验证库，再在保存的私有快照上执行 `context.query(task)` 与 `context.verify(certificate,task)`。这种会话的准备成本必须计入总成本。它与每次重新加载完整证明的冷重放是不同工作方式，比较报告分别记账。

## 单独核验证书

对 `--operation verify`，每行输入 `{"task":原任务,"certificate":完整证书}`，并选对应版本。返回的 `certificate_valid` 与 `target_met` 分别表示证明有效、原精度达标；完整原题条件必须一致。

B 的 `--operation verify-family` 接收库本身。检查结果不能只根据证书自行填写的成功状态判断。

## 结果的范围

本轮为公开开发对照，没有重新揭示历史随机种子，也不把已有 B001 题目当作新题。没有调用真实模型 API，因此没有实际模型 token 或跨模型优势结论。证书软件重放还不是形式化证明助手中的完整数学证明。
