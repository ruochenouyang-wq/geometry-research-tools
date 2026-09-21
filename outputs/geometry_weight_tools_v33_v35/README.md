# 几何权重与证明工具 V33–V35

本轮增加：全局矩阵系数下界、自动权重优化、针对原题二次代价的引理选择。实测见 [REPORT.md](REPORT.md)，完整论证见 [PROOF_V33_V35.md](PROOF_V33_V35.md)。

Python 3.9+，只需标准库，无需 API 密钥。程序不会自动理解任意自然语言题目，不会调用模型。

## 快速使用

在本目录运行：

```sh
python3 run.py solve examples/anisotropic_goal.json --output result.json
python3 run.py verify result.json
python3 run.py optimize examples/two_directions.json --output weights.json
python3 run.py verify weights.json
python3 run.py optimize examples/four_directions.json --output four.json
python3 run.py solve examples/weak_ansatz.json --max-leaves 16 --output weak.json
```

优化输出中的 `certificate.template` 是可复用引理。`--library` 可以直接读取优化结果、单个引理，或引理组成的列表；先核验输入再使用：

```sh
python3 run.py solve examples/reuse_goal.json --library weights.json --no-synthesis --output reused.json
```

引理是否适用于目标还要重新检查；同一目标的 trace 最优矩阵不一定在每个方向上都小于原题代价。优化结果有未闭合的全局误差时，其中已核验的可行引理仍可以使用。

复现实验：

```sh
python3 -m unittest discover -p 'test*.py' -v
python3 benchmark_weights.py --repeats 3
python3 run.py verify results/proofs
```

## 两种输入

原题使用结构化表达式。`t=cos(theta)`，谱编号 1 包括常数零模态，指最低特征值：

```json
{
  "geometry": "unit_sphere",
  "function_space": "full_sphere",
  "eigenvalue_index": 1,
  "parameters": ["a", "b"],
  "potential": "a*t+b*(t*t-1/3)",
  "penalty": "a*a/4+b*b/9",
  "domain": {"kind": "all_real"},
  "threshold": "0"
}
```

这要求证明完整参数域上的 `inf [λ₁(−Δ+potential)+penalty]≥threshold`。盒域改写为 `{"kind":"box","box":[[-2,2],[-3,3]]}`。

系数优化输入直接列出方向多项式的升幂系数：

```json
{
  "directions": [[0, 1], ["-1/3", 0, 1]],
  "cost": [[1, 0], [0, 1]]
}
```

其中 `cost` 可省略，默认为单位矩阵。优化的目标是 `trace(cost × G)`；不是原问题的特征值。每个方向次数不超过 6，方向数 1–4，模常数必须线性独立，cost 必须正定。

`optimize` 支持 `--tolerance 1/100000`、`--max-rounds 8`、`--max-range-leaves 128`。`solve` 支持 `--no-penalty-weights` 关闭本轮的目标权重步骤，`--no-synthesis` 关闭全部新引理生成；已有库仍可使用。

## 如何阅读结果

| 状态 | 含义 |
| --- | --- |
| `global_gap_closed` | 所定义的矩阵系数目标已有全局上下界，差距达到指定容限 |
| `global_gap_open` | 上下界有效，但预算内没有达到该容限 |
| `proved` | 原始声明目标已证明 |
| `refuted` | 有实际参数和试探函数构成原题反例 |
| `unresolved` | 当前方法与预算未决定原题 |
| `poisson_envelope_fits` | 这套 Poisson 构造的矩阵损失可由代价 Hessian/2 覆盖；原题的其他项仍需检查 |
| `poisson_envelope_incompatible` | 这套构造无法用这样的矩阵界覆盖代价；不代表原题错误 |
| `budget_test_unresolved` | 连这种矩阵覆盖检查也暂未决定 |
| `unsupported_or_invalid_input` | 输入不合法或没有适用方法，没有数学结论 |

`diagnostics` 有单独的证明对象，核验命令也会检查它们。数值搜索轨迹和耗时是实验记录，不是证明条件。

## 范围

固定单位二球面、轴对称多项式势、最低特征值。完整球面函数空间的适用性由原有约化或正函数证明支撑；这不支持任意曲面、非轴对称势和激发态。

矩阵目标的全局最优仍限于固定 Poisson 指数形式，不表示原谱不等式的锐利常数最优。对偶接受条件与上界重放使用有理数，但解析规则尚未进入形式化证明助手；部分代数实现共享。

V23–V32 的数学后端与测试在本包内保留，不修改之前的成果。模型对照仍未运行，状态见 `MODEL_EVALUATION_STATUS.json`。
