# 几何证明工具 V23–V32

本轮把已有谱计算工具接入“数学表达式 → 检查条件与约化 → 合成或复用引理 → 证明/反例/未解决”的流程。先读 `REPORT.md` 看实测结果，逐版论证见 `PROOF_V23_V32.md`。

这套程序不调用模型，不会自行读懂任意自然语言数学题。人或模型需要提供明确的势函数、参数域、函数空间、代价和目标阈值。工具会选择已实现的数学路线并输出可重放的证明链。

## 直接使用

Python 3.9+，只需标准库，不需要 API 密钥。在本目录运行：

```sh
python3 run.py solve examples/uniform_linear.json --output result.json
python3 run.py verify result.json
python3 run.py solve examples/counterexample.json --output counterexample.json
python3 run.py solve examples/unresolved.json --max-leaves 16 --output unresolved.json
python3 run.py lemma examples/lemma_spec.json --output lemma.json
python3 -m unittest discover -p 'test*.py' -v
python3 benchmark_proofs.py --repeats 3
python3 run.py verify results/proofs
```

独立编译可用 `run.py compile`。`--no-synthesis` 关闭新引理合成；`--library library.json` 接收已有引理证书组成的 JSON 列表。库中的每份引理仍会核验，散列值不代替证明。

## 输入和数学含义

```json
{
  "geometry": "unit_sphere",
  "function_space": "full_sphere",
  "eigenvalue_index": 1,
  "parameters": ["a"],
  "potential": "a*t",
  "penalty": "a*a/4",
  "domain": {"kind": "all_real"},
  "threshold": "0"
}
```

这里 `t=cos(theta)`，要求证明 `inf_a [lambda_ground(-Delta+a*t)+a*a/4] >= 0`，且函数空间是完整球面能量型空间。谱编号 1 指最低特征值，包含常数零模态；不自动表示第一非零特征值。

盒域写成 `{"kind":"box","box":[[-2,2]]}`，顺序对应 parameters。可用整数和有理数表达式，如 `1/3`；势对参数必须仿射、对 t 次数不超过 6，代价为二次函数。名称之外的符号、函数调用、浮点常数、参数分母和更一般几何会被拒绝。

## 十版分别做什么

| 版本 | 新能力 | 代码 |
| --- | --- | --- |
| V23 | 精确展开表达式，提取算子、参数与二次代价 | `compiler_v23.py` |
| V24 | 为轴对称势，把完整球面基态精确约化到轴对称计算 | `reductions.sphere_ground` |
| V25 | 吸收常数势，消去只影响代价的变量 | `reductions.absorb_and_eliminate` |
| V26 | 检查反射对称性并把参数域折半 | `reductions.reflect` |
| V27 | 检测势的实际参数秩，并精确压缩参数 | `reductions.compress` |
| V28 | 围绕指定不等式阈值搜索证明或反例 | `decision_v28.py` |
| V29 | 解多项式 Poisson 方程，生成全实参数的正函数谱引理 | `family_v29.py` |
| V30 | 把全参数矩阵界转换成一维多项式最值并优化系数 | `matrix_v30.py` |
| V31 | 按仿射参数匹配和势的单调性复用引理 | `transport_v31.py` |
| V32 | 搜索约化路线、选择引理或数值方法，核验完整证明链 | `planner_v32.py` |

V23–V27 也可作为独立的转换工具；V29–V30 生成可单独保存的参数族引理；V32 是总入口。原来的 V13–V22 模块保留为数值后端，原成果目录和 ZIP 未改写。

## 结果

- `proved`：原目标在完整声明域上成立。
- `refuted`：给出原域内的实际参数、非零试探函数及违背目标的精确能量。
- `unresolved`：当前方法和预算未能决定目标。
- `scale_enclosed` / `scale_gap_open`：引理内部放大系数的最值误差是否达标。较松的系数仍能生成有效引理；这两者不是原问题是否已解决的状态。
- `unsupported_or_invalid_input`：命令行拒绝了输入或当前没有适用后端，不作数学结论。

任何较小的数值优化差距，都不能代替原题的不等式判断。对 `E_*>=0`，只要下界仍为负数就不能标成 proved。

## 范围与验证

固定单位二球面、轴对称多项式势、最低特征值。V24 和参数族引理可覆盖完整球面上的所有函数；没有支持任意曲面、非轴对称势或激发态的同类约化。全实数值后端和秩压缩通常要求正定二次代价；参数族引理可处理部分半正定情形。

证明链依靠本文的解析规则、精确有理算术和可重放检查。谱核验另用稠密惯性，Bernstein 另用直接展开；部分编译及代数实现仍共享。`formal_assistant_checked=false`，尚未进入形式化证明助手。

`MODEL_EVALUATION_STATUS.json` 明确记录模型对照尚未运行。`MODEL_TOOL_GUIDE.md` 可交给模型作为调用说明。现有对照都是开发及构造性迁移测试，不是保留题库上的模型能力成绩。
