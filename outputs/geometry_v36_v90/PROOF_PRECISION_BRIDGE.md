# V88：精确极值与既有谱证书的独立桥接

`precision_bridge.py` 新增三个接口，把精确多项式极值接到已有正函数与矩阵
Poisson 包络上。没有修改 `exact_extrema`、`positive_search`、`family_v29`、
`matrix_refine`、`dual_v33` 或其他旧模块。新范围证据保留原生格式，分别调用
对应验证器；不会把 Sturm 证据伪装成旧 Bernstein 证书。

## 1. ExactRangeBackend：严格反号和区间绑定

```python
from fractions import Fraction as F
from precision_bridge import ExactRangeBackend

backend = ExactRangeBackend()
proof = backend.minimum([0, 1, 0, -1], F(1, 10**14), 33)
lower, upper = backend.verify_minimum([0, 1, 0, -1], proof)
```

稳定名称为 `exact_sturm_minimum_v1`。接口 `minimum(p,tolerance,max_leaves)`
将原多项式精确取负，调用

```text
exact_extrema.maximize(-p, interval=(-1,1), tolerance=tolerance,
                      max_nodes=max_leaves-1).
```

由于 `min p = -max(-p)`，若最大值证书给出 `[L,U]`，最小值区间为
`[-U,-L]`。验证器同时要求：

1. `exact_extrema.verify(proof)` 通过；
2. 原始绑定多项式恰为规范化的 `-p`；
3. 原始区间恰为 `[-1,1]`。

特别地，即使某个其他区间的独立 Sturm 证书自身有效，也不能用于这里的
全区间下界。平方坐标约化属于 `exact_extrema` 内部；桥接始终绑定其
**原多项式、原区间**，不拿工作坐标代替原变量。

叶预算为 1 到 4097 的整数，与 0 到 4096 个新增分区边界对应。只有系数、
区间和容差精确验证后才接受。一个合法 `open` 范围仍给出保守最小值区间，
但不声称达到所请求精度。

`positive_search` 的可插拔后端协议直接接受此对象：生成及后续验证都必须
显式提供同一受信任后端。默认旧验证器不会接受新名称。测试已覆盖下界记录
和完整的 `optimization_certificate`，不只是孤立地调用后端方法。

## 2. precise_family：同一 Poisson 恒等式，新证书格式

```python
from precision_bridge import precise_family, verify_family

certificate = precise_family([[0,1], [0,0,1]], [[2,"1/2"],["1/2",1]],
                             tolerance="1/100000000000000", max_nodes=64)
assert verify_family(certificate)
```

`precise_family(directions,weight,tolerance,max_nodes)` 调用原
`family_v29.assemble`，因而保留同一势输入限制、Poisson 解、平均值和正定
权重检查。记各 Poisson 解为 `r_i`，导数向量为 `v(t)`，对给定 SPD 权重 `W`
精确构造

```text
g_W(t) = (1-t²) v(t)^T W^-1 v(t).
```

求其最大值区间 `[δ_lower,δ_upper]`，并严格绑定 `g_W` 和 `[-1,1]`。
令 `G=δ_upper W`。由加权 Cauchy–Schwarz，对任意幅度向量 `a`，

```text
(1-t²)(a^T v(t))² ≤ g_W(t) (a^T W a) ≤ a^T G a.
```

所以 `G` 是完整区间上的 Poisson 矩阵包络。代入既有指数正函数的局部能量
恒等式，得到完整单位球面、所有实幅度向量的谱下界
`λ₁(-Δ+sum a_i p_i) ≥ sum a_i mean(p_i)-a^T G a`。

新证书 `method='precise_poisson_matrix_v88'` 包含：

- 原方向、权重、Poisson 解和平均值；
- `envelope_polynomial` 与原区间；
- `scale_lower`、`scale_upper`、`scale_error_bound`；
- `scale=scale_upper` 和 `quadratic_bound=scale_upper*W`；
- 原生 `range_proof`、`range_status` 和预算；
- 固定 Poisson ansatz 的明确范围及未形式化声明。

`verify_family` 从输入重算数学式、验证新范围并重构规范字段。它不调用
`precise_family` 或 `exact_extrema.maximize` 搜索。旧 `family.verify` 和
`family.verify_range` 按设计拒绝这种新格式。

常数方向允许得到 `g=0,δ=0`；单独 family 证书无需方向模常数线性独立。
矩阵优化的 seed 仍受旧 `dual.setup` 的独立性条件约束。

## 3. precise_matrix：保留旧 dual，重新包围同一最佳权重

```python
from precision_bridge import precise_matrix, verify_matrix

# seed_result 为 matrix_refine 的结果，或其中的旧 global_matrix_envelope_v33。
certificate = precise_matrix(seed_result, tolerance="1/1000000000000", max_nodes=64)
assert verify_matrix(certificate)
```

只读取 wrapper 的数学 `certificate` 并用旧 `dual_v33.verify` 验证。
搜索统计、迭代标志、耗时、warm samples 和外层声称的方向/成本都没有证明
权威。输入在生成时深复制，不因返回证书被修改而改变原 seed。

记旧 seed 中已验证的保留权重为 `W`，成本为 SPD 矩阵 `C`，旧 dual 下界为
`D`，旧上界为 `U_old`。该接口不再搜索权重，不更改 dual atoms，保持原
`D` 及整个 dual 证书。

先计算严格正的 `h=trace(CW)`，再把目标精度换算为范围容差

```text
epsilon_range = tolerance/h.
```

调用 `precise_family` 得到新候选 `G_new=δ_upper W`，从而
`U_new=trace(C G_new)=h δ_upper`。候选的固定权重范围误差不超过
`h(δ_upper-δ_lower)`，记录为 `range_objective_error_bound`。
这只是固定权重的范围误差，**不是**矩阵优化的全局间隙。

有限预算可能使新上界比旧上界宽，因此最终保留

```text
U = min(U_old,U_new),
G = 与所保留上界对应的已验证矩阵,
gap = U-D.
```

两者都为同一问题的有效上界，取更小者仍然有效；旧 dual 完整保留，所以
`D≤P*≤U`。其中

```text
P* = inf trace(C G),  G >= (1-t²)v(t)v(t)^T 对全部 -1≤t≤1。
```

全局状态仅由实际 `gap<=tolerance` 决定。新范围已达到精度，不代表原权重
或旧 dual 已把优化间隙闭合。反过来，一个尚未完成隔离的有效范围上界若已
足够证明小间隙，数学上的间隙结论仍然成立；范围自身状态继续独立披露。

## 4. 矩阵格式、绑定及负收益

新方法名是 `precise_global_matrix_v88`，主要字段为：

| 字段 | 含义 |
|---|---|
| `seed_certificate`, `dual` | 原已验证 primal/dual 证书和原 dual；lower 不变 |
| `precise_family` | 对 seed 同一权重的新范围候选 |
| `directions`, `cost`, `weight` | 与两边证据完全相同的原问题数据 |
| `quadratic_bound`, `upper_source` | 实际保留的可行矩阵及其来源 |
| `previous_upper_bound`, `candidate_upper_bound`, `upper_bound` | 旧、新候选、保留上界 |
| `candidate_upper_improvement` | `U_old-U_new`，允许为负 |
| `upper_improvement` | `U_old-min(U_old,U_new)`，一定非负 |
| `precision_outcome` | `improved`、`tied` 或 `negative_gain_seed_retained` |
| `range_objective_error_bound` | 固定权重范围误差的成本单位转换 |
| `lower_bound`, `gap`, `status` | 同一 Poisson 矩阵优化的下界、间隙及状态 |

`verify_matrix` 独立重放旧证书和新 family，逐项核对方向、成本和**同一权重**，
核对范围容差换算、实际所选矩阵、上下界、收益和 gap。它不重新运行范围或
矩阵搜索。顶层改成本、偷换另一个有效 family、替换权重、只改状态或下界
都不能通过验证。

严格规范类型比较也用于所有新字段，`0` 不能冒充 `false`，布尔值不能冒充
节点预算。旧数学证据仍使用旧验证器，其已定义的信任边界没有被暗中替换。

## 5. 真实改进、负收益和平局

所有结果均通过各自验证器。可重现入口是
`test_precision_bridge.comparison_rows()`，返回全部精确分数，没有计时或
模型调用。下表仅为便于阅读的近似展示；严格界以返回的有理证书为准。

前两行使用 `p(t)=t+t²/3+t⁴/5`、`W=C=[1]`，同一旧 dual 在 `t=0` 放置
质量矩阵 `[1]`。旧范围容差为 `10^-12`，V88 目标容差为 `10^-14`。
不同预算是有意展示增益和失败情形，**不是公平运行时间基准**。

| 场景 | 旧叶预算 | 新节点预算 / 实耗 | 旧上界约值 | 新候选约值 | 保留上界约值 | 候选收益约值 |
|---|---:|---:|---:|---:|---:|---:|
| 粗旧证书得到改进 | 1 | 64 / 9 | 0.31247372134039 | 0.27193749048181 | 0.27193749048181 | +0.040536230858575 |
| 紧旧证书，零新预算 | 128 | 0 / 0 | 0.27193749048203 | 0.31205720164609 | 0.27193749048203 | −0.040119711164060 |
| 已精确的 `p=t` | 1 | 0 / 0 | 0.25 | 0.25 | 0.25 | 0 |

第二行明确返回 `negative_gain_seed_retained`，原有有效上界保持不变，不把
候选失败包装成精度提升。第三行零 gap 只表示这一 Poisson 包络问题已精确
解决，不意味着原线性势的最佳谱二次常数为 `1/4`。

另一个真实 `matrix_refine.refine(...,budget=0)` wrapper 集成测试使用
方向 `[t+t²/3, t²+t³/5]` 和成本 `[[2,1/2],[1/2,1]]`，V88 在 64 节点预算、
`10^-12` 目标下取得严格正的上界改善，保留原 dual，并完成范围误差的成本
换算。此检查不调用浮点矩阵候选求解循环。

## 6. 验证记录和理论边界

`python3 -m unittest test_precision_bridge -v`：17 项测试通过。覆盖正函数
后端完整接入、反号/原区间绑定、open 范围、Poisson 公式重建、新旧格式
隔离、常数方向、正定权重检查、实际改善、负收益回退、真实 wrapper、成本
换算、旧 dual 不变、禁止重跑合成器的重放，以及方向/成本/权重/范围/收益/
scope 篡改。

理论基础是已有精确极值范围、Poisson 恒等式、加权 Cauchy–Schwarz 和旧
矩阵弱对偶。不新增 log-Sobolev 等外部分析定理。界宽更小只减少范围松弛；
不会自动优化权重、改进 dual、解除有限 Poisson ansatz 限制，或证明一般
几何谱不等式的锐常数。所有新证书均明确保留这些范围和限制。
