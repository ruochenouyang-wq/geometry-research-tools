# matrix_refine 独立数学审查

审查日期：2026-09-11。范围为 `matrix_refine.py`，对照 `dual_v33.py`、`optimize_v34.py`、`family_v29.py` 与 `matrix_v30.py`。本审查不修改实现或其测试，且未运行长计时基准。

## 当前结论

未发现把变换后问题的上下界冒充原问题上下界，或把固定 Poisson ansatz 的包络最优性冒充原谱问题全局最优的路径。最终证书使用**原始 directions 与 cost**重新构造，交给已有有理验证器重放。

审查复现了两项边界问题：有效有理原题的预条件坐标触发旧输入预算，导致合法 seed 无法零轮重放；对偶修复迭代耗尽后仍无条件返回，未必达到其声明的绝对误差。两者均已向主任务与实现者报告，**实现者修复后，独立的 13 项针对性复验全部通过**。后者原先仍返回有效下界，不能据此指称数学证明失效。

## 原问题与坐标恒等式

被证书界定的优化问题为

\[
P(C,v)=\inf_G\operatorname{tr}(CG),\qquad
G\succeq A(t)=(1-t^2)v(t)v(t)^T\quad\forall t\in[-1,1].
\]

这里 `v` 是给定方向的 Poisson 解的导数。方向含有常数项不影响其线性变换的正确性；原方向和成本仍由最终证书完整绑定。

对任意可逆有理矩阵 `T`（**无须对称**）：

| 量 | 变换 | 回到原坐标 |
|---|---|---|
| Poisson 导数 | `v_new = T v` | `v = T^-1 v_new` |
| 全区间约束 | `A_new(t) = T A(t) T^T` | `A(t) = T^-1 A_new(t) T^-T` |
| 原始变量 | `H = T G T^T` | `G = T^-1 H T^-T` |
| 成本 | `C_new = T^-T C T^-1` | `tr(C_new H) = tr(CG)` |
| 对偶原子 | `Z = T^-T Y T^-1` | `Y = T^T Z T` |

于是 `H-A_new(t) = T(G-A(t))T^T`，半正定性双向等价；且 `tr(ZA_new)=tr(YA)`，`sum Z <= C_new` 回拉后等价于 `sum Y <= C`。

代码中的 `congruence(transpose(inverse), c)`、`congruence(inverse, g)` 与 `congruence(transpose(transform), y)` 分别符合上述成本、原始变量和对偶原子的变换。已用非对称有理矩阵

```text
T = [[2, 1/3], [-1, 3]]
C = [[3, 1/2], [1/2, 2]]
directions = [[2/3, 1, 1], [-1/5, 2, -1]]
```

独立检验逆变换、trace 恒等式，以及 `t=-1,-2/3,0,1/7,1` 的约束和对偶目标恒等式，全部精确相等。

目标额外归一化为 `C_normalized=C_new/s` 时，数值求解器生成的对偶原子须乘 `s` 后再回拉。代码采用 `rational_mu=actual_mu*objective_scale`，符合该要求；原始上界始终由原成本 `C` 与回拉后的全区间模板计算，没有把归一化目标值当成原目标值。

## 对偶与原始证书的边界

- **对偶修复：** 每个有理原子必须 PSD；共同乘子最终通过精确惯性检验保证 `C-sum(Y_k)` PSD。原子合并仅对完全相同的 `A(t)` 进行，保持矩阵总和与目标值。PSD 松弛的精确 LDL 分解经过奇异与非奇异样例重建检验；补回松弛使原子和等于原成本，并增加非负目标贡献。
- **稀疏化：** 删除 PSD 原子保留对偶可行性；删除目标贡献用分数累计。旧较好下界只在新下界不低时替换，不会因稀疏化丢失已验证的最佳下界。
- **数值 proposal：** 仅提供搜索候选。有限样本可行性、障碍参数、停止标记或候选质量均不能独立成为全区间上界。候选权重回到原坐标后仍需 `family.build` 核验原多项式及整个 `[-1,1]` 的证据。
- **不可信区间 oracle：** 实际注入了返回另一多项式合法区间证书的 oracle。`family.build` 拒绝错绑；初始化回退到原有 evaluator，最终证书仍有效。
- **分离点和有限约束裁剪：** 只改变下一轮搜索，不删除保留的全区间原始证书或有限原子对偶证书。上下界的适用域因此仍是整个区间。

## seed 与结果绑定

`dual.verify(seed)` 单独只能验证 seed 内部自洽。续算还需要核对 seed 的原方向与成本等于当前请求；实现确实进行了这两项额外核对。独立检验了修改方向系数、修改方向常数项、修改成本以及直接篡改 seed 成本的情况，均被拒绝。

零轮续算在同题上保留 seed 的上下界，并以当前请求 tolerance 重新计算 gap 状态。proposal 失败的注入对照也保留原 seed 上下界。缓存采样点为不可信 proposal，越界或无效字符串不会成为证明证据。

`verify(result)` 明确只验证数学证书。`conditioning`、`warm_state`、迭代日志及统计元数据不是证书的一部分。外部使用者应从 `certificate.dual.directions`、`certificate.dual.cost` 读取证明绑定的题目，不能依靠可任意改写的 conditioning 元数据推断另一个目标。

## 复现的问题

### R1：合法 seed 被内部预条件坐标的输入预算拒绝

以下方向均满足现有原题输入约束，但初版默认预条件化后系数分母可能超过原题输入限额 `10^9`。`transform_problem` 再对内部坐标调用 `dual.setup`，在检查 seed 之前拒绝了它。

```python
from fractions import Fraction as F
import matrix_refine as m
import matrix_v30
import dual_v33 as dual

dirs = [[0, F(1,97), F(1,89), F(1,83)],
        [0, F(1,79), F(1,73), F(1,71)]]
cost = [[1, 0], [0, 1]]
template = matrix_v30.synthesize(dirs, max_leaves=4)
lower = dual.lower_certificate(dirs, cost, [
    {'t': '0', 'matrix': [['0', '0'], ['0', '0']]}])
seed = {'certificate': dual.certificate(template, lower, F(1,1000))}
assert dual.verify(seed['certificate'])
m.refine(dirs, cost, seed_result=seed, budget=0)
```

初版抛出 `ValueError: 势函数系数绝对值最多 1000，分母最多 10^9。`；无 seed 也触发。另一独立样例为 `[[0,1,1/999983],[0,1/999979,1]]`。

影响：合法输入/合法 seed 的可用性，尤其违背零轮安全重放的约定；没有产生无效数学结论。应将内部搜索坐标与原题输入预算分开，或安全退回恒等预条件。

### R2：对偶修复二分耗尽后没有验证目标误差

```python
from fractions import Fraction as F
import matrix_refine as m
import dual_v33 as dual

r = m.repair_dual([[0,1]], [[1]],
    [{'t': '0', 'matrix': [[10**100]]}],
    error=F(1,10**10), fill=False)
assert dual.verify_lower(r)
print(r['lower_bound'])
```

初版固定进行最多 256 次二分后返回 `0`；同一原子共同缩放到 `Y=1` 可得精确下界 `1/4`。因此它虽然保留对偶可行性，却没有达到所声明的 `1e-10` 目标误差。旧 `dual.verify_lower` 合理地只核验可行性与实际下界，无法发现这个搜索精度契约问题。

默认 `fill=True` 恰会在这个标量样例补满松弛，不能用它证明二分已达到误差；一般多维松弛补回也不保证复原原子最优比例。应在耗尽时明确失败/报告未达到目标，或以受控方式自适应建立足够窄的搜索区间。

## 检查规模与状态

独立执行的首批 **37 项轻量精确断言全部通过**，另执行不可信 oracle 回退、proposal 故障保留 seed、seed 成本篡改拒绝等小型对照。首轮实际数值 proposal 仅运行一次：非对称预条件、1 轮、8 Newton 步、4 区间叶片，约 **0.024 秒**，最终旧验证器通过且原题方向/成本保持不变。R1 与 R2 是另外的实际失败复现。

实现修复后的 **13 项针对性独立复验**：

1. R1 两组有理方向，分别带合法 seed 和不带 seed 的零轮调用，共 4 项通过。返回证书可由旧验证器重放，原方向/成本保持不变；带 seed 时上下界不变。
2. R2 原标量 `10**100` 原子分别在 `fill=False/True` 下均返回精确 `1/4`，共 2 项通过。新初始上界 `trace(C)/trace(sum(Y))` 是由 PSD 可行性的必要 trace 不等式得到的正确上界，且对原子的公共缩放不敏感。
3. R2 另用 `cost=diag(1,10**100)`、原子矩阵 `I`、方向 `[t,t*t]` 实际触发 256 步耗尽，新实现明确抛出 `Dual repair objective-error budget not reached`。此 1 项验证误差未达时不会静默返回。
4. 5 种损坏或不合类型的缓存 warm state 被安全忽略，合法 seed 仍可重放。
5. 修复后另执行一次非对称预条件的同规模小型实际 proposal，旧验证器通过，方向/成本仍绑定原题，共 1 项。

最终审查文件：

```text
matrix_refine.py SHA-256
6cfb0e7c25fa819a073564656ebc30675fa83992f0c631cf4c6a370bd2eaa1ef
```

此审查不宣称覆盖任意规模或任意恶意 Python 回调，也不评估数值性能。数学 scope 保持为固定 Poisson ansatz 的矩阵包络全局优化；并非原谱不等式的锐常数证明。
