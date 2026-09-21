# V89：把全实谱引理接回原始结构化目标

V89 是一个功能：接收 `compiler_v23` 的原始结构化目标，把经过验证的全实参数线性势谱引理用于该目标，再对**原参数域**上的导出二次式求全局最小值。代码位于 `original_goal.py`，测试位于 `test_original_goal.py`。没有把引理本身的证明状态直接当成原目标的证明状态。

## 1. 输入模型与明确的适用范围

原目标为

\[
\inf_{a\in D}\left[\lambda_1\left(-\Delta+q_0(t)+\sum_{i=1}^d a_i p_i(t)\right)
+Q(a)\right]\ge\tau,
\qquad Q(a)=c+\ell^Ta+\tfrac12a^THa.
\]

仅支持单位二球面、第一特征值（包含最低模态）、`full_sphere` 或 `axisymmetric` 两种声明的函数空间。`D` 是原始全实域或完整有理参数盒，参数数目 `0≤d≤4`。固定势 `q₀` 可以是次数至多六的多项式，每个参数方向必须为

\[
p_i(t)=\beta_i+\gamma_i t.
\]

编译器仍使用原来的受限表达式语言；不会推断任意自然语言或更一般几何。非仿射参数方向、非线性参数势、其他几何/谱编号、未知域以及未被接受的引理系数均明确报 `ValueError`。

证书完整保存 `compiled_goal`，包括原始 `source`、规范化的算子/参数域/代价和 `threshold`。复核首先重新编译原始表达式并比较规范结果。另有明确的 `operator_binding` 保存几何、非负 `−Δ`、函数空间、谱编号和测度归一化。

均匀球面概率测度与面积测度给出相同 Rayleigh 商。对坐标多项式，概率平均为 `∫₋₁¹q(t)dt/2`。以下常数试探的质量因此是 `1`，不把未归一化的面积积分混入谱值。

## 2. 两条全实参数谱引理分别验证

固定选择归一化方向 `p=t` 的全实引理

\[
\lambda_1(-\Delta+zt)\ge-Cz^2\qquad\forall z\in\mathbb R.
\]

`rule='log_sobolev'` 默认 `C=1/6`，调用并复核 `uniform_refine.synthesize_log_sobolev()`。它依赖明确披露的经典球面对数 Sobolev 定理：在单位 `S²` 的概率测度下 `Ent(|u|²)≤∫|∇u|²`；常数和归一化见 [Dolbeault、Esteban、Kowalczyk、Loss，Corollary 2，第 3 页](https://arxiv.org/pdf/1210.1853v1)。这一经典定理并未由本程序从公理证明或形式化。

`rule='barta'` 默认 `C=1/5`，调用并复核三阶正函数的完整 Bernstein 证书和远场拼接。该路径不增加对数 Sobolev 定理作为前提。

两者完整证明见同目录的 `PROOF_UNIFORM_REFINE.md`。V89 必须检查引理方法、`direction=['0','1']`、实际系数 `C` 和成功状态；Barta 路径还必须明确具有 `all_real=True`。已证明的有限盒引理也不能直接在此应用。`ansatz_obstruction`、`unresolved`、谱反例记录以及单独一个声称正确的系数均不能充当成功引理。

完整球面的谱下界也约束轴对称限制空间的谱底，因为后者的 Rayleigh 最小化集合更小。两种空间中均包含下文使用的常数试探函数。因此不需要在接入过程中错误改变原始函数空间，也不需要假定其他谱编号满足同一结论。

## 3. 固定势余项与原参数二次式

令 `b₀=[t]q₀`，删除固定势的一次项，保留其常数和全部高次项：

\[
r(t)=q_0(t)-b_0t,\qquad
z(a)=b_0+\gamma^Ta.
\]

原势于是准确分解为

\[
q(a,t)=r(t)+\beta^Ta+z(a)t.
\]

使用 `family_v29.maximum(-r)` 的完整范围证书，得到 `U≥maxₜ(−r)`，定义 `r_min=−U`。复核必须重新验证范围证书，并确认其多项式恰为**实际余项的负数**。范围预算未达到目标精度时，`U` 仍是有效全区间上界；该状态保留在证书中，不能把采样最大值当成 `U`。

点态势比较和全实谱引理给出

\[
\lambda_1(-\Delta+q(a,t))
\ge r_{\min}+\beta^Ta-C(b_0+\gamma^Ta)^2.
\]

加上原始 `Q` 得到参数二次下界

\[
F(a)\ge\widetilde c+\widetilde\ell^Ta+\tfrac12a^T\widetilde H a,
\]

其中

\[
\widetilde c=c+r_{\min}-Cb_0^2,\qquad
\widetilde\ell=\ell+\beta-2Cb_0\gamma,\qquad
\widetilde H=H-2C\gamma\gamma^T.
\]

证书保留 `amplitude_offset`、`amplitude_slopes`、`parameter_constants`、余项、范围证据及 `derived_quadratic`。验证器从原编译结果重算全部系数，再调用原有 `algebra.quadratic_min` 在**原域 D** 上求最小值；不放大域后又声称得到了原域的精确极值。

全实域中，二次式有有限最小值的充要条件是 `H̃` 半正定且 `−ℓ̃` 属于它的值域，此时任一精确驻点给出全局最小值。在参数盒中，程序枚举所有面上的端点和可解驻点，允许原二次式非凸。若一个最低值位于奇异自由 Hessian 的面内部，可沿零方向保持函数值不变地移至该面的边界，递推到更小的面；所以只枚举非奇异自由块不会遗漏最低值。这是已有有理二次最小化程序的数学依据。

如果全实域上的导出二次式向下无界，证书保存 `lower_bound=None` 和原因（负 Hessian 方向，或核空间中的非零线性项）。这只意味着**本下界松弛没有有限下界**，不意味着原始谱目标向下无界或被推翻。

## 4. 上界必须来自原始势的真实试探函数

为判断原目标是否被反例推翻，以及在可能时夹出精确下确界，V89 在原参数域中选择有限个有理参数候选。候选可来自上述松弛的最小点、零点/裁剪零点、盒中心和角点、或常数试探二次式的最小点。有限采样只负责提出上界见证，**不负责证明全域下界**。

对每个候选 `a*`，使用原问题中的 `u(t)=1`。在声明的两种函数空间中它均为合法试探，且

\[
\int u^2d\mu=1,\qquad \int|\nabla u|^2d\mu=0,
\]

从而

\[
\inf_{a\in D}F(a)
\le F(a_*)
\le \overline{q_0}+\beta^Ta_*+Q(a_*)=:U_*.
\]

验证器在原参数坐标重组完整的 `q(a*,t)`，精确计算每个偶次系数的平均 `qⱼ/(j+1)`，再计算原 `Q(a*)`。记录 `constant_trial` 保存质量、动能、原势、势平均、原代价以及最终目标上界。证书没有使用正函数局部能量上界、参数松弛的函数值或任意声称的数作为谱上界。

全部候选都逐点检查属于原域。最佳候选上界记为 `U`，有限二次下界存在时记为 `L`，于是 `L≤inf_D F≤U`。只有 `L=U` 时才标记 `infimum_exact=True`。

## 5. 原有 threshold 和状态语义

沿用旧目标的关系 `>=` 与原始阈值 `τ`：

- 有有限下界且 `L≥τ`：`proved`。
- 尚未证明，但某个真实原问题试探满足 `U<τ`：`refuted`。
- 其他情况：`unresolved`。

反例条件严格为 `<`。`U=τ` 不推翻 `>=τ`。导出二次下界无界时也仍可由一个合法试探给出真实反例；若没有这样的见证，则明确保持未解决。

例如 `q=at,Q=a²/7,τ=0` 的导出二次式向下无界，而常数试探在零点给出 `U=0`；V89 因此返回 `unresolved`。虽然另一条专门的非恒定试探已经可以推翻该目标，本次有界集成功能没有把那个试探暗中替换进来。相反 `Q=−a²` 时，在 `a=±1,u=1` 上原目标上界为 `−1`，可以真实返回 `refuted`。

## 6. 原 `a²/5` 目标现在得到精确下确界

原始结构化输入可为

```json
{
  "geometry": "unit_sphere",
  "function_space": "full_sphere",
  "eigenvalue_index": 1,
  "parameters": ["a"],
  "potential": "a*t",
  "penalty": "a**2/5",
  "domain": {"kind": "all_real"},
  "threshold": "0"
}
```

默认 `C=1/6` 给出导出二次式 `a²/30`，对应 Hessian `1/15`。其全实最小值为零；原问题在 `a=0,u=1` 上的 Rayleigh 商和原代价也都是零。因此

\[
\boxed{\inf_{a\in\mathbb R}\bigl[\lambda_1(-\Delta+at)+a^2/5\bigr]=0.}
\]

原盒 `[-2,2]` 同样得到下/上界均为零。独立的 Barta 路径 `C=1/5` 导出恒零下界，也由相同的原问题见证夹成精确零。两种路径的经典定理依赖分别披露。

## 7. 接口、核验和实际范围

```python
import original_goal

certificate = original_goal.solve(raw, coefficient=None, rule='log_sobolev')
assert original_goal.verify(certificate)
assert certificate['status'] == 'proved'

independent = original_goal.solve(raw, rule='barta')
assert original_goal.verify(independent)
```

返回对象自身就是 `method='original_spectral_goal_v89'` 的证书，不额外包裹 `certificate` 字段。`verify()` 返回真表示记录及其 `proved/refuted/unresolved` 状态诚实；接受原不等式还必须检查 `status=='proved'`。所有数值接受条件都使用有理数；用户可指定另一有理 `coefficient`，但该系数必须在所选规则下有成功的全实引理。

`test_original_goal.py` 的 22 项功能测试全部通过，覆盖原题精确零、两条独立引理、完整盒域、函数空间、固定高次势余项、仿射多参数混合、秩亏与被动参数、零参数、非凸盒二次式、真正反例、诚实松弛失败，以及原目标/算子/域/阈值/引理/余项/二次式/试探函数的篡改拒绝。

复核会调用经过验证的既有编译器、范围证书、有理矩阵和二次最小化实现，因而仍依赖这些代码与文中分析论证的正确性；共享基础运算不是完全独立内核。证书始终披露 `formal_assistant_checked=False`，对数 Sobolev 经典定理也明确标注 `proved_by_this_program=False`。
