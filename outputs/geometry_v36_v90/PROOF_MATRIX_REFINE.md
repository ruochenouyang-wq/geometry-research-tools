# 矩阵包络精度续算：八项组合改进与证据边界

本模块保持 V33 的数学对象与验证器不变。给定 1–4 个模常数线性独立的多项式方向，精确求解 Poisson 方程并记其导数向量为 \(v(t)\)。目标是

\[
P_C=\inf_G\operatorname{tr}(CG),\qquad
G\succeq A(t)=(1-t^2)v(t)v(t)^T\quad(-1\le t\le1),\qquad C\succ0.
\]

每个返回结果中的 `certificate` 都是原有 `dual_v33.certificate`，可以直接交给 `dual_v33.verify`。本文的浮点搜索、预条件信息、轮数、warm cache、误差分配和收敛标记不参与数学证据接受。只有经 `family_v29` 完整区间范围证书重放的上界，以及经 `dual_v33` 精确 PSD 检查的下界可以更新最终结果。

以下 MR 标签是内部研发标签，不是新的既定版本编号。它们对应不同算法步骤，可以在其它入口复用；没有把同一证明重复计成多个精度成果。

## MR1：精确坐标预条件与成本协变

对任意可逆有理矩阵 \(T\)，定义

\[
\widetilde v=Tv,\quad H=TGT^T,\quad
\widetilde C=T^{-T}CT^{-1}.
\]

于是 \(H\succeq (1-t^2)\widetilde v\widetilde v^T\) 当且仅当 \(G\succeq A(t)\)，且

\[
\operatorname{tr}(\widetilde C H)=\operatorname{tr}(CG).
\]

若内部对偶原子为 \(Z_k\)，拉回原坐标时必须使用 \(Y_k=T^TZ_kT\)。这样 \(\sum Z_k\preceq\widetilde C\) 推出 \(\sum Y_k\preceq C\)，且对偶目标也不变。不能将成本错误地变成 \(TCT^T\)，也不能在更换单位后继续沿用默认单位成本。

实现把各 \(v_i\) 的系数组成矩阵 \(B\)，对正定 \(BB^T=LDL^T\) 做精确 LDL 消元，再使用 \(T=SL^{-1}\)，其中 \(S\) 的对角元为有理的二次幂，近似 \(D^{-1/2}\)。这消除了系数之间的精确线性相关尺度，并避免把平方根写进证书。内部变换后的多项式只是搜索数据，不重新套用原题输入解析器的系数分母限制；最终证书始终绑定原始合法方向。

浮点阶段进一步使用 \(\widetilde C/\rho\)，\(\rho=\operatorname{tr}(\widetilde C)>0\)（极小值有数值尺度下限）。所有 barrier 误差都除以同一个 \(\rho\)，生成对偶原子时再乘回 \(\rho\)。成本归一化不改变原题目标。

## MR2：有效有限约束与满预算后的分离点替换

有限提议阶段，若精确检查得到 \(A(s)-A(t)\succeq0\)，则 \(H\succeq A(s)\) 已蕴含 \(H\succeq A(t)\)，后者可以删除。完全相同的 \(A(t)\) 先合并，再做 PSD 支配检查，避免互相删除所有相同约束。端点的零矩阵常常因此被移除；标量方向 \(p=t\) 的九个初始点精确压缩为 \(t=0\) 一个有效约束。

全域上界产生新的最差位置后，有限样本池即使已满，也必须接纳新分离点。实现保留新提议，再用现有对偶目标贡献排序旧点、移除贡献较小的点。保留的上界/下界证书不会随有限样本删除而消失；它们的有效性不依赖下一轮搜索用哪些点。

分离点来自完整范围证书中的 witness，以及多项式导数的浮点符号变化括区。\(t^2\) 坐标的平方根和浮点根均量化为有理点并加入正负两侧。这些根只是新有限约束，不是精确极值证明。偶重根或非常靠近的驻点可能被提议网格遗漏；最终上界仍由覆盖整个区间的旧范围证书保证。

## MR3：可续算的 barrier 候选搜索

新搜索器仅复用 `float_sdp` 的小矩阵 Cholesky/线性求解核，不修改旧搜索器。每轮把当前最好全域可行 \(G\) 变换到内部坐标作为起点，再添加正定 ridge。一个已接近最优的继承点通常在半正定锥边界；仅在最终 barrier 的量级上添加 ridge 会使 Newton Hessian 严重病态。实现使用较大的内部起始 ridge，并在进入下一轮时继承 barrier 的尺度，再逐级减小。

候选浮点矩阵量化后，以精确惯性检查全部有效有限约束；必要时按二次幂增加有理对角修正。失败时不接受该矩阵，不把浮点 solver 的 `reached_final_barrier` 当成全局误差证据。弱候选、失败候选与有限 Newton 预算都可正常返回 `global_gap_open`。

## MR4：经验证的上界、下界独立继承

`seed_result` 可以是已有结果字典，也可以是 V33 证书本身。首先重放 `dual_v33.verify`，再精确比较原方向和 \(C\)。不同问题或篡改证书会被拒绝。

各轮独立维护

\[
U_{k+1}=\min(U_k,U_{\rm new}),\qquad
L_{k+1}=\max(L_k,L_{\rm new}).
\]

上、下界不必来自同一轮。下一次续算即使使用不同容差、预算耗尽、候选失败，仍保持 \(L_{\rm new}\ge L_{\rm seed}\)、\(U_{\rm new}\le U_{\rm seed}\)。`budget=0` 允许只重放并返回原有界。缓存样本不可信；格式损坏或超出区间的缓存会被忽略，缓存没有证书地位。

## MR5：按原目标尺度分配误差

设目标绝对容差为 \(\varepsilon\)。实现分别分配有限搜索、对偶修复、量化、稀疏化和连续范围误差预算。初期 finite barrier 允许稍大的目标误差，随现有 gap 和轮数缩紧，以免把大量工作花在尚未加入关键约束的有限问题上。

若候选 \(W\) 的原目标为 \(s=\operatorname{tr}(CW)\)，而比例函数

\[
q(t)=(1-t^2)v(t)^TW^{-1}v(t)
\]

的范围上界多出了 \(\eta\)，相应矩阵目标上界最多额外增加 \(s\eta\)。因此使用与 \(\varepsilon/s\) 成比例的范围误差，而不是与坐标尺度无关的固定误差。

这些预算是算法工作目标，并非“只要输入该容差就一定达到”的承诺。旧范围器有 \(10^{-12}\) 的容差下限与最多 256 个叶子的预算；有理分母、Newton 轮数和样本数也有限。实际达到的唯一总误差是证书中的精确 `upper_bound-lower_bound`。任何预算失配都不能通过伪造状态变成闭合 gap。

## MR6：目标误差控制的双向对偶修复与余量分配

给定 PSD 提议 \(\widetilde Y_k\)，记 \(S=\sum\widetilde Y_k\)、\(\ell=\sum\operatorname{tr}(\widetilde Y_kA(t_k))\ge0\)。寻找最大公共倍数 \(\beta\ge0\)，使 \(C-\beta S\succeq0\)。相比旧修复，本实现既允许缩小，也允许放大原子。

若 \(S\ne0\)，迹给出精确上界

\[
\beta\le\frac{\operatorname{tr}C}{\operatorname{tr}S}.
\]

这个初始括区在全部提议原子共同放大 \(10^{100}\) 时也自动缩小，避免从固定区间盲目二分。对每个二分点直接检查有理矩阵惯性。停止条件是

\[
(\beta_{\rm hi}-\beta_{\rm lo})\ell\le\varepsilon_{\rm repair},
\]

因而控制的是原目标损失。若 256 次上限后仍不满足，就显式报提议失败，不能暗称达到误差预算。

修复后余量 \(R=C-\sum Y_k\succeq0\)。用精确 LDL 写成 \(R=\sum_j d_j u_ju_j^T\)，把每个 PSD 分量分配到候选样本中使其 \(\operatorname{tr}(d_ju_ju_j^TA(t))\) 最大的位置。增加项均非负，且新原子和恰好为 \(C\)，所以此步骤不会降低对偶下界。

标量例子：\(p=t,C=1\)，提议原子 \(Y(0)=1/100\)。旧修复保留很小的下界 \(1/400\)；双向修复加余量分配得到精确 \(1/4\)。这是一项实际数值改进，不需要浮点最优性断言。

## MR7：保持 PSD 的有理重构

对一个接近秩亏的 PSD 矩阵逐元素取整可能把它变成不定矩阵。本实现先精确分解

\[
Y=\sum_j d_j u_ju_j^T,
\]

再把 \(\sqrt{d_j}u_j\) 量化为有理向量 \(z_j\)，返回 \(\widehat Y=\sum_jz_jz_j^T\)。无论平方根和取整的浮点误差有多大，最终 \(\widehat Y\) 都由有理外积构成，PSD 是精确的。重构可能影响质量，之后仍必须做 MR6 的精确成本修复。

候选的共同分母由原目标量化预算决定，并受显式上限限制。分母预算是搜索参数；不因设定了分母就声称误差已达到，最终仍比较 V33 的确切 gap。

## MR8：保目标的对偶合并与可核验的稀疏化

若 \(A(t_i)=A(t_j)\)，则把两个原子合成 \(Y_i+Y_j\) 精确保留矩阵和与目标。它适用于同点，以及部分奇偶对称情形中不同点得到同一个 \(A\) 的情况。

进一步稀疏化时只删除 PSD 原子，删除的目标贡献逐项用有理数累计，并要求不超过给定误差。删除 PSD 原子保持 \(\sum Y_k\preceq C\)，因此删除后的证明仍有效；随后将产生的余量重新分配。主循环只有在重分配后的下界不低于原下界时才采用稀疏结果，不会为了更少原子牺牲当前最好界。

## API、运行与轻量实测

入口：

```python
result = matrix_refine.refine(
    directions, cost=None, tolerance=Fraction(1, 10**6),
    seed_result=None, budget={'rounds': 8, 'range_leaves': 128},
    range_evaluator=family_v29.maximum,
)
assert matrix_refine.verify(result)
assert dual_v33.verify(result['certificate'])
```

`budget` 也可以直接传轮数。字典支持 `rounds`、`newton_steps`、`range_leaves`、`max_samples`、`max_denominator`，以及 `precondition`、`dominance`、`sparsify` 开关；非法键/范围会被拒绝。命令行支持 `--directions`、`--cost`、`--tolerance`、`--rounds`、`--seed`，输出 JSON。

可替换范围器的契约是 `evaluator(q, epsilon, max_leaves)` 返回原 `family_v29` 格式的范围证书。`family.build` 必须验证该证书及其多项式绑定。当前模块不接受新的独立极值证书格式；如果另一个模块使用不同证明格式，需要外层统一验证适配，不能直接跳过旧 verifier。

执行：

```text
python3 -m unittest test_matrix_refine -v
python3 matrix_refine.py --tolerance 1/1000000 --rounds 8
```

这些是功能和回归测试，未开展长时间计时基准。测试覆盖非对称坐标矩阵、错误成本变换、非对角成本、四方向、同矩阵约束压缩、PSD 取整、极大对偶尺度、弱/失败候选、非法范围证书、seed 错绑、精确界单调继承和输入分母边界。实测代表结果如下；数字只说明这些有限案例，没有形成普遍精度/速度保证。

| 案例 | 预算与容差 | 已验证结果 |
|---|---|---|
| \(p=t,C=1\) | 零轮 | 精确 gap 为 0 |
| \(p=(t,t^2),C=I\) | 最多八轮，\(\varepsilon=10^{-6}\) | gap 闭合 |
| \(p=(t,t^2),C=\left[\begin{smallmatrix}2&1/2\\1/2&1\end{smallmatrix}\right]\) | 最多八轮，\(\varepsilon=10^{-6}\) | gap 闭合 |
| \(p=(t,t^2,t^3,t^4),C=I\) | 最多八轮，\(\varepsilon=10^{-5}\) | gap 闭合 |
| \(p=(1000t,t^2/1000),C=\operatorname{diag}(10^{-6},10^6)\) | 两个方法同为最多五轮，\(\varepsilon=10^{-5}\) | 新方法 gap \(\le10^{-5}\)，且小于 V34 gap 的百分之一 |

最后一个例子与 \((t,t^2),C=I\) 表示同一变换后的目标。比较的是相同原目标、相同绝对容差和相同交换轮数；它不声称两种算法的每轮工作量相同。

## 结论范围

最终 `global_gap_closed` 只表示 \(U-L\le\varepsilon\) 对上述固定 Poisson 指数构造中的矩阵目标成立。即使得到 \(U=L\)，也不能推出原谱不等式锐利、原问题已证明，或所有可能正函数构造都已最优。连续范围证明和 V33 的弱对偶论证没有改变；这些改进只提高提出好候选、保留有效证据和利用有限预算的能力。

验证使用标准库有理数与旧精确代数实现，部分代数代码为生成器与验证器共享；本模块没有宣称 Lean/Isabelle 内核形式化验证。
