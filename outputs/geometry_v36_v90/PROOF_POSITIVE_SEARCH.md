# 固定势的正函数凹优化：证明、证书与已复现实验

新增实现为 `positive_search.py`，验收为 `test_positive_search.py`。目标是固定单位二球面、固定有理多项式势下的研究工具，不涉及一般几何问题或开放问题的解决声明。已有 `family_v29.py` 使用固定 Poisson 对数试探；已有 `error_bounds.py` 使用数值 collocation 提议指数函数。本增量允许自由优化指定的对数多项式基，并为该优化本身提供全局误差证书。

## 六个实质增量

| 增量 | 验收对象 | 实际意义 |
|---|---|---|
| 1. 自由对数基与完整区间下界 | `lower_certificate` / `verify_lower` | 系数不再由 Poisson 方程固定；每个接受结果是完整球面的谱下界 |
| 2. 概率原子的精确二次上界 | `ansatz_upper_certificate` / `verify_ansatz_upper` | 给出指定指数族最佳 Barta 下界的全局上界，处理奇异矩阵与无界情况 |
| 3. 驻点和接触误差分解 | `optimization_certificate` / `verify` | 把可验证总误差分为平均梯度未平衡部分与接触/区间部分 |
| 4. 有预算的凹优化提议器 | `search` | 标准库 soft-min Newton、极小值见证交换、有理接受；保留弱提议与失败 |
| 5. Poisson 初始化、奇偶与能量平移 | `poisson_initialization`、`parity_certificate`、`shift_certificate` | 复用旧结果；偶势的对数偶化有精确不减下界证明；能量零点变化不要求重优化 |
| 6. 独立 Rayleigh 谱上界及实证区分 | `rayleigh_certificate` / `verify_rayleigh` | 真实谱区间与函数族优化区间分别给出；含严格改进和精确零误差例 |

所有输入/证书中的系数是整数、分数字符串或 `Fraction`；公共接口拒绝浮点系数。势次数和对数次数均不超过 6，独立基维数不超过 6。局部能量次数最多为 12，不能误截断到 6。概率原子最多 128 个，Bernstein 子区间最多 256 个。没有第三方数值库或模型调用。数值提议输出的候选会重新有理化，之后不信任数值停止标记。

## 1. 完整球面上的下界，而非仅限轴对称试验空间

令

\[
H=-\Delta_{S^2}+q(t),\qquad t=x_3\in[-1,1],\qquad
s_\theta(t)=\sum_{i=1}^d\theta_i b_i(t),\qquad \psi_\theta=e^{s_\theta}.
\]

多项式 $s(t)$ 在球面上光滑，指数严格为正。坐标端点是两极，不是带额外边界条件的区间边界；无需强加 $s'(\pm1)=0$。常数项只改变归一化，程序将每个基的常数项消去，并检查余下各基线性无关。

对轴对称函数，

\[
H_0s=2ts'-(1-t^2)s'',\qquad
E_\theta(t)=\frac{H\psi_\theta}{\psi_\theta}
=q(t)+H_0s_\theta(t)-(1-t^2)(s_\theta'(t))^2.
\]

这里的局部能量下界属于 Barta 型方法；正函数局部能量与本征值的关系可参见 Mouchet 的 Theorem 1 及其参数优化式 $4$。[原论文](https://arxiv.org/pdf/math/0505541)

本模型所需证明也可直接写出。对任意完整球面上的光滑函数 $u=\psi v$，分部积分得到

\[
\int_{S^2}(|\nabla u|^2+q|u|^2)
=\int_{S^2}\psi^2|\nabla v|^2+
  \int_{S^2}\frac{H\psi}{\psi}|u|^2.
\]

若 $E_\theta\ge L$，则右端至少为 $L\|u\|^2$。按变分原理及稠密性，完整球面的基态能量满足 $\lambda_1(H)\ge L$。这一论证允许任意 $u$，无需先证明基态具有轴对称性。

默认区间后端对 $-E_\theta$ 调用已有 `family_v29.maximum`。对于偶多项式使用 $x=t^2\in[0,1]$ 降维。二次多项式直接检查端点和内部驻点；更高次使用有理 Bernstein 系数及完整有序区间覆盖。Bernstein 基非负且和为 1，因此每段系数的最大值控制该段的多项式最大值。验证时重新进行有理仿射展开、检查完整覆盖和见证值，并把范围证明绑定到重建出的局部能量多项式。

证书中的 `spectral_lower` 是谱下界。`candidate_infimum_upper` 只控制该候选的局部能量下确界，二者的差 `range_gap` 是区间极值计算误差；后一个字段不是谱上界。

## 2. 凹性与有限概率原子的全局弱对偶

记 $w(t)=1-t^2$、$v_i(t)=b_i'(t)$、$h_i(t)=H_0b_i(t)$。则

\[
E_\theta(t)=q(t)+h(t)^T\theta-w(t)(v(t)^T\theta)^2,
\qquad \nabla^2_\theta E_\theta(t)=-2w(t)v(t)v(t)^T\preceq0.
\]

因此 $F(\theta)=\inf_t E_\theta(t)$ 凹。可直接验证：对任意 $0\le\alpha\le1$，先用每个 $t$ 处的凹性，再取下确界，可得

\[
F(\alpha\theta+(1-\alpha)\eta)
\ge\alpha F(\theta)+(1-\alpha)F(\eta).
\]

点态下确界及 log-sum-exp 的标准凸分析背景可见作者版 Boyd–Vandenberghe 讲义第 3 章。[作者讲义](https://web.stanford.edu/~boyd/cvxbook/bv_cvxslides.pdf)

定义本模块的优化目标

\[
M(q,b)=\sup_{\theta\in\mathbb R^d}F(\theta).
\]

任取有限个 $t_k\in[-1,1]$、有理权重 $\mu_k\ge0$，且 **精确** 满足 $\sum_k\mu_k=1$。令

\[
c=\sum_k\mu_kq(t_k),\quad
\ell=\sum_k\mu_kh(t_k),\quad
A=\sum_k\mu_kw(t_k)v(t_k)v(t_k)^T\succeq0.
\]

则对每个实参数向量，

\[
F(\theta)\le Q_\mu(\theta)
=\sum_k\mu_kE_\theta(t_k)
=c+\ell^T\theta-\theta^TA\theta.
\]

若精确线性方程 $2Az=\ell$ 有解，完成平方给出

\[
Q_\mu(\theta)=c+\frac12\ell^Tz-(\theta-z)^TA(\theta-z),
\quad U_\mu=c+\frac12\ell^Tz,
\quad M(q,b)\le U_\mu.
\]

程序允许奇异半正定 $A$，通过精确消元检查一致性，并储存一个有理最大化点 $z$。若不一致，存在核方向 $n\in\ker A$ 满足 $\ell^Tn\ne0$，沿一个符号的 $\theta+rn$ 有 $Q_\mu\to+\infty$。这种原子集不能提供有限上界，程序明确拒绝，不使用小量替代零主元。

只需弱对偶，无须宣称所有概率测度上的强对偶、最优测度存在或数值提议已找到最优原子。

最关键的范围区分是

\[
L\le M(q,b)\le U_\mu,\qquad L\le\lambda_1(H).
\]

一般不存在 $\lambda_1(H)\le U_\mu$ 的结论。本模块在每个相关证书中保存该限制，也不把 `ansatz_upper` 自动填入 `spectral_upper`。

## 3. 可验证的一阶条件和误差分解

对接受的候选 $\theta$，程序重建

\[
g_\mu=\ell-2A\theta,\qquad
U_\mu-L=
\underbrace{U_\mu-Q_\mu(\theta)}_{(\theta-z)^TA(\theta-z)\ge0}
+\underbrace{Q_\mu(\theta)-L}_{\ge0}.
\]

它们分别储存为 `averaged_gradient`、`quadratic_stationarity_gap` 与 `average_contact_and_range_gap`。第一部分刻画平均一阶条件的不足；第二部分同时含有限原子的接触不足及范围误差。不能把第二部分全称为 Bernstein 误差。

当 $g_\mu=0$、所有正权重原子都接触候选的真正极小值、区间下界精确时，总误差为零，便有全局最优性。一般情形只需 `ansatz_upper-ansatz_lower <= tolerance` 才标记 `ansatz_gap_closed`；这表示规定精度内的整个参数空间最优值包围，不表示获得精确最优参数，也不表示谱误差已达到同样精度。

## 4. 提议器、初始化、对称性和范围后端

数值提议优化有限点上的光滑凹函数

\[
G_\tau(\theta)=-\tau\log\sum_k\exp(-E_\theta(t_k)/\tau).
\]

在有限网格上，

\[
\min_k E_\theta(t_k)-\tau\log N
\le G_\tau(\theta)\le\min_k E_\theta(t_k).
\]

这些关系只涉及有限网格，不是完整区间下界。实现使用逐步降低温度的阻尼 Newton；soft-min 权重只提议概率原子。候选参数转为分母 $2^{32}$ 的有理数；概率权重取非负整数分子后精确归一化，并给每个网格点至少一个单位分子。后者避免网格 Gram 矩阵因舍入丢失全部必要内点。最后仍重新精确检查。

每轮候选经完整区间验收，并把极小值见证附近的有理点加入下一轮网格。偶多项式在 $x=t^2$ 的见证经平方根数值提议后重新有理化；该平方根没有进入证明。外部后端若没有默认见证结构，仍可产生可信下界，但当前搜索不从它自动交换新节点。

Poisson 初始化解 $H_0r=-q+\langle q\rangle$。若 $r$ 在指定基中，返回精确坐标；否则是多项式系数的最小二乘投影。因此结果字段刻意称为 `projected_poisson`，不是把截断后的解误称为原 Poisson 解。默认完整单项式基、或偶势的完整偶次基，均完整包含相应 Poisson 解。搜索也保留零对数候选，并保留一个均匀概率的有界二次上界，防止数值原子近退化时最终上界大幅恶化。

对偶势 $q$，写 $s=s_e+s_o$ 为偶奇分解。精确恒等式是

\[
E_{s_e}(t)-\frac{E_s(t)+E_s(-t)}2
=(1-t^2)(s_o'(t))^2\ge0.
\]

于是 $\inf E_{s_e}\ge\inf E_s$。因此在**全部次数不超过给定偶数的多项式对数**中，限制到偶多项式不降低最优值。对任意不对反射封闭的自定义基，不能据此宣称同样的无损降维；程序不会擅自改变自定义基。

能量平移 $q\mapsto q+c_0$ 使每个局部能量、每个原子二次上界和每个 Rayleigh 商恰好增加 $c_0$，参数及理论误差不变。默认后端的 `shift_certificate` 复用同一区间划分和见证，精确保持证书误差；自定义后端重新生成其范围证明，可能因新划分而得到不同的保守范围误差。

范围扩展接口 `range_backend` 是显式传入的可信对象/模块，需满足：

```python
name: str
minimum(polynomial, tolerance, max_leaves) -> proof
verify_minimum(polynomial, proof) -> (Fraction lower, Fraction upper)
```

这里 `lower <= inf$polynomial$ <= upper`。验证函数必须独立重建并检查证明，不能只读取两个声称的数字；无效证明抛出 `ValueError`。`verify_lower`、`optimization_certificate`、`verify`、`search` 均接受同一后端参数。非默认后端证书在未显式提供对应可信验证器时会被拒绝。这个接口可接入后续代数极值后端；本模块的默认结果不依赖该扩展。

## 5. 真正的 Rayleigh 谱上界

任取非零有理多项式 $p(t)$，无需它处处为正，完整球面变分原理给出

\[
\lambda_1(H)\le R[p]
=\frac{\int_{-1}^1[(1-t^2)p'(t)^2+q(t)p(t)^2]\,dt}
       {\int_{-1}^1p(t)^2\,dt}.
\]

圆周积分因子相消。程序直接使用有理单项式积分：奇次积分为零，偶次归一化均值为 $1/(j+1)$。保存动能、势能、范数及最终商，并重新计算验证。用于提议 $p$ 的小型 Jacobi 求本征向量不是证明的一部分。

只有这个独立 Rayleigh 证书被绑定到同一 $q$ 并通过验证后，联合结果才含 `spectral_lower`、`spectral_upper`、`spectral_gap`。它们与 `ansatz_lower`、`ansatz_upper`、`ansatz_gap` 并存，意义不同。

## 6. 已验证结果和保留的失败

以下小数是有理证书的便读近似，严格判断均使用证书中的分数。默认 `tolerance=1/10000`、`range_tolerance=1/1000000`。

| 固定势与对数基 | 固定 Poisson 谱下界 | 优化后的谱下界 $L$ | 函数族上界 $U_\mu$ | 函数族误差 $U_\mu-L$ | 独立 Rayleigh 谱上界 |
|---|---:|---:|---:|---:|---:|
| $6t^2$，基 $t^2$ | 1 | 1.438089574827 | 1.438182128513 | $9.25537\cdot10^{-5}$ | 1.571155842617 |
| $6t^2$，基 $t^2,t^4,t^6$ | 1 | 1.567438377485 | 1.567521294341 | $8.29169\cdot10^{-5}$ | 1.571155842617 |
| $20t^2$，基 $t^2,t^4,t^6$ | $-40/9$ | 3.639860741482 | 3.639943147345 | $8.24059\cdot10^{-5}$ | 3.656730749264 |
| $t+3t^2+2t^4$，基 $t,t^2,t^3,t^4$ | 0.061553686326 | 1.037460592636 | 1.037531080765 | $7.04881\cdot10^{-5}$ | 1.061203714989 |

第一与第二行给出一个直接证书反例：第二行已经证明真实基态能量至少为 1.567438…，严格高于第一行的函数族上界 1.438182…。因此把第一行 `ansatz_upper` 当作谱上界会立即得到错误结论。

第一行接受的精确参数为

\[
s=-\frac{3257485643}{4294967296}t^2.
\]

第二行接受的精确参数为

\[
s=-\frac{3370734683}{4294967296}t^2
  -\frac{397723263}{4294967296}t^4
  -\frac{152227129}{4294967296}t^6.
\]

一个更简单、无需任何数值优化的严格改进是 $q=6t^2$、$s=-4t^2/5$。其局部能量是

\[
E=\frac85-\frac{34}{25}t^2+\frac{64}{25}t^4,
\qquad \inf E=\frac{2271}{1600}=1.419375>1.
\]

另有一个**精确全局函数族最优**校准例：取 $q=\frac{20}{9}t^2$、基 $b=t^2$、参数 $\theta=-1/3$，则

\[
E_\theta=\frac23-\frac29t^2+\frac49t^4
=\frac{23}{36}+\frac49(t^2-\tfrac14)^2.
\]

故下界恰为 $23/36$。用单个概率原子 $t=1/2,\mu=1$，有

\[
c=5/9,\quad\ell=-1/2,\quad A=3/4,\quad z=-1/3,
\quad U_\mu=5/9+1/12=23/36.
\]

上下相等，严格证明整个实参数轴上最优值 $M=23/36$，并非声称谱值等于 $23/36$。

保留的失败/限制有明确实验对应：

* 对 $q=6t^2$、六次偶指数，仅给每温度 1 次 Newton、1 轮交换，得到下界约 1.559739960708。数值原子的二次上界约 $1.7812068894\cdot10^{11}$，暴露出权重近退化。它仍保存在 `ansatz_upper_attempts`。均匀概率保底上界约 2.409464864847 被选为最终上界；误差约 0.849724904139，状态严格保持 `ansatz_gap_open`。
* 只用 $b=t^2$、端点 $t=1$ 的单位概率原子时，$A=0$、$\ell=4$，二次目标无界。程序抛出明确异常，不输出伪有限上界。另一方面，基 $b=t$ 使用对称端点原子时 $A=\ell=0$，有限常值二次目标被正确接受。
* 测试通过注入一次提议器失败，验证 `failed_proposals` 保留原因、完整 Poisson 下界仍可接受，且不会误标优化完成。这是验收测试，不是一个成功运行的数值结果。
* 六次偶指数虽然在自己的函数族内误差很小，但 $q=6t^2$ 的真实谱区间宽度仍约 0.003717465132。进一步缩小该谱区间需要改进函数族、Rayleigh 试探或其他谱证书，不能靠把两个不同的误差字段混同。

## 复现与验收

在本目录运行：

```bash
python3 -m unittest -v test_positive_search
python3 positive_search.py
```

第二条命令输出四个小规模案例的完整 JSON 证书，包括所有接受的下界候选、概率上界提议及失败记录。当前 23 项测试全部通过，覆盖：精确标定、明显改进、非偶势、六次局部能量的完整十二次项、概率约束、奇异/无界、能量平移、奇偶恒等式、独立 Rayleigh、预算耗尽、提议失败保留、字段篡改、势/基绑定，以及外部范围后端的显式信任边界。

一个不运行提议器的最小验收示例：

```python
from fractions import Fraction as F
import positive_search as ps

q, basis = [0, 0, F(20, 9)], [[0, 0, 1]]
lower = ps.lower_certificate(q, basis, [F(-1, 3)])
upper = ps.ansatz_upper_certificate(q, basis, [{'t': '1/2', 'weight': '1'}])
certificate = ps.optimization_certificate(lower, upper)
assert ps.verify(certificate)
assert certificate['ansatz_gap'] == '0'
assert certificate['ansatz_lower'] == '23/36'
assert 'spectral_upper' not in certificate
```

这些是精确算术的软件证书及自足数学推导，不是形式化证明助理已检查的定理；字段 `formal_assistant_checked` 始终为 `False`。本模块没有改变旧 Poisson 全振幅结论的适用范围，而是对固定 $q$ 增加一条可验证的更强下界搜索路线。
