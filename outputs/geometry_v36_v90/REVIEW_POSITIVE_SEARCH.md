# positive_search 独立数学与信任边界审查

2026-09-11，Python 3.9.6。审查 `positive_search.py` 与 `test_positive_search.py`，新增独立验收文件 `test_positive_audit.py`。

**结果：未发现被验证接受的谱/函数族界混同、概率归一化缺失、半正定伪认证或异题证书拼接。发现并修复一项自定义区间后端平移的分派问题。新增 19 项独立测试与原有 23 项合计 42 项全部通过，实测 0.699 秒。**

唯一实现改动事先得到主任务明确授权：`shift_certificate` 的默认证明格式快速路径仅用于 `DEFAULT_RANGE_BACKEND` 对象，其他显式提供的可信后端通过自己的生成器与验证器处理。未改既有测试，未运行长计时基准或外部模型。

## 三类界的含义

设 `E_theta=H exp(s_theta)/exp(s_theta)`，`F(theta)=inf_t E_theta(t)`，以及指定基下的最优 Barta 值 `M=sup_theta F(theta)`。

| 字段 | 实际可证明的含义 | 不能据此声称 |
|---|---|---|
| `spectral_lower=L` | `L <= F(theta) <= M <= lambda_1(H)` | `L` 是准确谱值 |
| `candidate_infimum_upper` | 当前候选 `F(theta)` 的上界 | 谱上界或整个函数族的上界 |
| `ansatz_upper=U` | `M <= U` | `lambda_1(H) <= U` |
| `spectral_upper=R` | 独立有理多项式试探的 Rayleigh 商，`lambda_1(H) <= R` | `ansatz_gap_closed` 自动意味着谱误差关闭 |

实现只在独立 Rayleigh 证书验证通过、且其 `q_coefficients` 与下界的原势一致后，在联合证书中加入 `spectral_upper` 和 `spectral_gap`。没有 Rayleigh 证书的联合结果保留函数族上下界，不填入替代的谱上界。

谱下界适用于完整球面：`exp(s)` 是全局光滑正函数，令任意球面试探 `u=exp(s)*v`，分部积分恒等式给出

\[
\langle u,Hu\rangle
=\int e^{2s}|\nabla v|^2+\int E_s|u|^2
\ge L\|u\|^2.
\]

该论证不把任意 `u` 限制为轴对称函数；坐标端点是球面两极，也没有被误加为带边界条件的区间边界。

## 无数值优化的独立反例

采用固定势 `q=(20/9)t^2` 与标量对数基 `b=t^2`，取 `theta=-1/3`。独立展开得到

\[
E_\theta(t)=\frac{23}{36}+\frac49(t^2-\tfrac14)^2.
\]

因此候选谱下界恰为 `23/36`；单位概率原子 `t=1/2` 的函数族上界同样为 `23/36`，函数族误差严格等于零。

然而保持同一原势，将对数基扩大到 `t^2,t^4`，选取纯有理候选

\[
s=-\frac13t^2-\frac1{100}t^4,
\]

在最多 16 个区间叶片下，默认有理验证器给出已重放通过的谱下界

```text
7713232349/11796480000 > 23/36.
```

故真实谱值严格高于前述标量 `ansatz_upper`，直接否定把它当作谱上界的解释。此对照不运行数值搜索。

同一个精确标量例子用常数 Rayleigh 试探 `p=1` 得到 `spectral_upper=20/27`，因此 `ansatz_gap=0` 而 `spectral_gap=11/108>0`。独立测试同时验证两种 gap 并存而不被误标。

## 概率原子、PSD 与全参数上界

实现对每个原子检查 `t in [-1,1]`、非负有理权重，并要求权重和**精确等于一**。接近一但差 `±1e-12` 的原子质量也被拒绝。

由原子重建的

\[
A=\sum_k\mu_k(1-t_k^2)b'(t_k)b'(t_k)^T
\]

是半正定矩阵；代码额外用精确惯性检查并重建全部字段。若 `2Az=ell` 一致，完成平方得到

\[
Q_\mu(\theta)=c+\tfrac12\ell^Tz-(\theta-z)^TA(\theta-z),
\]

故 `U=c+ell^Tz/2` 确实控制整个实参数空间。奇异矩阵不使用浮点伪逆或虚拟正主元；不一致的核方向使二次式无界，明确拒绝有限上界。

独立测试补充了秩一而非满秩的合法例子：基 `t,t^3`、单位原子 `t=0`，得到 `A=diag(1,0)`、`ell=0`、有限上界 `0`。另检查原子拆分、重排和零质量原子不改变矩与目标值。

非正交有理基变换也经过对照：将基 `[t,t^2]` 变为 `[t+2t^2,3t-t^2]`，候选系数按逆转置变换，重建的局部能量及谱下界完全相等；同一概率测度的函数族上界亦完全相等。因此没有观察到基坐标变化时的目标错绑。

## Rayleigh 与证书篡改检查

对 `q=1+t+2t^2`、试探 `p=1+t`，独立单项式积分给出

```text
norm_squared         = 4/3
kinetic_numerator    = 2/3
potential_numerator  = 46/15
spectral_upper       = 14/5
```

程序结果逐项精确一致。试探 `p=t`、`p=-7t` 也给出相同的合法谱上界，确认 Rayleigh 路径不误要求试探处处为正，且满足缩放不变性。

以下独立篡改均被验证器拒绝：

- 修改概率 Gram 为负值、修改最大化点、线性系数或概率质量。
- 在函数族上界证书中额外加入伪 `spectral_upper`。
- 将 `ansatz_upper` 复制为联合结果的谱上界，改写 scope、限制文字或结论状态。
- 将另一原势的 Rayleigh 证书拼接到当前结果，或拼接不同势/不同基的函数族证书。
- 删除区间覆盖的一段、修改当前候选下确界上界，或移植另一局部能量的有效区间证明。

验证以重建后的完整证书字典为准，未只信任最终两个数字。`formal_assistant_checked=False` 保持为否，未声称形式化证明助手验证。

## range_backend 的显式信任边界与修复

默认后端从 `q,basis,theta` 重建整个局部能量多项式，再核验 `-E` 的有理区间证明。非默认后端必须由调用者显式提供可信的 `verify_minimum`；该对象的数学正确性是扩展接口的信任前提，模块不能凭名字替它证明正确。

独立测试确认：

- 外部后端的有效证书只有在提供对应验证器时通过。
- 将其名字改成默认名字，不会使默认验证器接受另一种证明格式。
- 返回浮点数、整数而非要求的 `Fraction`，或返回颠倒的界，会被拒绝。
- 包装后的范围证明仍需绑定当前局部能量；包装层不会给予异题证明权威。

**初版实际错误：** 一个合法 `WrappedSubclass(BernsteinBackend)` 重写 `minimum` 与 `verify_minimum`，采用 `{'wrapped': proof}` 的序列化形式与独立稳定名字。其下界、联合证书与显式验证全部通过；但 `shift_certificate(..., shift=1/7, backend)` 用 `isinstance(backend,BernsteinBackend)` 推断它仍采用默认 schema，随后读取 `proof['intervals']` 而抛出 `KeyError`。

**修复：** 默认格式的区间重用只对 `backend is DEFAULT_RANGE_BACKEND` 启用；任何其他后端，包括重写证明格式的子类，通过自身 `minimum` 重新生成并经自身显式验证器核验。没有扩大任何后端的隐式信任。

修复后的回归同时检查包装子类平移仍需显式验证器、四个相关上下界均精确平移，以及默认路径继续保留原完整分割和证明 gap。任意外部后端重新划分时，一般不承诺其保守 range gap 不变。

## 复现及文件指纹

在本目录运行：

```sh
python3 -m unittest -v test_positive_audit test_positive_search
```

输出为 `Ran 42 tests in 0.699s — OK`。新增 19 项审计测试只构造精确证书，不调用搜索提议器；原有 23 项包含小型真实搜索和失败注入。没有运行长基准、网络请求或模型 API。

```text
positive_search.py
08667022ba4e55f79638d31d67eb7418c17e50affe5fc309e26e338b27a4b404

test_positive_audit.py
0c27d4f0b22c69322daa5f7057409fe60aeb6b7b50ba6206ae2ecefe258cfbb3

test_positive_search.py
ccbf28840ca059a56d5bd4ec93593f73aaf3c1132d938ed5e05cb2d35c286c94
```

结论限于上述实现、信任前提及覆盖案例，未将软件证书重放提升为形式化内核验证，也未将某个有限对数基的优化精度当成真实谱误差。
