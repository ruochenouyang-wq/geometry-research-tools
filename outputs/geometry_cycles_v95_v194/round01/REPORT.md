# 第1轮：从谱存在性反驳到显式有理球谐违反者（V95–V104）

原题固定为单位球面、$q(t)=2t^2$、阈值 $12/5$，要求函数的球面平均值为零。冻结的 V94 baseline 已给出完整零均值空间基态谱区间约 $[2.382655406155,2.382655406157]$，因此能反驳阈值，但其证书明确标记没有包含显式函数。本轮针对这一缺口实现了十个不同的可调用功能，旧目录和 `common.py` 均未修改。

## 已交付的数学结果

一个便于手算复核的违反者是

$$u(x)=x_1(5-x_3^2).$$

它光滑且因 $x_1$ 的奇性而平均值为零。采用标准关联 Legendre 约定，它与

$$\left(P_1^1(t)-\frac1{36}P_3^1(t)\right)\cos\phi$$

相差非零常数 $-24/5$。其 Rayleigh 商和严格余量为

$$R[u]=\frac{722}{303},\qquad \frac{12}{5}-R[u]=\frac{26}{1515}>0.$$

这个两项球谐展开的角因子是 $\pi$，除去角因子后的精确质量为 $505/378$、总能量为 $1805/567$。完整球面的非轴对称 $m=1$ 函数是合法违反者；不能只在轴对称 $m=0$ 扇区寻找它。

V104 自动流程从 baseline 选择 $m=1$，经过量化和残差改进后返回明确的有理球谐展开。此次运行的精确商为

$$R=\frac{6350759013303120704}{2665412294574361111}\approx2.382655406156319,$$

严格余量为

$$\frac{12}{5}-R=\frac{231152468376729812}{13327061472871805555}>0.$$

具体所有 degree 与系数保存在 `certs/V104_counterexample.json` 的 `certificate` 和 `explicit_real_harmonic_expansion` 中，不以浮点本征向量代替函数。

## 十个真实顺序功能增量

`STAGES.json` 精确包含 V95–V104 十项；每项都有函数名、测试名、实际案例、结果文件与原始证书路径。测试和文档没有被单列成版本。

| 版本 | 可调用函数 | 本次案例中的实际作用 |
|---|---|---|
| V95 | `rayleigh_certificate` | 对任意非零有理球谐系数计算精确质量、能量和商；解析两项违反者得到 $722/303$ |
| V96 | `select_sector` | 验证完整谱证书后按扇区上界选择 $m=1$，保留仅选扇区尚未展示函数的限制 |
| V97 | `mass_normalized_proposal` | 将有限广义本征问题转成单位质量坐标，提议四模最低本征向量；仅量化后的证书有证明效力 |
| V98 | `quantize_proposal` | 比较 2、4、8、16、24 bits 的有理候选；4 bits 无收益时保留旧候选，其余较精细候选产生改善 |
| V99 | `sparsify` | 对含 $P_7^1$ 噪声的三项展开按质量贡献截断；一项截断使商升回 2.4 被拒绝，两项截断严格改善 |
| V100 | `residual_certificate` | 计算完整投影残差，包含势乘法创造的所有高阶；单独保存移除的 $l=0$ 项 |
| V101 | `residual_enrich` | 从单个 $P_1^1$ 的残差发现 degree 3，扩充支持并严格降低商；精确本征函数无新支持时诚实返回 |
| V102 | `two_trial_ritz` | 两个试探张成空间的精确 2×2 惯性下界，以及由有理组合给出的实际试探上界，误差小于 $10^{-8}$ |
| V103 | `residual_line_search` | 沿负残差解商导数的二次方程，接受精确/量化驻点及有理备选；$q=0,u=P_1+P_2$ 的有理驻点精确获得商 2 |
| V104 | `find_counterexample` | 自动选扇区、提议、量化、稀疏化和残差改进，给明确的 $R<\text{threshold}$ 证书或 `not_found` |

## 正确性和范围

实函数约定为 $\sum_l c_lP_l^m(t)\cos(m\phi)$。当 $m=0$ 时角向质量因子为 $2\pi$，否则为 $\pi$；在 Rayleigh 商中它相消。$m\ge1$ 的函数自然零均值；$m=0$ 的零均值约束必须排除非零 $l=0$ 系数。所有证书明确绑定 $q,m,\text{mean_zero}$、degrees、coefficients、阈值和 scope。系数可使用很小或很大的有理数，不沿用旧候选的小分母预算；势依旧限次数不超过 6。

对零均值 $m=0$ 空间，算子是投影后的 $\Pi(-\Delta+q)\Pi$。V100 先完成整个势的乘法，再移除 $l=0$；它不是先截断中间乘积。例如 $q=t,u=P_1$ 时，未投影残差为

$$\frac13P_0+\frac23P_2,$$

投影残差仅为 $\frac23P_2$。归一化残差平方正确值为 $4/15$，不是未投影的 $3/5$。`outside_support_residual_squared` 指非零试探支持之外的全部残差贡献，包括高阶和可能的支持空缺；不能把它误认成某个未声明截断窗口外的尾部。

V102 的 `ritz_lower` **只控制两个试探张成空间的最小商**。它不是完整球面或完整投影扇区的谱下界。相应 `actual_trial_upper` 来自一个已显式构造的有理线性组合，因此确实可作为同一完整函数空间的谱上界。证书和验证器拒绝将其 scope 改写为完整球面谱下界。

残差范数本身也不会被包装成基态下界；若后续使用 Temple 型估计，仍需要独立的谱分离证书。

## 保留的无收益与审查修复

本轮保存了粗量化无收益、过强稀疏截断被拒绝、零残差无支持增长、维数预算满时不扩充，以及 $q=0$、阈值 1 的 `not_found`。`not_found` 只表示当前搜索没有找到违反者，不能充当所求不等式的证明。

独立审查发现一个实际数值弱点：把合法试探乘以 $10^{60}$ 后，线搜索曾在转浮点判别式时溢出。已通过对导数二次多项式作精确有理归一化修复；若根提议仍有数值失败，只保留失败原因并继续精确有理备选。独立测试确认缩放 $1,10^{60},10^{-100}$ 得到相同商。另已补上支持扩充的总维数 64 上限，达到上限时诚实保留原候选。

V104 的 `verify` 保证最终显式证书与输入和选中扇区一致；搜索过程的所有日志不因此自动获得证明效力。每阶段单独保存并验证原始数学证书，登记在 `STAGES.json`，避免只验证最终结果就声称中间十项均已验收。

## 运行记录与复用

```text
python3 witness.py
python3 -m unittest -q test_witness round01.test_witness_review
```

真实运行已经产生十个阶段结果，所有登记的例子证书均通过。实现的 12 项测试与独立数学审查的 15 项测试合计 **27 项通过**；测试包含直接 Rodrigues 多项式积分交叉核验，不只镜像同一组装算法。独立审查材料见同目录 `INDEPENDENT_REVIEW.md`。

供后续轮次复用的核心 API：

```python
trial = rayleigh_certificate(q, coefficients, m=0, mean_zero=True,
                             degrees=None, threshold=None)
proposal = mass_normalized_proposal(q, m=0, mean_zero=True, modes=6)
trial = quantize_proposal(proposal, bits=(2, 4, 8, 16, 24))['certificate']
residual = residual_certificate(trial)
assert verify(trial, expected_q=q, expected_m=0, expected_mean_zero=True)
```

主要字段是 `rayleigh_quotient`、`norm_squared_divided_by_azimuth_factor`、`energy_divided_by_azimuth_factor`、`degrees`、`coefficients`；残差方差为 `normalized_residual_squared`，完整残差展开位于 `residual_coefficients`。

这些版本是数学工具的计算与证书能力增量，不是十项新数学定理，也没有声明形式化证明助理已经检查；`formal_assistant_checked` 始终为 `False`。
