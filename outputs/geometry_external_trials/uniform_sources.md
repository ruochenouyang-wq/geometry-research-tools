# 真实开放问题的有界文献筛选：球面谱与插值

## 检索结论

本次检索及主任务补充来源后，选出一个确有作者原文的球面 Gagliardo–Nirenberg 最优常数候选；有限后续检索未发现解决，但这不是对当前开放状态的穷尽确认。冻结 V90 可研究其固定 q 的多项式试探子空间，不能直接认证原问题涉及的全部函数和 q→∞ 极限。另保留一个有后续进展的历史方向，并排除一个已经解决的历史猜想。

本轮先做六轮批量网页查询/读取（最后一轮专门排除已解决问题），再按主任务提供的确定来源直接核对 2204.12414，没有继续广泛检索。作者网页的部分 PDF 直接读取失败，搜索索引仍返回相关节的正文。没有穷尽后续论文或引用链。没有改动冻结的 V90 文件。

| 条目 | 原文状态与后续检查 | 本轮处理 |
|---|---|---|
| Ilyin–Zelik：S² 零均值 GN 最佳 c_q 的渐近常数 | 2022 原文明确提出猜测式；有限 2025/2026 后续检索未找到解决 | 选作有来源候选，运行固定 q=4 的有限子空间试验 |
| Dolbeault 等：p>2 的谱正交约束改进 | 2014 明确开放问题，2015/2017 及更晚论文已有推进 | 历史线索；没有确认剩余范围的当前状态 |
| Chang–Yang：约束 Onofri 的 1/2 改进 | 2017 作者报告明确说明确认最佳常数 | 已解决，排除出当前开放题 |

## 0. 选定候选：S² 零均值 GN 最佳常数的渐近

Ilyin、Zelik，*On a class of interpolation inequalities on the 2D sphere*，[arXiv:2204.12414](https://arxiv.org/abs/2204.12414)，首次提交 2022-04-26；后发表于 2023 年 *Sbornik: Mathematics*。已直接核对 [原 PDF](https://arxiv.org/pdf/2204.12414) 的 §1、式 (1.8)、Corollary 2.1 的式 (2.4)、第 5 页的猜测式及第 9 页 Remark 2.2。

全部范数使用**面积测度 dσ**，而非概率测度。定义最佳常数

\[
c_q^*=\sup_{\substack{0\ne\varphi\in H^1(S^2)\\\int\varphi\,d\sigma=0}}
\frac{\|\varphi\|_{L^q(d\sigma)}}{
\|\varphi\|_{L^2(d\sigma)}^{2/q}
\|\nabla\varphi\|_{L^2(d\sigma)}^{1-2/q}},\qquad q\ge2.
\]

论文证明的显式上界是

\[
c_q^*\le(4\pi)^{-(q-2)/(2q)}\sqrt{q/2}.
\]

作者在第 5 页 §1 用 “one can suggest” 提出

\[
\boxed{c_q^*/\sqrt q\ \longrightarrow\ 1/\sqrt{8\pi}\quad(q\to\infty).}
\]

这是作者提出的猜测式，没有编号为某个 Conjecture。Remark 2.2 只论证 √q 的**增长幂次**最优，不能把它当成前导常数已经得到证明。也不能把 c_q^* 换成只含梯度范数的另一种 Sobolev 嵌入常数，再直接套用涉及 e 的渐近式；两个优化问题必须先严格比较。

主任务已进行标题、sharp 及 2025/2026 的有限后续查询，没有发现解决报告。本笔记对此只作“有限检索未找到”的陈述，不保证不存在未检出的解决或修正。

### 冻结 V90 能真正研究的 q=4 子问题

固定 q=4 时，定义四次比值

\[
K_{4,\sigma}(\varphi)=
\frac{\int\varphi^4\,d\sigma}{
\left(\int\varphi^2\,d\sigma\right)
\left(\int|\nabla\varphi|^2\,d\sigma\right)};
\qquad \sup K_{4,\sigma}=(c_4^*)^4.
\]

可选择整个子空间 `φ=aP₁+bP₃`，对所有 `(a,b)≠(0,0)` 最大化该比值，而非只采样系数。以下是本轮独立积分推导，使用 `dμ=dσ/(4π)`：

\[
M=a^2/3+b^2/7,\quad E=2a^2/3+12b^2/7,
\]

\[
D=ME=2a^4/9+2a^2b^2/3+12b^4/49,
\]

\[
N=\int\varphi^4d\mu=
a^4/5+(8/35)a^3b+(46/105)a^2b^2
+(48/385)ab^3+(241/5005)b^4.
\]

因此 `K₄,μ=N/D`，且 `K₄,σ=K₄,μ/(4π)`。由于 D 在非零参数处严格正，候选界可精确化为齐次四次多项式 `KD−N≥0`。使用比值 `b/a` 时还须覆盖射影边界 `a=0`。

对子空间的完整最优值给出**原全函数类最优常数的下界**；对子空间成立的上界不能直接成为全函数类的上界。q=4 的任何结果也不决定 q→∞ 的猜测。这是相关有限模型的真实算法试用，而非宣称原开放问题已经解决。

## 1. 确有作者原文的历史候选：p>2 下的谱正交约束改进

**出处和时间。** Jean Dolbeault、Maria J. Esteban、Michal Kowalczyk、Michael Loss，*Improved interpolation inequalities on the sphere*。预印本首次提交 2013-09-30；v2 为 2014-01-29；发表于 *DCDS-S* 7(4), 2014。原文 §6.3 的标题就是 “An open problem”。作者版本的排版页码为 21–22（不同版本页码稍有变化，应优先用节号定位）。

- [arXiv 原始论文及版本记录](https://arxiv.org/abs/1309.7931)
- [arXiv PDF](https://arxiv.org/pdf/1309.7931)
- [作者网页 PDF，§6.3](https://www.ceremade.dauphine.fr/~dolbeaul/Preprints/Fichiers/Super28.pdf)

**准确的问题类型。** 这不是一个已经指定唯一常数和唯一约束的猜想，而是作者提出的研究问题：如何把谱估计与恰当选择的正交约束结合，在 d>2 时的 2<p<2*=2d/(d−2)，以及 d≤2 时的 p>2，得到改进插值不等式。单位球面概率测度下的基准式为

\[
\int_{S^d}|\nabla u|^2\,d\mu
\ge\frac{d}{p-2}
\left(\|u\|_{L^p(d\mu)}^2-\|u\|_{L^2(d\mu)}^2\right).
\]

因此在 S² 上，其核心是寻找合适的约束类 A，使右侧的系数 2 可以被一个严格更大的、具有明确谱意义的常数代替，或得到明确的强化项。**把任意自行选择的 A 或某个具体更大常数写成作者原猜想是不准确的。**

**后续状态核查。** 作者后续论文 *Interpolation inequalities on the sphere: linear vs. nonlinear flows*（预印本 2015-09-30，期刊 2017）已经推进了带约束的改进，故不能再把 2014 年的宽泛提问原样标成没有进展的当前开放题。

- [后续作者论文 arXiv:1509.09099](https://arxiv.org/abs/1509.09099)
- [期刊原文，2017，26(2), 351 起](https://www.numdam.org/item/AFST_2017_6_26_2_351_0.pdf)
- [作者保存的后续论文版本](https://www.ceremade.dauphine.fr/~dolbeaul/Preprints/Fichiers/BDEL-14.pdf)

后续文献线索把限制指数写为 2#=(2d²+1)/(d−1)²；S² 对应 2#=9。结尾提出在 p>2# 时显式估计带积分约束最优常数的问题。但本轮没有完成该结尾与 2023 年以后结果的逐条覆盖比对，**此细化只列为下一步应核对的候选，不作为已经确认当前开放的条目**。

另检索到 Brigati、Dolbeault、Simonov 的后续 *Logarithmic Sobolev and interpolation inequalities on the sphere*，其中有专门的约束改进部分；尚未核对它的精确 p 范围能否排除上述候选。

- [作者后续稿 BDS2023-Sphere](https://www.ceremade.dauphine.fr/~dolbeaul/Preprints/Fichiers/BDS2023-Sphere.pdf)

**V90 的真实可用范围。** 冻结算法能对指定多项式势 q(t) 生成完整算子的有理谱界，也能处理其受支持的有限参数族。但是上述问题的原未知量是整个函数类 u，且约束可能是非线性矩条件。它没有现成的“该正交约束插值问题 → 原 V90 结构化谱目标”的已验证归约，故不能直接解原题。

可以做的相关子试验是：固定一个明确多项式势，计算其普通球面 Schrödinger 基态；或对指定参数盒上的多项式势家族检验一条单独写清的谱不等式。这些只能标记为与插值问题相关的有限模型。不能把有限阶 Legendre 采样、一个势族的结果或未验证的拉格朗日乘子推导提升成整个约束函数类的结论。

## 2. 必须排除的历史猜想：Chang–Yang 约束 Onofri 常数 1/2

同一篇 2014 论文 §6.2 写到 Chang–Yang 提出的问题：在 S² 上若

\[
\int_{S^2} e^u x_j\,d\mu=0\quad(j=1,2,3),
\]

是否有

\[
\log\int_{S^2}e^u\,d\mu-\int_{S^2}u\,d\mu
\le\frac18\int_{S^2}|\nabla u|^2\,d\mu?
\]

这里 dμ 是概率测度；在面积测度 dω 下右边等价于 (1/(32π))∫|∇u|²dω。没有约束时 Onofri 的对应概率测度系数是 1/4。

**它不应作为当前开放题。** Gui–Moradifam 的 *The Sphere Covering Inequality and Its Applications* 已给出相关解决；作者 Moradifam 在 UCI 于 2017-05-30 的报告页明确声明确认了 Chang–Yang 于 1987 年猜测的最佳常数。

- [原研究预印本 arXiv:1605.06481，2016-05-20](https://arxiv.org/abs/1605.06481)
- [作者报告原始网页，明确说明解决，2017-05-30](https://www.math.uci.edu/node/35707)

此外本轮检索找到 Gui、Li、Wei、Ye 于 2026-08-13 的预印本 *A positive answer to the generalized Chang-Yang conjecture on S^N*，其摘要宣称证明 N≥3 的推广。这是新的研究预印本声明，本笔记没有独立核验其证明，也没有把它当成已经完成同行审稿的结论。

- [arXiv:2608.13497](https://arxiv.org/abs/2608.13497)

V90 目前也不能直接编码原题的指数重心约束和指数积分，因此即使作为已知题复现，也需要额外的、经过证明的转换与非多项式积分范围控制。

## 3. 查过但未冒充开放题的相关谱文献

**球面基态最佳范数估计。** Dolbeault、Esteban、Laptev，*Spectral estimates on the sphere*，arXiv:1301.1210，2013；期刊 2014。作者给出第一特征值与势范数的最佳关系。本轮在打开的正文中检索 open/unknown/conject 没有命中；不能仅因某个最佳常数函数不是闭式，就替作者创造一个开放猜想。

- [原文](https://arxiv.org/pdf/1301.1210)

**球面 Lieb–Thirring 常数。** Ilyin、Laptev、Zelik，*Lieb–Thirring constant on the sphere and on the torus*，2020-09-01。Theorem 2.1（第 3 页）对均值为零的完整球面正交族给出上界 3π/32；Remark 1.2 同页给出半经典下界 1/(2π)。这是一个未闭合的常数区间，但本轮没有发现作者把某个具体等号常数明确标成猜想。

- [2020 原论文](https://arxiv.org/pdf/2009.00527)

Kowacs、Ruzhansky 的 *Lieb-Thirring inequalities on the spheres and SO(3)*（2024-06-21，PDF 版本日期 2024-07-16）给出较高维的新界，其 Table 1 仍引用 S² 的 3π/32；没有据此确认截至本轮日期的最优常数已被穷尽核查。

- [2024 原论文](https://arxiv.org/pdf/2406.15134)

**重要模型差别。** 上述零均值族对应带投影的算子，不能把它简单当成普通 −Δ+q 的“第二个特征值”。一般势会把常数模态与零均值空间耦合；此外 Lieb–Thirring 涉及完整球面全部角动量和负谱和，而冻结 V90 的普通基态/轴对称多项式工具没有直接认证这些投影算子与多重性的总和。不能为了跑出一个数字就悄悄换问题。

## 建议的诚实汇报方式

可报告：选择 Ilyin–Zelik 的有来源渐近猜测作为研究候选，有限后续检索未找到解决，并在固定 q=4 的完整多项式子空间执行真实试验；同时记录一个有后续进展的历史方向，排除一个早已解决的猜想。应把“已知题复现”“有限子模型探索”“当前状态未穷尽确认的开放候选”分开记录，不把子空间最优或单个 q 的结果宣传为原猜测已解。
