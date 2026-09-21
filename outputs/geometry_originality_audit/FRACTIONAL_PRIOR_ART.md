# 分数幂适配与 Temple 认证：原创性核对

检查日期：2026-09-18。只读检查 `outputs/geometry_precision_followup/enriched_trial.py` 及 `REPORT.md`；未改算法。本文件是针对一个核心机制的文献对照，不是穷尽检索或专利结论。

## 结论

**不能把“按奇点选择分数幂基，再用 Rayleigh/强残差给谱界”作为新数学机制。** 非多项式奇点适配谱基、Müntz 空间、广义 Jacobi 空间和 Temple 认证都有明确先例。项目可保留的贡献候选，是针对球面赤道势 `q(z)=-|z|^(-1/4)` 的具体结构推导，以及将精确矩、原算子完整残差、单实方位角分量的间隙界、其他分量与无穷尾部覆盖合成可重放证书的实现。**“有自己的实现”与“已证明学术新颖性”是两回事。** 尚未找到与此具体完整组合逐项相同的论文，但未找到不构成首创证明。

## 六项相关主源

1. **Huiyuan Li、Zhimin Zhang，2016 预印本／2017 正式论文，*Efficient Spectral and Spectral Element Methods for Eigenvalue Problems of Schrödinger Equations with an Inverse Square Potential*.** [作者预印本全文](https://arxiv.org/pdf/1606.06388)，[SIAM 正式版](https://epubs.siam.org/doi/10.1137/16M1069596)。

   相似：引言直接提出按特征函数奇异性修改非多项式谱基；径向分数幂乘 Jacobi 多项式，并结合球谐。不同：主要是球体、扇形和多边形上的反平方势与 Dirichlet 问题；本项目是球面上的赤道奇异势、零均值空间与非正交分数幂单项式。其预印本明确区分只加入少数奇异项和拟合整个奇异结构；不能据本项目的普通多项式对照认定超越该类成熟适配谱法。本次读了全文引言、基函数构造及附录，未把其数值误差表当成有理端点证书。

2. **Xiu Yang、Li-Lian Wang、Huiyuan Li、Changtao Sheng，2023，*Müntz ball polynomials and Müntz spectral-Galerkin methods for singular eigenvalue problems*.** [作者预印本](https://arxiv.org/abs/2303.05020)，[全文 §4.2](https://arxiv.org/pdf/2303.05020)。

   这是最接近的先例：§4.2 明确研究球体上带有有理分数幂势的 Schrödinger 特征问题，利用变量替换选择 Müntz 基，球谐分解后精确推导质量与刚度矩阵元素。不同：径向点奇异性与本项目的赤道曲线奇异性不同；其正交基设计、稀疏矩阵与本项目的稠密有理矩阵也不同。本次检查的 §4.2 未见本项目这种完整强残差加 Temple、逐分量排除后的有理证书流程；该观察不意味着整条组合必然新颖。

3. **Dianming Hou、Chuanju Xu，2017，*A fractional spectral method with applications to some singular problems*.** [出版方页面，Advances in Computational Mathematics 43, 911–944](https://link.springer.com/article/10.1007/s10444-016-9511-y)。

   相似：以经典 Müntz 分数幂多项式构成逼近空间，适配有限正则解；文章给出变量变换后光滑时的谱收敛分析。不同：对象主要是积分微分和分数阶微分方程，不是本项目的球面二阶整数阶算子谱认证。此次只取得出版方摘要，因此不对其全部实现细节下结论；“用分数幂提升奇异函数逼近精度”的大思路已经存在则证据明确。

4. **Sheng Chen、Jie Shen、Li-Lian Wang，2014 预印本／2016 正式论文，*Generalized Jacobi Functions and Their Applications to Fractional Differential Equations*.** [作者预印本与全文](https://arxiv.org/abs/1407.8303)。

   相似：让基函数的非整数参数适配奇异行为，并在合适加权空间中分析精度。不同：Riemann–Liouville／Caputo 分数阶微分算子、Petrov–Galerkin 设计和相关误差理论与本项目不同。**“分数幂基”不等于“分数阶算子”**；这篇是适配空间思想的背景先例，不能称为本项目的同一算法。

5. **Gerald Teschl，2009，*Mathematical Methods in Quantum Mechanics: With Applications to Schrödinger Operators*，Theorem 4.13，印刷页 120。** [作者提供的书稿](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf#page=132)。

   这是报告已引用的直接来源。该定理给出以完整强残差及相邻谱点界控制特征值的 Temple 不等式；书中还明确说明地面态的 Rayleigh 上界和第二特征值下界如何搭配。因此代码 `mu - residual/(beta-mu)` 是标准公式的实例化，不是新不等式。项目特定工作是可靠取得 beta、处理 cos/sin 重数、核验定义域并覆盖其他方位角，而不是创造 Temple 公式。

6. **Xuefeng Liu，2026-07-05 预印本，*Guaranteed Lower Eigenvalue Bounds for Spectral Galerkin Methods with Application to Schrödinger Operators*.** [作者预印本元数据与摘要](https://arxiv.org/abs/2607.04247)。

   作为近期边界检查：该文提出投影常数与辅助投影器路线，为谱 Galerkin 计算下界；其摘要明确区分与需要相邻特征值信息的 Kato／Weinstein–Temple 界。项目仍采用 Temple 并另外证明谱间隙，不是该论文的新投影定理。此次只核对作者摘要，未独立审查其全部证明，不能将摘要里的优先性声明当成定论；但发表项目时不能笼统声称“首次把谱法与严格下界结合”。

## 项目中的具体推导与已知机制

以 `x=sqrt(1-z^2) cos(phi)` 表示一个实 m=1 分量，当前试探为 `x P(|z|)`。代码第 45–50 行作用式对应

\[
H(x|z|^s)=x\{(s+1)(s+2)|z|^s-s(s-1)|z|^{s-2}-|z|^{s-1/4}\}.
\]

从常数首项的 `-|z|^(-1/4)` 与二阶导数项平衡，得到首个非光滑幂 `s=7/4`，系数为 `-1/[(7/4)(3/4)]=-16/21`。以 7/4 和 2 加法生成幂次，是当前方程指数平衡的自然闭包；本次没找到相同赤道方程的相同系数发表记录，但该计算属于具体方程的局部展开，不能单凭数值不同宣布新方法。

代码第 37–42 行的积分 `4/[(r+1)(r+3)]` 来自直接积分 `(1-z^2)|z|^r`。第 59–99 行构造质量、算子和算子像的 Gram 矩阵，从而精确得到 Rayleigh 商和强残差。这些是标准 Galerkin／后验残差量；值得展示的是使用有理数对原函数做完整积分，不依赖采样误差或只算投影残差。

第 102–173 行把单一实分量的 Temple 下界和其他分量、无穷角向尾部界合并，是具体完整认证流程。该流程可能具有工程组合价值，但要把它升级成论文中的原创算法贡献，仍需说明相对于现有验证谱法究竟新增了哪项可证明能力、复杂度收益或可复用理论。

## GitHub 表述建议

可用：**“基于奇点适配谱逼近与经典谱夹逼的可重放球面谱证书原型；包含赤道幂奇异势的分数幂试探与精确残差实现。”**

暂不宜用：全新分数幂算法、首个奇异谱精确算法、已超越现有谱法、证明 AI 数学研究获得若干亿倍加速。当前普通偶多项式基与分数幂基的巨大证书宽度比，只能支持这两条具体已测流程之间的精度对照。

下一步最有判别力的原创性工作是把该固定势推广为一个明确的参数族，证明自动选幂与残差收敛的条件，并与成熟的 Müntz／广义 Jacobi 适配空间进行同题、同认证规则、明确预算的对比。它比继续提高同一算例的小数位数更能证明方法贡献。

## 检索边界

检索覆盖 Müntz spectral method、fractional potential eigenvalue、singularity-enriched Schrödinger spectral basis、generalized Jacobi、Temple residual bounds 等词，以及 7/4、16/21、equatorial singular potential 的特定组合。优先采用作者 arXiv、作者大学文档及正式出版方页面。没有进行付费全文库穷尽查重、全部软件仓库逐文件相似度审计，也没有证明任何数学优先权。代码独立创作和许可来源需结合主审计另行核查。
