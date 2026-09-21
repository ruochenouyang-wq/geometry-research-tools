# 核心算法的已有工作对照

审查日期：2026-09-18。只读审查当前实现与已有公开主源；本文件不改动历史版本，不执行发布。本次覆盖五类核心方法，不是逐行代码溯源、专利查新或全部数学文献的穷尽检索。

**结论：核心数学方法有明确且接近的已有工作。当前最有依据的定位是“面向单位球面谱问题的可复核计算工具及其工程组合”，尚无充分证据支持“发明了新的通用数学研究算法”或“新的谱估计理论”。** 这个判断来自下面的方法对照，不来自代码风格、版本数量或者精度改善倍数。采用已知方法并不等于抄袭代码；独立编写实现也不等于发明数学方法。

## 1. 有限谱块、Schur 尾部、原函数精确矩

本项目位置：[早期谱证书器](/PROJECT/outputs/geometry_v36_v90/spectral_certifier.py:69)、[直接矩装配](/PROJECT/outputs/geometry_precision_followup/direct_moments.py:104)、[尾部比较矩阵](/PROJECT/outputs/geometry_precision_followup/direct_moments.py:142)。

已有工作的直接对应是 Dusson、Sigal、Stamm 的 *Analysis of the Feshbach–Schur method for the Fourier spectral discretizations of Schrödinger operators*（2020 预印本，2021 修订）。其第 2.2 节把 Schrödinger 算子分成有限投影及其补空间，证明尾部下界，再用依赖谱参数的 Schur 补实现有限维约化；该文也明确研究较低正则性的势。它使用周期 Fourier 基，本项目使用球谐基及有理不等式包围，不能称二者实现完全相同，也不能把“有限维计算控制无限维谱尾”作为本项目首次提出的思想。[原论文，式 (5)–(7)](https://arxiv.org/pdf/2008.10871)

同三位作者的 *The Feshbach-Schur map and perturbation theory*（2021）给出更一般的自伴算子框架和带显式常数的谱估计，进一步说明 Schur 约化与严格误差控制已是既有路线。[原论文](https://arxiv.org/abs/2105.02058)

代码的 `matrix_assembly` 对原势 q 和 q² 求矩，并扣除有限投影和显式近尾部投影，构造剩余耦合的 Gram 矩阵。在正交归一记号下，其核心恒等式是

\[
(QqP)^*(QqP)=Pq^2P-(PqP)^2,\qquad Q=I-P.
\]

这是由正交投影直接展开的恒等式；代码里的质量矩阵、零均值删除和近尾部扣除是其具体坐标实现。Legendre/Rodrigues 公式同样是标准特殊函数公式。[NIST DLMF，Legendre 函数](https://dlmf.nist.gov/14)

**分类：已知数学方法的专门实现与工程组合。** 值得保留的工作量包括：q/q² 精确积分、零均值空间绑定、所有方位角及未算方位角覆盖、有理惯性检查、可重放证书。这些足以构成有用软件，但还不是新理论的证据。若要争取算法论文，需要与上述框架逐式比较，证明当前 Gram/近尾部策略具有此前未有的误差阶、适用范围或复杂度收益，而不能只与自身旧版本比较。

## 2. V235–V237：逼近误差、相对二次型和固定参考

本项目位置：[L² 投影](/PROJECT/outputs/geometry_general_v235_v237/function_models.py:228)、[固定参考与误差传递](/PROJECT/outputs/geometry_general_v235_v237/spectral_transfer.py:156)。

一致范数误差 ε 把谱区间扩成 [L−ε,U+ε]，是 Rayleigh 商比较及极小极大原理的直接推论。L² 路线先用 Hölder 不等式，再用概率球面上的

\[
\|u\|_4^2\le \|\nabla u\|_2^2+\|u\|_2^2.
\]

这正是 Dolbeault、Esteban、Kowalczyk、Loss 论文式 (1) 的 d=2、p=4 特例；该论文也给出此类插值不等式的更早来源。因此该常数和不等式本身不是本项目发现的。[原论文，式 (1)、Theorem 1](https://arxiv.org/pdf/1210.1853)

“选取一个固定参考形式，证明扰动相对于它的界小于 1，再比较谱”属于标准相对二次型技术。Teschl 的第 6.5 节、Theorem 6.24（KLMN）给出相对界小于 1 时的闭形式及自伴算子构造；第 4 章的变分原理用于谱比较。**本项目的具体传递式是将这些已知工具用于当前逼近器的推导，未检索到它作为独立新理论的依据。** 这也不等于声称 Teschl 已逐字给出当前证书格式。[作者提供的书稿，第 4 章、第 6.5 节](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf)

**分类：已知分析工具的可认证适配。** 原始表达式绑定、相消前定义域检查、精确 L² 余项、保持固定参考量不随多项式点态下界恶化，都是具体实现贡献。若以后提出新参考选择算法，应分别证明误差包围的有效性和参考选择的性能；当前有限案例及 θ/ε 恒等式不提供全局最优参考的证明。

## 3. Poisson 矩阵包络、半正定对偶与精确修复

本项目位置：[矩阵包络对偶](/PROJECT/outputs/geometry_v36_v90/dual_v33.py:43)、[精度改进](/PROJECT/outputs/geometry_v36_v90/matrix_refine.py:163)、[证明及范围](/PROJECT/outputs/geometry_v36_v90/PROOF_MATRIX_REFINE.md)。

问题是最小化 tr(CG)，同时对所有 t∈[−1,1] 满足 G≽A(t)。有限 PSD 原子 Yₖ、ΣYₖ≼C 给出 Σtr(YₖA(tₖ)) 下界，属于半正定锥的弱对偶。用有限样本求候选、寻找最坏违反点、再加约束，是半无限规划的交换方法。弱对偶、内点法和坐标预条件有成熟基础。[Boyd–Vandenberghe，*Convex Optimization*，第 5、11 章](https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf)

连续区间约束也有通过多项式及 SDP 求解的既有路线，例如 Papp 的 *Semi-infinite programming using high-degree polynomial interpolants and semidefinite programming*。它并非本项目的同一个 Poisson 包络求解器，但足以排除“用多项式/SDP 处理全区间约束”这一宽泛首创说法。[原论文](https://arxiv.org/abs/1512.06796)

“先浮点寻找，再重构精确有理证据”同样有直接先例。Peyrl–Parrilo（2008）专门研究从数值 SOS 近似获得可精确核验的有理分解。当前程序的外积量化、缩放修复和余量分配不等于完整实现该论文，但这种数值—符号认证思路不是新发现。[作者公开原文](https://www.mit.edu/~parrilo/pubs/files/PeyrlParrilo-ComputingSumOfSquaresDecompositionsWithRationalCoefficients.pdf)

**分类：专用半无限 SDP 的工程组合与启发式改进。** 目前有证据支持的是更好的预条件、候选继承、双向缩放修复、成本协变和证书稀疏化；尚无新收敛率、复杂度或最优性定理。这里检查的模块主要是矩阵包络及范围证书，不应包装为已经实现通用 SOS/Lasserre 层级求解器。其“全局最优”只针对固定 Poisson 函数族目标，不等于原谱常数全局最优。

## 4. 全局搜索、Sturm 极值与区间证书

本项目位置：[全局区域核验](/PROJECT/outputs/geometry_v36_v90/global_core.py:68)、[Sturm 与极值工具](/PROJECT/outputs/geometry_v36_v90/exact_extrema.py:90)、[既有来源声明](/PROJECT/outputs/geometry_v36_v90/PROOF_GLOBAL.md)。

最低特征值是仿射 Rayleigh 商的下确界，从而对势参数为凹函数；顶点插值给出下界后做分支定界，是标准变分与全局优化的组合。Boyd 的课程材料明确描述保留上下界、细分下界最小区域、以 gap 停止的框架。[作者课程材料](https://see.stanford.edu/materials/lsocoee364b/17-bb_slides.pdf)

更直接的谱优化对照是 Mengi、Yıldırım、Kılıç 的 *Numerical Optimization of Eigenvalues of Hermitian Matrix Functions*：其建立分段二次下估计并研究全局收敛。当前代码不用它的全部曲率理论，而采用自己的低维凹性区域界；这种差异值得明确记录，不能据此跳到“新全局谱优化方法”的结论。[原论文](https://arxiv.org/abs/1109.2080)

一元多项式的端点和导数根穷举、平方自由化、Sturm 根计数、Bernstein 区间包围、区间 Newton 都属于已有计算代数技术。SymPy 官方接口提供 Sturm 序列、实根计数与隔离，说明这些基础能力已有成熟公开实现；**本审查没有声称项目复制了 SymPy 代码**。[SymPy 官方文档](https://docs.sympy.org/latest/modules/polys/reference.html)

**分类：已知全局方法的精确证书实现。** 驻点 Taylor 恒等式、临界余式与多种包围取交集可以是本项目的推导和组合，但当前没有新颖性对照或复杂度证明。合理主张是“为所支持目标实现完整覆盖、诚实 open 状态及独立于搜索历史的验证”，而不是发现 Sturm 定理或一般全局优化。

## 5. GPT 工具分工、token 回执、内容哈希缓存

本项目位置：[回执、Store、JSON Pointer](/PROJECT/outputs/geometry_v36_v90/protocol.py:60)、[模型评测边界](/PROJECT/outputs/geometry_v36_v90/MODEL_EVALUATION_STATUS.md)、[公开编码器工程关口](/PROJECT/outputs/geometry_v36_v90/results/token_gate_result.json)。

PAL（2022 预印本，ICML 2023）已明确把语言模型的问题分解与程序执行分开，用外部运行时负责求解。其任务和本项目研究型谱证书不同，但“让模型调程序处理可靠计算”不是新的通用思路。[原论文与正式发表页](https://proceedings.mlr.press/v202/gao23f.html)

本地代码将规范化 JSON 按哈希保存、通过引用取回、按字段和页取证据。这些分别对应成熟的内容寻址和结构化局部访问设计；Git 官方文档介绍内容寻址对象存储，RFC 6901 定义 JSON Pointer。这里是技术原理对照，不是认定当前序列化格式与 Git 对象格式相同。[Git 官方文档](https://git-scm.com/book/en/v2/Git-Internals-Git-Objects)、[RFC 6901](https://www.rfc-editor.org/rfc/rfc6901)

生成者提供证据、接受者由检查器核验的架构也有长期先例；Necula 的 proof-carrying code 是其中一种。该论文的应用对象是程序安全，不是当前数学证书，引用只支持架构思想已有先例。[原始会议论文](https://doi.org/10.1145/263699.263712)

**分类：面向数学证书的接口工程。** 范围与目标绑定、保留向外舍入标记、按需读取大证据、缓存重新核验，值得作为产品特点。当前 token 关口测的是指定公开词表下可见文本大小，项目自己的状态文件明确没有真实 GPT-5.6/Astra 任务对照；不能据此宣传模型推理 token 已降低多少，或已达到“半小时胜过八小时”。

## 来源、代码原创与可公开主张

当前代码及报告已经承认多项外部依赖：球面插值公式在 `spectral_transfer.py` 的 `EMBEDDING` 中直接标出论文；全局证明文档引用 Boyd 和 Mengi 等；精度报告引用 Temple 定理。应将这些分散引用汇总到 GitHub 的方法来源表，避免读者将自写实现理解为数学原理首创。

可选编码器已有 [THIRD_PARTY_NOTICES.md](/PROJECT/outputs/geometry_v36_v90/THIRD_PARTY_NOTICES.md)，说明使用 tiktoken 0.11.0，并记录随包词表的来源和哈希。上游代码采用 MIT 许可证，可直接对照其官方许可证；具体发布文件仍应按各自已记录来源保留通知。[上游许可证](https://github.com/openai/tiktoken/blob/main/LICENSE)

本次方法审查未进行全网逐行匹配，没有检出代码相似性的结论可供发布使用。数学方法相似不能证明复制；局部文件没有版权头也不能证明独立创作或无许可要求。代码溯源、依赖目录筛选和项目自身许可证应作为单独发布准备核验。

建议对外表述：

> 本项目实现了一套针对单位球面 Schrödinger 谱及相关不等式的可复核计算工具，组合精确积分、Schur 谱尾估计、相对二次型比较、全局优化和内容寻址证书接口。数学基础沿用所引文献，当前贡献主要是具体实现、适用问题的认证流程及可复现精度案例。尚未声称通用定理证明能力、核心数学理论首创或实测模型加速。

若希望主张更高层次原创性，需补充明确的新命题或新算法、逐项最接近工作对照、比现有通用实现更强的适用范围/复杂度/收敛分析，以及固定未参与开发案例上的公平实验。版本累计数和与自身弱基线的巨大精度比都不能代替这些证据。

本次检索主要使用：Feshbach–Schur spectral discretization / irregular potentials；relative form boundedness / KLMN；rational SOS / semi-infinite SDP；Sturm / eigenvalue global optimization；PAL / proof-carrying code / content-addressable storage。来源取作者论文、作者书稿、出版页面、官方项目文档和 RFC；未把搜索摘要中的二手网站作为技术结论依据。没有打开的 DOI 页面也未用来证明定理内容。
