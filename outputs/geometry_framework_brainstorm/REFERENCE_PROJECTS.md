# 相似项目与框架参考

检索日期：2026-09-18。用途：继续 brainstorm，比较已有架构，不在本轮选型或实现。阅读对象为作者论文、官方仓库、官方文档及部分搜索代码；没有安装运行这些项目，没有复测其成绩。以下“可借鉴”均为本次分析，不代表对方已经实现我们的几何研究目标。

## 1. AlphaGeometry / AlphaGeometry2：模型提出构造，符号引擎推进证明

原始 AlphaGeometry 将题目表示为几何对象和关系图。先运行 DD+AR 推导，未完成时由模型提出辅助点等构造；过滤非法构造后重新推导，按 beam search 保留候选分支。实际入口代码保存图、模型上下文和原题表达，并输出带依赖的证明步骤。[官方仓库](https://github.com/google-deepmind/alphageometry)、[实际搜索代码](https://github.com/google-deepmind/alphageometry/blob/main/alphageometry.py)

可借鉴：为模型定义有数学意义的研究动作，例如提出辅助函数、分解或引理，让确定性工具处理其后的常规推导。其几何对象语言也说明，搜索空间的设计本身是核心工作。

范围：主要是相应语言覆盖的欧氏几何题，不是一般微分几何或无限维谱研究。原版公开实现移除了部分内部并行优化。AlphaGeometry2 当前官方仓库公开的是 DDAR 符号核心及示例；部分难题示例人工提供辅助点，不能称为完整神经搜索系统已公开。[AlphaGeometry2 仓库说明](https://github.com/google-deepmind/alphageometry2)

## 2. LeanDojo / ReProver：把证明状态变成可搜索、可检索的对象

ReProver 从当前 Lean 状态检索可用前提，模型生成 tactic，best-first search 选择分支，Lean 执行后返回状态。搜索代码维护已见状态，避免重复节点，并分别记录模型与环境执行时间，设定节点和时间预算。[作者仓库](https://github.com/lean-dojo/ReProver)、[搜索代码](https://github.com/lean-dojo/ReProver/blob/main/prover/proof_search.py)、[论文](https://arxiv.org/abs/2306.15626)

可借鉴：将研究状态中的已知界、未证条件、可用引理和历史动作明确保存，以便检索、去重、回退。长期资产是可用数学前提及状态关系。LeanDojo-v2 另提供 tracing、题目数据库、训练器和证明器的模块化组织，属于后续可继续查阅的参考。[官方 v2 仓库](https://github.com/lean-dojo/LeanDojo-v2)

范围：主要处理已形式化命题；Lean 对形式输入的检查不会自动保证该输入忠实表达自然语言原题，也不会自动建立我们需要的奇异算子误差理论。

## 3. DeepSeek-Prover-V2：把提出引理与完成引理分开

论文的方法将自然语言分析与 Lean 子目标对应起来，用大型模型拆分问题，再由较小 prover 递归处理子目标；成功子目标成为后续前提，组合完成的证明用于训练。[论文 §2.1](https://arxiv.org/html/2504.21801v1#S2.SS1)

可借鉴：规划和执行使用不同成本；明确引理间依赖；把成功分解积累成后续经验。这里描述的递归分解尤其是训练数据构造流程，不能直接当作公开即用的同款研究 agent。

公开程度：本次看到官方仓库的论文、说明、权重/数据链接及样例证明，没有看到完整递归搜索和 RL 训练流水线源码。[官方仓库](https://github.com/deepseek-ai/DeepSeek-Prover-V2)

## 4. AlphaProof：遇到难题时研究附近的题

AlphaProof 结合策略网络、树搜索和 Lean。其题目专属训练会生成简化、推广、辅助引理等变体，筛选并形成课程，在目标题附近继续学习。[Nature 原文，Methods](https://www.nature.com/articles/s41586-025-09833-y)

可借鉴：把“卡住后研究哪些邻近问题”作为明确能力，保存变体与原题的关系。对我们而言，可以先考虑积累局部结论和有效策略，是否更新模型权重留待后续判断。

范围与公开程度：依赖较大训练投入；论文说明提供伪代码、部分基准及远程服务申请，不是完整系统源码。解决变体也不自动解决原题，必须确认结论可以传回。[同文 Code availability](https://www.nature.com/articles/s41586-025-09833-y)

## 5. FunSearch：搜索产生答案的程序

固定问题骨架和评估器，让模型修改关键函数；合格程序进入多个群体，按成绩和程序特征采样形成下一轮提示，并通过群体更新维持不同探索方向。[Nature 论文](https://www.nature.com/articles/s41586-023-06924-6)

可借鉴：搜索生成函数基、候选构造或误差证书的小程序，使一次发现可以用于多个输入。多个群体可承载机制不同的路线。评估器的正确性和覆盖范围决定了结果能说明什么，进化搜索本身不保证全局最优。

官方仓库提供进化与代码操作的核心和单线程流程，不包含原模型、执行沙箱及完整分布式基础设施。[官方代码](https://github.com/google-deepmind/funsearch)

## 6. AlphaEvolve：同时搜索构造、搜索算法与中间解

白皮书允许不同抽象层级的代码参与进化；用户提供可修改区域和评估函数，模型提交修改。系统使用多指标档案、群体、逐级评估及异步执行组织候选。[白皮书 §2](https://arxiv.org/html/2506.13131v1)

可借鉴：研究对象可以包含算法本身；廉价检查先筛选，昂贵核验留给有希望的候选。数值评估与精确恒等式验证在不同任务中的强度不同，需要分别看待。

作者结果仓库公开数学结果及部分验证代码，明确不包含运行 AlphaEvolve 的系统代码。第三方复现不能视为官方原实现。[官方结果仓库](https://github.com/google-deepmind/alphaevolve_results)

## 7. ShinkaEvolve：将候选档案、成本与经验总结结合

ShinkaEvolve 从分群档案选择父程序和启发程序，使用修改、重写、交叉等方式产生候选，过滤重复后评估。论文描述按相对收益选择模型、周期性总结经验等机制。[论文 §3](https://arxiv.org/html/2509.19349v1)

可借鉴：保存多条有差异的成功路线与产生它们的策略，按实际反馈分配模型调用。公开仓库包含可恢复的运行框架及本地/集群评估支持，适合后续阅读实现。[Sakana AI 官方仓库](https://github.com/SakanaAI/ShinkaEvolve)

范围：任务仍需要外部评估器。“novelty”主要针对运行档案中的代码相似度，不是学术原创性检查；候选获得高分也不等于得到全局最优性证明。

## 8. Chebfun / Spherefun：函数是由系统维护表示的对象

Chebfun 自动调整函数的数值表示，并支持分段及部分奇点处理。官方文档给出已知指数的结构化表示，也允许尝试自动估计某些未知奇异指数。这与我们此前“自动选择函数”的设想直接相关：必须进一步区分函数奇点识别、从方程推导未知解的结构，以及严格认证。[奇点处理，Guide 9](https://www.chebfun.org/docs/guide/guide09.html)

Spherefun 面向球面函数，采用双 Fourier 扩展和保持对称性的低秩构造，使表示适合球极结构及后续微分积分。[Spherefun，Guide 17.3](https://www.chebfun.org/docs/guide/guide17.html)、[公开代码](https://github.com/chebfun/chebfun)

可借鉴：让函数对象保存几何结构，表示由系统管理，后续运算通过统一接口使用。范围：数值解析度及接近机器精度的结果不等于严格区间证书；文档明确说明部分奇点功能可靠性较弱。

## 9. pyMOR：把昂贵计算变成可复用的小模型

pyMOR 用 VectorArray、Operator、Model 抽象接入不同求解器，通过降阶器生成小模型；误差估计驱动贪心选择新的参数点或基向量。适当结构下，离线阶段完成高维工作，在线求解与误差估计只需较小规模的计算。[技术架构](https://docs.pymor.org/latest/technical_overview.html)、[贪心算法](https://docs.pymor.org/latest/autoapi/pymor/algorithms/greedy/index.html)

可借鉴：将跨问题族复用作为框架能力，明确保存模型适用的参数范围和误差条件。它与“研究很久，后续问题很快”的目标直接相关。

范围：不同问题需要不同估计器；降阶误差不自动覆盖原连续问题的离散化误差。在线优势需单列离线成本及其摊销。代码公开。[项目仓库](https://github.com/pymor/pymor)、[架构论文](https://arxiv.org/abs/1506.07094)

## 10. IntervalRootFinding.jl：搜索顺序与证明判定分别负责

该项目保存区间区域构成的搜索树；Newton、Krawczyk 等方法收缩或判定区域，遍历策略决定下一步处理哪里。已证无根区域可以排除，已证唯一根区域保存，无法判定则细分或返回 unknown。[内部架构](https://juliaintervals.github.io/IntervalRootFinding.jl/stable/internals/)

可借鉴：模型可以建议研究顺序，但已证、排除、未定的状态由确定性检查维护。保证标记还反映数值来源和运算条件，不把未知当成失败或成功。[保证机制](https://juliaintervals.github.io/IntervalRootFinding.jl/stable/decorations/)

范围：主要针对所给有限区域的求根；用于无限维谱或几何全局最优前，还必须建立相应覆盖、尾部和紧性论证。[公开代码](https://github.com/JuliaIntervals/IntervalRootFinding.jl)

## 11. LeanEuclid：原题表达是否忠实也需要单独检查

LeanEuclid 在 Lean 中实现欧氏几何的形式系统 E，使用模型翻译文本步骤，SMT 辅助补足图形推理；还尝试检查自动生成的命题与参考命题的逻辑关系。[ICML 论文](https://proceedings.mlr.press/v235/murphy24a.html)

可借鉴：为问题表达设置单独的语义检查，避免证明了修改过的目标或遗漏几何假设。公开仓库自己指出，使用外部 SMT 的做法存在 soundness 代价；不能仅因项目使用 Lean 就宣称整条链都由 Lean 内核独立验证。[官方 README，Acknowledgements](https://github.com/loganrjmurphy/LeanEuclid)

## 这些参考如何改变我们的 brainstorm

以下是本次综合推论，暂不作为确定的架构决策：

- **把模型动作设计成有意义的数学操作。** AlphaGeometry 的辅助构造提示我们，可以研究哪些操作能真正开启新的推导路径。
- **维护研究状态与依赖关系。** ReProver 的状态图、AlphaProof 的邻近题和 LeanEuclid 的语义检查，支持把原题、子题和转换关系显式组织起来。
- **同时保存多个候选及其适用条件。** 程序进化项目中的档案可扩展为带误差证据、证明缺口和迁移范围的研究档案。
- **把表示与参数复用交给数学软件。** Chebfun 与 pyMOR 已经实现相当多这方面能力，值得作为强基线，也值得研究与模型策略之间的连接。
- **明确谁有权宣布成功。** 求解器分数、形式证明闭合、区间认证和语义一致性回答不同问题；原题验收需要组合适用的证据。

目前最值得保留的框架假设，是将“模型提出数学动作、程序推进并检查、档案保存可迁移关系”连接起来。要判断它是否比直接调用现有项目更有价值，仍需同题、同验收标准、计入全部成本的比较。单纯组合这些已有机制并不证明算法原创性。
