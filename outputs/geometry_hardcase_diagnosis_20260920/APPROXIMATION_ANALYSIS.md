# H11–H13 函数逼近失败：支持边界与近零误差

**三题失败的直接原因是现有分段多项式接口只接受指数 `-1/4`；三题本身都可以在概率 L² 意义下逼近。** 它们没有进入精度规划，也没有耗尽十秒预算。这里区分原始运行事实、精确数学推导和用于解释尺度的理想化计算，不把诊断当成已经实现的修复。

本文只读冻结源码及已经公开的最终结果；未读取随机种子，未修改求解器、证书、协议或原始评测。详细参数、九十条已有运行记录、精确分数及来源哈希见 [approximation_diagnostics.json](/PROJECT/outputs/geometry_hardcase_diagnosis_20260920/approximation_diagnostics.json)。

## 1. 原题与实际失败位置

统一记源函数为 `f(z) = A |z|^α + b = A |z|⁻ᵖ + b`，其中 `p = -α > 0`，`a = |A|`。**源码中的指数 α 为负数；下面的 p 为正的奇异阶数。** 三题均为函数逼近，没有平均零谱空间约束；目标 τ 是 L² 范数，而非其平方。

| 题目 | 原函数 | 目标 τ | 六方法的实际结果 |
|---|---|---:|---|
| H11 | `-(1/3)|z|⁻²⁄⁹` | `1e-6` | 每方法 0/5，均无有效证书 |
| H12 | `-(3/5)|z|⁻²⁄¹¹` | `1e-8` | 每方法 0/5，均无有效证书 |
| H13 | `-(3/5)|z|⁻²⁄⁹` | `1e-10` | 每方法 0/5，均无有效证书 |

原题见 [公开实例 H11](/PROJECT/outputs/geometry_strategy_v317_runtime/evaluation/holdout_instances.json:229)、[H12](/PROJECT/outputs/geometry_strategy_v317_runtime/evaluation/holdout_instances.json:250)、[H13](/PROJECT/outputs/geometry_strategy_v317_runtime/evaluation/holdout_instances.json:271)。完整最终报告为 [holdout_20260920T042408501655Z/report.json](/PROJECT/outputs/geometry_strategy_v317_runtime/evaluation/holdout_20260920T042408501655Z/report.json)。

已有三十个 session 中与这三题有关的九十条记录均保留；各方法的所有相关求解分支都报告 `ValueError: This structured model supports only exponent -1/4`，没有这些题的超时记录。例如新版第一次运行的 [H11 拒绝原因](/PROJECT/outputs/geometry_strategy_v317_runtime/evaluation/holdout_20260920T042408501655Z/release_r0_output.json:8974)、[H12](/PROJECT/outputs/geometry_strategy_v317_runtime/evaluation/holdout_20260920T042408501655Z/release_r0_output.json:9027)、[H13](/PROJECT/outputs/geometry_strategy_v317_runtime/evaluation/holdout_20260920T042408501655Z/release_r0_output.json:9080)。这一次评估器计入每题的在线时间约 48–67 微秒，远小于十秒期限；该数值不包含整个 session 的冷启动。外层 `execution_ok: true` 表示评估器完成了记录流程；内层求解分支明确为 `execution_ok: false`，不能解释为求解成功。

实际调用链如下：新版 `evaluate_task` 调用其持有的原服务方法，原服务的函数逼近分支调用 `adaptive_approx.solve`；该函数第一步即调用 `pieces.normalize_function`。后者直接要求指数等于 `-1/4`，所以规划器尚未搜索分区或次数就退出。相关入口依次是 [新版入口](/PROJECT/outputs/geometry_strategy_v317_runtime/release_service.py:78)、[原服务逼近分支](/PROJECT/outputs/geometry_strategy_v238_v317/research_service.py:139)、[规划器入口](/PROJECT/outputs/geometry_strategy_v238_v317/adaptive_approx.py:149)、[严格指数检查](/PROJECT/outputs/geometry_precision_followup/piecewise_models.py:45)。私有运行时仍加载冻结的原服务源码，见 [request_runtime.py:107](/PROJECT/outputs/geometry_strategy_v317_runtime/request_runtime.py:107)。

## 2. 为什么不能只删掉指数检查

`-1/4` 不只是接口白名单。下列数学公式和证据格式都利用了 `p=1/4`：

| 位置 | 当前写死的关系 | 对一般 p 的含义 |
|---|---|---|
| [幂矩与分区](/PROJECT/outputs/geometry_precision_followup/piecewise_models.py:103) | `t = R⁴`；幂矩系数 `4/(4k+3)`，根的幂为 `4k+3` | 它精确计算 `∫t^(k−1/4)dt`，不能直接用于 `k−2/9` 或 `k−2/11` |
| [单区源函数平方范数](/PROJECT/outputs/geometry_precision_followup/piecewise_models.py:114) | 平方项为 `2 A²(R_right²−R_left²)` | 这是 `∫A²t^(-1/2)dt` 的特例 |
| [缩放模板](/PROJECT/outputs/geometry_precision_followup/piecewise_models.py:147) | 系数按 `1/R` 缩放，平方残差按 `R²` 缩放 | 一般 `t=B x` 应分别按 `B⁻ᵖ`、`B^(1−2p)` 缩放 |
| [全域误差合成](/PROJECT/outputs/geometry_precision_followup/piecewise_models.py:216) | 核心系数 `2/9`，层因子 `r^(2L)`，环带级数分母 `1−r²` | 这些由四次根网格及 `p=1/4` 的投影误差共同决定 |
| [全域源范数](/PROJECT/outputs/geometry_precision_followup/piecewise_models.py:220) | `b² + (8/3)bA + 2A²` | 一般为 `b² + 2bA/(1−p) + A²/(1−2p)` |
| [自动规划误差与最低层数](/PROJECT/outputs/geometry_strategy_v238_v317/adaptive_approx.py:68) | 重用同一 `2/9`、`r²` 公式；固定次数误差下限也据此计算 | 即使构造器被扩展，旧规划公式仍会估错是否达标 |

优化版也重复保留了这些专用关系，见 [fast_math.py:39](/PROJECT/outputs/geometry_strategy_v238_v317/fast_math.py:39)、[fast_math.py:105](/PROJECT/outputs/geometry_strategy_v238_v317/fast_math.py:105)。证书还写有 `left_fourth_root`、`right_fourth_root` 和固定缩放说明。独立验证器重新归一化原函数并重建整份旧格式证书，见 [piecewise_models.py:247](/PROJECT/outputs/geometry_precision_followup/piecewise_models.py:247)；规划结果也必须通过该重建检查，见 [adaptive_approx.py:132](/PROJECT/outputs/geometry_strategy_v238_v317/adaptive_approx.py:132)。因此只放宽入口会使源函数与被计算的积分不一致，不能形成有效修复。

另一个潜在共享条件是形状缓存当前只以 `(degree, ratio)` 为键，见 [adaptive_approx.py:25](/PROJECT/outputs/geometry_strategy_v238_v317/adaptive_approx.py:25)。这在只有一个奇异指数时合理；若数学族扩展到多个 p，形状误差也依赖 p，不能继续跨指数复用同一缓存项。这里仅指出假设边界，没有修改缓存。

## 3. 球面概率 L² 的归一化与可积性

原代码使用单位球面上的概率测度 `dμ = dσ/(4π)`，见 [proof_transport.py:11](/PROJECT/outputs/geometry_strategy_v238_v317/proof_transport.py:11)，函数证书声明同一概率 L² 范数，见 [piecewise_models.py:226](/PROJECT/outputs/geometry_precision_followup/piecewise_models.py:226)。

可以直接推导所需测度因子：球坐标下 `z=cos θ`，`dσ=sin θ dθ dφ=|dz|dφ`。积分掉 `φ∈[0,2π)` 后，z 在 `[-1,1]` 的概率密度恰为 `1/2`。因此对只依赖 z 的误差，

\[
\|e\|_{L^2(S^2,d\mu)}^2=\frac12\int_{-1}^1|e(z)|^2\,dz.
\]

若误差关于 z 为偶函数，则这又等于 `∫₀¹|e(t)|²dt`。**双侧近零带的 2 与概率密度 1/2 抵消；不能再多乘或少乘一个 2。** 代码中的一般轴函数矩 `1/(α+k+1)` 和平方范数 `1/(2α+1)` 与这一归一化一致，见 [function_models.py:217](/PROJECT/outputs/geometry_general_v235_v237/function_models.py:217) 及 [function_models.py:235](/PROJECT/outputs/geometry_general_v235_v237/function_models.py:235)。

对非零奇异幅度，积分 `∫₀¹t⁻²ᵖdt` 有限当且仅当 `p<1/2`，也就是 `α>-1/2`。H11/H13 的 `p=2/9`、H12 的 `p=2/11` 都满足条件；它们甚至都弱于当前支持的 `p=1/4` 奇异性。`p=1/2` 是对数发散边界，`p>1/2` 也发散；这个结论针对非零奇异项，不能用零幅度常数函数当例外来解释当前三题。

z=0 对应赤道，虽然该带附近函数无界，赤道本身的球面面积为零；在该处单独赋任意有限值不改变 L² 类。当前源码也明确采取了几乎处处代表元约定，见 [piecewise_models.py:233](/PROJECT/outputs/geometry_precision_followup/piecewise_models.py:233)。无界不妨碍 L² 逼近，但不能据此要求有界多项式在包含赤道的整个区域上一致逼近无界源函数。

## 4. “直接扔掉近零区”需要多小的带宽

以下是假设明确的理想模型：在 `|z|<ε` 内把奇异项置零；常数偏移若存在则原样保留；在外部 `|z|≥ε` 精确保留原函数。三题实际偏移均为零。于是只有近零区产生误差，且对 `0<p<1/2`：

\[
E_0^2(\varepsilon)=\frac12\int_{-\varepsilon}^{\varepsilon}a^2|z|^{-2p}\,dz
=\frac{a^2}{1-2p}\varepsilon^{1-2p}.
\]

要求 L² 范数 `E₀≤τ` 时，要先把目标平方：

\[
\varepsilon\le\left[\frac{(1-2p)\tau^2}{a^2}\right]^{1/(1-2p)}.
\]

| 题目 | `1−2p` | ε 上限的精确表达式 | ε 上限近似 | 若按 z 坐标二分，到这一尺度的层数 |
|---|---:|---|---:|---:|
| H11 | `5/9` | `(5×10⁻¹²)^(9/5)` | `4.55141051×10⁻²¹` | 68 |
| H12 | `7/11` | `[7/(39600000000000000)]^(11/7)` | `1.76164652×10⁻²⁵` | 83 |
| H13 | `5/9` | `[1/(64800000000000000000)]^(9/5)` | `2.18355929×10⁻³⁶` | 119 |

层数是 `ceil(log₂(1/ε))`，只描述一种几何尺度，**不是运行中实际使用的配置、所需均匀点数或所有算法的复杂度下界**。例如尺度 `2.18e-36` 仍可表示为普通浮点数，并非浮点下溢；困难在于用粗而均匀的分区解析相应近零变化。几何分区能用逐层缩小的带宽描述这些尺度。

表格还把全部误差预算交给了近零区，而外部被假设为精确。在真实分段多项式中，外部各区误差和向外取整余量也必须计入；不能拿表中的 ε 直接宣布整个算法达标。若仅给核心一部分平方误差预算，例如 `τ²/2`，需把公式中的 `τ²` 替换为 `τ²/2`。

## 5. 当前代码的核心区其实比“丢掉”更好

当前 `-1/4` 模型在最靠近零的区间使用**最佳常数 L² 投影**，并非零近区。构造器明确以 `degree=0` 调用核心投影，见 [piecewise_models.py:210](/PROJECT/outputs/geometry_precision_followup/piecewise_models.py:210)。一般 p 时，在 `[0,ε]` 的最佳常数为

\[
c_\varepsilon=b+\frac{A}{1-p}\varepsilon^{-p}.
\]

由“源函数平方积分减去常数投影平方积分”，核心误差精确为

\[
E_{\rm const}^2(\varepsilon)
=A^2\left(\frac1{1-2p}-\frac1{(1-p)^2}\right)\varepsilon^{1-2p}
=\frac{a^2p^2}{(1-2p)(1-p)^2}\varepsilon^{1-2p}.
\]

当 `p=1/4`，系数正好为 `2/9`；旧网格 `ε=r^(4L)` 又把幂变成 `r^(2L)`，与源码逐项吻合。这也解释了为什么把专用 `2/9` 直接沿用到其他指数不正确。

在仍假设外部精确、全部预算给核心的前提下，最佳常数核心允许的 ε 比零近区更大：H11 约 `4.13800458e-19`，H12 约 `1.99008848e-23`，H13 约 `1.98522598e-34`，分别约放大 90.9、113.0、90.9 倍。它改善系数，却仍保留 `ε^(1−2p)` 的衰减幂。这里是数学推导，不是当前代码已支持这些新指数的证据。

## 6. 已有一般指数工具与表示方式的差别

不能把当前入口限制说成“整个历史代码从未支持一般指数”。较早的 [function_models.l2_model](/PROJECT/outputs/geometry_general_v235_v237/function_models.py:228) 已支持 `α>-1/2` 的全局 Legendre 投影，最高次数为 24。但当前自动逼近路线选择的是 `-1/4` 分段模型，当前 [结论格式登记](/PROJECT/outputs/geometry_strategy_v238_v317/verification.py:13) 也没有登记旧的 `sphere_axis_profile_projection_v1` 格式。

本次只读调用那份未改动的旧全局投影，对 degree=24、80 位向外平方根上界做了一次诊断，并用其原有验证器复查；没有替换正式评测结果：

| 题目 | 旧全局 24 次投影的 L² 误差约值 | 目标 |
|---|---:|---:|
| H11 | `0.0590711` | `1e-6` |
| H12 | `0.0690383` | `1e-8` |
| H13 | `0.106328` | `1e-10` |

旧证书保存精确投影余量 `error_squared`，且三题均严格大于目标平方；这不是仅凭一个偏大的上界就宣布失败。利用正交投影最优性，可以说该一维全局多项式子空间内次数不超过 24 的模型达不到这三道目标；不能推成更高次数或其他表示也不可行。精确分数及复查结果已写入附带 JSON，计算依据见 [function_models.py:237](/PROJECT/outputs/geometry_general_v235_v237/function_models.py:237)。

分段表示在概念上更匹配这里的误差分布：把近赤道的奇异带单独计误差，在远离零的每个区间对光滑函数做局部逼近。对一般 p，需要的基本积分是

\[
\int_\ell^r t^{k-p}\,dt=\frac{r^{k+1-p}-\ell^{k+1-p}}{k+1-p}.
\]

若 `p=m/d` 为既约有理数，把端点表示成有理数的 d 次幂，则上述幂对应整数指数 `d(k+1)−m`；平方源项对应 `d−2m>0`。这说明保持有理精确积分在数学上有自然途径。单元缩放应使用 `B^(1−2p)`，而非固定旧 `R²`。这些关系解释可扩展性，不保证某个尚未实现的分区、次数、证书大小或十秒预算一定够用。

另一种表示是显式保留奇异项 `|z|⁻ᵖ`，只近似剩余平滑部分。对于这三道纯幂输入，若允许该奇异函数本身成为输出基函数，加上常数项便在代数上精确表示原函数；**这改变了输出表示的允许范围**，不能当成现有分段多项式接口已经解决了任务。仍需独立定义并验证相应表示的误差语义，也不能把源函数表示精确误认作谱问题已解。

因此，本次诊断可确认的是：三题的正式失败属于已知实现路线的支持缺口；严格精度还暴露了单一低次全局多项式的表示局限。L² 可积性与精确积分公式没有构成数学上的不可逼近障碍。是否能在既定预算内扩展现有认证流程，仍需另一个明确授权的实现与验收任务；本次没有开展该迭代。
