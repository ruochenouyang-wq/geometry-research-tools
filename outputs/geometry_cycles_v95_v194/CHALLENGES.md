# 十轮候选压力题库：冻结 V94 的起点

固定清单 SHA-256：`e583a5bf8d028b968b832c3b71c8fd8580ad09a27e5d3994b01f0ab3ec1c3e05`。共 10 类、每类 2 个输入；机器可读原题见 [challenges.json](challenges.json)，实际初始探测及完整证书见 [baseline_probes.json](baseline_probes.json)。本题库不修改 round01 的已定基线，也不预判 V95–V194 的任何一版成功。

所有谱问题在单位 S² 上，t=z=cosθ，默认概率测度 dμ=dσ/(4π)、完整实 H¹ 零均值空间；非恒定势下是压缩二次型 Π(−Δ+q)Π。谱 Rayleigh 商不受面积归一化影响，GN 常数会改变。m>0 在完整实谱中有两份重数；固定扇区第 k 值不能直接当成全空间第 k 值。

这些是已知理论校准和自行构造的具体压力题，不是开放问题清单。成功标准要求真实全域或全谱证明；拒绝、缺接口、有效宽区间和资源用尽都保留。初始探测仅用于定位瓶颈，不作新旧速度优势声称。

## 1. 显式违反者

**C01 · 零势显式违反者**（精确校准）。

输入：`{"q": ["0"], "threshold": "21/10", "trial": "u=z"}`。

- Return an exact nonzero mean-zero trial and directly replay its moments: mass=1/3, energy=2/3, Rayleigh=2.
- Verify threshold*mass-energy=1/30>0; spectral existence alone is insufficient.

冻结 V94 起点：weighted_poincare is applicable but records existence-only refutation.

**C02 · 椭球势的严格有理违反者**（构造压力题）。

输入：`{"q": ["0", "0", "2"], "threshold": "12/5", "trial": "u=x*(1-z^2/10)"}`。

- Independently integrate the displayed polynomial trial on S2 and prove mass>0, mean=0 and Rayleigh<12/5.
- Report the actual positive violation margin with original potential and full-sphere scope.

冻结 V94 起点：Existence refutation is supported; no explicit function in V94 output.

## 2. 固定势极高精度

**C03 · 固定m=1势的二十四位证书宽度**（构造压力题）。

输入：`{"q": ["0", "0", "2"], "m": 1, "mean_zero": true, "k": 1, "target_width": "1/1000000000000000000000000", "initial_modes": 8, "bits": 100}`。

- Replay an exact interval of width <=1e-24, or report the remaining open width and resource budget.
- The independent Usov rounding band [2.382655406155,2.382655406165] is a compatibility check only; it cannot certify 24 digits.

冻结 V94 起点：Fixed sector supported; probe fixed N=8,bits=100 to distinguish truncation error from arithmetic precision.

**C04 · 强势全角向二十位闭合**（构造压力题）。

输入：`{"q": ["0", "0", "100"], "target_width": "1/100000000000000000000", "initial_modes": 8, "initial_max_modes": 8, "initial_max_m": 4, "bits": 100}`。

- Replay a full-space interval of width <=1e-20 including an explicit omitted-angular-sector lower bound, or keep certified_bound_open_gap.
- Do not infer high precision merely from bits=100 or reuse an m=0-only interval as a full-space result.

冻结 V94 起点：full_ground supports the input; limited truncation is deliberately probed.

## 3. 高于六次的势

**C05 · 八次非负势**（构造压力题）。

输入：`{"q": ["0", "0", "0", "0", "0", "0", "0", "0", "1"], "target_width": "1/100000000"}`。

- Accept degree 8 with exact coefficient binding and global potential range.
- Coarse independent guard: 2 <= full mean-zero ground <=67/33, using u=x and E(z^(2r))=1/(2r+1); requested narrow result must independently replay.

冻结 V94 起点：V94 degree cap 6 rejects this input.

**C06 · GN六次试探平方产生十二次势**（构造压力题）。

输入：`{"q": ["-1/49", "0", "0", "0", "0", "0", "2/7", "0", "0", "0", "0", "0", "-1"], "construction": "q=-(t^6-1/7)^2", "target_width": "1/100000000"}`。

- Bind the expanded degree-12 polynomial exactly to the factored trial square.
- Independent guard q>=-36/49 gives full mean-zero ground >=62/49; a successful spectral certificate must cover the whole interval, not samples.

冻结 V94 起点：V94 degree cap 6 rejects this natural trial-square input.

## 4. 参数区间统一谱界

**C07 · 整个线性参数区间的统一正下界**（构造压力题）。

输入：`{"q(a,t)": "a*t", "parameter_box": {"a": ["-4", "4"]}, "claimed_uniform_threshold": "1"}`。

- Prove lambda_mean_zero(a*t)>=1 for every a in [-4,4].
- One admissible proof uses concavity of the constrained ground energy in a plus verified endpoint lower bounds; a finite grid without a remainder is insufficient.

冻结 V94 起点：Pointwise endpoints supported; no parameter-box certificate API.

**C08 · 非负参数盒与端点精确达到**（构造压力题）。

输入：`{"q(a,t)": "a*t^2", "parameter_box": {"a": ["0", "100"]}, "claimed_uniform_threshold": "2"}`。

- Prove the uniform threshold 2 from q>=0 and Poincare, with equality at a=0 and u=x.
- Bind the full parameter box; sampling the endpoints must not be presented as the general proof.

冻结 V94 起点：Fixed coefficients supported; uniform quantifier and exact minimizer are outside the API.

## 5. 参数全局最优与分支竞争

**C09 · 带惩罚的参数全局最优**（构造压力题）。

输入：`{"q(a,t)": "a*t+t^2", "parameter_box": {"a": ["-4", "4"]}, "objective": "lambda_full_mean_zero(q(a,t))+a^2/5", "target_global_gap": "1/100000000"}`。

- Return rigorous lower and attained/trial upper bounds for the global objective with gap <=1e-8, or retain an explicit open gap.
- Cover the whole parameter interval and all relevant angular branches; attach proof to any reported minimizer enclosure. No optimal parameter/value is assumed in advance.

冻结 V94 起点：A fixed a can be evaluated; no global parameter optimizer/certificate API.

**C10 · 正负椭球势的角分支竞争**（构造压力题）。

输入：`{"q(a,t)": "a*t^2", "parameters_to_compare": ["-2", "0", "2"], "parameter_box_for_optional_global_search": {"a": ["-2", "2"]}}`。

- At a=-2 and a=2 establish competing m=0/m=1 intervals and identify a winning branch only after separation; retain undecided labels if they overlap.
- At a=0 prove full mean-zero eigenvalue 2 with real multiplicity 3 (one m=0 plus two m=1). Do not infer uniqueness or a whole-box branch ordering from only these three points.

冻结 V94 起点：Pointwise sectors and full ground supported; branch-family claims require new machinery.

## 6. 额外正交约束

**C11 · 删除完整l=0和l=1空间**（精确校准）。

输入：`{"q": ["0"], "orthogonality": ["1", "x", "y", "z"]}`。

- Ground energy is exactly 6; prove all constraints and accept an l=2 trial such as 3z^2-1.
- Excluding only m=0,l=1 is insufficient because the x/y modes remain.

冻结 V94 起点：Only the constant constraint is available in V94.

**C12 · 只删除常数和z仍保留横向模**（精确校准）。

输入：`{"q": ["0"], "orthogonality": ["1", "z"]}`。

- In the full sphere the ground is exactly 2, with admissible x and y.
- For the separately declared axisymmetric space the same constraints give ground 6; scope binding must distinguish the two.

冻结 V94 起点：No extra-constraint API; full and axisymmetric meanings must not be conflated.

## 7. 完整高阶谱、重数与负谱和

**C13 · 完整高阶谱与重数**（精确校准）。

输入：`{"q": ["0"], "mean_zero": true, "ordered_full_eigenvalues_k": "1..9"}`。

- Return [2,2,2,6,6,6,6,6,12], including real angular multiplicities.
- Each sector may be certified independently, but a full-spectrum merge must also exclude omitted levels/sectors.

冻结 V94 起点：Fixed-sector kth eigenvalue is supported; full_ground only returns full-space k=1.

**C14 · 负谱个数及负部迹**（精确校准）。

输入：`{"q": ["-7"], "spaces": ["full_mean_zero", "full_unprojected"], "quantity": "sum_j max(0,-lambda_j), with multiplicity"}`。

- Mean-zero result: negative count 8 and negative-part trace 20 (l=1:3*5; l=2:5*1).
- Unprojected result: negative count 9 and trace 27; prove no omitted negative eigenvalues. These are exact test identities, not a new Lieb-Thirring bound.

冻结 V94 起点：Sector energies are available; full counting and negative-trace sum API are absent.

## 8. GN 变分与集中试探

**C15 · GN候选常数的直接反例**（精确校准）。

输入：`{"inequality": "E_mu(u^4)<=C*E_mu(u^2)*E_mu(|grad u|^2)", "C": "4/5", "trial": "u=z", "mean_zero": true}`。

- Direct moments are E(u^2)=1/3, E(|grad u|^2)=2/3, E(u^4)=1/5; quotient=9/10>4/5.
- If using area d_sigma instead, the same trial quotient is 9/(40*pi), not 9/10. A weighted linear spectral threshold is not a GN certificate.

冻结 V94 起点：V94 certifies linear quadratic forms, not this quartic functional.

**C16 · 明确集中序列的有理矩与GN试探**（构造压力题）。

输入：`{"trial_family": "u_n(t)=(1+t)^n-2^n/(n+1)", "n_values": [4, 16, 64], "normalization": "mu=d_sigma/(4*pi)", "mean_zero": true}`。

- For each prescribed n, independently integrate mean, M2, M4 and Dirichlet energy and return the exact rational GN quotient.
- M2=2^(2n)*(1/(2n+1)-1/(n+1)^2); energy=n*2^(2n-1)/(2n+1). M4 follows binomial expansion with E((1+t)^j)=2^j/(j+1).
- All three cases stay in the record even if a resource limit is reached. No claim that these trials are globally optimal or attain the sharp constant.

冻结 V94 起点：No GN moment/search interface; potential squares would exceed degree 6 for these trials.

## 9. 非轴对称势

**C17 · 旋转等价的非轴输入q=x**（构造压力题）。

输入：`{"q_xyz": [{"coefficient": "1", "powers": [1, 0, 0]}], "mean_zero": true, "target_width": "1/1000000000000"}`。

- Certify the full mean-zero ground and reproduce the q=z result under an explicit rotation of S2.
- Preserve the mean-zero constraint and full-space scope; comparing only the old m=0 sector is not a rotation test.

冻结 V94 起点：The axisymmetric q(t) API cannot directly represent x; it can supply the rotated calibration q=z.

**C18 · 真正三轴二次势**（构造压力题）。

输入：`{"q_xyz": [{"coefficient": "1", "powers": [2, 0, 0]}, {"coefficient": "2", "powers": [0, 2, 0]}, {"coefficient": "3", "powers": [0, 0, 2]}], "mean_zero": true, "target_width": "1/100000000"}`。

- Replay a full-sphere certificate for x^2+2y^2+3z^2, which has no rotation making it a function of one axis alone.
- Independent coarse guard is 3 <= ground <=18/5: q>=1 and the trial u=x has potential quotient 8/5 and Dirichlet quotient 2.

冻结 V94 起点：Nonaxisymmetric coefficients/couplings are outside V94 input scope.

## 10. 跨工具串联与续算

**C19 · 同一谱证书支持两个后续判定**（构造压力题）。

输入：`{"q": ["0", "0", "2"], "thresholds": ["119/50", "12/5"], "workflow": ["certify", "persist_by_content_hash", "reload_and_verify", "derive_two_threshold_decisions"]}`。

- One replayed full-space certificate with L>=119/50 and U<12/5 yields proved and refuted-by-spectral-existence respectively.
- Preserve certificate identity and avoid recomputing the spectrum in both downstream steps; any explicit violating-function extension must be separately checked.
- Reject changed q or scope bindings and corrupt content hashes.

冻结 V94 起点：V94 already supports content-addressed storage, verify and fetch. Its weighted_poincare operation recomputes the spectrum and takes no existing certificate id. This challenge extends reuse, not basic persistence.

**C20 · 开放证书的预算续算**（构造压力题）。

输入：`{"q": ["0", "0", "100"], "first_budget": {"modes": 4, "max_modes": 4, "max_m": 0, "bits": 60, "tolerance": "1/100000000000000000000"}, "target_width": "1/100000000000000000000", "continuation_budget": "explicitly recorded larger mode/angular budget"}`。

- First preserve and replay the honest open certificate. A resumed result must keep best lower=max(old lower,new lower) and best upper=min(old upper,new upper) for the same operator.
- Record actual reused work and certificate ancestry; a fresh full recomputation may be a correct fallback but is not evidence of warm-start reuse.
- If the larger budget still cannot close 1e-20, retain open status. Reject an inherited certificate with changed coefficients or projection.

冻结 V94 起点：V94 has persisted certificates but full_ground has no seed/resume option.

## 实际初始探测

| 探测 | 对应题 | 结果 | 解释 |
|---|---|---|---|
| P01 | C01 | refuted_by_spectral_existence | 独立重放 True；宽度≈0 |
| P02 | C02, C19 | refuted_by_spectral_existence | 独立重放 True；宽度≈7.529e-11 |
| P03 | C03 | certificate_returned | 独立重放 True；目标宽度满足=False；宽度≈2.338e-15 |
| P04 | C04 | certified_bound_open_gap | 独立重放 True；宽度≈0.01084 |
| P05 | C05 | rejected_or_computation_failed | 势函数 q 的系数按升幂排列，次数最多为 6。 |
| P06 | C06 | rejected_or_computation_failed | 势函数 q 的系数按升幂排列，次数最多为 6。 |
| P07 | C07 | certified_target_met | 独立重放 True；宽度≈9.007e-12 |
| P08 | C07 | certified_target_met | 独立重放 True；宽度≈9.007e-12 |
| P09 | C10 | certified_target_met | 独立重放 True；宽度≈1.008e-13 |
| P10 | C08, C10 | certified_target_met | 独立重放 True；宽度≈0 |
| P11 | C10 | certified_target_met | 独立重放 True；宽度≈7.529e-11 |
| P12 | C13 | certificate_returned | 独立重放 True；宽度≈0 |
| P13 | C14 | certified_target_met | 独立重放 True；宽度≈0 |
| P14 | C17 | certified_target_met | 独立重放 True；宽度≈1.099e-13 |
| P15 | C20 | certified_bound_open_gap | 独立重放 True；宽度≈34.36 |

C02 的独立有理试探检查：Rayleigh=`24076/10089`，阈值差=`688/50445`，严格为正。该检查说明题目存在可验收的具体答案，不代表冻结 V94 已能输出该函数。

参数盒、额外约束、GN、非轴势等缺接口结论来自已保存的方法签名和 JSON 服务字段检查；未伪造不存在的调用结果。C19 特别保留 V94 已有的内容寻址存储与 verify/fetch 能力，其压力点是已有谱证书的后续复用；C20 压力点是继续预算与保留旧界。

## 来源和参考边界

- [DLMF 14.30: Spherical and Spheroidal Harmonics](https://dlmf.nist.gov/14.30)：Unit sphere harmonic eigenvalues and angular basis; exact constant-potential calibrations also derived directly.
- [Usov, Spheroidal quantum well, arXiv:2307.04124v1, Eq.7a and Table1c](https://arxiv.org/pdf/2307.04124v1)：For q=+c^2 t^2, m=1, c^2=2, l=1: Lambda=2.38265540616. Only compatibility with the fixed rounding interval plus/minus 0.5e-11; not a rigorous high precision enclosure.
- [DLMF 30.2: Differential Equations](https://dlmf.nist.gov/30.2)：Convention conversion Lambda=lambda_DLMF+c^2 for q=+c^2 t^2. Positive c^2 here is prolate.
- [Ilyin and Zelik, On a class of interpolation inequalities on the 2D sphere, arXiv:2204.12414v1](https://arxiv.org/pdf/2204.12414)：Eq.1.8 gives a nonsharp scalar zero-mean interpolation upper bound. At exponent r=4 its probability-measure form is E(u^4) <= 4 E(u^2) E(|grad u|^2). It does not identify the sharp GN constant. The paper also motivates orthonormal-density / spectral trace work.

Usov 的既定 11 位十进制带不能用于证明 20–24 位目标；DLMF 的特征参数换算明确写入清单。Ilyin–Zelik 的归一化 GN 上界只作理论背景，不声称其为 sharp 常数；负谱和例子来自常数势球谐精确谱，未套用未经核对的 sphereLT 常数。

后续每轮应先选择并记录题目及资源，再实施变更。需要改变输入、容差或验收条件时保留本版本和失败结果，明确说明原因；不能事后筛掉失败题。
