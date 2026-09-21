# 第 2 轮：V105–V114，精度驱动的全球谱区间

本轮仅新增 `precision.py`、`test_precision.py` 和本目录的阶段、证书与结果。V90、V94、external_trials 全部作为只读依赖。`STAGES.json` 给出十个不同可执行入口；`python3 -B precision.py` 重放并生成逐阶段证据。

## 实际问题和结果

原始 `precision_candidate.json` 的问题是：对零均值全 S² 上 `q=2t²`，目标绝对区间宽度 `1e-18`，V94 默认固定 44 位二分，两个已算扇区都扩到 N32，仍为 open，宽约 `1.5718972078458572e-13`。本轮保持数学算子不变，按容差选择位数、粗界排除不可能胜出的扇区，再由径向残差决定工作。

| 输入 | 本轮验证后宽度 | 最终径向工作 | 状态 |
|---|---:|---|---|
| `q=2t²` | `4.216156861368222e-20` | m0:N4；m1:N9 | 达到 `1e-18` |
| `q=t+t²` | `3.796275121728574e-20` | m0:N4；m1:N9 | 达到 `1e-18` |
| `q=100t²` | `7.625418624961757e-25` | m0:N4；m1:N28；m2:N4；m3、m4仅解析界 | 达到 `1e-18` |
| `q=2t²,max_steps=0` | `6/5` | 解析试函数与角尾界 | 有效 open |

所有完整证书均经 `verify_full` 重放。时间仅为本机运行记录，包含生成和验证；没有模型速度或跨进程性能比较。结论不依赖这些时间。

Temple 并非每例都有收益：`q=2t²`、`q=t+t²` 的最终区间由原 Schur 证据主导，Temple 单独宽度分别约 `6.47e-17`、`2.94e-17`，相交没有进一步缩小。`q=100t²` 的最后 Schur 宽约 `8.85e-21`，Temple 将它缩到 `7.63e-25`，宽度改善约 11606 倍。代码始终取精确交集，保留无收益路径的真实结果。

## 十个增量与证据

| 版本 | 独立入口 | 实际能力与结果 |
|---|---|---|
| V105 | `diagnose_precision` | 重放旧证据后，指出固定 44 位瓶颈，给出增加位数的下一动作；也区分角尾缺口、径向缺口。 |
| V106 | `choose_bits` | 精确比较绝对容差与尺度，自动选择位数；同 N12，44 位宽 `1.57e-13`，68 位宽 `1.41e-20`。 |
| V107 | `temple_certificate` | 完整投影残差加同 m 的第二特征值下界，生成独立 Temple 证书。 |
| V108 | `recover_gap` | 故意给定弱的合法分隔下界 0，重新认证后恢复到 6；对 μ=12 的高阶试函数拒绝 Temple 下界。 |
| V109 | `residual_schedule` | 按径向窗口分割全部残差，尾主导时实际 N4→N6，完整方差从约 `3.79e-4` 降为 `1.12e-7`。 |
| V110 | `refine_coefficients` | 增加有理系数位数，使用精确 LDL 逆迭代；N8 时方差从约 `1.29e-3` 降为 `9.18e-12`。 |
| V111 | `intersect_bounds` | 精确交叉 Schur、Temple 和历史同扇区证据；任一项无收益也不会放宽结果。 |
| V112 | `checkpoint` / `full_ground(saved=...)` | 真正续算混合势案例，验证旧证据、输入与投影绑定，并单调继承上下界。 |
| V113 | `winning_sectors` | 以已证下界剪掉不可能获胜的扇区；q2 的 m0 停在 N4，m1 精化到 N9。 |
| V114 | `full_ground` | 全 S² 零均值绝对精度入口，包含无限角尾；工作不足时交付有效 open。 |

每项证据位于 `results/v105.json` 到 `results/v114.json`，完整映射在 `STAGES.json`。`results/summary.json` 保留三例、旧结果和零预算结果；`certificates/` 中包含完整重放材料。验收测试为 `python3 -B -m unittest -v test_precision.py`，14 项通过。

## 接口与接受条件

```python
full_ground(q, mean_zero=True, tolerance=Fraction(1, 10**18),
            modes=4, max_modes=32, max_m=16, max_steps=16,
            use_temple=True, saved=None)
verify_full(cert, expected_q=None, expected_mean_zero=None)
checkpoint(cert)
verify_checkpoint(saved, expected_q=None, expected_mean_zero=None)
```

新格式分别为 `same_sector_temple_v107`、`precision_sector_intersection_v111`、`precision_checkpoint_v112`、`precision_full_ground_v114`。不能把它们交给旧格式验证器并假称成功。验证器只进行证据重放，不调用浮点提议或求解搜索。`status=certified_target_met` 当且仅当完整全局证据宽度不超过所声明的有理容差；`verify_full=True` 也可以对应诚实的 open 状态。

`max_steps` 是精确谱精化轮次数，不是墙钟时间、CPU 工作量或递归验证次数上限。初始解析界与验证也会消耗有限工作。参数沿用 q 次数至多 6、系数范围和分母限制，m≤64、N≤64，目标容差在 `1e-30..1`。续算链深度最多 8，扇区相交链最多 32，超过即拒绝；更长历史需要外部重新压缩证明，本轮未实现。

## Temple 的严格范围

固定 `q,m,mean_zero` 后令 H 为该扇区的自伴投影算子。m0 且零均值时使用 `P H P`，原算子产生的 l0 项必须在完整作用后投影移除；m≥1 没有此投影。由 R1 `witness.residual_certificate` 精确重放一个非零有理试函数 u，得

`μ=<u,Hu>/<u,u>`，`v=||(H-μ)u||²/<u,u>`。

v 包括 q 作用产生的所有高阶系数；不把有限矩阵残差当作全残差。输入光滑有限球谐组合在算子定义域内。另一条独立的 V94 Schur 证书必须证明同 q、同 m、同零均值投影的 `g<=λ₂`，并且必须严格满足 `μ<g`。

由于 `λ₁<=μ<g<=λ₂`，谱上的 `(H-λ₁)(H-g)` 非负。因此

`0 <= v + (μ-λ₁)(μ-g)`，所以 `λ₁ >= μ-v/(g-μ)`。

这是同一扇区的第一特征值下界，绝不使用全 S² 的第二特征值；后者可能具有 m=±1 的重数，不能提供这里需要的分隔。μ≥g（包括相等）时输出 `gap_accepted=False,lower=None`，仅保留 Rayleigh 上界。坏间隙不会通过增加标签或忽略分母符号变成成功证书。

该推导引入的理论前提是标准自伴算子的谱定理和 Temple 二次多项式论证，并沿用 V94 已披露的球谐完备性、尾部比较和惯性定理；`formal_assistant_checked=False`，没有宣称形式化证明助理已经检查。

## 诊断与续算的边界

自动位数是调度建议，不能单凭一个估计的二分尺度声称误差达标，最终始终检查实际有理区间。残差分割按完整径向窗口进行；稀疏试函数支持集以外的项可能仍在该窗口内，不能全部叫作径向尾。窗口内残差较大不证明它完全来自量化；增加系数位数与精确逆迭代是实际实验，其 Rayleigh 商和完整方差都需精确改善才接受。

checkpoint 的摘要用于检测存储变化，并非数学信任来源。验证器重放嵌套完整证书，逐扇区核对 q、m、零均值绑定，核对 scope 与角尾覆盖，再精确相交新旧全局区间。q 或投影不同、漏掉中间 m、篡改上下界、状态或残差都会拒绝。零预算结果与故意无效 gap 已保留，不会补造达标结果。
