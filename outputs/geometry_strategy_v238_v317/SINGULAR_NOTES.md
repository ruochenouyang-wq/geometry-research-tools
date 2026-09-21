# 自动分数幂候选：V268–V277

本轮把旧算例中的人工幂表改为从算子结构生成的候选，并增加失败诊断和参数点迁移。策略仍由人设计；没有证明通用自动研究能力、跨模型优势或新数学定理。十个版本的真实源码、验证和比较按执行时序保存在 `iterations/v268` 至 `iterations/v277`，历史算法目录只读。

## 原题与适用范围

处理单位球面、球面概率测度、全部实平均零 $H^1$ 函数中的最低谱值，势为

\[
q$t$=a|t|^\alpha+b,\qquad -\tfrac12<\alpha<0,
\]

其中 $a,b,\alpha$ 为有理数，$t$ 是指定坐标轴的高度。当前完整尾比较还要求

\[
\eta\ge |a|\sqrt{\frac1{2\alpha+1}-\frac1{$\alpha+1$^2}},\qquad\eta<1.
\]

即使满足此条件，Temple 的谱隙也可能不足。此时 `solve` 返回经过重放的原题源区间和 `certified_open`；不会把局部试探的精度当成全局完成。

试探函数为一个实 $m=1$ 分量中的偶径向函数：

\[
u$t,\phi$=\sqrt{1-t^2}\sum_s c_s|t|^s\cos\phi.
\]

所有幂满足 $s=0$ 或 $s>3/2$，从而没有赤道导数跳跃，完整强算子作用属于 $L^2$。极点处的平方根和角向因子合并为光滑坐标因子。sin 分量与 cos 分量等谱；两者不会合并排序后误用第二特征值。

## 核心数学

令 $\psi_s=\sqrt{1-t^2}|t|^s$，则

\[
H\psi_s=\sqrt{1-t^2}\big[((s+1)(s+2)+b)|t|^s-s(s-1)|t|^{s-2}+a|t|^{s+\alpha}\big].
\]

先删零系数项，再用

\[
\int_{-1}^{1}(1-t^2)|t|^r\,dt=\frac4{(r+1)(r+3)},\quad r>-1
\]

构造精确有理矩阵 $M,A,R$，分别表示质量、算子形式和完整强作用的 Gram 矩阵。设

\[
\mu=\frac{c^TAc}{c^TMc},\quad
\sigma^2=\frac{c^TRc}{c^TMc}-\mu^2,\quad d=Ac-\mu Mc.
\]

投影内残差为 $d^TM^{-1}d/(c^TMc)$，其余部分为完整残差减去该值。若候选使用降维列变换 $C$，投影矩阵必须改为 $C^TMC$，残差坐标改为 $C^Td$。该分解只指导搜索；最终认证仍使用完整 $\sigma^2$。

偶势保持径向奇偶子空间。令中心 (\bar q=a/$\alpha+1$+b)，已有嵌入估计给出形式下界 $$1-\eta$E+$\bar q-\eta$M$。偶 $m=1$ 子空间的第二自由谱值为 12，奇子空间的第一自由谱值为 6。因此

\[
\beta_e=12$1-\eta$+\bar q-\eta,\quad
L_o=6$1-\eta$+\bar q-\eta,
\]

\[
L_{m=1}=\max\left\{L_{m=1,\mathrm{source}},
\min\left$\mu-\frac{\sigma^2}{\beta_e-\mu},L_o\right$\right\},\quad\mu<\beta_e.
\]

再与源证书中其他角向分量及无限角向尾下界取最小值，获得完整原题下界。上界取源上界与试探 Rayleigh 商的较小者。源证书本身重新核验，不假设最低分量始终是 $m=1$。数学依据沿用旧方法的 [球面嵌入式 (1)，取 d=2,p=4](https://arxiv.org/pdf/1210.1853) 和 [Temple 不等式，Theorem 4.13，印刷页 120](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf)。

## 十步实际结果

| 版本 | 实际改变 | 对照和限制 |
|---|---|---|
| 268 | 自动枚举生成元 $2,2+\alpha$ 的幂半群 | 固定原题恢复原幂表，未声称固定题精度增益；不宣称级数收敛定理。 |
| 269 | 合并完整残差，从幂 $r$ 提议修正幂 $r+2$，限额比较 | 四项小实验比固定顺序差：宽度约 $1.3422\cdot10^{-3}$，对照 $8.9420\cdot10^{-4}$。保留为可选策略。 |
| 270 | 精确固定首项比例 (c_{2+\alpha}/c_0=a/[$2+\alpha$$1+\alpha$]) | 同四项宽度约 $4.9486\cdot10^{-4}$；同时少了一个自由系数，不能当作同自由度优势。 |
| 271 | 内外残差分解决定系数迭代停止 | 1/3/6 项实验由各 14 次降至 1/1/3 次；宽度从约 $3.4321\cdot10^{-6}$ 变为 $3.4337\cdot10^{-6}$，诊断本身也耗时。 |
| 272 | 以认证下界生成严格低于谱底的逆迭代移位 | 同两次迭代，六项宽度从约 $4.0759\cdot10^{-6}$ 降至 $3.4321\cdot10^{-6}$。移位不改变被认证的算子。 |
| 273 | 精确质量正交化及合同变换 | 已测 32 bit、十项例没有改善，且略慢；不默认开启。 |
| 274 | 对 $R-2\mu A+\mu^2M$ 做限额正则化逆步，按完整 Temple 宽度保留 | 四项例从约 $8.9421\cdot10^{-4}$ 改善到 $3.2652\cdot10^{-4}$。不声称有限维全局最优。 |
| 275 | 偶分量 Temple 加奇分量完整排除 | 同六项候选宽度从约 $3.4321\cdot10^{-6}$ 降至 $7.1303\cdot10^{-7}$。 |
| 276 | 幅度与常数平移，原函数的 $H,H^2$ 全绑定 | 常数平移的源端点精确平移，无需重新二分；残差与宽度不变。正负幅度均有点测试。 |
| 277 | 推广有理奇异指数及适用性判断 | 支持新参数点，并保留谱隙不足和其他分量阻塞的实际失败。不是整个连续参数盒的统一证书。 |

## 精选最终对照

原题 $a=-1,\alpha=-1/4,b=0$，80 bit、4 次逆迭代：

| 配置 | 全局区间宽度，约 | 本次完整运行秒数 |
|---|---:|---:|
| 8 项，默认奇偶证明 | $3.88892\cdot10^{-9}$ | 0.092 |
| 相同 8 项候选，旧完整径向谱隙 | $1.87189\cdot10^{-8}$ | 仅证明对照，未另算搜索时间 |
| 8 项，加一次完整残差优化 | $3.26250\cdot10^{-10}$ | 0.106 |
| 10 项，默认奇偶证明 | $2.09496\cdot10^{-11}$ | 0.090 |

这些是单次本机观测，源缓存有复用；不能换算为跨模型、token 或通用速度优势。完整有理证书与精确成本记录位于 `iterations/v277/FINAL_VALIDATION.json`。

新参数点的十项、四次迭代检验：

- ($\alpha,a$=(-1/8,-1))：宽度约 $4.18995\cdot10^{-13}$，达到 $10^{-8}$。
- ($\alpha,a$=(-1/3,-1/2))：宽度约 $2.50430\cdot10^{-11}$，达到 $10^{-8}$。
- ($\alpha,a$=(-3/8,-1/2))：完整残差平方已约 $2.77933\cdot10^{-10}$，全球宽度仍约 0.04388。奇分量粗界使新的偶分量 Temple 不能完成原题，最终保留原 $m=1$ 源下界。此时应加强证明，不应继续把增项当成主要解法。
- ($\alpha,a$=(-1/3,-1))：虽然 $\eta<1$，现有 Temple 谱隙不足，返回原题的已认证宽区间。

以上参数在本轮已参与开发和检查，不能再当作未接触的最终测试集。

## 接口与验证

```python
import adaptive_singular as s
run = s.solve(s.ORIGINAL, tolerance='1/100000000',
              max_terms=8, precision_bits=80, iterations=4)
assert s.verify(run['certificate'], expected_function=s.ORIGINAL,
                expected_mean_zero=True, expected_tolerance='1/100000000')
```

可选预算为 `source_modes/source_bits`、`basis_policy='nested'|'residual'`、`candidate_budget`、`leading_cancellation`、`adaptive_iterations`、`shift_policy='certified'|'none'`、`orthogonalize`、`residual_steps=0..4`、`gap_policy='even_odd'|'full_radial'`。默认采用自动半群、可靠移位和奇偶证明；其他候选策略按需启用。

公开数学接口包括 `operator_terms`、`trial_matrices`、`statistics`、`residual_terms`、`residual_diagnosis`、`cancellation_transform`、`mass_orthogonal_transform`、`gap_bound`、`certificate`。`verify` 从原函数、试探和源证据重新构造证明，不调用候选搜索。`iterations=0` 只使用初始候选，不伪造收敛。只支持完整平均零原题，不接受轴对称范围冒充全局范围。

本地 25 项测试覆盖弱点、解析校准、定义域、真实残差、参数绑定、变换、谱隙和篡改。独立审查另见 `test_independent_math.py` 与 `INDEPENDENT_MATH_REVIEW.md`。本模块没有形式化证明助手背书。
