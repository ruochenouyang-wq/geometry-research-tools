# 分数幂试探与完整残差证书

本模块只处理单位球面上原势 (q(z)=-|z|^{-1/4}) 的全局平均零谱问题。使用贴合赤道奇点的试探函数后，16 项给出约 (1.14555\times10^{-18}) 的严格区间宽度。所有接受的矩、系数、残差和区间端点均为有理数；不是采样误差估计。

## 公式与定义域

在一个实数 (m=1) 角向分量中，取

\[
u(z,\phi)=\sqrt{1-z^2}\sum_s a_s|z|^s\cos\phi.
\]

允许 (s=0) 或 (s>3/2)。后一条件保证赤道处没有导数跳跃产生的 delta 项，且完整强算子作用属于 (L^2)。(s=0) 时先删去零系数项，不能计算无意义的负次幂积分。极点处的平方根与角向因子组合为光滑坐标因子。

令 (w=1-z^2)、(\psi_s=\sqrt w|z|^s)，则

\[
H\psi_s=\sqrt w\big[(s+1)(s+2)|z|^s-s(s-1)|z|^{s-2}-|z|^{s-1/4}\big],
\qquad
\int_{-1}^{1}w|z|^r\,dz=\frac4{(r+1)(r+3)},\quad r>-1.
\]

据此精确构造 (M_{ij}=\langle\psi_i,\psi_j\rangle)、(A_{ij}=\langle\psi_i,H\psi_j\rangle)、(R_{ij}=\langle H\psi_i,H\psi_j\rangle)。角向积分的共同因子消去。候选系数通过有理逆迭代和二进制有理舍入生成；证明重放不运行搜索。

\[
\mu=\frac{a^TAa}{a^TMa},\qquad
\sigma^2=\frac{a^TRa}{a^TMa}-\mu^2.
\]

这里的 (\sigma^2) 是完整强残差，包含试探空间外的部分。首个奇异修正 (a_{7/4}/a_0=-16/21) 恰好抵消算子作用中的 (|z|^{-1/4}) 项，解释了分数幂基的用途。

## Temple 下界与全局范围

球面概率测度下已知 (\|v\|_4^2\le E(v)+M(v))。原势中心为 (-4/3)，中心化势的 (L^2) 范数为 (sqrt{2/9})。取其有理上界 (eta<1)，最小最大原理给出

\[
\beta=(1-\eta)6-\frac43-\eta
\]

作为**单一 cos 分量、包含两种径向奇偶性的 (m=1) 算子第二特征值下界**。sin 分量与其等谱。不能将 cos 与 sin 合在一起排序，再将这个数称作第二特征值下界。嵌入常数见 [Dolbeault 等，式 (1)，取 (d=2,p=4)](https://arxiv.org/pdf/1210.1853)。

当 (\mu<\beta) 时，Temple 不等式给出 (\lambda_{1,m=1}\ge\mu-\sigma^2/(\beta-\mu))。这里使用第二特征值的下界即可，见 [Teschl，Theorem 4.13 及其后说明，印刷页 120](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf)。

模块重放 `certificates/singular_n4.json`，核验原势、球面、平均零约束和所有角向范围。将 Temple 的 (m=1) 下界与源证书的 (m=0)、其他已算分量、未算角向尾下界取最小值，得到全局下界；上界取原上界与 (\mu) 的较小者。有限试探空间并未替代全局函数空间。

## 实测结果

同一嵌套基、有理系数精度 112 bit、14 次逆迭代，保留粗阶未达标结果：

| 项数 | 全局区间宽度，约 | 完整残差平方，约 | 目标 (10^{-8}) |
|---:|---:|---:|---|
| 1 | (1.58277\times10^{-2}) | (2.84298\times10^{-1}) | 未达到，沿用原证书区间 |
| 3 | (1.06788\times10^{-2}) | (8.88116\times10^{-3}) | 未达到 |
| 6 | (3.43211\times10^{-6}) | (2.85446\times10^{-6}) | 未达到 |
| 10 | (1.00839\times10^{-10}) | (8.38670\times10^{-11}) | 达到 |
| 16 | (1.14555\times10^{-18}) | (9.52748\times10^{-19}) | 达到 |

16 项有理证书的一个向外舍入十进制包围为

\[
0.535142576744105895552
\le\lambda_1\le
0.535142576744105896699.
\]

本机三次完整搜索（清空本模块矩缓存，含最终证书重放）的中位耗时：目标 (10^{-8}) 约 0.056 秒，目标 (10^{-16}) 约 0.124 秒。这只反映本题和本次运行，不是模型能力或通用复杂度结论。17 项测试覆盖独立弱能量校准、完整残差积分、定义域、证据篡改、全局范围绑定及禁止搜索的 JSON 回放。

## 复现入口

在本目录运行：

```sh
python3 -B enriched_trial.py
python3 -B -m unittest test_enriched_trial -v
```

默认在达到 (10^{-8}) 后停止，通常使用 10 项。取得 16 项结果：

```python
import json
import enriched_trial as e

with open('certificates/singular_n4.json') as stream:
    source = json.load(stream)
result = e.enrich(source, tolerance='1/10000000000000000')
assert e.verify(result['certificate'], expected_source=source)
```

`trial_certificate` 认证给定有理试探；`verify` 重建矩阵、残差、Temple 界并独立重放源证书。当前成果是特定奇异势的认证工具，不是新不等式定理，也没有证明这组试探基在所有可选基中最优；未经过形式化证明助手核验。
