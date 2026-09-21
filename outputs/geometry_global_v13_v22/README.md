# 几何研究工具：V13–V22

这十版继续研发能给出全局误差证明的研究工具。它们适合检验候选函数、判断局部极小值是否值得相信、处理更多连续参数及约束，而不是直接宣称解决某个开放几何问题。

先读 `REPORT.md` 看实测结论；`ITERATIONS.md` 列出每版改动；完整公式与适用条件见 `PROOF_V13_V22.md`。所有结果和未成功算例位于 `results/`。

## 能做什么

| 任务 | 入口 | 版本 |
| --- | --- | --- |
| 1–2 参数盒上的函数能量全局优化 | `global_core.search` | V13–V15 |
| 全实参数与 t、t² 矩关系排除 | `global_core.all_real_problem` + `search(domain='all_real')` | V16 |
| 3–4 连续参数与代数下界细化 | `global_core.search` | V17–V18 |
| 线性等式约束下的全实参数优化 | `equality_v19.search` | V19 |
| 连续不确定参数下的鲁棒谱设计 | `robust_v20.search` | V20 |
| 对整个函数空间施加一个期望值约束 | `constrained_v21.search` | V21 |
| 按误差来源选择谱精度、代数细化或参数划分 | `global_core.search(version=22)` | V22 |

V19–V21 改变了问题类型。各入口不应被视为同一目标函数的不同速度档位。底层公共代码可复用；历史 V11–V12 文件保留以便对照。

## 直接运行

Python 3.9 或更高版本，仅标准库，不需要安装依赖、不需要 API 密钥。进入本文件所在目录后：

```sh
python3 run.py demo --version 22 --output result.json
python3 run.py verify result.json
python3 -m unittest discover -p 'test*.py' -v
python3 benchmark.py --repeats 3
python3 run.py verify results/certificates
```

`--version` 可选 13 至 22；每版都有对应的示例。这里没有联网调用 GPT。GPT 可负责提出问题、选择参数化与候选函数，数值工具负责核算和拒绝不成立的全局结论。

一个自定义主问题：

```python
from fractions import Fraction as F
from global_core import problem, search

p = problem(
    q0=[0], directions=[[0, 1]], box=[[-6, 6]],
    linear=[F(1, 50)], hessian=[[F(1, 5)]],
)
r = search(p, tolerance=F(1, 1000), version=22)
c = r['certificate']
print(c['status'], c['global_lower'], c['candidate_upper'])
```

它搜索 `λ₁(−Δ+a t)+a²/10+a/50` 在完整实区间上的全局下确界。函数空间是固定球面的全部轴对称能量型函数，有限 Legendre 计算之外的尾部也在证书中。

## 怎样读结果

- `epsilon_global`：候选能量减去完整域的已证明下界，不超过所给容限。
- `global_gap_open`：上下界有效，但当前预算不足以证明所需精度。
- `feasibility_not_found`：V21 尚未找到满足约束的函数，没有证明不可行。
- `certified_infeasible`：V21 用全域矩范围证明目标不可行。

十版算法均有适用范围。主问题最多 4 个连续参数；鲁棒版最多 2 个设计参数与 2 个不确定参数；约束版当前处理一个期望值不等式。全部基于固定球面轴对称模型，尚未扩展到任意几何。`formal_assistant_checked` 为 false；这些结果经过精确有理算术核验，但尚未进入形式化证明助手。

## 文件

`point_oracle.py` 实现 V13–V14；`simplex_v15.py` 实现单纯形二次最小化；`bernstein_v17.py` 实现 V17–V18；`global_core.py` 接入 V15、V16 和 V22 的区域证明与调度。V19–V21 分别在同名模块。`test_new_versions.py` 包含本轮针对性测试；其余测试回归底层谱和范围算法。

证书文件可以复制到独立目录并单独复核；复核不需要重新做全局搜索。`MANIFEST.sha256` 记录发布文件的内容哈希。
