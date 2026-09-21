# 几何误差估计工具：第二轮

先读 [误差分析与函数选择](ERROR_ANALYSIS.md) 和 [实测比较](EXPERIMENTS.md)。本包独立运行，只需 Python 3.9+ 标准库。

在本文件所在目录执行：

```sh
python3 run.py bound --q "0,1" --tolerance 1e-10 --out results/example.json
python3 run.py verify results/example.json
python3 -m unittest -v
```

所有计算仍针对单位二球面上 `-Δ+q(cos θ)` 的轴对称子空间。新方法计算最低特征值，包含 q=0 时的零模态。q 系数从常数项起排列，只接受整数或分数，次数最多六；负号开头使用 `--q="-1,2"` 写法。

## 方法

| 参数 | 用途 |
| --- | --- |
| `--method auto` | 先检查简单指数形状，随后使用残差与对称性，必要时尝试 Schur 方法 |
| `--method temple` | 完整残差平方与经过认证的相邻特征值下界，不做对称性约化 |
| `--method symmetry` | 对偶 q 使用偶函数空间的残差认证；其他 q 使用普通 Temple |
| `--method exponential` | 搜索 `exp(s(t))`，对其局部能量做严格全域包络 |
| `--method banded` | 保留原误差公式，使用更快的带状精确消元 |
| `--method baseline` | 第一轮的原始算法，用于比较 |

`--tolerance` 指认证区间的**宽度**，使用精确十进制或分数输入，如 `1e-10`。若取区间中点，其绝对误差不超过宽度的一半。达到目标返回 0；预算内未达到目标返回 2，保留当前最佳证书和尝试记录。有效但较宽的区间不能伪装成达到精度目标。

普通与对称方法的 `--max-modes` 均限制原始 Legendre 索引 `0,...,N-1`，默认 24、最大 32。对称方法会跳过奇索引。指数方法单独搜索 log 函数次数 1 至 16，不使用此模式数上限；此差异在基准中单独标注。

```sh
python3 run.py bound --q "0,0,-100" --method symmetry --out results/even.json
python3 run.py bound --q "0,20,-100" --method exponential --out results/exponential.json
python3 benchmark.py
```

`benchmark.py` 会重新生成 `results/` 内的测量及代表性证书。它比较相同 `1e-10` 区间宽度目标，记录到达目标前的全部搜索开销，包括失败尝试。毫秒计时不包括 Python 启动或文件写入，因此不能直接等同于点击运行的总等待时间。

显示小数使用 12 位向外舍入；当有理证书比 `1e-12` 更窄时，请查看 JSON 中的 `lower`、`upper` 和 `exact_width`，不能用显示小数反推出全部认证精度。

## 主要文件

- `error_bounds.py`：残差、对称性、正试探函数、全域多项式包络、带状惯性和复核器。
- `spectral_certifier.py`：第一轮的原始实现，原样保留作基线与交叉检查。
- `run.py`：有明确误差目标的搜索与证书复核。
- `benchmark.py`：全部性能实验；原始记录在 `results/benchmark.json`。
- `test_error_bounds.py`：30 项正确性及错误证书测试。
- `SCHUR_PROOF.md`：第一轮的模型定义与尾部证明。

本轮未调用在线 GPT API。GPT 可以提出候选函数和研究方向，程序按数学条件检查候选。浮点求解器只负责提出候选；数学接受条件使用有理数重算，检查器仍依赖本地实现和文档中的解析证明，并非形式化证明助手。
