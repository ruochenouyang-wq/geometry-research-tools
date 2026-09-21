# V11–V12：从局部候选到全局最优差距

本轮实现了两个新迭代：**V11 对整个参数域建立下界并进行分支定界；V12 精确消去可解参数，把问题转成带全局证书的非线性函数优化。** 不以梯度小、数值稳定或多次重启一致作为全局成功标准。

先读 [结果与失败例](REPORT.md)，数学论证见 [全局证明](PROOF_GLOBAL.md)。核心算法、测试与复核只需 Python 3.9+ 标准库。

在此目录运行：

```sh
# 默认 V12：覆盖所有实参数与完整的轴对称函数空间
python3 run.py optimize --problem examples/tilted.json --epsilon 1e-3 --out results/example.json
python3 run.py verify results/example.json

# V11：只在输入的参数框内做全局优化
python3 run.py optimize --problem examples/tilted.json --domain box --out results/box_example.json

# V11 的另一条全实参数路线：证明无穷远处不可能更优
python3 run.py optimize --problem examples/tilted.json --all-real-method growth --out results/growth_example.json

# 两个实参数，仍对完整函数空间给出全局界
python3 run.py optimize --problem examples/two_parameters.json --epsilon 1e-2 --out results/two_example.json

python3 reproduce.py
```

`--domain all-real` 为默认值，会由数学证明确定一个充分搜索框，输入文件中的 `box` 只供显式选择 `--domain box` 时使用。`--bound lipschitz` 切换到较粗的区域下界作对照。`--max-leaves` 和 `--modes` 分别限制参数分区与谱空间，默认 256 和 10。`--epsilon` 控制目标值距全局最优的差距。

## 输入与输出

输入的势是 `q(a,t)=q0(t)+Σ aᵢ pᵢ(t)`。JSON 中 `directions` 的每行是一个 pᵢ 的多项式系数，从常数项起排列。

代价约定为 `constant + linearᵀ a + 1/2 aᵀ hessian a`，不是省略 1/2 的另一种约定。全实参数的 V12 路线要求 Hessian 严格正定；框内路线可以处理非凸二次代价。未满足全实参数路线的充分条件时会报错，不会宣称问题无界。

输出包含：

- `candidate_upper`：一个具体有理试探函数所实现的目标值。
- `global_lower`：对全部允许参数、全部允许函数都有效的下界。
- `global_gap`：两者差值。
- `epsilon_global`：差距不超过输入 epsilon；`global_gap_open`：预算内尚未达到。
- 完整参数覆盖树、原始谱证书、候选函数，以及扩展到所有实参数时使用的附加证明。

终端小数仅为近似展示；JSON 有理数是正式数值。计算达标返回 0，预算内未达标返回 2，输入或证明失败返回 1。复核通过表示声明与证据一致，因此一个有效的 `global_gap_open` 结果也可以通过复核。

“全局”指明确定义的数学目标与函数空间。当前模型是单位球面上的轴对称 `−Δ+q(cos θ)`；不是所有曲面、所有函数族、所有研究策略或最快算法的全局最优。相关固定模型的旧证明保留在 `BASE_MODEL.md`；旧版测试和源码也保留用于回归。

`experiments.py` 生成性能数据和局部陷阱证据，`audit.py` 重新核验全部保存声明。`reproduce.py` 顺序运行测试、实验和审查，避免并行计时相互干扰。
