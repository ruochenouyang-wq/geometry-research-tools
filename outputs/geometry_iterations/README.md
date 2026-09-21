# 八个新增版本的几何研究工具：V3–V10

这是可运行的研究原型，目的是把“近似算出来的对象”变成带有明确条件和误差界的证据。已完成八个实质迭代，之前的 V1、V2 不计入这八个版本。只需 Python 3.9+ 标准库，无需联网、付费软件或 GPT API 密钥。

先看 [实测结果与选择建议](REPORT.md)，再按 [版本台账](ITERATIONS.md) 查看每项数学依据。问题来源见 [研究背景](RESEARCH_CONTEXT.md)。所有命令在本目录执行。

## 可以直接做什么

```sh
# 一般谱编号；这里 k=2 表示从最低值开始的第二个特征值
python3 run.py spectrum --q "0,1/1000000,-100" --index 2 --tolerance 1e-10 --out results/example_spectrum.json

# 修正试探函数；目标是由函数导出的特征值区间宽度
python3 run.py function --q "0,1" --modes 4 --max-modes 24 --tolerance 1e-40 --out results/example_function.json

# 首两个特征函数张成的空间，输出投影误差的平方界
python3 run.py cluster --q "0,1/1000000,-100" --size 2 --modes 20 --out results/example_cluster.json

# 证明在 a∈[-1,1] 上，q_a(t)=a*t 的最低谱值统一大于等于 -0.3
python3 run.py family --q0 "0" --direction "0,1" --interval="-1,1" --threshold=-3/10 --out results/example_family.json

# 从原始输入重新检查，不信任保存的搜索过程
python3 run.py verify results/example_spectrum.json

# 所有测试、八版实验、跨代对照、全部保存结果的复核
python3 reproduce.py
```

`spectrum --version 3/4/5/6/10` 可选各代谱算法。V3 用 `--modes` 和 `--bits` 固定空间及二分步数；其余版本用 `--max-modes` 限制自适应搜索。`refine` 对应固定空间的 V7；`range` 对应 V5 多项式全域范围；V8、V9 分别由 `cluster`、`family` 调用。每个命令均有 `--help`。

保留的 V2 也可独立调用：`python3 previous_run.py bound --q "0,1" --out results/example_v2.json`。旧版本在部分简单问题上仍然合适，具体见报告中的同目标对照。

## 结果是什么意思

- `target_met`：有理区间的**宽度**达到输入容差。取中点时，特征值绝对误差不超过该宽度的一半。
- `target_not_met` / `unresolved`：预算内没有完成请求，保存的较宽区间仍可能有效。
- 参数族的 `proved` / `disproved`：分别有完整区间覆盖证书或一个严格反例参数。
- 谱簇的 `projector_hilbert_schmidt_squared_upper`：投影差 Hilbert–Schmidt 范数的平方上界；`nontrivial_bound=false` 表示只得到平凡界。

计算成功或得到严格反例返回码为 0；预算内未完成返回 2；输入或证明失败返回 1。`verify` 返回 0 只表示已保存的数学声明与证书一致，并不把预算失败变为成功。实际数学端点是 JSON 中的有理数 `lower`、`upper`；12 位展示小数可能明显宽于高精度证书。

## 当前支持的数学对象

单位二球面上 `H_q=-Δ+q(cos θ)` 的**轴对称函数子空间**。q 是次数≤6 的有理系数多项式，系数绝对值≤1000、分母≤10⁹。系数从常数项起排列，使用整数或分数字符串；容差可使用 `1e-40`。负号开头的命令行值建议使用 `--q="-1,2"` 的形式。

算子采用光滑球面轴对称函数能量型闭包，不能另换区间端点条件。谱编号从 1 开始，q=0 时是 0、2、6、12……，包含零模态。[模型定义与初始证明](BASE_MODEL.md) 给出完整约定；其中 32 模式上限属于 V1，新的 Schur 路线支持至 64，函数和谱簇路线仍为 32。

这是无限维算子的解析尾部控制与有理证书。独立检查使用另一种稠密惯性消元；矩阵装配和模型证明仍有共享依赖，没有导出到 Lean 等形式化内核。当前不接受任意曲面，不覆盖几何离散、数值积分、网格误差或任意中间谱段。

## 文件与复现

`v03_tail.py` 至 `v10_adaptive.py` 是八版实现，`V03.md` 至 `V10.md` 是对应推导与限制。`bench_vXX.py` 和 `results/vXX/` 保存各版对照；`bench_final.py` 是 V2 与 V10 的同目标比较。`TEST_RESULTS.txt`、`VALIDATION.json` 分别记录测试和证书复核。`reproduce.py` 会更新这些记录及实验结果，计时因机器与系统负载而异。

GPT 可以用这些接口提出 q、试探函数、参数范围和预算，再读取证书与失败原因决定下一步。这里没有连接在线模型，也没有用模型的自我评价作为证明接受条件。
