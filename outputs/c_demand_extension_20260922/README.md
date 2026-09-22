# C2：按证明需求分配计算量

这是旧 C 的独立扩展。完整结果与失败分析见 `REPORT.md`；本目录保留原始运行、阶段轨迹和冻结快照。旧算法与旧验证器的文件没有修改。

## 能做什么

输入单位球面上的轴对称势函数及目标误差，输出最低谱值的可靠上下界。最低谱值可理解为给定能量表达式能够取得的最小能量；“平均零”版本额外要求试函数在球面上的平均值为零。支持范围仍然受现有证书方法限制，不能据此宣称适用于任意几何问题。

新增加两条搜索策略：阶跃势逐级提高试函数的多项式次数；较强负幂势先构造可验证的粗谱间隔，再把实际残差与目标误差相比较。其余已支持的负幂题复用 C1，函数逼近题沿用基线。这些来源在每次回执里分别标记。

证书是包含确切有理数、证明所需界限和重放信息的数据。只有原来的验证器按照原题的函数、空间、精度重新检查通过，才算成功。候选函数、较宽的有效区间、超时和不支持都不算解题成功。

## 使用

需要保留本目录及原项目的 `outputs` 依赖目录；仅复制 `controller.py` 无法单独运行。Python 代码使用标准库，本轮统一测试使用 Python 3.12.14。

在项目根目录运行：

```sh
python3 -B outputs/c_demand_extension_20260922/controller.py
```

每行输入一份 JSON，例如：

```json
{"kind":"spectrum","function":{"kind":"axis_profile","profile":"step","axis":0,"amplitude":"-1/3","offset":"0"},"mean_zero":false,"tolerance":"1/1000000","budget":{"wall_seconds":10}}
```

程序输出一行 JSON。`target_met` 表示验证后达到原精度；`certificate_valid` 表示证书通过验证，二者不同。`route` 取 `native_c2`、`reused_c1` 或 `baseline_fallback`；`fallback_reason` 和 `attempts` 解释搜索过程。`unverified_candidates` 仅用于恢复中断工作，不能当作数学证书。

程序内部预算为协作检查：一个正在进行的代数运算可能越过检查点。严格的总时间限制由独立评测器从外部终止整个进程组来实现。

## 重现比较

`cases.json` 固定 24 道公开测试题；`protocol.json` 固定三次冷启动重复和每次十秒上限。“冷启动”表示每次求解从新的 Python 进程开始。两个超出现有方法支持范围的控制题单独统计。轴坐标变化只是球面旋转检查，不是推广到新空间维度。

```sh
python3 -B outputs/c_demand_extension_20260922/evaluate.py check
python3 -B outputs/c_demand_extension_20260922/evaluate.py run
```

重现要求当前文件与已冻结快照一致；不要为了通过检查重写旧冻结记录。完整时间包含加载、求解、输出和评测器独立重放证书。并行开发期间的功能检查不用于提速比较。

数值评测程序没有调用模型 API，开发过程中各代理的 token 未单独计量。输出字节数只表示数据量；它不是 token 数，也不能证明某个较小模型已经胜过较大模型。
