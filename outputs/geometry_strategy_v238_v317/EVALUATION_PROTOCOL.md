# 独立评测协议（不占算法迭代编号）

协议 v1 在最终留出题揭示之前固定。机器可读规则为 `evaluation/protocol.json`，由 `evaluation_protocol.public_rules()` 生成。依据《研究指导报告》1.1 的第 12–14 节，七个维度分别报告，没有总分。

1. **认证精度**：函数逼近只报告原函数的概率 L² 误差上界；谱任务只报告指定全空间或零均值空间的完整谱区间宽度。两者不可互相替代。成功要求完整证书独立重放通过、原函数/结论类型/空间绑定正确、非负误差界达到任务容差，且没有超出统一在线预算。常数势另以全空间 `c`、零均值 `2+c` 检查区间包含性。
2. **可靠性**：重放不读取 solver 的 `verified`、`status`、`target_met` 作为证据。针对每系统及结论类型各取一份有效证书，测试错原函数、错空间、负界/反向区间、只有成功标志的伪证书。任何误接受阻止可靠性提升声明。重放与搜索职责独立，但共享数学后端公式和实现；不是完全独立数学证明，也不是形式化证明。
3. **自动选择**：同题比较旧默认、固定递增、验证集选择最佳固定配置，以及新适配器。固定配置依次为 small/medium/large，具体参数事先写入协议。验证集以达标数、有效证书数、小配置优先的字典序只选一个全局配置，不作逐题事后最优。全部 3×5 条选择实验及失败都保留，完整选择成本计入该基线。程序策略是比较对象，未调用外部模型。
4. **新参数迁移**：开发 6 题、验证 5 题；最终留出 6 题在冻结后生成。固定 `-|z|^(-1/4)` 的函数逼近属于回归而非留出。留出只覆盖已知常数/step/abs_power 家族的新参数，两个 abs_power 指数不同且未在公开题使用，并满足旧后端的稳定性条件。旋转与偏移不称作新数学族；本协议不声称跨结构泛化。
5. **运行资源**：每系统每题在线求解墙钟上限 20 秒；递增策略共享这 20 秒。每份证书重放上限 20 秒，重放成本另计并进入总成本。记录墙钟和本进程 CPU、所有尝试、异常、超时及验证选择成本。主 Unix 线程以信号中止超时；其他调用环境须在报告中标明没有强制超时。新服务内部完整返回数据全部保留。多 agent 并行开发会干扰墙钟；当前单进程交错执行不提供严格冷启动或缓存加速结论。内存峰值和未取得的研发/模型成本记未知。
6. **可见 token**：本评测台没有接入文本编码计数器，token 数为未知，不以字符或字节冒充。交互组件另有独立的离线文本编码实验，不能直接套用到本评测返回值。完整请求/返回的字符数与 UTF-8 字节数另存。实际模型输入、输出、推理与各 agent 用量均未知。本评测测完整证据适配器，不代表精简交互视图的文本成本。
7. **可用性**：记录接口异常、缺失证书、无效证书、超时；所有原始响应与证据保存。没有每题人工调参。可重放与可追查属于已测可用性；真实用户理解或模型使用能力尚未测量。

## 验收与比较目标

- 正确性优先；出现错误证书被接受时不宣称改进。
- 留出比最强基线至少多完成 1 题，且可靠性不退步，可称为本小题组的认证成功率收益。
- 若成功数相同，包含验证选择准备的总 CPU 至少减少 20%，才记为初步资源收益；并发干扰、样本量和单次测量限制必须随结果说明。后续独占环境复测才能提出稳定速度结论。
- 7 维各自保留负面、未知与失败；不把宽度比写成速度比，不把文本长度写成真实模型成本，不声称已完成跨模型或半小时/八小时挑战。
- 三种基线全部展示，尤其不能只与较弱的默认配置比较。旧默认保留旧工具的默认计算参数，仅绑定当前原函数、空间和容差；旧 `approximate` 入口没有 `tolerance`，容差只在独立验收时应用。

## 留出规则与隔离

`evaluation/private_seed.json` 保存 32 字节随机种子，只允许评测侧查看。公开 `seed_commitment.json` **只包含 SHA-256 hash**。种子不进入报告、不发给开发者或根 agent。公开生成规则、候选参数和验收目标已预先固定；具体实例及答案不在冻结前构造或披露。共享目录是纪律隔离而非权限隔离。

冻结流程为完成开发、保存验证选择，再显式调用 `write_freeze_manifest()`。manifest 哈希覆盖新目录全部 Python 模块、协议和基线选择。`phase='holdout'` 要求 manifest 存在、状态 frozen 且所有文件哈希吻合，然后才生成 `holdout_instances.json`。题名为无答案的 H01…H06。已揭示后的源代码变化不允许重新冻结原留出；若用其修复，必须另立协议与新独立留出。

## 调用

```python
from evaluation_protocol import run_evaluation
from research_service import ResearchService
from verification import verify_any

report = run_evaluation(ResearchService(cache=False).evaluate_task,
                        phase='development', verifier=verify_any)
```

`service(task)` 接收 `{kind, function, tolerance, mean_zero?, budget?}` 并返回完整证书或含 `certificate` 的完整字典。当前固定题不额外传内部预算关键字，以便各方法接受相同数学任务；统一预算由评测器执行。新服务可直接用 `ResearchService.evaluate_task`。

省略 verifier 时，仅接受冻结旧后端格式。`verifier` 要接受 `expected_function`、`expected_mean_zero`、`expected_tolerance` 参数并返回严格布尔值。`run_evaluation(None)` 仅运行三种旧基线，便于集成完成前建立开发参照。

```sh
python3 -B outputs/geometry_strategy_v238_v317/evaluation_protocol.py
python3 -B outputs/geometry_strategy_v238_v317/evaluation_protocol.py --new-service
```

报告保留全部 solver 响应与失败分支，评测文件都在 `evaluation/` 下。不得把 `private_seed.json` 加入面向开发侧的批量阅读或发布材料。
