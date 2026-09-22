# C2 独立评估器交接

本轮只使用根预先锁定的 24 道公开开发题。20 道 new_public（其中 2 道函数逼近）、2 道旧失败回归、2 道 expected_unsupported_control 分别报告。22 道数学任务与 2 道方法边界控制分开；明确不支持不是数学求解成功，其他失败和超时也不能自动称为正确拒绝。

比较 baseline/C1/C2，各题 3 次新进程，总计 216 条。题序与系统顺序按重复轮换，始终串行。10 秒从进程启动前开始，覆盖导入、剩余预算求解、完整 JSON 编解码、独立完整证书重放、原任务精确端点/误差核对以及原始输出刷新。父进程独立杀死整个进程组；阶段事件和吞异常不能延长截止。实际执行统一使用启动评估器的 Python（子进程继承 sys.executable）。

C2 的 progress 回调原样保存为每次运行的 progress.jsonl，行级刷新。无论正常失败还是被硬杀，都保留最后完整内部动作及原因。原始 response、certificate、全部 attempts、逐阶段事件和验证结果保留；硬杀期间尚未返回的内部对象不会虚构为已取得证据。

路线明确区分 native_c2、reused_c1、baseline_fallback。C1 也按其真实 baseline_fallback_calls 标记基线回退；缺少显式回退原因时，仅依据保存的原始尝试记录注明推断来源。正常复用 C1 本身不是错误回退，不要求编造失败原因。C2 真正回退基线时须保留 fallback_reason。数学成功只由原函数/空间/容差绑定、完整重放和精确误差达标决定，不靠 solver 状态字符串。

每阶段分开记录 SELF 与已回收 CHILDREN 的 user/system CPU 及其总和。整体 CPU 直接取外层 os.wait4，不能再加一次 CHILDREN。Darwin 合成校准已验证：外层 CPU 0.281255 秒，worker 自身 0.110281 秒与已回收子进程 0.157896 秒合计 0.268177 秒，差额包含最后写出和退出。峰值记录取内核外层高水位与观测到的已回收子进程高水位之最大值；单位为 bytes。合成子进程峰值 84,475,904 bytes，高于 worker 的 32,112,640 bytes，外层同样覆盖前者。这是可观测的最大单进程高水位，不是进程树同时内存之和。

若外层硬杀时仍有未回收后代，其最终 CPU/RSS 可能无法完整取得。对应行以及失败成本总账明确标为可能下界；不能以这类失败行为作为加速证据。稳定成功配对才比较 CPU、wall、峰值 RSS、输出与证书的规范 UTF-8 JSON bytes。新题、回归题、控制题与不同精度目标保持分层；实际误差/容差比例原样给出。字节体积不是 token，实际模型与缓存 token 未测量；测量中的数值程序不调用模型 API。

冻结通过上一轮清单的精确路径核对并原文保全 604 文件，不运行任何旧目录扫描器或旧初始化。新目录所有 Python（含分析程序）、cases.json 与 protocol.json 一并加入源快照，运行结束再次严格核验。冻结后新增或修改任何本目录 Python 都会使比较失效；普通 Markdown 与结果可另写。

由根统一执行：

```sh
/PYTHON/bin/python3 -B outputs/c_demand_extension_20260922/evaluate.py freeze
/PYTHON/bin/python3 -B outputs/c_demand_extension_20260922/evaluate.py run
```

每完成一条向 stdout 刷新简短进度，并轻量追加 completed_index.jsonl。checkpoint.json 只记已完成数量及最新条目；完整证据始终在各运行目录的 result.json/output.json 等文件中，避免反复改写巨量响应。全部完成后一次写 report.json，顶层为 rows、summary、freeze、identity_before/after、complete_fixed_denominator、valid_frozen_comparison。

summary 按系统给每题 3/3 稳定成功、精确误差、实际路线、回退原因、失败阶段与资源原值；by_classification/by_dimension/by_kind 分层。数学任务主分母在 mathematical_tasks_excluding_method_boundary_controls，方法控制在 expected_unsupported_control_outcomes。pairwise_comparisons 仅采用各系统与 baseline/C1 都三次成功的同题交集，不能跨不同交集直接排名。

最终 62 项纯合成自检通过，涵盖吞 BaseException、硬杀时阶段/进度保留、完整验证超时、任务绑定、子进程 CPU/RSS 校准、进程组终止、源文件/源快照篡改拒绝、固定分母和控制题分类。未运行 C2 数学调参或全量性能；正式运行由根代理在其数学集成完成后启动。
