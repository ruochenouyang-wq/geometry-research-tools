# 独立评估器交接

`benchmark.py` 与 `selftest.py` 只写新目录，不调用旧评估器初始化；全部题目直接读根固定的三份 JSON。阶段预算固定为冷单题整个进程 10 秒、批量准备 30 秒、每查询 10 秒、整会话 200 秒、每张完整证书额外冷重放 10 秒。批量准备、查询和完整验证均计真实费用，失败不减少分母。这里的“维度”是精度、势强度、空间等评估轴，数学几何仍为单位 S2。

最终只由根依次执行：

```sh
/PYTHON/bin/python3 -B outputs/parallel_cross_dimension_20260921/benchmark.py freeze
/PYTHON/bin/python3 -B outputs/parallel_cross_dimension_20260921/benchmark.py run
```

开始前核对原 567 文件、上一轮 572 源码/12 静态资料及其源快照；新快照按路径并集合保全，包括原清单中 5 份 math_research 依赖的原文。新目录全部 Python、cases/families/protocol、environment.json 和 notes/ANALYTIC_CONTROLS.md 一并冻结，运行结束再次严格核验。其余 Markdown 和结果可后写；新建或修改本目录任何 Python 都会使冻结失效。需做分析时宜用一次性解释器代码，或在本目录以外保存新分析程序。

最终 `evaluation/final_*/report.json` 顶层：

- `common_rows`：360 条冷单题原始记录。每条含 case、原 task/hash、repetition、preparation、完整 response、verifier_result、assessment、timings、process、JSON bytes、回退、失败分类和原始目录。
- `common_summary`：按系统给每题 5/5 稳定成功、真实精确误差与误差/容差比例、资源五次原值/中位数/范围；`by_kind` 分谱/函数逼近，`by_dimension` 按评估轴。四系统共同稳定成功交集是主比较集；另给各系统相对基线的稳定成功交集，不能拿不同交集跨算法排总名次。全 18 题的计数仅用于完整性核对，不是混合能力排名。
- `batch_sessions`：40 会话，每个保留 16 查询与 16 个冷重放槽位；无证书则明确未执行重放。`checkpoints` 是 N=1/4/16 最后一题结束时的含准备前缀费用与附加冷重放费用；`complete_resource_totals` 是整个会话真正退出后再加全部冷重放的完整总账，两者不能混称。
- `batch_summary`：按 family→systems 汇总固定分母、每题稳定成功、准备/证书体积、会话 RSS、完整费用及前缀配对。
- `freeze`、`identity_before/after`、`complete_fixed_denominators`、`valid_frozen_comparison`：冻结与分母状态。冻结有效不表示所有数学题都成功。

进程 wall 从启动前到回收/输出管道结束；每个进程直接用 `os.wait4` 得 user/system CPU 和 peak RSS，含被杀进程。Darwin 的 ru_maxrss 原始单位就是 bytes；每阶段的 RSS 只标记累计高水位，不做虚假差分。冷任务的 full verify wall/CPU 可单列，但其峰值仍属于同一进程；批量额外冷重放有独立进程 RSS。序列化体积使用规范 UTF-8 JSON 的字节长度，准备证书和序列化证明库体积另列，不表示 Python 堆大小或 token。

所有返回证书都经过原任务绑定与精确端点/误差再次评分；solver 的 target_met 不作为裁决。unsupported、能力未完成、精度不足、验证失败、超时和进程故障区别保留。同一失败耗时仍计入所有尝试的费用。求解、验证、JSON 与原始输出刷新都在对应截止内；父进程的独立全局截止不会因子进程阶段事件重置。

开发检查：首轮 42 项通过，含一次原基线阶跃轻冒烟；最后一轮 45 项纯合成检查通过，覆盖吞 BaseException、验证超时保留证书、超时 RSS、16 项固定分母及不同共同成功子集。没有运行 A/B/C 数学调参或全量比较。真实模型 API 调用为 0；实际模型/研究代理 token 与缓存 token 均未知，JSON bytes 不能代替它们。
