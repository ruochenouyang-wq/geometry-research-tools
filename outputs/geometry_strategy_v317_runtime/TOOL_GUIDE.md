# 验收修订版：怎样使用

这是已经完成的 v238–v317 的发布修订。保留原数学方法和所有证书格式，改进重复检查和初始化的开销。**证书**是保存原题、上下界和复查数据的完整记录；**重放**是重新检查这些记录，不必重新寻找候选。

调用方式和字段继承[原完整使用指南](../geometry_strategy_v238_v317/TOOL_GUIDE.md)，只将程序入口替换为本目录的 `release_service.py`。默认仍返回数学摘要；`view="full"` 返回完整证据、研究过程和单次请求重放复用统计。

## 一次真实请求

在当前机器运行下列命令。输入要求保持一行；它计算单位球面上势函数 `q(z)=-|z|^(-1/4)` 在平均值为零的函数空间中的最低谱值，要求认证区间宽度不超过 `1e-8`。这里 `z` 是球面高度，**谱区间**是已经通过数学检查的下界和上界组成的区间。

```sh
cd /PROJECT
/PYTHON/bin/python3 -B outputs/geometry_strategy_v317_runtime/release_service.py <<'JSONL'
{"op":"research","target":"spectrum","function":{"kind":"axis_profile","profile":"abs_power","exponent":"-1/4","amplitude":"-1","offset":"0","axis":2},"mean_zero":true,"tolerance":"1/100000000"}
JSONL
```

只做函数逼近时把 `target` 改成 `approximation`，并删除 `mean_zero`。**L² 误差**是函数差的平方在球面上的平均值再开平方；它衡量整体逼近误差，不能代替谱结论。当前自动逼近仍只支持阶跃函数和指数为 `-1/4` 的绝对幂。

先查看 `certificate_valid`、`target_met`、`within_budget`：它们分别表示证据有效、达到本次精度、在预算内返回。有效证据可以尚未达标；正确结果也可能超时。谱请求的 `budget.wall_seconds` 是完整服务调用的软时间预算；**软预算**不能中断已开始的同步计算，程序会在返回时标明超时。独立评测另有能够中断运行的硬期限。

## 保存和复查

摘要将完整证据保存在本目录 `runtime_evidence/`，并返回 `evidence_ref.sha256`。把这个 64 位字符串作为 `fetch` 请求的 `reference` 值即可取回；读取时会再次检查。也可使用 `GEOMETRY_EVIDENCE_ROOT` 环境变量指定保存目录。

```json
{"op":"fetch","reference":"将上次返回的64位证据标识放在这里"}
```

标识用于核对内容身份，不是数学证明。`fetch` 不附带新的原题目标，因此 `target_met` 为 `null`；要针对具体原题复查，使用原指南的 `verify` 请求，明确函数、空间、结论类型和精度。批处理 `batch` 继续逐项处理，一个错误不会取消后续请求。

## 复用为什么不改变数学结果

同一次请求中，只有完整证书内容、全部目标绑定和验证器版本均相同，才能复用已经成功的检查。不同空间、不同函数、不同精度或不同数学字段都重新检查。非标准 JSON 对象按原类型检查，不把元组变成列表后冒充原验证器接受。请求结束后关闭并清空该缓存；通过复制执行上下文也不能延长其有效期。

普通请求先使用现有数学模块，奇异谱计算需要额外的私有模块时才构造它们；这些首次构造的时间仍计入请求。完整证据存储和交互组件也在首次 `handle` 调用时初始化，目录有效性在构造与首次使用时均检查。用于生成候选的模块和验证依赖在服务私有实例中接合，不替换共享数学后端的全局函数。每个新证据的首次检查继续使用冻结的验证方法；发布评测还会在缓存之外用冻结验证器独立重放。

复用统计仅附在完整视图及程序评测接口，默认摘要不增加这组计量字段。`replay_reuse.hits` 是本次省去的完全相同检查次数，不是省下的模型 token，也不是新的数学定理。内部嵌套计时存在重叠，总时间应读取完整服务或评测记录，不能相加。

## 复现实验

进入本目录后，运行 `verify_release.py --tests` 检查全部新测试、冻结 v317 回归和文件保留情况。运行 `release_evaluation.py run-validation` 或 `run-regression` 可以重复公开数据比较。最终留出已经放行后，可运行 `release_evaluation.py run-holdout` 重复冻结实现的计时；重复计时不会产生新的独立数学样本。

程序接口为 `ResearchService(evidence_root=None, cache=True, replay_reuse=True)`，提供 `handle(request)` 和 `evaluate_task(task)`。关闭 `replay_reuse` 是比较缓存效果的**消融实验**，即只关闭这一项机制而保持其他行为。模型真实消耗、模型理解率和跨模型能力实验仍未测量。

## 可见文本优化的实际入口

默认摘要保留成本和执行状态，没有自动套用旧交互实验的请求绑定压缩。该压缩是可单独调用的 `interaction.choose_response`：它比较保留全部语义的两种表示，仅在提供的各编码方案都不增长时选压缩表示；没有编码器时保持独立摘要。接收方必须保有原请求，必要时用 `expand_bound_summary` 恢复并核对。

旧[完整文本实验](../geometry_strategy_v238_v317/iterations/v317/TOKEN_COMPARISON.json)报告的是这个可选组件的固定脚本计数。默认完整服务包含动态计时及其他使用数据，不能直接套用该组件的百分比；真实模型输入、输出和推理用量均尚未测量。
