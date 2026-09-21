# 本地数学研究工具使用指南

本指南对应 `research_service.py` 的集成入口。先选择要解决的任务，再读取它给出的数学证据与完成状态。示例只使用本地程序，不调用 GPT 或其他在线模型。

## 1. 先确定目的

目前主要支持单位球面上的两类任务。**单位球面**指半径为 1 的球的表面。

| 你的目的 | `target` | 完成标准 |
|---|---|---|
| 用较简单的分段函数逼近给定函数 | `approximation` | 原函数的 L² 误差上界不超过目标 |
| 估计给定势函数下的最低谱值 | `spectrum` | 完整原题的谱区间宽度不超过目标 |

**L² 误差**是把两函数之差平方，在球面按面积取平均，再开平方；它衡量整体差距，不保证每一点的误差都小。**势函数**是能量式中给定的函数；工具中的谱任务估计该能量在规定函数空间中的最低值。**谱区间**由已证明的下界和上界组成，其宽度为上界减去下界。

函数逼近与谱估计是两种不同任务。阶跃函数可以被精确表示，L² 误差为零；这不意味着该势函数的谱问题也已经精确求解。

## 2. 输入怎样写

程序使用 **JSONL**：每一行是一个完整的 JSON 对象，每行输入对应一行输出。JSON 是用字段名和值表示数据的文本格式。复制示例时，保持每个请求占一整行；不要把行内字符串拆开。

建议从 `op="research"` 开始，它自动选择现有方法及部分计算参数。

| 字段 | 含义 |
|---|---|
| `op` | 操作名称，自动研究使用 `research` |
| `target` | 必须明确为 `approximation` 或 `spectrum` |
| `function` | 原函数的精确描述，示例见下文 |
| `tolerance` | 允许的误差或区间宽度，默认 `"1/100000000"` |
| `mean_zero` | 仅用于谱任务；默认 `true`，要求参与求解的函数在球面上平均为零 |
| `budget` | 计算预算，两个任务接受的字段不同 |
| `policy` | 可选谱路线策略，入门使用时可以省略 |
| `view` | `summary` 返回摘要，`full` 返回完整过程与证书；默认 `summary` |
| `request_id` | 可选的请求标签，便于对应输入与输出，不是数学证明 |

**有理数**指整数或两个整数之比。函数参数、指数和精度应写成精确字符串，例如 `"-1/4"`、`"3/2"`、`"1/100000000"`；不要用 `0.00000001` 代替精确目标。精度范围为 `1e-30` 到 `1`。计数预算则使用 JSON 整数，例如 `24`。布尔值为小写 `true` 或 `false`。

### 函数的描述

`kind="axis_profile"` 表示函数只随一个坐标轴上的高度变化。`axis` 可以是 `0`、`1`、`2`，分别对应 x、y、z 轴，默认是 z 轴。高度在 `[-1,1]` 内。

- `profile="step"`：高度不大于零时取 `offset`，大于零时取 `offset+amplitude`。
- `profile="abs_power"`：取 `offset + amplitude * |高度|^exponent`；指数必须明确给出。

`amplitude` 表示幅度，默认 `"1"`；`offset` 表示整体加上的常数，默认 `"0"`。例如以下原函数为 `q(z)=-|z|^(-1/4)`：

```json
{"kind":"axis_profile","axis":2,"profile":"abs_power","exponent":"-1/4","amplitude":"-1","offset":"0"}
```

自动函数逼近目前支持阶跃及指数恰为 `-1/4` 的绝对幂。谱路线支持的输入更广，但还须通过相应数学条件检查，不能据此推断任意指数、任意幅度都可解。

### 两类预算

自动函数逼近的常用 `budget` 字段如下。**次数**指每一段所用多项式的最高幂次；**层数**指向奇异区域逐步缩小的分区数量。

| 字段 | 默认值 | 含义 |
|---|---:|---|
| `max_degree` | 24 | 次数上限，允许 0–32 |
| `max_levels` | 64 | 层数上限，允许 1–64 |
| `max_shape_evaluations` | 264 | 计划阶段新计算标准形状的次数上限，可以为 0 |
| `max_certificate_slots` | 2113 | 展开后系数槽数量上限，至少 2 |
| `sqrt_bits` | 40 | 向上取整开平方的起始二进制位数 |
| `max_sqrt_bits` | 256 | 上述位数的上限 |

还可给 `root_ratios`，例如 `["1/2","2/3","3/4"]`，指定尝试的分区比例；或用 `root_ratio` 只指定一个比例，两者不能同时出现。服务自行管理标准形状的复用，JSON 请求中不能传 `shape_cache`。这些计数限制不构成整个请求的硬时间上限。

自动谱研究的 `budget` 字段为：

| 字段 | 默认值 | 含义 |
|---|---:|---|
| `max_attempts` | 8 | 最多启动多少次外层求解调用 |
| `work_units` | 2000000 | 控制器分配工作的估算单位，不等于秒数、实际运算次数或 token |
| `wall_seconds` | 60 | 完整服务调用的软墙钟预算，即现实时间预算 |

**软预算**表示程序在调用之间检查时间，无法中断已经开始的同步计算。一个调用可能超过预算才返回；其有效证书仍会保留，同时标记 `within_budget=false`、`status="budget_exceeded"`。最终状态还会把证书检查、摘要和存储时间纳入完整调用预算；它仍不是操作系统级硬超时，不包含进程启动和输出传输。工具不会把迟到的正确结果算作预算内完成。

可选谱策略包括 `compare_routes`、`stagnation_patience` 和 `min_relative_gain`；它们影响路线比较与无收益停止。先使用默认值即可。未达标或停止搜索不等于证明原题无解。

## 3. 怎样解释结果

**证书**是保存原题、界限及检查所需计算的完整数据。**验证器**是按既定数学规则重新检查这些数据的程序。**重放**指读取证书后再做这些检查；它不必重新寻找候选解。

先看以下字段，再看小数：

| 字段 | 可以得出的结论 |
|---|---|
| `execution_ok` | 请求的执行流程是否正常返回；正常返回不等于数学目标达成 |
| `certificate_valid` | 证书是否通过实际检查；`true` 本身不表示精度足够 |
| `target_met` | 有效证据是否达到这次请求的数学精度；不单独保证在预算内完成 |
| `within_budget` | 完整谱服务调用是否在软墙钟预算内返回；没有该字段不能当作 `true` |
| `status` | 当前综合状态；常见值见下表 |
| `stop_reason` | 如有此字段，说明为什么停止，例如时间、尝试次数或达到目标 |
| `cost` | 本次数学工具的实际计数或调度估计，各字段的单位不同 |

常见状态：

- `certified_met`：摘要中的数学精度达标；完整视图通常对应 `target_met`。
- `certified_open`：已有有效证书，但目标精度仍未达到。
- `budget_exceeded`：超出软时间预算，可能同时保留达到数学精度的有效证书。
- `verification_failed`：证据未通过检查，不能使用它声称原题成立。
- `request_failed`：输入或执行出现错误，应查看 `error` 与 `reason`。
- `open` 或其他停止原因：可能尚无证书，例如计划预算为零；查看 `certificate` 和具体说明。

不能只凭 `ok=true` 判断完成。在一批请求中，最外层的 `ok=true` 仅表示批处理返回了结果；必须逐项检查。

### 摘要与完整视图

`summary` 保留原题范围、精确界、检查状态、成本和证据引用，适合日常阅读。逼近误差在 `bounds.error_upper` 中；谱的精确宽度在 `bounds.exact_width` 中。`decimal_...` 字段是便于阅读的小数展示；最终精度判断以精确分数为准。

摘要中的 `evidence_ref.sha256` 是完整证书内容对应的 64 位十六进制标识。它用于定位和检查内容有没有变化，**标识本身不是数学证明**。摘要会把证书保存到本项目的 `outputs/geometry_strategy_v238_v317/runtime_evidence/`。

`full` 把完整证书放在 `certificate` 中，并保留尝试记录、诊断等细节；输出可能很长。完整视图不自动产生摘要引用。希望稍后按引用取回证书时，先请求摘要，再使用 `fetch`。

相同进程中的相同请求可以复用计算。另起一个命令进程后，内存中的请求缓存与形状缓存会重新开始；磁盘上的证书仍可取回。缓存命中与证书通过验证是不同的事情。相同 verify 请求可能复用同一代码版本下的检查结果，并在 interaction.cache_hit 中注明；fetch 每次都会读取完整证据并重放，补读与失败的成本也计入 usage。

## 4. 完整示例

### 先进入项目目录

以下命令适用于当前机器，Python 的确切路径为 `/PYTHON/bin/python3`。若项目搬到其他目录，只需相应修改目录路径。

```sh
cd /PROJECT
```

### 示例 A：自动选择函数逼近方案

目标是逼近 `-|z|^(-1/4)`，要求球面 L² 误差上界不超过 `1e-8`。不用手工指定最终次数和层数。

```sh
/PYTHON/bin/python3 -B outputs/geometry_strategy_v238_v317/research_service.py <<'JSONL'
{"op":"research","target":"approximation","function":{"kind":"axis_profile","axis":2,"profile":"abs_power","exponent":"-1/4","amplitude":"-1","offset":"0"},"tolerance":"1/100000000","budget":{"max_degree":24,"max_levels":64},"view":"summary","request_id":"approx-example"}
JSONL
```

查看 `certificate_valid`、`target_met` 和 `bounds.error_upper`。这个请求不应添加 `mean_zero`，因为它不是谱任务。

### 示例 B：自动研究同一势函数的最低谱值

这里 `mean_zero=true` 明确固定平均为零的原题。工具选择并检查现有谱方法，目标为完整谱区间宽度不超过 `1e-8`。

```sh
/PYTHON/bin/python3 -B outputs/geometry_strategy_v238_v317/research_service.py <<'JSONL'
{"op":"research","target":"spectrum","function":{"kind":"axis_profile","axis":2,"profile":"abs_power","exponent":"-1/4","amplitude":"-1","offset":"0"},"mean_zero":true,"tolerance":"1/100000000","budget":{"max_attempts":8,"work_units":2000000,"wall_seconds":60},"view":"summary","request_id":"spectrum-example"}
JSONL
```

查看 `bounds.lower`、`bounds.upper`、`bounds.exact_width`，并同时检查 `target_met` 和 `within_budget`。把 `view` 改为 `full` 可以读取每次尝试、停止原因和完整证书；这两个视图对应同一数学任务。

### 示例 C：一次提交两个任务

**批处理**指一次发送多个请求并逐项得到结果；本入口按顺序处理，不等于并行。最多 64 项，不允许在批处理中再嵌套批处理。

下面第一项精确表示一个阶跃函数；第二项求常数势 `3/2` 在零均值空间中的最低谱值。第二项有解析答案 `2+3/2=7/2`，适合检查输入和返回格式。

```sh
/PYTHON/bin/python3 -B outputs/geometry_strategy_v238_v317/research_service.py <<'JSONL'
{"op":"batch","requests":[{"op":"research","target":"approximation","function":{"kind":"axis_profile","profile":"step","amplitude":"3","offset":"-2"},"tolerance":"1/100000000","view":"summary"},{"op":"research","target":"spectrum","function":{"kind":"axis_profile","profile":"step","amplitude":"0","offset":"3/2"},"mean_zero":true,"tolerance":"1/100000000","view":"summary"}]}
JSONL
```

返回的 `responses` 数组与输入顺序一致。某一项失败不会让其他项自动变成失败或成功。

### 示例 D：生成摘要，再自动取回完整证书

`fetch` 的 `reference` 必须是 `evidence_ref.sha256` 中的 **64 位字符串**，不是整个 `evidence_ref` 对象，也不是文件路径。以下完整例子会先保存一行摘要，然后自动提取该字符串并生成正确的 `fetch` 请求，无需手工复制长标识。

```sh
/PYTHON/bin/python3 -B outputs/geometry_strategy_v238_v317/research_service.py <<'JSONL' > outputs/geometry_strategy_v238_v317/example_summary.jsonl
{"op":"research","target":"approximation","function":{"kind":"axis_profile","profile":"step","amplitude":"3","offset":"-2"},"tolerance":"1/100000000","view":"summary"}
JSONL

/PYTHON/bin/python3 -B -c 'import json,sys; result=json.load(sys.stdin); print(json.dumps({"op":"fetch","reference":result["evidence_ref"]["sha256"]}))' < outputs/geometry_strategy_v238_v317/example_summary.jsonl | /PYTHON/bin/python3 -B outputs/geometry_strategy_v238_v317/research_service.py
```

`fetch` 检查保存内容并重放证书，返回完整 `certificate`。它没有收到新的原题及精度要求，所以返回 `binding="unbound_evidence_replay"`、`target_met=null`。`null` 在这里表示“未对新的目标作判断”，不能当成完成新的问题。

### 运行测试

下面命令运行本目录中的测试。**测试**用于检查实现行为和已覆盖案例；通过测试不等于证明工具对所有数学问题正确。

```sh
/PYTHON/bin/python3 -B -m unittest discover -s outputs/geometry_strategy_v238_v317 -p 'test_*.py' -v
```

只检查自动函数逼近时：

```sh
/PYTHON/bin/python3 -B -m unittest discover -s outputs/geometry_strategy_v238_v317 -p test_adaptive_approx.py -v
```

## 5. 范围、成本与下一步

当前几何固定为单位球面；函数须符合已支持的沿轴变化形式。奇异谱方法还需要可积性、定义域和用于证明的谱隙等条件。**谱隙**指相关特征值之间的间隔；候选很准确，但间隔的可靠估计不足时，原题仍可能未达标。工具不会用局部候选的好结果替代遗漏方向的检查。

自动逼近中“有限配置内最少系数槽”的结论只针对当前列出的比例、次数和层数范围。它不是所有函数、所有网格中的全局最优，也不是最快运行时间的证明。

程序使用精确有理数计算和自有验证器，尚未接入 Lean 等形式化证明助手。**形式化证明**是把命题与推导编码为严格逻辑对象，再由证明内核检查；当前的数值证书与它不同。保留整个项目目录及所依赖的历史代码，不能只复制本文件夹中的某一个求解器。

`interaction.usage`、顶层 `usage` 和 `cost` 记录程序运行信息。阶段计时可能相互包含，不能简单相加当作墙钟时间；失败尝试和读取详情的成本也应保留。`work_units` 只作调度估算。

**token** 是模型处理文本时采用的编码单位。本项目另有本地文本编码计数功能，它只估计明确文本在指定编码下的长度，不能测得模型内部推理消耗、服务额外开销或账单。服务本身不调用模型；`actual_api_usage=null` 表示没有真实 API 用量数据，不表示整个研发项目的模型消耗为零。摘要更短也不能单独证明 GPT 解题更快。

进一步阅读：[自动逼近](APPROX_NOTES.md)、[分数幂谱方法](SINGULAR_NOTES.md)、[研究控制器](CONTROLLER_NOTES.md)、[证据保存](EVIDENCE_NOTES.md)、[结论复用](TRANSPORT_NOTES.md)、[请求与状态](CONTRACT_NOTES.md)。
