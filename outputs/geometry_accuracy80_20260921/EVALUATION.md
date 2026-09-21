# 困难工具任务的独立评测协议

本协议衡量已有数学任务族中新参数的可验证工具能力，不能把通过有限子任务称为解决其关联的开放问题，也不能称为对新数学结构的泛化。80% 的最终口径仍服从用户选择；当前默认提案是 20 道独立新困难题至少 16 道稳定成功，并完整报告旧 H05–H14 十道困难回归。旧 H05 是成功对照，仍保留。旧结果和源码只读。

当前状态：主任务已于 2026-09-21 批准下述固定参数池、20 题布局和抽样前纯身份去重；困难强度不得因结果不佳而改动。公开回归和评估器边界测试已准备。没有访问任何种子，没有生成正式留出实例，没有冻结验收代码。`evaluation/evaluate.py` 明确拒绝种子路径，命令行也不提供正式留出运行入口。若用户更改验收口径，必须明示修订，并在承诺/冻结前完成。

## 固定的 20 题分层与参数

| 任务 | 数量 | 分层 |
|---|---:|---|
| 负幂势谱 | 12 | 3 个精度 × full/meanzero × 中心化 L² 范数 η<1/η≥1 |
| 一般指数 L² 函数逼近 | 6 | 每个精度 2 题，指数不同于固定 −1/4 |
| 负阶跃势谱 | 2 | full、meanzero 各一题，均为 1e−8 |

三个精度固定为 1e−6、1e−8、1e−10。所有函数轴为 z，offset=0；不加入常数精确题或低精度易题冲高比例。任务包括当前实现不支持或已知失败的结构，不按支持性、成功率或耗时删题。

负幂函数为 `q(z)=a|z|^alpha`，a<0，−1/2<alpha<0。球面概率测度下，`eta_squared = a² alpha² / ((1+2alpha)(1+alpha)²)`。强度分层只用这条精确公式，不能据求解结果选题。

| 参数池 | 指数 alpha | 幅度 a |
|---|---|---|
| 弱谱 | −1/7, −3/13, −4/15, −5/17, −5/14, −7/18 | −k/23，k=6,…,13 |
| 强谱 | −5/12, −7/16, −9/20, −11/24, −12/25, −13/27 | −k/23，k=14,…,22 |
| 一般指数逼近 | 上述全部 12 个指数 | −k/23，k=6,…,22 |
| 阶跃谱 `a 1_{z>0}` | 无 | −k/23，k=6,…,22 |

弱谱池全部满足 η<1，强谱池全部满足 η≥1。有限参数池和强度范围只代表这次提案覆盖的任务，不能外推到所有势函数、几何或开放问题。高强度、接近 −1/2 的参数保留，即使最终全部失败。

已获批准的采样规则：先冻结全部公开已见实例目录，以**规范化的原函数 + 结论种类 + 谱空间**定义身份，忽略 tolerance。在求解前从笛卡尔参数池机械排除公开身份；从余下池用 `SHA256(seed || 8-byte big-endian counter)` 和无偏拒绝采样抽取索引，不放回抽样。拒绝的仅是导致模偏差的随机整数，不是数学任务。每个身份一旦选中，要从相同种类、空间的所有精度池移除。排除目录、原因及哈希全部记录。只有确定得到 20 个新身份才允许正式验收；禁止求解后替换失败题。

固定抽样顺序：按 1e−6、1e−8、1e−10 依次生成 full 弱/强、meanzero 弱/强谱题；随后按同精度顺序各生成 2 道逼近题；最后生成 full、meanzero 阶跃谱。此处只定义规则，不生成任何实例或种子。承诺、冻结、授权揭示步骤尚待主任务批准执行，规则不得默改。

## 独立评估和成本

当前方法入口是 `solver.solve(task) -> dict`，返回证书本体或含 `certificate` 的完整响应。独立重放入口是 `solver.verify_certificate(cert, *, expected_function, expected_mean_zero, expected_tolerance, expected_kind) -> bool`。验证器必须重建数学证据，不能调用 solve 或依赖 solve 的成功标志。旧比较对象是 `geometry_strategy_v317_runtime/release_service.py`，用冻结旧注册表 `verification.verify_any` 重放。不同方法始终装载于不同冷进程，防止同名模块污染。

评估器自行规范化原函数，再绑定结论类型、原函数、空间和任何嵌入的目标误差。谱成功要求 `upper−lower <= tolerance`；逼近成功要求 L² 范数上界 `error_upper <= tolerance`，不是平方误差小于 tolerance。存在 `error_squared` 时还检查平方误差非负且被 `error_upper²` 包含。只接受验证器明确返回布尔 `True`。

每个方法、每道题、每次重复均启动独立冷进程。初始化有共同 30 秒上限，单列 wall/CPU；随后共同 10 秒期限覆盖求解、首次使用时的延迟导入、**完整响应的 JSON 编码和解码**、独立重放，以及包含完整响应的原始记录序列化和落盘。父进程接收工作开始的单调时钟时间戳后，使用同一绝对截止时刻。工作线程必须是 Unix 主线程，SIGALRM 抛出 `BaseException` 子类；即使求解器捕获该异常，父进程仍终止整个工作进程组。没有合作式超时降级。

完整冷成本由父进程测量，覆盖解释器启动、导入/初始化、在线工作、完整输出、异常与退出。CPU 使用操作系统对子进程的实际统计。逐阶段计时、完整冷 wall/CPU、每次重复逐题成本及总和全部保留；进程被强杀后无法得到的在线阶段 CPU 标为未知，绝不填成实际零。该任务总冷成本仍在分母和总成本中。当前实现串行运行；共享主机背景活动仍可能影响时间，因此不作统计显著性承诺。

所有方法题目、期限与重复次数相同。方法顺序按重复轮转，题序也轮转。默认正式提案为 5 次冷重复，稳定成功必须 5/5；3 次可用于稳定开发回归；单次仅作开发冒烟，另列观察成功题数，不能宣称稳定通过或正式达到 80%。重复不是新的数学实例。报告按任务族、空间、强度、精度分层，并保留原始失败、每次 wall/CPU、配对增益与丢失题目。本轮没有附加旧任务的 20% 提速门槛。

不支持、异常、超时、不可序列化、原题不匹配、证书重放失败、证明正确但精度不够、初始化失败，都计入原来的题目分母。初始化失败或工作进程缺失不会跳过题目。开发报告明确 `formal_acceptance=false`。运行前后代码身份变化会标记为不能用于正式比较，原始数据仍保存。

记录实际模型身份、输入/输出/推理 token 与费用的字段目前为未知，不用字符数或字节数替代 token。数学计算成本与模型资源不得混报。

## 开发数据和冻结边界

开发只使用旧公开 H05–H14 和旧 13 道验证题，或明确标为开发的用户授权实例；不能根据新增参数的成功结果反过来决定正式题池。旧十题的 function、space、kind、tolerance 均保持原样。评估器自动把每次已运行题目加入 `seen_development_cases.json`。任何已揭示实例，只要参与调参就永久成为开发或回归；同函数换精度也不能再次标为独立新实例。

正式验收前须冻结全部求解与验证代码、评估器、参数规则、已见目录、回归验收指针及其文件哈希，并保存种子承诺；只有根任务明确授权后才能读取新种子和揭示新参数。已揭示但失败的正式批次完整保留，继续优化后必须使用后续预先承诺的新批次，不把旧批次重新包装为独立验收。

连续 5 次真实实现迭代没有提高稳定开发成功率时，必须记录路线复盘及继续或更换方法的理由。更窄区间、较小残差可以帮助诊断，但不算成功率提升，也不能靠版本名递增重置无进展次数。

## 可运行命令

使用 Python：`/PYTHON/bin/python3 -B`。下列路径均相对工作区根。

```text
<Python> outputs/geometry_accuracy80_20260921/evaluation/evaluate.py prepare-draft
<Python> outputs/geometry_accuracy80_20260921/evaluation/evaluate.py selftest
<Python> outputs/geometry_accuracy80_20260921/evaluation/evaluate.py run --phase regression --repetitions 1 --label smoke
<Python> outputs/geometry_accuracy80_20260921/evaluation/evaluate.py run --phase regression --repetitions 3 --label all_routes
<Python> outputs/geometry_accuracy80_20260921/evaluation/evaluate.py run --phase regression --repetitions 5 --label final_regression
<Python> outputs/geometry_accuracy80_20260921/evaluation/evaluate.py run --phase development --cases outputs/geometry_strategy_v317_runtime/evaluation/validation_cases.json --repetitions 3 --label old_validation
```

缺省同时跑 `current frozen_runtime`，也可用 `--methods` 限定开发测量对象。单次旧十题、两方法最坏在线约 200 秒；三次约 600 秒；五次约 1000 秒，另加单列初始化、启动和记录开销。20 道正式新题、两方法、五次的在线上限为约 2000 秒；当前草稿没有该运行入口。

`evaluation/test_boundaries.py` 的合成证书只用于评估器接线和超时边界测试，绝不进入数学准确率分母。真实新证书仍须通过各自数学验证器及其反例/篡改测试；合成测试通过不能替代新证明格式的正确性审查。

## 已实现的一次性正式工作流

以下是锁定参数协议之后新增的执行机制，没有改变 `protocol.json` 的任何字节、参数池、布局、预算或成功门槛。协议文件中的 `formal_release_state` 记录的是锁定当时尚无正式入口的历史状态；现在的执行入口为 `evaluation/workflow.py`。真实承诺、冻结、揭示均须由根任务在审查后显式调用，本节写作时未执行。

公开目录还包括 `evaluation/public_test_cases.json`，它静态登记模块测试和独立审查使用或公开定义的函数，并引用源码行号。该清单涵盖一般逼近的四个额外参数及接近可积边界的失败、三组独立逼近审查参数、full/mean-zero 谱与阶跃测试、四个强势压力实例。它不运行求解器，也不生成新题。`evaluate.public_catalog()` 合并该清单、历史公开题、所有开发运行及所有已揭示批次；身份始终忽略精度。仅做局部数学计算的公开函数也保守地登记，清单不会被用来挑选成功任务。

正式工作顺序固定如下：

1. 所有 `.py` 作者停写；旧十题和公开四道强势压力各跑五次、两方法。根任务写 `evaluation/acceptance_regression.json`，格式为 `{"reports":[{"report_path":"evaluation/<旧十题目录>/report.json","sha256":"<原文件SHA256>"},{"report_path":"evaluation/<四强势目录>/report.json","sha256":"<原文件SHA256>"}]}`。路径相对本项目目录；不用参考题成功来删题。冻结要求两份报告各有完整的 case × method × 5 原始行、源码运行前后相同，并与冻结时所有源码字节一致。
2. `commit --batch B001` 先保存已见数学身份快照，再仅用 `os.urandom(32)` 产生随机字节、以独占创建方式写入 0600 私有文件，只公开 SHA-256 承诺；此步骤绝不读取或输出种子。先有公开目录，后有承诺。提交后如出现新的公开身份，旧批次必须保留并放弃，不能偷偷扩目录或换种子；同一身份重复测量、改变容差不产生新身份。
3. `freeze --batch B001` 在读取种子之前冻结所有新旧几何 Python 源码、测试、审查代码、协议、公开目录、承诺、回归指针和两份报告。它额外执行并保存 `runtime.release_evaluation.checked_manifest()`、`verify_old_preservation()`、旧 `audit_campaign.audit_frozen()` 三项只读历史保全结果。这些历史审计明确不读取旧私有种子，覆盖旧数学依赖及既有清单中的非 Python 资料。全部源文件及关键验收资料还归档为 `source_snapshot.zip`，内部 `.blob` 加索引保存原始路径、字节和 SHA-256；后续代码优化不会销毁该批次的可复现源文件，也不会让归档递归进入 `.py` 快照。
4. `authorization-template --batch B001` 只输出默认 `authorized:false` 的模板。根任务审查后手工写该批次的 `authorization.json`，必须明确 `authorized:true`、`authorized_by:"root"`、`action:"reveal_once"`，并绑定 batch、冻结文件、承诺文件及协议的准确哈希。工具没有自动授权功能。
5. `reveal --batch B001` 先验证授权和全部冻结资料，独占创建不可重试的揭示开始标记，才允许一次性读取种子、核验承诺和抽样。计数器从 0 起、32 字节随机字与大端 8 字节计数器拼接；候选池按原题身份的规范 JSON 字符串排序；256 位 SHA-256 值达到 `2^256−(2^256 mod 池大小)` 时只丢弃该随机整数并递增计数器。每个选中身份和每项纯去重排除均保存审计记录，不按成功筛选。即使揭示失败也不能重冻或重复读取同批种子。
6. `run --batch B001` 严格运行 20 题 × 2 方法 × 5 次，共 200 个独立冷 worker；每题仍使用同一个 10 秒求解/完整输出/独立验证协议。需要所有 200 行完整、源码与历史资料复查通过，且 current 至少 16 题全部 5 次通过，才记录本批正式 80% 验收通过。单个 worker 或后置保全异常会保留已有行和原 20 题分母、记录失败，不允许选择性重跑覆盖。

命令均使用前文指定 Python：

```text
<Python> outputs/geometry_accuracy80_20260921/evaluation/workflow.py selftest
<Python> outputs/geometry_accuracy80_20260921/evaluation/workflow.py commit --batch B001
<Python> outputs/geometry_accuracy80_20260921/evaluation/workflow.py freeze --batch B001
<Python> outputs/geometry_accuracy80_20260921/evaluation/workflow.py authorization-template --batch B001
# 根任务显式授权并写入本批 authorization.json 后：
<Python> outputs/geometry_accuracy80_20260921/evaluation/workflow.py reveal --batch B001
<Python> outputs/geometry_accuracy80_20260921/evaluation/workflow.py run --batch B001
```

所有批次的题目和失败永久保留。后续新承诺自动排除所有已揭示身份，包括评测中断的批次；完成失败报告后才能新建批次，或由根任务用 `abandon --batch B001 --reason "具体原因"` 明示放弃且保留全部文件。新批仍使用完全相同的锁定参数池和规则，不能因困难分布未达 80% 而改池。

生命周期测试只在临时目录使用固定合成字节和分母为 101 的替代幅度池，刻意与正式分母 23 的参数池分离；不读取真实种子，不从真实正式参数池抽样，也不调用真实生命周期命令。真实流程可用性与数学求解性能分别评估。
