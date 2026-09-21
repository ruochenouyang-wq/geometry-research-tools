# v298–307：可复用证据与研究关系

本组件负责保存、核验和连接已有数学证据。它不搜索新的谱值，也不把研究图中的关系名称当作证明。默认数学检查复用冻结后端；源文件没有修改。

## 十步与真实检查

每一步在实现和测试后调用 `campaign.record`，保存当时的 Python 源文件。`iterations/v298` 至 `v307` 是十份逐步变化的快照，最终源码的 13 项测试含多组功能与攻击断言。

| 版本 | 实际能力 | 主要对照或攻击 |
|---|---|---|
| 298 | 完整规范内容的 SHA-256 引用与去重 | 字段顺序不同复用一个对象；嵌套下界不同得到不同引用 |
| 299 | 严格 JSON 与小型完整内容摘要 | 非字符串键、元组和非有限浮点数拒绝；末尾隐藏数据改变摘要哈希 |
| 300 | 原子发布、读后检查与哈希读取 | 模拟发布失败无半文件；修改已存界后读取和重复保存都拒绝 |
| 301 | 引用、文件和格式边界 | 路径穿越、短哈希、符号链接、大文件和过深输入拒绝；未知格式明确未认证 |
| 302 | 默认数学重放及显式检查回调 | 关闭谱搜索/试探搜索仍可复核分段证据；篡改误差与仅有 `verified:true` 的数据不能过关 |
| 303 | 按完整内容与检查器版本复用核验 | 重复请求只调用一次检查器；版本变化重新检查；缓存命中前仍查文件哈希 |
| 304 | 证书与完整原题绑定 | 函数、空间、几何、测度、谱编号、目标精度任一变化均拒绝错绑 |
| 305 | 六种有类型的研究关系 | 等价、蕴含、上下界、近似和启发均可记录；仅有名称和假设不获得认证 |
| 306 | 循环依赖检测与拓扑顺序 | 自依赖、反向闭环、多前提闭环均拒绝；失败时保留原图 |
| 307 | 有检查规则的关系及只读重放 | 合法弱化界通过；过强结论、错范围、未证条件、错近似函数均拒绝；重放时禁用全部写入仍通过 |

开发中有两次测试先失败再修正，已写入对应版本记录：v302 最初误将冻结验证器共用的精确重建函数也禁用；v305 测试方法插入位置造成名称错误。没有删除这些失败说明，也没有把修测试额外算作版本。

## 接口

```python
from evidence_store import EvidenceStore, ResearchGraph, problem_claim

store = EvidenceStore("path/to/local/evidence")
ref = store.put(certificate)
value = store.get(ref)
receipt = store.summary(ref)
check = store.verify(ref, expected_problem=problem_claim(certificate))
assert check["verified"]
```

`put` 接受有限、严格 JSON：字典键必须是字符串，拒绝会静默变形的 Python 对象。它允许保存候选、失败、未知格式和图节点；保存成功只表示完整数据能够按引用取回。

`get` 接受完整的 64 位小写十六进制内容引用。它对文件大小、普通文件类型、完整字节哈希与规范 JSON 逐次核对。默认文件上限为 32 MiB，可由 `max_bytes` 设置为 64 字节至 256 MiB。没有任意路径读取接口。原子发布使用同一目录的临时文件和独占硬链接，不能覆盖已有不同内容；不声称这是恶意进程并发替换祖先目录的沙箱，也不声称断电持久性已验证。

`summary` 哈希覆盖整个对象；其显示字段只是预览，始终带 `summary_is_certificate:false`。数字类型不混同：`True`、`1`、`1.0` 的内容引用不同。JSON 层允许有限浮点数用于诊断；数学验证器仍决定证据数值类型是否合法。

`verify` 返回字典，其中 `verified` 只有实际数学检查或合法缓存命中才为真。默认识别三种冻结格式：完整直接矩谱、分段 L² 逼近、奇异分数幂 Temple 谱。未知格式返回 `unknown_format`，回调错误返回 `verifier_error`。默认验证共享既有精确装配代码，独立于搜索历史，并不是另一个独立实现或形式化内核。

自定义检查器用 `EvidenceStore(root, verifier=callback, verifier_version="my-rule-v1")` 注入。只有严格布尔值 `True` 表示通过。自定义回调本身是信任依赖；无版本回调不缓存。默认版本在构造时记录冻结 Python 源码指纹与 Python 版本。缓存仅驻留当前进程，绑定完整引用、版本和回调身份，每次仍重读并核对完整文件。修改检查器逻辑或其环境应更换版本或新建存储实例。

`problem_claim(certificate)` 仅提取已支持格式的原题身份，不检查其真伪。谱身份包括函数、单位球面、空间、零均值要求、谱编号、测度和声明容差；近似身份包括原函数、几何、范围和误差范数。传入 `expected_problem` 必须完整且精确匹配，检查在缓存之前进行。有效的 `certified_open` 证据可以通过重放，这不表示精度目标达成。

## 研究图中的证明与想法

```python
problem = problem_claim(certificate)
graph = ResearchGraph(store)
premise = graph.add_node("certificate", problem, evidence=ref)
claim = graph.add_node("claim", {"problem": problem, "value": certificate["upper"]})
edge = graph.add_relation(
    "upper", [premise], claim,
    scope=problem, conditions={"mean_zero": True},
    rule="certified_bound_v1",
)
assert graph.relation_status(edge)["verified"]
manifest = graph.export()
restored = ResearchGraph.from_export(store, manifest)
replay = restored.replay()
```

节点类型为 `problem`、`claim`、`candidate`、`certificate`。关系类型为 `equivalent`、`implies`、`upper`、`lower`、`approximation`、`heuristic`。这些是研究语义，不能由名称直接授予数学权威。图边表示有方向的论证依赖；等价边不会自动产生反向论证依赖。节点最多 512 个，关系最多 2048 条，每条最多 32 个不同前提。

内置两条有限的数学关系规则：

1. `certified_bound_v1`：检查原证书和原题；允许把已证上界调高或下界调低。目标语句必须恰为 `{"problem": ..., "value": ...}`。它不会通过关系图收紧未经证明的界。
2. `certified_approximation_v1`：检查分段 L² 证据，并绑定 `digest(certificate["cells"])` 标识的具体近似。目标语句为 `{"problem": ..., "error_upper": ..., "approximant_digest": ...}`，只能保留或放宽已证误差。零近似误差不成为谱精度结论。

两条规则都要求 `scope` 完整等于目标问题身份，前提证书绑定同一身份。`conditions` 只能引用该已核验身份中已成立的精确字段；不存在或不一致的条件不能作为免费假设。当前是保守的完全匹配，不求解一般范围包含或等价变换。

任意 `equivalent`、`implies` 没有默认检查规则。可以用 `ResearchGraph(store, relation_rules={name: {"version": ..., "kinds": [...], "check": callback}})` 注入有版本的可信规则。检查器接收 `(relation, sources, target, context)`，其中上下文包含已验证前提证书。回调仍须证明它声称的转换；框架不会使任意回调自动可靠。条件、范围或前提失败时不会调用它；`heuristic` 即使被指定检查器也始终不认证。关系检查不持久缓存，重放重新运行规则。

`from_export` 是只读的：重新读取全部节点和关系，检查哈希、字段模式、端点、证据引用和依赖环。`replay` 分别报告结构合法性、各关系状态、认证关系数及是否全部关系已认证。结构合法但含未证明想法的图很正常；它不是完整证明。没有关系的空图也不会报告全部关系已认证。

## 验证与当前边界

在本目录运行 `python3 -B -m unittest test_evidence_store -v`。测试临时文件都置于本目录并自动清理，不写冻结历史目录，不联网，不调用模型。测试证明的是列出能力及反例防护；不提供模型 token 或跨模型解题收益的实验。

本组件没有实现一般定理证明、数学等价搜索、跨原题自动迁移或密码学身份认证。哈希保证引用与内容一致，不能单独证明数学正确性或作者身份；认证仍依赖所选检查规则和原题绑定。
