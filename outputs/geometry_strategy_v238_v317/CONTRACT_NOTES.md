# 请求与证据边界：V238–V247

本组件只改进接口可靠性，不改变冻结后端的数学算法，不构成十项新数学方法。历史文件保持只读。

## 十步实际变更

每一步在修改实现和相关测试、运行累计测试之后调用 `campaign.record`，保存当时的两份 Python 源码。十份 `request_contract.py` 快照的 SHA-256 均不同，最终源码与 V247 快照一致。

| 版本 | 独立假设与变更 | 累计测试 | 实际对照 |
|---|---|---:|---|
| V238 | 重复键不能覆盖原题绑定 | 2 | 旧 JSON 解析保留最后一个 function；新入口拒绝，包含嵌套对象 |
| V239 | 非有限数和非对象根不能进入数学计算 | 4 | 拒绝 NaN、±Infinity、1e999，以及直接调用中的非有限浮点 |
| V240 | 显式 null 不能关闭期望条件 | 5 | 旧工具的 function、mean_zero、tolerance 三种 null 都可返回 verified；新入口拒绝 |
| V241 | 函数逼近证书不能满足谱结论请求 | 7 | 同一合法证书在 expected_kind=approximation 时通过，在 spectrum 时绑定失败 |
| V242 | 奇异搜索无效预算应先拒绝 | 9 | 四种无效输入旧路径调用来源计算 4 次，新路径 0 次；零次搜索仍合法 |
| V243 | 其他路由统一检查预算与类型 | 12 | 无效谱容差旧路径执行一个扇区，新路径 0 个；布尔值不能冒充整数预算 |
| V244 | 并行请求需要身份和语义关联 | 15 | 等价默认值及约分有理数保持同摘要；改容差改变摘要；错误期望摘要不派发 |
| V245 | 成功调用不等于完成数学目标 | 20 | 有效但区间过宽仍 certified_open；伪造 solver 状态无法绕过真实证书重放 |
| V246 | 失败类别应能指导后续处理 | 22 | 区分坏输入、缺证书、未知格式、坏证据、错目标及尚未达标 |
| V247 | 单次故障不能中止后续请求或变成数学反例 | 26 | 请求、执行器与验证器故障后下一请求继续；中断信号不吞掉 |

V242 开发时发生一次补丁位置错误导致的 SyntaxError，修复后才记录通过；失败已写入该步记录，没有把失败运行算成通过。另实测 JSONL 进程连续处理坏 JSON、合法但未达精度的谱请求、无效预算三行，全部输出对应记录且进程正常结束。

## 对外入口

- `parse_request(text)`：严格解析并校验一个 JSON 对象；错误抛出 `ContractError`，通过 `.code` 和 `.field` 读取机器可识别原因。
- `validate_request(request)`：检查操作、字段、原函数、预算、绑定条件及元数据，不执行求解。
- `execute(request, executor=None, verifier=None)`：默认调用冻结工具，数学证据另行重放；执行器或验证器故障返回结构化未知状态。请求自身无效时仍抛出 `ContractError`。
- `process_lines(lines, executor=None, verifier=None)`：逐行生成字典响应；跳过空行，保留原始 `line_number`，隔离普通异常。脚本入口将这些字典逐行编码输出。
- `evidence_state(request, certificate, verifier=None)`：可供未来新路由独立复用的证据/绑定/目标分类。
- `request_digest(request)`：对已校验请求的数学语义生成摘要；调用者应先校验。

新增元数据为 `expected_kind`、`request_id` 和 `expected_request_digest`，均在转交旧工具前移除。`expected_kind` 只有 `spectrum` 与 `approximation` 两值。请求 ID 是调用方关联标签，不是签名；本组件不保证不同请求 ID 的全局唯一性。

语义摘要会补全当前路由默认参数、归一化原函数和有理数；它包含期望结论类型和证书内容，但排除 request_id 与期望摘要本身。verify 请求没有提供的绑定仍保持未提供，不能用证书里的条件冒充调用方条件。

## 状态的解释

- `execution_ok`：求解或检查流程正常返回；不是证明。
- `certificate_valid`：实际检查了证书是否成立。故障时为 null，表示没有完成判断。
- `bindings_match`：提供的原函数、空间、容差和结论类型与有效证书匹配。
- `binding_complete`：谱验收有原函数、空间、目标容差和结论类型；函数逼近验收有原函数与结论类型。生成请求的操作与默认值可以提供这些条件，独立 verify 不能猜调用方意图。
- `verified`：证书有效且已提供的条件匹配。仍要看 binding_complete，不能单独把它用作原题完成声明。
- `target_met`：针对当前请求精度目标的判断，依据重放后的端点差精确计算。不看 solver 或证书自述的 status。没有当前目标或流程故障时为 null。
- `certificate_target_met`：有效谱证书相对于其内部声明容差的判断；与调用方是否提供完整目标分开。

因此，有效但未满足精度的谱结果为 `certified_open`；不带外部原题条件的重放为 `verified_unbound`；没有目标误差要求的函数逼近为 `certified_no_target`。不能把它们都算作完成了原题。

`reason_code` 区分：duplicate_key、malformed_json、nonfinite_number、invalid_root、null_binding、invalid_budget、invalid_rational、unknown_operation、unknown_field、invalid_type、invalid_kind、missing_function、request_digest_mismatch、invalid_request、missing_certificate、malformed_certificate、unsupported_certificate_format、invalid_evidence、conclusion_kind_mismatch、binding_mismatch、unbound_verification、target_not_met、no_target_requested，以及运行中的 backend_rejected、backend_failure、resource_limit、internal_error。一般错误文本只供解释，不是稳定的分类接口。

## 集成边界

验证器回调签名为 `verifier(cert, expected_function=None, expected_mean_zero=None, expected_tolerance=None) -> bool`，必须返回真正的布尔值。验证器是受信任组件，不是模型给出的函数；错误输入应返回 False，异常会被记录为验证器故障，不会当成数学反证。

新证书格式需要在 `CERTIFICATE_KINDS` 中明确登记，并配套正确验证器。现有 `validate_request` 只接受四个冻结操作；未来 research 操作应由集成层独立校验，然后复用 `evidence_state`。本轮没有冒然把未知操作或未知格式当成已支持。

资源计量和简短回执由后续阶段负责，这里没有自行保存数学结果或估算模型 token。整数/有理数预算检查不等于强制墙钟超时。同进程异常处理无法恢复操作系统终止、解释器崩溃或失控进程，也不是操作系统沙箱。通过现有精确验证器仍不是 Lean 等形式化证明。

## 复现

在项目根目录运行：

```text
python3 -B -m unittest discover -s outputs/geometry_strategy_v238_v317 -p test_request_contract.py -v
```

目前累计 26 个针对真实误用与故障的测试通过。各步完整对照和当时源码保存在 `iterations/v238` 至 `iterations/v247`。这些结果只支持接口可靠性及无效请求的零求解浪费，尚未测量有效任务速度、真实模型 token 或跨模型优势。
