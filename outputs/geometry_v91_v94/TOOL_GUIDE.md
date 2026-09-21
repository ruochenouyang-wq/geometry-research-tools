# V91–V94 独立研究接口

在本目录运行 `python3 research_tool.py`，从标准输入逐行发送 JSON 对象。每个非空输入行恰好返回一个 JSON 行，无附加日志。只使用本地 Python 数学程序，不执行客户字符串，不联网，也没有模型或 token 节省评测声明。

| op | 必需字段 | 可选字段 |
|---|---|---|
| `sector` | `q` | `m=0, mean_zero=true, k=1, modes=12, bits=44` |
| `mean_zero_ground` | `q` | `modes=8, max_modes=32, bits=44, max_m=16, tolerance="1/10000000000"` |
| `weighted_poincare` | `q, threshold` | 与 `mean_zero_ground` 相同 |
| `verify` | `certificate_id` | 无 |
| `fetch` | `certificate_id` | 无 |

字段在请求顶层，不嵌套于 `options`。未知字段全部拒绝。后两种数学操作固定在完整球面的均值零函数空间，**不接受 `mean_zero` 覆盖字段**，即使其值为 true。`sector` 只能证明一个方位角扇区，不能把它当作完整球面结论。

生成证书还必须绑定本次请求的容差及截断规模：扇区证书的 `modes` 与请求相同；完整谱证书的 tolerance 必须精确相等，已有扇区的 `m` 不得超过 `max_m`，每个扇区 modes 必须处于请求的 `[modes,max_modes]`。同题但宽松容差的合法证书不能替代严格请求。`bits` 会做输入校验，但证书不记录执行步数，接口不声称独立验证实际迭代次数。

`q` 是升幂有理系数列表，次数最多 6；系数绝对值最多 1000、分母最多 `10^9`。JSON 中使用整数或 `"-3/7"` 形式字符串，拒绝浮点数、布尔值、表达式和小数字符串。`threshold` 和 `tolerance` 同样要求精确有理数。`m/max_m` 在 0..64，`modes/max_modes` 在 1..64，`1<=k<=modes`，`modes<=max_modes`，`bits` 在 8..160；`tolerance` 在 `[1/10^30,1]`。整数预算不接受布尔或浮点值。

示例：

```json
{"op":"sector","q":[0],"modes":2}
{"op":"mean_zero_ground","q":[0],"modes":2,"max_modes":2,"max_m":1}
{"op":"weighted_poincare","q":[0],"threshold":3,"modes":2,"max_modes":2,"max_m":1}
```

完整球面的非恒定势演示（已实际经 JSONL 运行并跨进程重放）：

```sh
python3 research_tool.py <<'JSONL'
{"op":"weighted_poincare","q":[0,0,2],"threshold":"238265540615/100000000000"}
{"op":"weighted_poincare","q":[0,0,2],"threshold":"12/5"}
JSONL
```

第一行返回 `proved`，第二行返回 `refuted_by_spectral_existence`。量词都是**完整单位球面上的所有均值零 H¹ 函数**，不是选定的二维试探空间。第二行只认证违反阈值的函数存在，没有构造显式反例函数；使用返回的精确有理上下界或 `fetch` 中的完整证书核查结论。

默认返回 `status, scope, lower, upper, width, certificate_id`。界和宽度使用精确分数字符串，避免十进制近似改变含义。scope 明确记录原势、均值零约束、谱序号；扇区结果另外记录 `m` 和投影。`weighted_poincare` 的 scope 同时记录阈值。

- `certified_sector_bound`：只认证所请求扇区的指定本征值包围。
- `certified_target_met` / `certified_bound_open_gap`：完整均值零谱界达到/尚未达到目标宽度；两者的实际界都经过认证。
- `proved`：完整均值零谱下界不低于阈值。
- `refuted_by_spectral_existence`：完整均值零谱上界严格低于阈值，证明存在违反不等式的函数；**没有给出显式反例函数**。回执保留 `refutation_does_not_include_an_explicit_function=true`。
- `undetermined`：当前有效谱包围尚不能决定阈值问题，不等于真或假。

完整证书存放在本新目录的 `certificates/service`，ID 为规范 JSON 内容的完整 64 位小写 SHA-256。复制返回的真实 ID 即可重放或提取：

```json
{"op":"verify","certificate_id":"替换为实际返回的64位小写hex ID"}
{"op":"fetch","certificate_id":"替换为实际返回的64位小写hex ID"}
```

`verify` 重新校验内容哈希、证明格式及数学证据，返回原结论和 `verified=true`。`fetch` 做相同校验并额外返回完整 `certificate`；不能用它跳过证明验证。两者只接受 ID，不接受任意路径、外部证书或客户指定作用域，也拒绝证书目录中的符号链接。证书跨进程保留，但依赖原证书目录；本接口没有宣称持久化搜索缓存。

输入错误、未知操作、证据缺失或篡改返回 `status="invalid_request"`；底层算术异常返回 `status="computation_failed"`。两者均附错误类型、说明及 `mathematical_verdict=null`，后续请求继续处理。`verified` 表示本程序重放数学证据；证书中的 `formal_assistant_checked` 仍为 false。

运行接口回归测试：`python3 -m unittest -v test_research_tool`。测试使用隔离的临时证书目录，不修改 V90，也不改正式服务证书。
