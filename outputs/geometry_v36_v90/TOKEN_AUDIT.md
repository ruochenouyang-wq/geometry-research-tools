# V36–V55 token 协议独立审计

2026-09-11，Python 3.9.6。新增 `test_token_protocol.py`，共 **60 项测试，全部通过，实测 1.840 秒**。审计工作未修改 `protocol.py`、`token_meter.py` 或版本台账；下述核心修复由主任务完成，再由独立测试核验。

## 复现

在 `outputs/geometry_v36_v90` 目录运行：

```sh
python3 -m unittest -v test_token_protocol
```

测试的证据库及 CLI 请求文件全部使用自动清理的临时目录；不读取或修改 `results` 基准目录，不运行耗时数值基准，不调用模型。真实词表测试强制使用本地 `token_cache`，并封锁该测试中的 socket 连接。其余 token 计数和 CLI 使用同一现有本地词表，无网络请求。

固定数学样本由现有 V35 算法现场生成，再用独立 verifier 重放：

| 原题 | 参数域、阈值 | 实际结论 |
|---|---|---|
| `a*t`，罚项 `a*a/4` | 全实数，`0` | `proved` |
| `a*t`，罚项 `a*a/5` | `[-2,2]`，`0` | `unresolved` |
| `a*t`，罚项 `a*a/10` | `[-4,4]`，`-1/10` | `refuted` |

均显式绑定单位球面、完整球面函数空间和第一本征值。小型求解上限为 `max_leaves=1, max_modes=4, max_states=4`。测试同时包含精确标量族证书与矩阵包络证书，不依赖优化基准结果。

## 覆盖范围

| 版本 | 检验内容 |
|---|---|
| V36 | `cl100k_base` 与 `o200k_base` 真实已知词条、规范 JSON、显式会话项逐项计数、不支持的编码拒绝 |
| V37–V38 | 三种结论保持不变，去除未受信任的包装层声明，原题作用域可见，原证书可取回；正负分数、数量级边界及小/大数的精确向外舍入 |
| V39–V41 | 内容寻址稳定性、输入别名隔离、修改/缺失证据检测、路径注入拒绝、JSON Pointer、完整分页重组及结束标记 |
| V42–V45 | 显式 profile 与完整目标等价，禁止隐式补充题目条件，打包多项式/三角 Hessian 与独立手写表达式精确等价，零参数目标，补丁保留未声明条件并重新验证 |
| V46–V50 | 批量完整 inventory、异常结论筛选、同 Service 缓存命中与证据重验、受检 lemma 复用、索引隔离、诊断合法性及原题归属 |
| V51–V55 | 实际文本 token 预算、诊断缩略/省略视图、不足时显式报错且不截断、轮询、检查点和库引用恢复、上下文、单调用提交以及实际 JSONL/CLI 链路 |

篡改检验独立修改状态、下界、原题阈值、函数空间、Poisson 解、诊断尺度、原题绑定和证据文件；这些篡改均被拒绝。将另一原题的**合法**诊断移植到当前证明的对照也被拒绝，因此该检查不只是侦测格式错误。

## 审计发现与修复复验

以下问题均在最初实现上实际复现，修复后的上述运行全部通过对应回归测试。

| 问题 | 最小触发方式 | 修复后行为 |
|---|---|---|
| 合法省略视图无法进入后续流程 | 对 `receipt(..., diagnostics='omit')` 调用 `budget`、`poll`、`checkpoint` | 接受并验证省略视图；`available_in_evidence` 缩略视图亦可复用 |
| lemma 索引与不可变证据失配 | `add_lemma([lemma])` 后修改原始 `lemma['directions']`，或修改 `lemma_index()` 返回值 | 输入和输出均隔离；索引仍对应实际存储内容 |
| 原题引用类型未验证 | 把 `store.put('context', 原题)` 当成回执或检查点的 goal 引用 | 即使内容相同，也拒绝错误引用类型 |
| 无关原题诊断可移植 | 将 `a*t**3`、`a*a/100`、盒域的合法诊断贴到 `a*t`、`a*a/4`、全实数原题证明 | 诊断必须对应当前原题可验证的约化表示 |
| JSON Pointer 接受非法转义 | `project({'~2': 1}, '/~2')`，以及孤立 `~` | 拒绝非法转义，保留合法 `~0`/`~1` 选择 |
| 恢复未检查证据及库引用类型 | 在检查点中用内容合法的 `context:` 引用冒充 proof 或 lemma | 恢复拒绝错误类型；另覆盖 goal/proof/lemma 任一缺失情况 |
| 无效表达式终止 JSONL 会话 | 先提交 potential `(` 或 `1/0`，再提交合法原题 | 返回结构化错误，`mathematical_verdict=null`，后续请求正常运行 |

回归测试名集中包含 `omitted_diagnostic`、`lemma_index`、`reference_kind`、`unrelated_problem_transplant`、`invalid_escape`、`restore_rejects_wrong` 与 `invalid_expression_is_nonfatal`，可单独重跑定位。

## 结果边界

- `hello world` 的本地词条对照为：`cl100k_base=[15339,1917]`，`o200k_base=[24912,2375]`。预算检验比较的是相同编码下完整输出文本的实际编码长度。
- 测试未发生模型/API 调用。文本 token 计数不是计费量、隐藏推理 token、模型效果或模型成本的测量；会话报告明确保留 `actual_api_usage=null` 与 `model_mapping_asserted=false`。
- `proved` 表示现有有理证书重放通过，未宣称形式化证明助手内核验证。族与矩阵包络回执明确保留 Poisson 指数 ansatz 的范围。
- 缓存测试证明相同 Service/JSONL 会话中复用已验证结果；跨 Service 测试证明 lemma 库可重新加载，未据此宣称独立 CLI 进程间共享求解缓存。
- 此审计验证协议保真性及轻量端到端行为，不替代 V56–V90 数学精度、全算法覆盖或性能基准评估。

被测文件 SHA-256（便于后续改动判断是否需复验）：

```text
protocol.py             d70680868d7c4170403d733db2ee85498f2438caeb32973f65d632ab8dbe0a28
token_meter.py          f11f46e925810db7cd7173ca71b3384a13c1fea26072fafa3d8093142e9ab9b2
run.py                  8b1098e68b519a597266c0887b699162f995cc135e1c8a0bf1475a9c417ea222
test_token_protocol.py  e62889a535e928b58c2b185828856b911fb62dd092af8f97bc41add1d44cc2b5
```

## 最终集成复验补记

上面的时间与哈希是V55阶段审计快照。V90随后增加完整返回包装计量、新数学格式和逐任务上下文范围；最终489项测试包括这60项全部通过。当前文件指纹与完整包预算复验见 `REVIEW_INTEGRATION.md` 和 `results/final_validation.json`。`verified`表示证据重放，原不等式还需成功状态；它不把unresolved变成proved。
