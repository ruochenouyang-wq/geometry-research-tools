# V90 研究接口整合独立审查

审查日期：2026-09-11，Python 3.9.6。审查 `research_registry.py`、`protocol.py` 新研究操作、`run.py`、`test_research_service.py`，并检查其与 V88/V89 当前接口的绑定。

新增 `test_integration_audit.py` 共 **21 项独立测试**，与已有研究接口的 22 项合跑，**43 项全部通过，2.517 秒**。本次未修改被审整合代码，未重跑历史全套测试或长基准；发现的问题由主任务修复，再独立复验。

## 结论与测试范围

所检查的新证书摘要保持实际的数学作用域与结论，未发现把函数族上界当成谱上界、把 ansatz 障碍当成原题反例，或将异题证据绑定到当前任务的接受路径。

| 证书类别 | 摘要保留的意义 |
|---|---|
| Sturm 极值 | 原多项式与精确原区间上的全局最大值；无球面几何或谱结论标签 |
| 自由正函数优化 | 指定对数基的最佳 Barta 值；`ansatz_upper_is_spectral_upper=False`，独立谱界分开保存 |
| V88 Poisson 族 | 全实振幅、完整球面及 Poisson 指数 ansatz；上界来自完整区间梯度商证明 |
| V88 矩阵优化 | 原 directions/cost 的 `trace(C G)` 包络目标，范围限于 Poisson ansatz |
| 一致谱不等式 | 全实数或明确的有限振幅域；成功、实际反例及 ansatz 障碍分别保留 |
| V89 原题 | 原势、罚项、参数域和阈值绑定，界控制 `inf_parameters(lambda_1(-Delta+q)+penalty)` |

测试中构造了上述类别的 **10 个证书/状态样本**，包括原题 `proved`、`refuted`、`unresolved`，一致不等式 `proved`、`disproved`、`ansatz_obstruction`。所有公开 bound 字段逐项与精确证书对照。

## 舍入、作用域与经典定理信任

- 名为 `lower` 或以 `_lower` 结尾的数值均向下舍入；upper、gap 和误差等字段向上舍入。`None` 保留未知界，不转换为零。
- 正函数精确标定例具有 `ansatz_gap=0` 但 `spectral_gap>0`，经过摘要后仍如此。原题未获有限下界时，摘要保留 `lower=None`。
- `counterexample_objective_upper` 是实际反例目标的上界，向上舍入后仍保持严格负号；`ansatz_obstruction` 不输出伪反例界。
- 采用 log-Sobolev 的一致不等式与 V89 结论都保留 `analytic_dependency`，明确写出经典 S2 log-Sobolev 定理**不是本程序证明的**。`formal_kernel=False` 保持不变。
- 该依赖声明经过完整预算包装、证据回取和检查点恢复仍保留。试图删除回执中的声明会被回执重建校验拒绝。
- 因此 `verified=True` 应理解为**在所示依赖条件下证书重放通过**，不是宣称程序单独证明了经典解析定理。Barta 路径的摘要没有加入 log-Sobolev 依赖。

## 注册表、证据与任务绑定

`research_registry.backend` 只接受内置 Bernstein 与固定的 `exact_sturm_minimum_v1` 名称，后者加载仓库中预定的 `ExactRangeBackend`。任意 JSON 名称、模块路径、路径穿越字符串或对象不能成为被信任的验证器。

独立对照将一个**自身有效**的 Sturm 区间证明改为 `[0,1]`，再拼接到要求 `[-1,1]` 的正函数证书。注册表拒绝该跨域证据，说明不只是检查底层证明内部自洽。

所有新回执都能在新 Service 实例中恢复精确任务引用与受检证据。错误任务内容绑定、证据文件被改写、移除信任说明等情况均被拒绝。投影前先重验 proof，不能通过只取一个字段绕过证据完整性。

V88 精确族有自己的证书格式；它虽能由新注册表验证，但不能未经转换直接进入 V35 仅理解旧 Poisson 格式的 lemma 库。独立测试确认该格式混淆被拒绝。

## JSONL、完整文本预算与错误恢复

实际子进程运行了五项新操作：`extrema`、`positive`、`uniform`、`refine`、`research`。使用有界设置：极值最多 4 个细化节点，矩阵零优化轮，正函数 9 点网格、1 次交换、每温度 1 次 Newton。每个成功结果都能从其 proof 引用重验。

预算测试独立覆盖 `o200k_base` 与 `cl100k_base`：

1. `text_tokens` 等于内层回执文本的实际编码长度。
2. `wire_tokens` 等于包含 JSON 转义、编码名、fits、token 字段等**完整外层结果**的实际编码长度。
3. 对真实 JSONL 输出行也直接计数，包括本次实际输出的结尾换行，结果与 `wire_tokens` 一致。
4. 设置预算恰等于完整包长度时成功；不足时返回 `fits=False`，无截断 text、无伪数学 status。
5. 两种合法展示形态间选择较小的完整包需求值，避免把更长的诊断替代说明当成最低需求。

以上数字只度量明确文本编码，不是模型 API 计费、隐藏推理或外部传输框架开销。`fits=False` 明确表示无法满足预算，不声称报错消息本身一定能装入极小预算。

## 实际发现并修复的问题

### I1：错误类型引用令 JSONL 会话退出

最小输入：

```json
{"op":"solve","goal_ref":null}
{"op":"extrema","spec":{"polynomial":[0,1]}}
```

初版第一行触发 `.startswith` 的 `AttributeError`；JSONL 异常处理未覆盖此类型，进程退出 1，后续合法请求不执行。`restore(ref=null)` 及整数引用有同类问题。

主任务修复后，solve/restore 严格检查引用为字符串，JSONL 同时补充 `AttributeError` 与 `ArithmeticError` 的错误处理。独立回归确认三种坏引用均返回结构化错误，`mathematical_verdict=null`，随后合法极值请求继续成功。

### I2：不足预算时高报可行完整包的最低需求

初版 log-Sobolev 回执在该复现中完整包为 **208 tokens**；以 **207** 为预算时，错误返回 `required_wire_tokens=211`，原因是把空 diagnostics 换为 `available_in_evidence` 反而增加长度，但仍只报告后者需求。实际上 208 已能容纳完整包。

主任务修复为按 `wire_tokens` 选择 full/short 中较小的完整包，再报告相应需求。两个编码的独立边界回归均通过。

### I3：研究混合会话头部沿用统一球面 profile

旧 `context` 在只含多项式极值或混合研究回执时，仍输出顶层 `profile='sphere_ground'`。单个回执的 scope 正确，但会话头部可能误导使用者将多项式极值也理解为球面谱问题。

主任务将其改为 `available_input_profile='sphere_ground'` 与 `scope='per_task'`，明确 profile 只是可用输入缩写，数学作用域依各任务而定。新增回归分别检查纯 extrema 与混合 extrema/原题/正函数会话，恢复后的回执也保持原 scope。

## 复现与文件指纹

在本目录运行：

```sh
python3 -m unittest -v test_integration_audit test_research_service
```

本次结果：`Ran 43 tests in 2.517s — OK`。全部可写证据与 CLI 状态均在自动清理的临时目录中。未接触结果基准目录。

```text
research_registry.py
780cb1b17f14dce579308c8ec3602f67e96fa0c8890199182e7cb66a1fcf5574
protocol.py
73eb1841cb2b88341eab0d0bdc404da1c1f094f3ff37b6618ab2aa49913ae171
run.py
ca643a40a0e9c9b5f846da4be42e8d03433318afbe289031082b92676307c680
test_research_service.py
9c7213608be8bfef1fd042107ecf20e0da680124586771aa93d77abca68d7b8f
test_integration_audit.py
dcef4907a86a1fe85903d912c2e46449fe6f3e96511146d1030b22c91dce87e9
precision_bridge.py
875b0fb51f83395089cb7c1a072734ebee61eacd3adde9407e0cbdb5dcd6e744
original_goal.py
0e0718b0dc2aef155ac09a968cc233e111fa03ad6ac56dafeb7677bc4c4a2274
```

本审查检查的是整合层的接受与展示边界，不能替代各数学模块的独立证明审查，也未声称任意输入规模的性能覆盖。
