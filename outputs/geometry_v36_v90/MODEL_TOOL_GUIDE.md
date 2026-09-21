# 给研究模型的工具说明

先用 `capabilities` 查看入口；每个请求是一行 JSON。数学字段用整数或有理字符串，避免把浮点数当作精确输入。几何 profile `sphere_ground` 明确指定单位 S²、包含零模的第一本征值、完整球面函数空间。坐标 `t=x_3` 在 `[-1,1]`。原始势必须对最多四个参数仿射、对 `t` 至多六次；罚项是参数二次式。

## 先选所需结论

原题是 `inf_a [lambda_1(-Delta+q(a,t))+penalty(a)] >= threshold`，且每个参数方向仿射于 `t` 时：

```json
{"op":"research","goal":{"profile":"sphere_ground","parameters":["a"],"potential":"a*t","penalty":"a*a/5","domain":{"kind":"all_real"},"threshold":"0"}}
```

默认采用有经典定理依赖的 log-Sobolev 规则；加 `"rule":"barta"` 使用已验证的正函数全参数规则。`research` 对不支持的输入报错，不会悄悄更改题意。其他多项式参数方向使用旧入口 `submit`，其状态和失败诊断仍绑定原题。

需要研究统一二次系数时：

```json
{"op":"uniform","spec":{"direction":[0,1]}}
{"op":"uniform","rule":"barta","spec":{"direction":[0,1],"coefficient":"1/5","all_real":true}}
```

第一条返回最优系数 `1/6` 的解析规则证据；第二条生成全实参数的张量覆盖。省略 `all_real` 的 Barta 请求只覆盖其显示的有限幅度区间。

固定势的正函数搜索与全局矩阵优化：

```json
{"op":"positive","spec":{"potential":[0,0,6]},"options":{"degree":6,"range_backend":"exact_sturm_minimum_v1"}}
{"op":"refine","spec":{"directions":[[0,1],[0,0,1]]},"tolerance":"1/1000000","budget":8}
{"op":"extrema","spec":{"polynomial":[0,1,0,-1],"interval":[-1,1]},"tolerance":"1/1000000000000","max_nodes":32}
```

默认 `positive` 用有理 Bernstein 范围证书；示例显式选择登记的 Sturm 后端。后端是程序内预先登记的验证器，JSON 不能提供自定义执行代码。`refine` 的 `budget` 可以是轮数，或包含 `rounds,newton_steps,range_leaves,max_samples,max_denominator,precondition,dominance,sparsify` 的受检对象，详见其实现/证明文档。更多预算只增加搜索机会，不保证任一指标单调改善。

## 解读回执

返回外层的 `text` 是决策 JSON，`text_tokens` 是该文本计数，`wire_tokens` 是整个已知外层 JSON 的计数。默认预算 500，可指定 `token_limit`。`fits=false` 意味着必需信息容纳不下，不是数学反例。

| 字段或状态 | 正确理解 |
|---|---|
| `ground_goal` 的 `proved` | 已证明原目标下界不低于 threshold |
| `ground_goal` 的 `refuted` | 具体参数和真实试探函数给出严格反例 |
| `unresolved` / `open` | 当前证据或预算不足；已有界仍可能有效 |
| `ansatz_gap_closed` | 指定对数多项式函数族的全参数最优值已在容差内包围 |
| `global_gap_closed` | Poisson 矩阵包络的全局原始/对偶间隙满足容差 |
| `ansatz_obstruction` | 这个正函数候选失败，不能推断原谱不等式为假 |
| `uniform` 的 `disproved` | 此格式包含原谱不等式的实际反例 |
| `spectral_upper` | 由独立 Rayleigh 试探给出；缺失/null 时没有谱上界 |
| `ansatz_upper` | 仅界定所选函数族的最佳 Barta 下界；不能替代谱上界 |
| `analytic_dependency` | 这条结果依赖列明的经典解析定理，程序未证明它 |

展示下界向下取整，上界和误差向上取整。`gap` 包围的是精确证书中两界的差；显示取整后的两端之差可能更大。最终判断使用精确有理数，不使用小数显示。

## 只读取下一步所需证据

`task` 是精确问题的内容引用，`evidence` 是完整证据的引用。将实际返回的引用填入后续请求，不要猜测引用：

```json
{"op":"inspect","ref":"<返回的 task 引用>"}
{"op":"inspect","ref":"<返回的 evidence 引用>","pointer":"/certificate"}
{"op":"inspect","ref":"<返回的 evidence 引用>","pointer":"/certificate/lower_certificate/theta"}
```

JSON Pointer 路径取决于证书类型；选择字段减少读取成本。数组可以加 `offset` 和 `limit` 分页；必须检查 `total`、`next`。读完整证明时，所有取回内容都会耗费文本 token。

对旧结构化目标，可 `register` 一次，再 `solve`/`revise` 或批量 `batch`。仅完全相同的目标、求解选项和引理库复用缓存。`exceptions_only` 保留每题状态清单，`deduplicate` 以 `unique+order` 无损编码重复结果。引理库只接收原 V35 Poisson 格式；新格式不能冒充旧格式加入库。

`checkpoint`/`restore` 恢复的是精确任务、证据引用及已验证旧引理。它们不会恢复一个并不存在的数值迭代中间态。保留 `--store` 指定目录才能恢复；跨机器可以导出完整证书再核验。

## 计量边界

o200k_base、cl100k_base 只用于公开词表的真实可见文本编码。没有声称它们对应 GPT-5.6 或 Astra 的实际模型映射；未计未知 API framing、隐藏推理或模型理解错误。基准计入明确给出的能力说明、操作列表、请求、完整已知返回包装和证据取回；这份额外说明若放进模型上下文，也需另计。程序本身不调用模型。
