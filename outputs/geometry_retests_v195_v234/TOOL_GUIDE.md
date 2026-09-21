# 数学研究工具使用说明

入口 `research_tool.py` 从标准输入逐行读取 JSON 请求，每行返回一个 JSON 结果。它只用 Python 标准库与工作区已有的冻结数学模块，不需要外部账号或模型 API。请保留同一 `outputs` 目录下的 V36–V90、V91–V94 和 V95–V194 依赖。

从工作区根目录运行：

```sh
python3 -B outputs/geometry_retests_v195_v234/research_tool.py
```

随后每次输入下面的一行。数学数值使用整数或精确有理字符串，例如 `"1/100000000"`；次数 n、L 等使用 JSON 整数。势的 `"a,b,c"` 表示 xᵃyᵇzᶜ 的系数。

## 原题的四种调用

完整球面函数同时正交于 1、z：

```json
{"op":"constraints","q":{},"constraints":[{"0,0,0":"1"},{"0,0,1":"1"}],"L":2}
```

返回精确区间 [2,2]。如需判断不等式，把 `op` 换为 `constrained_inequality` 并加入 `"threshold":"3"`；结果会包含反例决定，其完整证明中的试探函数可通过 `fetch` 读取。

两种空间下的负谱计数和完整负部迹：

```json
{"op":"spectrum","q":["-7"],"spaces":["full_mean_zero","full_unprojected"],"k":9,"quantities":["ordered","count","trace"]}
```

这里 q 的数组表示升幂的一元轴对称多项式 q(z)。结果按空间分别返回，不能混用。

原集中函数的精确矩与 GN 商：

```json
{"op":"cap","n":64,"original":true}
```

`original:true` 表示 `(1+t)^n-2^n/(n+1)`；`false` 表示 `y^n-1/(n+1)`。二者 GN 商相同，但质量、能量和四次矩不同。需要逐项核验已有原多项式时，调用 `cap_original` 并传入 `raw_t_coefficients`：完整的 n+1 个升幂系数。原题的 65 项输入与精确复核已保存在 `acceptance/C16.json` 和对应证明记录中。

三轴势的完整球面最低谱值：

```json
{"op":"parity_ground","q":{"2,0,0":"1","0,2,0":"2","0,0,2":"3"},"target":"1/100000000","start_L":3,"max_L":9,"bits":40,"wall_budget_seconds":90}
```

结果中 `status` 和 `exact_width` 决定目标是否完成。内部时间预算在步骤间检查；需要严格墙钟限时的调用方应另加外层进程限时，本轮验收即采用该做法。

## 让模型复用证明

每次成功计算会把完整证明保存在本目录 `service_certificates/`，只向模型返回结论、适用范围和内容哈希 `certificate_id`。复制返回的真实哈希调用：

```json
{"op":"verify","certificate_id":"把这里替换为返回的64字符哈希"}
```

`verify` 复核已有证明及数学输入绑定，不重新做优化搜索；需要展开详细证明时用 `fetch`。请求空间、势、原函数及目标宽度都参与绑定。运行耗时与缓存工作量属于观测信息，不属于数学证明；不要将证明复核误读为严格时间保证。

输入项数、次数或计算预算超过支持范围会被明确拒绝。一个有效但较宽的区间仍可能未达标。试探函数的 GN 商仅给出试探下界，不能作为全函数空间的最优上界。

## 本地复现

只读发布与证明核验：

```sh
python3 -B outputs/geometry_retests_v195_v234/verify_release.py
```

完整测试：

```sh
python3 -B -m unittest discover -s outputs/geometry_retests_v195_v234 -p 'test_*.py'
```

重新实际计算四个原题并复核旧题证明：

```sh
python3 -B outputs/geometry_retests_v195_v234/acceptance.py
```

最后这条会刷新本次验收记录和计时，从而使冻结的发布哈希发生变化；如需保留当前发布，请先复制整个项目或在副本中重测。
