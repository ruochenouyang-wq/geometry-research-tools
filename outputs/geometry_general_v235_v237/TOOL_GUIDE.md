# 三版一般函数接口

在工作区根目录运行：

```sh
python3 -B outputs/geometry_general_v235_v237/research_tool.py
```

每行输入一个 JSON 请求。返回原函数、完整函数空间、谱值区间、目标是否完成、误差分项和内容哈希。完整证书留在本目录 `service_certificates/`；不需要外部服务或模型 API。

## V235：光滑复合函数

下面的原函数是 exp(z/2)+sin(z/3)+log(1+z/4)：

```json
{"op":"solve","method":"linf","order":12,"function":{"kind":"analytic_sum","terms":[{"function":"exp","argument":{"0,0,1":"1/2"}},{"function":"sin","argument":{"0,0,1":"1/3"}},{"function":"log1p","argument":{"0,0,1":"1/4"}}]}}
```

`analytic_sum` 可含 `polynomial` 和最多 12 个 `terms`。每项支持 `exp`、`sin` 或 `log1p`，`argument` 是 xyz 多项式，`coefficient` 默认为 1。指数键 `"a,b,c"` 表示 xᵃyᵇzᶜ。`log1p(t)` 表示 log(1+t)，当前要求参数多项式系数绝对值之和小于 1；这是可验证的充分条件，可能拒绝部分实际上有定义的表达式。

## V236：跳跃或无界势

```json
{"op":"solve","method":"l2","degree":8,"function":{"kind":"axis_profile","axis":2,"profile":"step"}}
```

这里原势为 1_{z>0}。`axis` 为 0、1、2，分别对应 x、y、z；`amplitude` 默认为 1，`offset` 默认为 0。

负奇异势 −|z|⁻¹ᐟ⁴：

```json
{"op":"solve","method":"l2","degree":8,"function":{"kind":"axis_profile","axis":2,"profile":"abs_power","exponent":"-1/4","amplitude":"-1"}}
```

`exponent` 必须大于 −1/2，以保证概率 L² 可积。奇点本身按几乎处处等价的代表处理，不能据此使用点态有界势公式。

## V237：固定参考能量

```json
{"op":"solve","method":"fixed","degree":12,"reference":{"0,0,0":"-4/3"},"function":{"kind":"axis_profile","axis":2,"profile":"abs_power","exponent":"-1/4","amplitude":"-1"}}
```

`reference` 可以是已支持的有理多项式。必须证明它到当前逼近的距离 θ<1，且真实余项除以 1−θ 后仍小于 1。固定参考并不总是优于普通 L² 传递，应比较最终区间。报告中的参考为事先固定的 4 次投影，其完整系数保存在相应证书内；上例用较简单的常数参考，数值结果不同。

## 通用设置与证据复用

`mean_zero` 默认为 true，改为 false 会计算不同的完整函数空间。`target_width` 默认为 `"1/100000000"`，仅控制完成状态；函数逼近次数不会自动增长。`source_options` 可设置源谱求解的 `modes`、`max_modes`、`max_m`、`bits`、`tolerance`，多轴后端可用 `L`。这些是求解设置，不是经证明的运行耗时保证。

所有数学系数使用整数或有理字符串；浮点数不会被当作精确数据。`order`、`degree` 等次数使用 JSON 整数。

把返回的真实 `certificate_id` 放入以下请求：

```json
{"op":"verify","certificate_id":"替换为返回的64字符哈希"}
```

`verify` 不重新搜索谱值。`fetch` 则在复核后返回完整证明。证明验证会检查原函数、近似多项式、范数与测度、完整空间源谱、传递参数和端点；只传一个没有证明的误差数字不能生成证书。

`target_met` 表示总宽达到目标；`certified_open` 表示有有效界但仍未达目标；`not_certified` 表示本次输入、预算或认证条件没有通过。后两者不能视为高精度结果。

## 发布检查与重算

```sh
python3 -B outputs/geometry_general_v235_v237/verify_release.py
python3 -B -m unittest discover -s outputs/geometry_general_v235_v237 -p 'test_*.py'
```

重算全部示例：

```sh
python3 -B outputs/geometry_general_v235_v237/run_examples.py
```

重算会更新本版的计时记录和证据文件，影响发布哈希。如需保留当前发布，请在完整项目副本中运行。保留相邻旧版本目录，它们是只读计算依赖。
