# 有条件的结论复用：v288–297

本组件不搜索新特征函数，而是复用已经验证的结论。范围固定为单位球面 S²，势函数沿坐标轴变化；旧证书只读。所有数学源证据嵌入新证书并重新检查，不以摘要、文件名或缓存命中代替证明。实现使用既有的谱事实和不变性，不主张这些数学规则原创。

## 数学规则与条件

1. **常数势。** 对 q=c，球面 Laplace 特征值为 l(l+1)。完整 H¹ 空间允许 l=0，最低值为 c；零均值空间剔除常数，最低值为 2+c。常数生成器要求原函数的 amplitude=0，不把非零幅度误当常数。`constant_for_function` 保留原 profile、exponent、axis 与 offset。
2. **轴旋转。** 把一个坐标轴旋转为另一个坐标轴可由球面等距变换实现。它保持归一化面积测度、梯度能量、H¹ 空间与零均值条件。因此最低谱值不变；L² 函数误差也不变。仅支持轴 0、1、2，不接受任意未证明映射。
3. **谱常数平移。** 在任一固定允许空间中，单位范数函数的 Rayleigh 商在 q→q+c 后增加 c，故上下界同步增加 c，宽度不变。不能把 q→a q 当成整个算子的谱缩放，因为梯度项未同步缩放。
4. **函数逼近仿射变换。** 若 ‖q−p‖₂≤ε，则 ‖(a q+b)−(a p+b)‖₂≤|a|ε。支持负 a、a=0，以及常数 b；证书明确保留 L² 范围，不推导谱误差。
5. **谱区间交集。** 对同一规范化函数、同一空间的全部已验证谱区间，取最大下界与最小上界。每个前提均须检查；不同原题、不同空间、空交集均拒绝。
6. **逼近证据选择。** 对同一目标函数的已验证 L² 证书，选误差上界最小者，同时保存其源索引。不同逼近函数不能只抽取一个较小数字而丢掉对应代表元。`evaluate` 沿被选证据评价实际代表元。
7. **链规范化。** 先验证原链，再按顺序合成外层轴变换和仿射变换。谱链只能 factor=1；L² 链允许有理 factor。规范化保留最内层数学源证据；不会删除区间交集或证据选择中仍需的前提。

## API

```python
constant_spectrum(value=0, mean_zero=True, axis=2)
constant_for_function(source_function, mean_zero=True)
spectral_axis(source, axis, verifier=None)
spectral_shift(source, shift, verifier=None)
approximation_axis(source, axis, verifier=None)
approximation_shift(source, shift, verifier=None)
approximation_scale(source, factor, verifier=None)
spectral_intersection(sources, verifier=None)
approximation_select(sources, verifier=None)
normalize_transport(certificate, verifier=None)
evaluate(certificate, t, verifier=None)
verify(certificate, expected_function=None, expected_mean_zero=None,
       expected_tolerance=None, verifier=None)
```

返回格式 `geometry_proof_transport_v1`；`FORMATS` 可用于注册。谱证书提供 function、mean_zero、scope、lower、upper、exact_width 与 full_infinite_space_covered；逼近证书提供 function、scope、norm、error_upper 与 spectral_transfer_claimed=False。

`expected_tolerance` 是调用方要求：检查谱区间宽度或 L² 误差上界是否不大于它，不在源证书上虚构不同目标。逼近证书不能绑定 mean_zero。`expected_function` 要求规范化表示严格相等；虽然不同轴上的常数函数数学等价，也不会静默忽略目标轴字段。

`verifier(source) -> bool` 可接入其他新格式的数学验证器。回调属于可信验证边界，必须检查数学，不应仅查摘要。默认回放旧完整谱、分数幂谱和分段 L² 格式。额外谱格式须具有规范化 function、geometry=unit_S2、measure=d_sigma/(4*pi)_on_unit_S2、明确 mean_zero 与对应 scope、eigenvalue_index=1、full_infinite_space_covered=True、lower/upper。额外逼近格式须提供 probability L² scope/norm/error_upper，且 spectral_transfer_claimed=False。

`evaluate` 的 t 是目标轴上的坐标，范围 [-1,1]，不接收三维位置。它支持旧分段代表元与本模块变换；对没有注册代表元评价器的外部格式，证明可以有效但评价明确拒绝。

## 10 轮记录与实际检查

| 版本 | 新能力 | 当时累计针对测试 |
|---|---|---:|
| 288 | 两种空间的常数势解析证书 | 2 |
| 289 | 谱结论的轴旋转 | 3 |
| 290 | 谱结论的常数平移 | 4 |
| 291 | 保留 L² 范围的轴旋转 | 5 |
| 292 | L² 常数平移与代表元评价 | 6 |
| 293 | L² 有符号幅度缩放 | 7 |
| 294 | 相同原题谱区间交集 | 8 |
| 295 | 相同目标的 L² 证据选择 | 9 |
| 296 | 有证明的多步变换合成 | 10 |
| 297 | 原题/空间/目标绑定、范围防护和严格常数接口 | 14 |

每一步测试通过后才调用 `campaign.record` 保存当时源代码。记录位于 `iterations/v288` 至 `iterations/v297`，不是将最终同一份代码追补十次。

开发用五步逼近变换链，规范化前后同一原题、同一误差与代表元相同。后续回放中的运输规则检查从 5 次减为 1 次；规范 JSON 从 3037 字节减为 1528 字节。两者都仍重放原数学源证据。这不是 token 计费实验或模型提速证据；规范化本身先验证原链，准备成本也存在。

## 限制与保留的拒绝案例

14 项最终检查覆盖真实旧谱证书、奇异函数分段证书、常数解析谱、不同空间、仿射顺序、代表元点值、篡改端点、错误轴、错误目标、未知参数、循环来源和外部格式回调。全部通过。被拒绝的数学运输是功能边界，不改标成成功：L²→谱、势幅度→谱整体缩放、不同原题证据交集、非零幅度直接走常数捷径。

证明图最多深 32 层、256 个运输节点、每步 32 个源证据；超过范围拒绝，不省略前提。区间与误差是已有证据的可靠复用，尚未自动发现新变换、尚未覆盖一般曲面、任意旋转矩阵或跨空间谱比较，也未接入形式化证明助手。

复现：`python3 -B -m unittest discover -s outputs/geometry_strategy_v238_v317 -p test_proof_transport.py -v`。
