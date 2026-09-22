# C3：只复用构造工作，保留原验证边界

本模块实现 `UnverifiedAssemblyContext`，即“尚未认证的组装上下文”。它仅生成旧格式的全空间证书候选；它既不验证谱间隙，也不宣告数学成功。控制器必须把候选序列化，再对原任务调用冻结的 `baseline.assess`，完整重建谱间隙、矩阵、统计量和最终区间。不得以候选的状态字段、摘要或上下文身份替代这一步。

## API 与适用范围

```python
ctx = prepare_fullspace(task, gap)
pending = ctx.assemble(powers, coefficients,
                       prepared=prepared_trial, statistics=exact_statistics)
# 控制器负责预算、序列化、原验证器冷核验及最终成功状态。
```

`prepared` 和 `statistics` 均可省略。`prepared` 接受 `function_key`、`powers`、`M/H/R` 属性或映射字段，与 `fast_candidate.PreparedTrial` 对接；矩阵必须为精确有理数、方阵且严格对称。`function_key` 为规范化原函数的标准 JSON。省略 `prepared` 时只在这里生成一次原始 Gram 矩阵。

还提供 `ctx.evaluate(...)` 与 `ctx.assemble_evaluation(receipt)`：前者冻结矩阵并计算统计，后者组装同一上下文中的结果。receipt 是防止误混用的局部算术记录，不是数学证明。对相同任务重新建立的另一个上下文也不能复用该 receipt。

当前仅适用于单位球面 S²、`kind='spectrum'`、`mean_zero=False`、负幅度且指数位于 `(-1/2,0)` 的 `abs_power` 势。谱间隙必须声明完整实径向 `m=0` 空间的第二特征值。没有把全空间结论移植给平均零空间，也没有修改阶跃路径；它们继续使用各自的旧构造和验证。

## 信任边界与输入绑定

上下文先把原任务和 gap 复制为私有 JSON 快照，再检查函数、轴、空间、几何、谱序号和端点结构。任务绑定包含 kind、规范原函数、平均零约束和容差；预算与标签不是数学原题身份，允许控制器传入合法的剩余时间。

候选矩阵也会复制成不可变的精确有理数元组。模块用它们重算质量、算子型和算子范数平方三个二次型，进而取得 Rayleigh 商与完整残差。若提供 `statistics`，必须与重算值完全一致。调用者标注的 `verified=True` 或 `certificate_valid=True` 没有认证效力。

这些结构检查仍不能证明候选矩阵确实来自原算子。伪造但对称的矩阵、结构相容的坏 gap，都可能形成候选；外层结果始终是：

```text
certificate_valid = None
target_met = False
verification_pending = True
```

`candidate_target_met` 只表示按候选算术得到的宽度符合目标。内层旧格式的 `status` 为维持完全一致而保留，也必须等原验证器通过后才有证据意义。输入/输出复制、冻结数据类及 owner 身份只防普通误用和别名突变，不构成对恶意同进程 Python 反射代码的安全隔离。

## 数学与可节省的工作

构造使用原有严格条件 `mu < beta` 和非负完整残差平方，按 Temple 不等式写出 `lower = mu - residual_squared/(beta-mu)`、`upper = mu`。矩阵摘要、gap 摘要、全空间 Fourier 支配说明、强算子域和其他所有证书字段与旧 `fullspace_singular.certificate` 完全一致。摘要用于完整性比较，不证明数学正确。

本轮只消除构造阶段的重复劳动：复用搜索已产生的 M/H/R 时不再建立 Gram 矩阵；不再重复核验搜索阶段已经用于控制搜索的 gap。三个精确二次型仍重算一次，最后原验证器也仍完整重建。若旧有理数子进程没有返回矩阵，组装器需重建一次 Gram，只省去本地 gap 的重复重放。不改变候选算法、不缩减最终证明义务；实际速度收益必须由统一计时实验评估，不能由调用次数直接声称。

## 已完成的功能检查

只使用已公开的 E09、6 项试函数和既有证书，没有运行新问题搜索或性能比较。Python 3.12、`-B` 的 15 项单元检查全部通过；该次工具返回约 0.7 秒的测试运行记录，仅作为功能测试结果，不是算法速度对照。

检查覆盖：有/无 prepared 两条路径与旧构造及已有证书完全一致；序列化后原核验成功；prepared 路径不重建 Gram 且不调用本地 gap 验证；无 prepared 路径恰好建立一次 Gram；原题/轴/空间/容差错配拒绝；输入、返回值与矩阵突变隔离；跨上下文 receipt 和错误统计拒绝；伪矩阵、伪 verified gap 只产生 pending 并被原验证器拒绝；最终端点、摘要、容差篡改被拒；非对称、维数不符、浮点矩阵被拒；`mu=beta` 不获放行。

独立审查者另以同一公开证书检查了整体伪 M/H/R、坏 gap、突变及跨任务控制，无阻断。审查结论与其他 C3 模块的边界说明记录于同目录 `C3_REVIEW.md`。本模块不维护跨任务数学缓存，不修改任何冻结源码。
