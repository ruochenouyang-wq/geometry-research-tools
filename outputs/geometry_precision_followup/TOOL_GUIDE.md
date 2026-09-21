# 精度研究工具

入口是 `precision_tool.py`：每行接收一个 JSON 请求，返回一份 JSON 结果。只使用本地数学程序，无模型调用或网络调用。所有输入系数使用整数或有理数字符串，禁止用浮点数冒充精确数据。

## 构造更贴近的分段函数

原函数为单位球面上的 `q(z) = -|z|^(-1/4)`。以下请求构造分段多项式，并给出整个球面的概率 L² 误差上界：

```json
{"op":"approximate","function":{"kind":"axis_profile","profile":"abs_power","exponent":"-1/4","amplitude":"-1"},"levels":48,"degree":20,"root_ratio":"2/3","sqrt_bits":40}
```

这是**函数逼近证书**。它不包含谱值结论，也没有把高精度分片多项式转换为旧后端能够消费的一条全局多项式。

`levels` 是正半轴几何环片数，另外保留靠近零点的核心片；`degree` 是每个环片的多项式次数，核心片使用常数。另一半轴由偶对称表示。高阶环片只计算一个形状模板，其余环片按精确缩放关系构造。

跳跃函数 `{"kind":"axis_profile","profile":"step"}` 使用两片常数精确表示，L² 误差为零。两种函数都支持有理 `amplitude`、`offset` 和 `axis: 0|1|2`。`abs_power` 在此分片工具中须明确写 `exponent:"-1/4"`；省略时沿用旧工具默认指数 1，并因当前分片方法不支持该指数而拒绝。

## 直接计算原函数的谱界

以下请求直接计算原跳跃函数的矩，认证完整零均值球面空间的最低谱值：

```json
{"op":"spectrum","function":{"kind":"axis_profile","profile":"step"},"modes":12,"near_tail":16,"bits":36,"tolerance":"1/100000000"}
```

将函数替换为上面的负奇异函数即可计算其原谱。`modes` 是每个方位角中的有限球谐个数，`near_tail` 是额外显式计算耦合的尾部个数；再远的无穷尾部仍有整体界。设置 `mean_zero:false` 可以计算包含常数的完整空间。

接口支持 `step` 和 `abs_power` 的精确矩。后者要求指数大于 −1/2，且当前中心化 L² 尾部界满足严格稳定性条件。谱算法可能拒绝某些确实存在谱界的函数，因为这条估计路线不足以认证；这不代表原数学问题不存在解。

默认精度目标始终为 `1e-8`。输出 `certified_open` 说明区间有效但未达目标，不能当成高精度结果。

## 用分数幂试探求高精度原谱

对于明确的原题 `q(z)=-|z|^(-1/4)`，以下请求会先生成覆盖整个零均值球面的粗源证书，再增加分数幂试探项，直到通过精度目标或达到预算：

```json
{"op":"precise_singular","function":{"kind":"axis_profile","profile":"abs_power","exponent":"-1/4","amplitude":"-1"},"tolerance":"1/100000000","max_terms":16}
```

默认目标下实际在 10 项停止。把 `tolerance` 改为 `"1/10000000000000000"`，同一预定基函数序列会计算到 16 项。每次尝试的未达标宽度保留在 `attempts` 中。

这一路径目前只对上述 z 轴、幅度 −1、零偏移、零均值原题认证。候选系数用有理逆迭代选取，再以完整强残差及已经认证的第二径向特征值下界形成 Temple 谱界。程序不会用有限投影残差代替完整残差，也不会把局部渐近展开当作答案。

## 独立于搜索的重放

把上一响应的 `certificate` 原样放入如下请求：

```json
{"op":"verify","certificate":{},"function":{"kind":"axis_profile","profile":"step"}}
```

实际调用时用真实证书替换空对象。谱证书还可以绑定 `mean_zero` 和 `tolerance`；函数逼近证书不能绑定谱空间或谱精度。重放会重建精确积分、误差和必要矩阵，并验证给定端点，不重新搜索最优谱值。

## 运行

在项目根目录中，将一个请求保存为文件后运行：

```sh
python3 -B outputs/geometry_precision_followup/precision_tool.py < request.json > response.json
```

完整证明适合保存在文件里，日常查看误差、谱区间和状态即可。高精度分片证明可能超过 1 MB；这是证据文件大小，不是模型 token 成本测量。
