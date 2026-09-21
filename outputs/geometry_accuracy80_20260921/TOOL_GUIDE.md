# 几何精度工具使用指南

本目录的统一入口是 `solver.py`。它接收一个明确的数学任务，尝试求解，并返回可以重新检查的完整证书。本文先说明怎样使用，再说明结果的边界。已完成的本轮评测结果见 [研发与验收报告](REPORT.md)。

## 1. 工具用来做什么

目前的几何对象是**单位球面**，即半径为 1 的球的表面。输入函数只随某个坐标轴上的高度变化；绕该轴旋转时函数不变，这称为**轴对称**。

| 目的 | 输入 `kind` | 要求达到的精度 |
|---|---|---|
| 用分段多项式表示给定函数 | `approximation` | 原函数与近似函数的 L² 误差上界 |
| 估计给定势函数下的最低谱值 | `spectrum` | 原题最低谱值所在区间的宽度 |

**分段多项式**指把范围分成若干段，每段使用一个有限次多项式。**L² 误差**是两函数之差平方在球面上的面积平均值，再开平方；它衡量整体误差，不保证每一点都同样准确。

**势函数**是能量公式中给定的函数。这里的**最低谱值**可理解为：在题目允许的函数中，先把函数大小统一，再寻找能量的最低水平。工具给出它的下界与上界，组成**谱区间**；区间宽度等于上界减去下界。

这两个任务的误差不同。函数逼近误差很小，不能直接推断谱区间也很窄。

## 2. 输入中的函数、精度和空间

每条任务是一个 **JSON** 对象。JSON 是用字段名和值记录数据的文本格式。统一入口接受以下字段：

| 字段 | 作用 |
|---|---|
| `kind` | 必须为 `spectrum` 或 `approximation` |
| `function` | 原函数的精确描述 |
| `tolerance` | 误差上界或区间宽度的目标，默认 `"1/100000000"` |
| `mean_zero` | 仅用于谱任务，选择是否要求参与求解的函数平均为零；默认 `true` |
| `budget` | 本入口只接受 `wall_seconds`，默认 10 秒 |
| `id` | 可选标签，不参与选择数学方法，也不改变验收条件 |

**有理数**是可以写成两个整数之比的数。数学参数应使用整数或精确分数字符串，例如 `"-2/5"`、`"-1/3"`、`"1/1000000"`。不要把目标写为 `0.000001` 或 `"1e-6"`。统一入口的 `tolerance` 范围为 `1e-30` 到 `1`。

### 原函数 `function`

`kind="axis_profile"` 表示沿坐标轴变化的函数。`axis` 为 `0`、`1`、`2` 时分别使用 x、y、z 坐标；默认是 `2`。把选定轴的高度记为 \(t\)，它在 `[-1,1]` 内。

- `profile="abs_power"` 表示 \(q(t)=a|t|^\alpha+b\)。`amplitude` 是 \(a\)，`exponent` 是 \(\alpha\)，`offset` 是 \(b\)。指数必须明确给出。
- `profile="step"` 表示阶跃函数：\(t\le0\) 时为 `offset`，\(t>0\) 时为 `offset+amplitude`。阶跃输入没有 `exponent` 字段。

幅度 `amplitude` 默认为 `"1"`，常数偏移 `offset` 默认为 `"0"`。

当前主要范围是单位球面上的轴对称负幂

\[
q(t)=a|t|^\alpha+b,\qquad -\tfrac12<\alpha<0,
\]

以及轴对称阶跃势。**负幂**指指数为负，函数在高度趋向零时可能无界。谱计算还受具体路线的数学条件与资源限制约束；例如完整空间的负幂谱路线要求 `amplitude<0`，较强奇异势的替代证明还有指数分母等限制。不能把范围条件理解为其中所有输入、所有精度都保证成功。

### `mean_zero` 为什么重要

**函数空间**指题目允许哪些函数参加求解。本工具的谱任务包含有限能量条件，也可增加平均值条件：

- `mean_zero=false`：完整空间，允许常数函数。
- `mean_zero=true`：平均为零的空间，排除非零常数函数。

这是两个不同的数学问题，不能为了得到较好结果而擅自更改。轴对称描述的是输入势函数，最终谱证书仍须覆盖所选的完整允许空间。函数逼近任务不需要这个谱约束，因此不得填写 `mean_zero`。

## 3. 一个标准输入示例

以下示例仅使用已经公开的开发题 **H06**：

\[
q(z)=-\tfrac13|z|^{-2/5},
\]

选择完整空间 `mean_zero=false`，要求谱区间宽度不超过 `1e-6`。它已经参与方法开发，不能作为新的独立测试题。

**标准输入**是程序从终端或文件读取文字的通道。这里使用 **JSONL**：每行写一个完整的 JSON 对象，程序逐行返回结果。下面命令可在终端直接复制；请求对象须保持一整行。

使用本次固定的 Python 3.12.14 运行时：

```sh
cd /PROJECT

/PYTHON/bin/python3 -B outputs/geometry_accuracy80_20260921/solver.py <<'JSONL' > outputs/geometry_accuracy80_20260921/h06_example_result.jsonl
{"kind":"spectrum","function":{"kind":"axis_profile","profile":"abs_power","axis":2,"exponent":"-2/5","amplitude":"-1/3","offset":"0"},"mean_zero":false,"tolerance":"1/1000000","budget":{"wall_seconds":10}}
JSONL
```

结果保存在 `outputs/geometry_accuracy80_20260921/h06_example_result.jsonl`，里面包含完整证书。把命令末尾的输出重定向去掉，就会直接在终端显示结果。项目搬到其他目录后，应修改前面的目录路径；固定 Python 路径也应指向实际安装的同版本运行时。

本入口直接接收任务对象并输出完整结果。不要额外添加本接口未列出的搜索参数；更底层的实验接口见文末方法文档。

## 4. 怎样判断结果

**证书**是保存原函数、数学界限及检查所需计算的完整数据。**验证器**是重新计算并核对这些数据的程序。工具不能仅凭自己写了“成功”就通过检查。

先看以下字段：

| 字段 | 解释 |
|---|---|
| `certificate` | 完整证书；若为空，当前没有可使用的证书 |
| `certificate_valid` | 证书是否通过原题绑定和数学检查 |
| `metric` | 谱任务为精确区间宽度，函数逼近任务为 L² 误差上界 |
| `target_met` | 本次统一求解结果是否通过精度验收；超过所报预算时会强制改为 `false` |
| `within_budget` | 本次 `solve` 调用是否在其协作时间预算内返回 |
| `status` | 完成、未达标、超时或失败的状态 |
| `attempts` | 已尝试路线及其细节；其中可能保留有效但未达标的阶段 |
| `total_elapsed_seconds` | 此次 `solve` 调用的实测时间，不等于完整进程或整个研究的总成本 |

常见状态如下：

- `target_met`：证书通过，精度达标，并未被统一入口判为预算超限。
- `certified_open`：已有有效数学界，但误差或区间宽度仍超过目标。
- `no_verified_certificate`：这次没有取得可验证证书。
- `verification_failed`：证书没有通过检查，不能用来声称原题成立。
- `execution_failed`：路线执行出现异常，详情在 `error` 中。
- `budget_exceeded`：时间预算已用尽或最终调用超时；可能仍保留一份有效证书，但本次不计作按预算完成。

输入本身非法时，命令行入口还可能只返回 `target_met=false` 和 `error`，此时不能假定其他字段一定存在。无论哪一种失败，都不意味着已经证明原数学问题无解。

谱证书中的 `lower`、`upper` 和 `exact_width` 使用精确分数。显示成普通小数后，两个很接近的端点可能看起来相同；应以精确区间宽度为准。

### 时间预算的边界

本入口的 `wall_seconds` 是**协作预算**：代码在约定位置检查时间，不能强行中断一个正在运行的同步计算。`wall_seconds=0` 会立即返回预算超限，不启动求解。

`total_elapsed_seconds` 包含路线执行、证书转换和统一入口的复核，但不包含 Python 进程启动及终端读写等全部成本。最终正式评测还需要由外部程序设置**硬期限**，即到时间便终止整个被测进程，并将证书序列化、独立重放等工作计入约定预算。因此，单看一次命令返回的时间，不能代替正式评测。

## 5. 不重新搜索，独立重放完整证书

**重放**指重新读取证书并检查数学计算，而不是重新寻找候选解。下面命令读取上一例的文件，重新固定 H06 原函数、完整空间和原精度，再检查证书。它不会调用 `solve` 去重新搜索。

```sh
/PYTHON/bin/python3 -B - <<'PY'
import json
import sys
from pathlib import Path

root = Path("outputs/geometry_accuracy80_20260921").resolve()
sys.path.insert(0, str(root))
import solver

task = solver.normalize_task({
    "kind": "spectrum",
    "function": {
        "kind": "axis_profile", "profile": "abs_power", "axis": 2,
        "exponent": "-2/5", "amplitude": "-1/3", "offset": "0"
    },
    "mean_zero": False,
    "tolerance": "1/1000000",
    "budget": {"wall_seconds": 10}
})
result = json.loads((root / "h06_example_result.jsonl").read_text())
certificate = result.get("certificate")
if certificate is None:
    raise SystemExit("上次运行没有返回证书，请先检查原始失败状态。")

check = solver.assess(certificate, task)
print(json.dumps({
    "certificate_replay": check,
    "original_run_status": result.get("status"),
    "original_within_budget": result.get("within_budget")
}, ensure_ascii=False, indent=2))
PY
```

这里 `solver.assess` 重新调用对应数学验证器，并单独计算精度是否达到。它不会重新判断上次运行的时限；例如超时后留下的有效证书，数学重放仍可能达标，但这不会把原先的超时改成预算内成功。

需要布尔型证书检查时，也可以直接使用：

```python
valid = solver.verify_certificate(
    certificate,
    expected_function=task["function"],
    expected_mean_zero=task["mean_zero"],
    expected_tolerance=task["tolerance"],
    expected_kind=task["kind"]
)
```

对函数逼近证书，`expected_mean_zero` 应为 `None`。返回 `True` 仅表示证书及其原题绑定有效；仍须比较目标误差。不能改变函数、空间或目标精度后，继续拿原证书的状态文字当作新题的证明。

“独立重放”在这里指不依赖候选搜索过程进行核验；各模块仍有共享的精确计算与数学前提。它不等于另一套完全独立实现，更不等于已由形式化证明助手检查。

## 6. 测试与方法说明

**测试**用于检查实现行为和已覆盖的算例，并不能证明对所有输入都正确。可在项目根目录运行本目录测试：

```sh
/PYTHON/bin/python3 -B -m unittest discover -s outputs/geometry_accuracy80_20260921/tests -p 'test_*.py' -v
```

本指南编写时只检查示例的文字与语法，没有重新运行上述数学求解，也没有生成新评测任务。

进一步阅读：

- [一般有理负幂的函数逼近](GENERAL_APPROX.md)：怎样保留近奇点区域的完整误差。
- [完整空间奇异谱方法](FULLSPACE_METHOD.md)：H06 所用路线及其条件。
- [平均零空间与谱分离](SPECTRAL_GAP.md)：怎样证明候选与其他谱值之间有足够间隔。
- [阶跃势方法](STEP_METHOD.md)：怎样处理赤道两侧的跳跃和完整残差。

保留本项目原有目录结构。新实现仍使用已有版本中的部分数学代码，仅复制 `solver.py` 不能组成可运行工具。

## 7. 80% 目标意味着什么

本轮的 80% 是**固定数学题族中新参数工具任务**的通过率目标。按当前计划，预先固定 20 个任务，至少 16 个满足既定精度、证据复核和重复预算条件才达到目标。代码冻结后才用于正式评测的新参数，与已经公开、用于改进方法的 H06 等开发题分开统计。

这不是“80% 的几何开放问题已经解决”，也不是对任意几何、任意函数的成功保证。**开放问题**指尚无已知完整解答的数学问题；本工具做的是与研究有关的、明确限定的计算和认证子任务。有限参数点的成功不能自动证明整个连续参数范围，更不能替代开放猜想的证明。

最终通过率与总成本应以冻结后正式评测记录为准。本指南不预先填入成功率。

**token** 是模型处理文字时的编码单位。当前没有实际模型调用用量或跨模型收益的实测证据；更短输出、本地运行变快或证书变小，都不能单独证明模型更省 token，也不能证明较低成本模型已经超过更强模型。整个工具目前仍是可复核的专用数学研究原型。
