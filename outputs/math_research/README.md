# GPT 数学研究原型：非线性系统的稳定性

这是一个可运行的研究算法原型，目标是为多项式微分方程寻找**全局渐近稳定性证明**。GPT 负责选择候选函数的表示，数值搜索负责寻找系数，精确验证器负责接受或拒绝证明。

核心能力是**在表示失败后换一种数学表示**：从二次到高次多项式，或从多项式切换到对数函数。支持 2–3 个变量，输入为有理系数多项式向量场。

当前完成的是算法与基准实验。没有声称解决公开未解问题，也没有测出 GPT 相比传统搜索的性能优势。现有的对数表示来自已知数学方法；它的意义是检验算法能否跳出多项式搜索的限制。

## 先运行

在本文件所在目录执行：

```sh
python3 run.py demo
```

会运行三个基准，结果写入 `results/`：

| 基准 | 结果 | 检验目的 |
| --- | --- | --- |
| 稳定但没有多项式 Lyapunov 函数的文献系统 | 切换到对数函数，证书通过 | 数学表示转换 |
| 无二次 Lyapunov 函数的三次系统 | 找到四次函数，证书通过 | 扩大搜索空间 |
| 原点不稳定的线性系统 | `not_solved` | 防止输出错误证明 |

本机启动器会自动找到已有的 NumPy 运行环境。其他机器需要 Python 3.9+ 和 NumPy：

```sh
python3 -m pip install -r requirements.txt
python3 run.py demo
```

启动器不会自动安装任何包。**精确证书复核不依赖 NumPy**：

```sh
python3 stability_research.py verify results/stability_no_polynomial.json
python3 stability_research.py verify results/stability_quartic.json
```

## 接入 GPT

自动模式调用 OpenAI Responses API，模型需支持结构化输出。API 密钥通过环境变量 `OPENAI_API_KEY` 提供，模型通过 `--model` 或 `OPENAI_MODEL` 指定。密钥不要写进问题文件或研究记录。

```sh
python3 run.py run examples/stability_no_polynomial.json \
  --provider gpt --model "你的可用模型ID" \
  --rounds 6 --out results/gpt_research.json
```

每轮 GPT 收到向量场、之前选择的表示、搜索失败原因及可用的反例，返回一个受约束的计划。程序求解并验证该计划对应的证书；失败会反馈给下一轮。程序不执行模型生成的代码。

当前已实现的计划空间只有两族：2/4/6 次多项式，以及按变量分离的对数与二次项组合。GPT 目前不能任意添加新函数族。其价值是否超过固定顺序搜索，需要后续在独立测试集上比较。

本次环境没有 API 密钥，因此**真实 API 联调未执行**；API 请求格式、响应解析、失败处理已通过模拟响应测试。结构化输出接口依据[官方文档](https://developers.openai.com/api/docs/guides/structured-outputs)实现，实际模型权限与可用性由你的账户决定。

## 不使用 API 密钥

可以先生成研究提示，将其粘贴给 GPT：

```sh
python3 stability_research.py prompt examples/stability_no_polynomial.json
```

把 GPT 返回的计划放进 JSON 列表，保存为一个文件，再运行：

```sh
python3 run.py run examples/stability_no_polynomial.json \
  --plans examples/representation_plan.json \
  --out results/manual_research.json
```

示例计划文件由本次助手编写，用于演示接入格式，并不是一次独立的 GPT 基准测试。需要把失败记录反馈给 GPT 时：

```sh
python3 stability_research.py prompt examples/stability_no_polynomial.json \
  --history results/manual_research.json
```

## 输入自己的系统

问题 JSON 包含 `title`、`dimension`、`vector_field`。第 i 个分量是一组单项式；`powers` 为各变量幂次，`coefficient` 是整数或分数字符串。例如 `{"powers":[1,1],"coefficient":"1/2"}` 表示 `xy/2`。

直接复制 `examples/stability_no_polynomial.json` 修改方程即可。约束如下：变量数为 2 或 3，总次数不超过 8，每个分量最多 150 项，原点必须是平衡点；研究定义域固定为整个实空间。

```sh
python3 run.py run my_problem.json --rounds 6 --iterations 5000 \
  --seconds 20 --out results/my_research.json
```

预算控制是每轮迭代上限和软时间上限。时间检查发生在迭代中，初始化的精确消元和数值分解不受硬超时约束；整个程序不是资源隔离沙箱。搜索变量上限为 600。

## 如何理解输出

- `certified_global_asymptotic_stability`：在内置数学规则及验证代码正确的前提下，证书证明原点全局渐近稳定。
- `not_solved`：当前表示和预算没有找到证书，不能解释为系统不稳定。
- `search_failed`：某一轮搜索失败，不是“不存在该函数”的证明。
- `counterexample`：只反驳记录中那个近似候选的规定裕量，不反驳原系统稳定性。
- `provider_error`：GPT 接入失败，保留已有记录并停止；不自动重复请求。

验证器逐项重算有理系数恒等式，并以精确分数运算检查 Gram 矩阵半正定。数值残差很小不构成证明。当前没有接入 Lean/Isabelle；对数函数的正性、无穷远增长以及 Lyapunov 定理属于内置的数学规则。

## 文件与测试

- `stability_research.py`：主要研究算法、GPT 计划接口、数值搜索和精确检查。
- `ALGORITHM.md`：算法定义、证明条件、复杂度与扩展方向。
- `EXPERIMENTS.md`：实际实验结果及其边界。
- `run.py`：运行入口。
- `math_research.py`：API 公共接口，以及保留的初等求和与固定除数验证器测试原型。
- `test_*.py`：测试；运行 `python3 run.py test`。

本次完整测试共 45 项通过。包含表示切换、三维输入、证书篡改、极小负特征值、反例反馈、API 模拟响应及预算停止。
