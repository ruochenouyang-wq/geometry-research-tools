# 几何研究工具 V95–V194：调用说明

这一轮把题目中暴露的弱点变成十组可执行算法。工具的主要对象是球面上的谱、不等式、参数优化和试探函数。它需要 Python 3.9 或更新版本及标准库；不需要 API 密钥、安装软件或联网。

在本目录运行：

```sh
python3 -B portal.py
```

逐行输入 JSON，逐行得到答复。系数按常数项、一次项、二次项排列；精确小数请改写为分数，例如 `"12/5"`。不要传入浮点数。返回的长证书存入本目录的 `service_certificates/`；默认答复只返回结论、范围和证书编号。

## 从一次计算得到多个结论

```json
{"op":"ground","q":[0,0,2],"tolerance":"1/1000000000000000000"}
```

这求的是完整零均值空间中的最低能量，包含全部角向和径向尾部。获得 `certificate_id` 后，把下面的占位符替换成实际编号：

```json
{"op":"threshold","certificate_id":"实际64位编号","q":[0,0,2],"threshold":"119/50"}
{"op":"threshold","certificate_id":"实际64位编号","q":[0,0,2],"threshold":"12/5"}
```

两次判定重放同一份证据，不重新求谱。若需要真正的反例函数，再运行：

```json
{"op":"counterexample","q":[0,0,2],"threshold":"12/5"}
```

`threshold` 的反驳可能仅证明反例存在；`counterexample` 会另外验证具体函数。`fetch` 可取得完整函数系数、积分和残差：

```json
{"op":"fetch","certificate_id":"实际64位编号"}
{"op":"verify","certificate_id":"实际64位编号"}
```

## 常用题型

| 入口 | 实际任务 | 范围 |
|---|---|---|
| `ground` | 自适应精度，Schur 与 Temple 区间相交 | 次数≤6的轴对称势，完整球面 |
| `wide_ground` | 较高次势的完整谱底 | 次数≤24；有限预算可返回宽区间 |
| `exponential` | 原始势 `exp(a*t)` 的谱底 | Taylor 统一余项也进入结论 |
| `family` | 整个参数盒上的不等式 | 轴对称仿射势族，≤3参数 |
| `global_minimum` | 参数全局最小值及所有可能最优点区间 | 一参数谱底加次数≤4的多项式惩罚 |
| `constraints` | 排除低次球谐或加入有限矩约束 | 完整次数排除可覆盖全球面；一般矩限单个实角分量 |
| `ordered_spectrum` | 前 k 阶谱、重数、严格计数及负谱和 | 完整零均值空间 |
| `gn_ratio`、`gn_refine` | GN 商的精确计算和试探改进 | 有限次数轴对称试探；不是无限维全局最优证明 |
| `gn_diagnostic` | 将 GN 试探转为高次势并诊断 | 只给固定势必要性诊断 |
| `anisotropic`、`rotation` | 非轴对称势的完整耦合或旋转转移 | xyz 多项式次数≤4；保留空间预算很小 |
| `radius` | 非单位半径的谱、不等式及半径区间搜索 | 精确区分归一化坐标和物理坐标 |

例如，对所有 `s∈[-1,1]` 检验 `q_s=2t²+s*t` 的统一阈值：

```json
{"op":"family","q0":[0,0,2],"directions":[[0,1]],"box":[[-1,1]],"threshold":"23/10"}
```

带惩罚的全局搜索：

```json
{"op":"global_minimum","objective":{"q0":[0,0,2],"direction":[0,1],"domain":[-2,2],"penalty":[0,0,"1/10"]},"max_leaves":32,"tolerance":"1/10000"}
```

零势的完整前九阶谱：

```json
{"op":"ordered_spectrum","q":[0],"k":9}
```

非轴对称势 `q=x*y`：

```json
{"op":"anisotropic","q":{"1,1,0":"1"},"L":3,"bits":28}
```

半径为 2 的球面，零势的精确 Poincaré 阈值为 1/2：

```json
{"op":"radius","q":[0],"radius":2,"threshold":"1/2"}
```

若势写成物理坐标 `q(x3)=x3²`，请明确使用 `"coordinates":"ambient_x3"`。默认 `normalized_t` 表示 `t=x3/R`。二者在改变半径时是不同问题。

GN 试探 `u=P1+(2/5)P3=(2/5)t+t³`：

```json
{"op":"gn_refine","raw":[0,"2/5",0,1],"iterations":2,"degree_budget":12}
```

预算用尽时，用原精度证书继续：

```json
{"op":"resume","certificate_id":"原ground返回的编号","q":[0,0,100],"tolerance":"1/100000000000000000000","max_modes":32,"max_m":16,"max_steps":24}
```

这会保留以前已经认证的上下界。当前实现重新安排扇区求解，并不复用所有特征向量或矩阵分解；不能据此宣称实现了热启动速度提升。

## 怎样解释结果

- `proved`：记录所声明的空间、势和参数范围内，证据支持不等式。
- `refuted...`：检查是否包含具体函数，还是只通过谱界证明反例存在。
- `open`、`undetermined`、`budget...`：仍有可量化的间隙。有效证书不意味着达到所请求精度。
- GN 试探商提高，意味着最佳常数的已知下界提高；有限平面上的上界不能提升为整个函数空间的上界。

这些是按明确分析规则核验的有理数软件证书，尚未导入 Lean 等形式化证明助手。接口不调用语言模型，因而没有真实模型 token 消耗或 5.6/Astra 对照结论。

## 复现和目录

`round01` 到 `round10` 各保存原题弱点、十项能力、运行证据及独立审查。`PROGRESS.json` 只在主任务验收后增加完成数。`challenges.json` 是冻结的二十题压力清单，包含精确校准和自构造题；不是一份声称仍未解决的开放问题清单。

运行整套新增检查：

```sh
python3 -B -m unittest discover -p 'test_*.py'
python3 -B -m unittest round01.test_witness_review
```

原 V90、V94 是只读依赖；请保持它们与本目录的相对位置。
