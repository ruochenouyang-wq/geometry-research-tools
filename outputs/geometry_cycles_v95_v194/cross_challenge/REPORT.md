# 冻结 20 题的独立交叉验收

题库 SHA-256：`e583a5bf8d028b968b832c3b71c8fd8580ad09a27e5d3994b01f0ab3ec1c3e05`。原题、坐标、阈值、规定 n 和首轮预算均保留；每题进程预算 90 秒。

所有保存的数学 API 证书均在落盘前调用相应验证器。独立校准/范围解释另列，未伪装成模块支持。

| 题号 | 原题 | 结果 | 实际区间宽度或说明 |
|---|---|---|---|
| C01 | 零势显式违反者 | solved | 见逐题记录 |
| C02 | 椭球势的严格有理违反者 | solved | 见逐题记录 |
| C03 | 固定m=1势的二十四位证书宽度 | narrow_target_met | 2.5657e-25 |
| C04 | 强势全角向二十位闭合 | narrow_target_met | 7.6254e-25 |
| C05 | 八次非负势 | narrow_target_met | 2.8096e-12 |
| C06 | GN六次试探平方产生十二次势 | narrow_target_met | 5.8467e-12 |
| C07 | 整个线性参数区间的统一正下界 | solved | 见逐题记录 |
| C08 | 非负参数盒与端点精确达到 | solved | 见逐题记录 |
| C09 | 带惩罚的参数全局最优 | narrow_target_met | 9.7505e-09 |
| C10 | 正负椭球势的角分支竞争 | solved | 见逐题记录 |
| C11 | 删除完整l=0和l=1空间 | solved | 见逐题记录 |
| C12 | 只删除常数和z仍保留横向模 | partially_supported | No full-space assembly API for arbitrary constraints distributed across m sectors. |
| C13 | 完整高阶谱与重数 | solved | 见逐题记录 |
| C14 | 负谱个数及负部迹 | partially_supported | 见逐题记录 |
| C15 | GN候选常数的直接反例 | solved | 见逐题记录 |
| C16 | 明确集中序列的有理矩与GN试探 | partially_supported | 见逐题记录 |
| C17 | 旋转等价的非轴输入q=x | narrow_target_met | 3.4355e-15 |
| C18 | 真正三轴二次势 | certified_open | 2.3617e-06 |
| C19 | 同一谱证书支持两个后续判定 | solved | 见逐题记录 |
| C20 | 开放证书的预算续算 | narrow_target_met | 7.6254e-25 |

分类计数：`{"solved": 9, "narrow_target_met": 7, "certified_open": 1, "partially_supported": 3, "unsupported": 0}`。

C02 严格使用原题 u=x(1-z²/10)：概率质量 1121/3500、能量 6019/7875、Rayleigh 商 24076/10089，阈值乘质量减能量为 172/39375>0。没有换用较优但不同的试探函数。

C09 覆盖完整 [-4,4] 后达到原 1e-8 全局 gap；全部极小点的已证包围为 [-1/512,1/512]，不声称唯一性或精确极小参数。C18 保留至当前支持上限 L=4，宽度约 2.362e-6，仍大于原 1e-8 目标。

C12 保留单扇区额外约束与完整球面约束的区别；C14 无投影负谱接口仍不支持；C16 的 n=64 原样提交并保留拒绝，未替换为较小 n。

C19 两个下游判定在禁止谱求解函数的环境中重放同一内容哈希证书。C20 复用的是已验证上下界和证明祖先，不是迭代向量；目标是否闭合独立记录。

状态含义：solved 完成无窄区间目标的原题；narrow_target_met 达到原数值宽度；certified_open 保持完整 scope 的有效宽界；partially_supported 原题仅部分接口/子项受支持；unsupported 无受支持输出。有效但宽的界没有计为闭合。

完整原输入、各次预算、负测试、异常、证书路径与哈希见 result.json 和 cases/*.json。本报告不是新数学定理或已穷尽性能评估。

新进程重放：原题绑定 20/20；证书文件哈希及数学验证 56/56 通过，见 replay.json。首次 C12 通用二分输出被 runner 误要求精确等于 6 的断言记录保留在 negative/；随后用已有有理惯性构造器直接验证端点 [6,6]，未修改题目或模块。
