# Geometry Research Tools

用于几何分析研究的可复核计算原型：函数逼近、球面谱界、误差估计和结构适配的试探函数。项目保存算法、完整证据、失败案例及历次研究记录。

本项目基于已有数学方法，当前贡献主要是针对具体问题的实现、认证流程与可复现实例；尚未证明核心数学方法的学术原创性。发布定位与逐项来源见 [原创性检查报告](outputs/geometry_originality_audit/REPORT.md)。

## 当前入口

- [最新精度研究报告](outputs/geometry_precision_followup/REPORT.md)
- [运行与输入说明](outputs/geometry_precision_followup/TOOL_GUIDE.md)
- [独立数学审查](outputs/geometry_precision_followup/DIRECT_MOMENT_REVIEW.md)
- [前一轮函数一般化](outputs/geometry_general_v235_v237/REPORT.md)

最新高精度示例针对单位球面上的负奇异势 `q(z) = -|z|^(-1/4)`。默认目标为完整零均值空间最低谱值的区间宽度不超过 `1e-8`，结构适配的 10 项分数幂试探已达到该目标。具体区间与限制以报告及证据为准。

## 快速运行

最新精度模块使用 Python 3 标准库。请保留整个 `outputs/` 目录：最新算法会以只读方式复用旧版本中的精确算术和验证程序。

将以下内容保存为 `request.json`：

```json
{"op":"precise_singular","function":{"kind":"axis_profile","profile":"abs_power","exponent":"-1/4","amplitude":"-1"},"tolerance":"1/100000000"}
```

然后运行：

```sh
python3 -B outputs/geometry_precision_followup/precision_tool.py < request.json > response.json
python3 -B -m unittest discover -s outputs/geometry_precision_followup -p 'test_*.py' -v
python3 -B outputs/geometry_precision_followup/verify_release.py
```

`response.json` 包含原函数谱证书和逐预算尝试。`certified_open` 表示界有效但未达到目标；`target_met` 表示证书通过所指定的宽度目标。

## 文件结构

- `outputs/geometry_precision_followup/`：最新精度工具、函数逼近、结构适配试探和全部复核证据。
- `outputs/geometry_general_v235_v237/`：非多项式函数误差向完整谱界的传递。
- 其他 `outputs/` 子目录：历史算法、报告与实验。
- `work/`：早期工作记录和本地实验。部分旧 token 实验使用额外依赖，按其各自文档运行。

历史版本保持不变，以便核对已保存证据。该项目是数学研究原型，不宣称解决一般未解问题、实现任意函数的高精度谱求解或证明不同大模型之间的速度优势。精确有理数软件复核也不等同于证明助手中的形式化证明。

## 方法来源

奇点适配谱基参见 [Li–Zhang](https://arxiv.org/abs/1606.06388) 和 [Yang 等的 Müntz 球多项式](https://arxiv.org/abs/2303.05020)；谱尾控制参见 [Dusson–Sigal–Stamm](https://arxiv.org/abs/2008.10871)；残差下界采用 [Teschl 书稿中的 Temple 不等式](https://www.mat.univie.ac.at/~gerald/ftp/book-schroe/schroe.pdf#page=132)。本项目与这些工作的具体区别、更多方法来源及尚未完成的比较，见原创性检查报告。
