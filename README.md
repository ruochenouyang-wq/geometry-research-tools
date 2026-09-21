# Geometry Research Tools

用于几何分析研究的可复核计算原型：函数逼近、球面谱界、误差估计，以及辅助模型组织研究的工具设计。

长期目标是让成本较低的模型借助可靠工具完成更困难的数学任务。当前已验证的是限定数学题族中的本地算法能力，尚未验证跨模型优势或实际模型 token 节省。

## 从哪里开始

- [最新研究方向报告](outputs/research_directions_20260921/研究方向探索报告.md)：九条研发方向、参考项目、最小实验和转向条件。
- [长期研究指导报告](outputs/geometry_research_strategy/研究指导报告.md)：项目目标、整体框架与评估思路。
- [最新算法与验收报告](outputs/geometry_accuracy80_20260921/REPORT.md)：三次实质优化及完整结果。
- [输入与证书说明](outputs/geometry_accuracy80_20260921/TOOL_GUIDE.md)：任务格式、适用条件、结果解释。
- [原创性检查](outputs/geometry_originality_audit/REPORT.md)：已有方法来源及尚未确认的原创贡献。
- [发行说明](PUBLISHING.md)：文件范围、历史路径与复现边界。

## 当前能力与证据

主要对象是单位球面上的轴对称负幂函数和阶跃函数。工具可计算带误差上界的函数逼近，或给出最低谱值的严格上下界；支持范围仍受原函数、函数空间及证明路线的条件约束。

2026 年 9 月 21 日的冻结评测中，新版在 20 道新参数任务上全部稳定达标，旧版为 3/20。每题要求五次独立运行均满足精度、证书复核和预算。任务与原始记录见[正式报告](outputs/geometry_accuracy80_20260921/evaluation/batches/B001/formal_report.json)。

这不表示解决了 80% 或 100% 的几何开放问题。当前没有 GPT-5.6 与 Astra 的实际对照数据。证书的软件复核也不等于证明助手中的完整形式化证明。研究方向报告中的下一阶段实验仍是提案。

## 快速运行

最新求解入口使用 Python 标准库；现有评测环境为 Python 3.12.14。建议使用 Python 3.12。早期实验可能另需依赖，见对应目录的说明。

在仓库根目录运行：

```sh
python3 -B outputs/geometry_accuracy80_20260921/solver.py < examples/h06.jsonl > response.jsonl
```

示例 H06 是已经用于开发的公开任务：完整函数空间上的势函数 `q(z) = -|z|^(-2/5)/3`，目标谱区间宽度为 `1/1000000`。它不属于新的独立评测题。

结果中的 `certificate` 是完整证书；`certificate_valid` 表示通过检查；`target_met` 表示本次达到精度与入口预算要求。返回失败不等于原数学问题无解。时间预算的具体含义见[使用指南](outputs/geometry_accuracy80_20260921/TOOL_GUIDE.md)。

检查当前接口的核心行为：

```sh
python3 -B -m unittest discover -s outputs/geometry_accuracy80_20260921/tests -p 'test_solver.py' -v
```

运行当前算法的完整测试集：

```sh
python3 -B -m unittest discover -s outputs/geometry_accuracy80_20260921/tests -p 'test_*.py' -v
```

## 仓库结构

- `outputs/geometry_accuracy80_20260921/`：当前求解器、数学方法、测试与冻结评测证据。
- `outputs/research_directions_20260921/`：多方向研究报告及十一份专题笔记。
- `outputs/geometry_research_strategy/`：长期指导报告。
- 其他 `outputs/` 目录：历次算法、失败记录、审查与完整对照证据。
- `output/pdf/`：已有长期指导报告的 PDF。
- `work/`：早期工作记录；下载的第三方依赖不在发行包内。

请保留整个目录结构。最新实现复用历史模块，只复制一个求解脚本不能运行。

## 方法来源与权利说明

项目基于已有数学方法，当前主要贡献是具体实现、认证流程与可复现实例，尚未证明核心方法的学术原创性。逐项方法来源见[原创性检查](outputs/geometry_originality_audit/REPORT.md)及各模块报告。

本项目按 [MIT License](LICENSE) 开源。第三方材料保留各自版权与许可，见 [第三方说明](THIRD_PARTY_NOTICES.md)。公开副本已将历史本机路径替换为中性占位符；转换与证据边界见 [发行说明](PUBLISHING.md)。
