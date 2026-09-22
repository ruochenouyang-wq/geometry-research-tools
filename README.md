# Geometry Research Tools

用于几何分析研究的可复核计算原型：函数逼近、球面谱界、误差估计，以及辅助模型组织研究的工具设计。

长期目标是让成本较低的模型借助可靠工具完成更困难的数学任务。当前已验证的是限定数学题族中的本地算法能力，尚未验证跨模型优势或实际模型 token 节省。

## 从哪里开始

- [最新研究方向报告](outputs/research_directions_20260921/研究方向探索报告.md)：九条研发方向、参考项目、最小实验和转向条件。
- [长期研究指导报告](outputs/geometry_research_strategy/研究指导报告.md)：项目目标、整体框架与评估思路。
- [最新 C3 / C3.1 迭代报告](outputs/c_targeted_iteration_20260922/REPORT.md)：180 次完整比较、48 次定向回归，以及精度退步和恢复成本。
- [C3.1 使用说明](outputs/c_targeted_guard_20260922/README.md)：最新入口、条件与预算边界。
- [与外部方法的比较](outputs/c_external_comparison_20260922/REPORT.md)：已有方法、组件实验与证据限制。
- [上一轮算法与验收报告](outputs/geometry_accuracy80_20260921/REPORT.md)：先前实质优化及完整结果。
- [输入与证书说明](outputs/geometry_accuracy80_20260921/TOOL_GUIDE.md)：任务格式、适用条件、结果解释。
- [原创性检查](outputs/geometry_originality_audit/REPORT.md)：已有方法来源及尚未确认的原创贡献。
- [发行说明](PUBLISHING.md)：文件范围、历史路径与复现边界。

## 当前能力与证据

主要对象是单位球面上的轴对称负幂函数和阶跃函数。工具可计算带误差上界的函数逼近，或给出最低谱值的严格上下界；支持范围仍受原函数、函数空间及证明路线的条件约束。

2026 年 9 月 22 日的 C2 / C3 完整比较有 28 道数学题、2 道方法边界控制，每版每题重复三次，共 180 次。稳定达标题由 25/28 增至 26/28；双方稳定成功的 25 题累计用时减少 17.9%，但较强负势路线内存中位数增加约 24%。这不是所有题都更快或更精确。

C3.1 随后修复了 N29 的有效区间退步。8 题、48 次定向回归中，N29 三次恢复旧 C2 精度，仍未满足原目标；用时中位数从 3.75 秒增至 7.07 秒。两版均稳定达标 5/7 道数学题，另一题为方法边界控制。此回归不能替代全 30 题比较，N30 仍无证书。完整范围与失败记录见[最新报告](outputs/c_targeted_iteration_20260922/REPORT.md)。

历史上，2026 年 9 月 21 日的另一组冻结评测在 20 道新参数任务上全部稳定达标，旧版为 3/20，每题五次运行；它使用不同题单和协议，不能与本轮分母合并。记录见[历史正式报告](outputs/geometry_accuracy80_20260921/evaluation/batches/B001/formal_report.json)。

这不表示解决了 80% 或 100% 的几何开放问题。当前没有 GPT-5.6 与 Astra 的实际对照数据。证书的软件复核也不等于证明助手中的完整形式化证明。研究方向报告中的下一阶段实验仍是提案。

## 快速运行

最新入口是 C3.1。建议 Python 3.12，要求 Unix 主线程顺序调用；暂不支持 Windows 或同进程多线程并发。NumPy 是可选的快速候选依赖，缺少时保留有理计算后备；正式测量使用 Python 3.12.14、NumPy 2.3.5。算法不需要模型 API 密钥。

在仓库根目录运行：

```sh
python3 -B outputs/c_targeted_guard_20260922/controller.py < examples/c31_e18.jsonl > response.jsonl
```

示例 E18 是已经用于开发的公开任务：平均零函数空间上的势函数 `q(z) = (2/3)|z|^(-2/7)`，目标谱区间宽度为 `1/1000000`。平均零表示函数在整个球面的平均值为零。它不属于新的独立评测题。保留整个 `outputs/` 结构，因为入口复用 C3、C2、C1 及旧数学模块。

结果中的 `certificate` 是完整证书；`certificate_valid` 表示通过检查；`target_met` 表示本次达到精度与入口预算要求。返回失败不等于原数学问题无解。时间预算的具体含义见[使用指南](outputs/geometry_accuracy80_20260921/TOOL_GUIDE.md)。

检查 C3.1 的后备、证书选择与截止保护：

```sh
python3 -B outputs/c_targeted_guard_20260922/test_controller.py GuardControls -v
```

安装 NumPy 后可运行 C3 的组合功能检查：

```sh
python3 -B -m unittest discover -s outputs/c_targeted_iteration_20260922 -p 'test_*.py' -v
```

## 仓库结构

- `outputs/c_targeted_guard_20260922/`：最新 C3.1 入口、48 次定向回归与独立审查。
- `outputs/c_targeted_iteration_20260922/`：冻结 C3、180 次完整比较与总报告。
- `outputs/c_demand_extension_20260922/`、`outputs/geometry_parallel_v1_20260921/`：保留的 C2、C1 与并行方案。
- `outputs/parallel_cross_dimension_20260921/`：早期跨维度比较；大体积历史结果以公开压缩归档保存，解压说明见[本次发行说明](docs/PUBLICATION_20260922.md)。
- `outputs/geometry_accuracy80_20260921/`：基础求解器、数学方法、测试与历史冻结评测证据。
- `outputs/research_directions_20260921/`：多方向研究报告及十一份专题笔记。
- `outputs/geometry_research_strategy/`：长期指导报告。
- 其他 `outputs/` 目录：历次算法、失败记录、审查与完整对照证据。
- `output/pdf/`：已有长期指导报告的 PDF。
- `work/`：早期工作记录；下载的第三方依赖不在发行包内。

请保留整个目录结构。最新实现复用历史模块，只复制一个求解脚本不能运行。

## 方法来源与权利说明

项目基于已有数学方法，当前主要贡献是具体实现、认证流程与可复现实例，尚未证明核心方法的学术原创性。逐项方法来源见[原创性检查](outputs/geometry_originality_audit/REPORT.md)及各模块报告。

本项目按 [MIT License](LICENSE) 开源。第三方材料保留各自版权与许可，见 [第三方说明](THIRD_PARTY_NOTICES.md)。公开副本已将历史本机路径替换为中性占位符；转换与证据边界见 [发行说明](PUBLISHING.md)。
