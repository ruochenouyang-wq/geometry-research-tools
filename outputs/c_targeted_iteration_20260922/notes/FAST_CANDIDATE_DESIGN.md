# C3 快速候选器接口与边界

`try_candidate(function, powers, gap, tolerance, *, prepared=None, full=None, progress=None, deadline=None)` 只在 6/10 项空间尝试 NumPy binary64 Cholesky、白化和最小 Ritz 向量；16 项不尝试，原因是上一固定实验的明确失败证据。`deadline` 是 perf_counter 的绝对截止，阶段边界协作检查，不能中断正在执行的库调用；根仍须保留外部硬期限。

返回 `accepted / fallback_required / reason / candidate / prepared / phases / elapsed_seconds / cpu_seconds`。candidate 含旧格式 powers、精确有理 coefficients、完整精确 statistics、precision_bits=53 和筛选宽度；即使宽度未达标仍保留 candidate 供根记录。只有精确 mu<调用方 beta 且 residual_squared/(beta-mu)<=原 tolerance 才 accepted。这里 beta 只是调用方提供的候选筛选上下文，不宣称已证明谱界；候选器绝不返回最终 certificate_valid。

NumPy 仅在实际尝试时懒加载，单列 import wall/CPU 和进入时是否已加载。导入失败、Cholesky 失败、非有限系数、异常、精确残差无效、mu>=beta、宽度未达标及预算停止均返回明确原因。根随后调用原 C2 有理候选器；本组件不隐藏额外调用、不加正则化、不降低目标、不对外部向量追加残差细化。阶段计时包含失败，完整组件 elapsed 是总费用，不能重复叠加嵌套计时。

`prepare_exact(function, powers, *, full=None, matrices=None, transform=None)` 提供不可变 PreparedTrial，字段为 function_key（`full.canonical(full.normalize(q))`，排序键、紧凑 JSON 编码）、powers、M/H/R/C、mr/hr，矩阵均为 Fraction 元组。默认调用冻结 full.trial_matrices 和 cancellation_transform；可接入同一请求中调用方已有的精确矩阵。复用时检查规范化原函数和 powers 一致，拒绝跨题/跨空间使用；无持久缓存。传入矩阵只作结构和精确类型检查，不把它的来源认证为数学证据，最终旧验收必须独立重建原题矩阵。

认证组件可取 `prepared.M/H/R` 与原系数直接构造旧格式候选证书，省去同次搜索重复准备，但最终 baseline.assess 原样独立重放。快速筛选与证书真实性之间保持这个明确边界。53 位 Ritz 优化与 C2 的 112 位残差策略不同，暖微基准结果不代表包含 import、准备和最终验证的完整求解提速。

测试限于 mock/解析控制，以及已有公开 E09/10 与 E19/6 两个固定功能样例；不跑封存新题，不因功能结果试调参数。

## 首轮功能检查

使用 bundled Python 3.12.14 执行 `python3 -B test_fast_candidate.py`，16 项检查全部通过，总测试 0.643 秒。检查包含：不可变精确准备、跨函数/幂绑定拒绝、错误 gap 空间绑定、16 项廉价跳过、缺失 NumPy、Cholesky 失败、非有限/零向量、负残差、mu=beta（即使残差为零也拒绝）、实际精度不足、零预算不启动、计算统计后截止保留候选，以及无最终证书有效性声明。

预定两个公开功能样例未调参：E09/10 accepted，精确统计与冻结模块重新准备 M/H/R 后的结果逐项相等，旧 full.certificate 和 full.verify 通过；E19/6 返回 exact_width_exceeds_tolerance，保留完整候选/统计，旧验证仍确认其为有效 open 证书。两者精确筛选宽度均与旧证书完全相等。

E09/10 首次 NumPy 导入在组件内记录 wall/CPU 为 0.036423/0.020959 秒，组件合计 0.041698/0.026218 秒；随后 E19/6 已导入，import 约 0.000005 秒，组件合计约 0.001779 秒。这些只是这次功能运行的阶段记录，不是完整求解比较；测试调用方预先加载冻结模块、取得历史公开已证 gap，相关费用不包含在这些组件时间中。控制器必须从完整 solve 入口计费，不能据此宣称端到端提速。

与认证作者已对齐：`prepare_trial(function,powers,*,full=None)` 是不导入 NumPy 的公开精确准备入口，可供旧有理路径复用。PreparedTrial 是 Python 对象，日志 JSON 应剥离它或记录其 function_key/powers；candidate 本身可 JSON 化。认证组件读取 M/H/R、重新做精确二次型并比对统计，最终 baseline.assess 继续独立重建。
