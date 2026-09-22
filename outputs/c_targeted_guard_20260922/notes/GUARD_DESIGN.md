# C3.1：保留已验证 open 证书的一次后备

入口保持 `solve(raw_task, progress=None)`。首先完整调用冻结 C3 的 controller.solve；从包装入口开始计算共同 wall/CPU，模块加载和后备启动也消耗原任务预算。仅原题为正幅 abs_power、mean_zero=True 的谱问题，且 C3 返回 certificate_valid=True、target_met=False 的证书时考虑后备。没有题号、特定幅度或指数分支，不改变原容差。

固定保留 1.1 秒供父进程最终原 assess 和输出；除这份余量外至少还有 0.2 秒才启动一次后备。子进程调用冻结 baseline.solve，收到完全相同原题、仅减少 wall_seconds；启动、输入输出 JSON、失败调用均在这个子期限内。旧响应先完整 JSON 序列化到父进程，再由父进程原 baseline.assess 绑定原题重放。只有证书有效且 exact_width 严格更小才替换；相等、较宽、损坏、异常、缺证据或超时都保留原 C3 已验证证书，绝不重复尝试。

加载器使用明确绝对源码路径和私有模块名；C3 的 runtime 及本地依赖仅在受控导入作用域绑定到冻结 C3 目录，并检查 baseline/solver 实际路径，防止误用同名 runtime。没有修改任何冻结源码或校验函数。

该导入作用域临时绑定 sys.modules/sys.path，仅承诺同进程顺序调用；未承诺多线程并发安全。数学验证函数从未替换。完整后备响应保留在 guard.fallback_worker.response，可能增加最终输出体积；其 JSON 传输/检查仍计费，不为了压缩日志丢弃审计证据。后备响应设置 1600 万字符解析上限，过大的响应作为失败保留旧证书。

本后备明确限定冻结的单进程正幅 mean-zero abs_power 路线。已只读核对 baseline.solve→solver._route→singular_solver.solve→spectral_gap、adaptive_singular、direct_moments/enriched_trial 的数值调用链，没有创建子孙进程。后备继承外层评估器进程组，不 start_new_session；局部超时 terminate/kill 并回收这个直接子进程，外层 killpg 仍覆盖它。若未来更换为会创建子孙进程的后备，本实现的局部终止假设须重新审查，不能静默扩展。

共同 total_elapsed_seconds 从入口计到序列化检查结束。CPU 汇总父进程 process_time 与已回收子进程 rusage 增量，涵盖本固定单进程后备；分别列出分量，不把嵌套阶段重复相加。父进程原 assess 用 Unix 主线程 SIGALRM 定时作用域保护，自定义 BaseException 避免被旧验证器的普通错误处理吞掉；截止设为总 deadline 前 0.25 秒，保留最终序列化/输出空间。超时不替换并返回原最佳证书，成功重放后也检查尚在该截止前；原数学函数未改。已有活动进程定时器时不夺用它，放弃后备验收并保留原证书。外层整组硬截止仍必需。超出总预算时仍保留已验证 best，但不得标记预算内达标。actual_model_calls=0，actual_model_tokens=None，未进行真实模型 token 计量。

测试主要使用 mock 验证无改善、伪证据、异常、超时、预算不足和已有成功时保留 best；另只运行已公开 N29 与 E18 各一次原容差、10 秒功能回归，保留成败，不按结果继续调参。

本修正和定向回归发生在 N29 原结果公开之后，标记为 post-reveal；不能计入原 180 次固定正式评估，也不是新的盲测。

## 开发验证记录（保留首次成败）

bundled Python 3.12.14 首次执行 14 项检查时，两个既定真实回归都通过：N29 保持有效 open，宽度缩至 0.0005088078846477835，wall/CPU 为 5.412358/5.324765 秒（父 2.824629、子 2.500136 秒），后备恰一次；E18 达标，宽度约 7.242073270893462e-7，wall/CPU 为 1.868641/1.855805 秒，没有后备。这些是单次 post-reveal 功能记录，不是正式提速结果。

该首次套件有一项计费断言失败：不变的子 CPU 经 `a+b-c-d` 差分产生约 -1.7e-18 的舍入误差，导致汇总比父 CPU 小这个微量；已改为分别求 user/system 增量并把负舍入钳到零。这不是证书或数学失败，原始失败未隐藏。

独审随后要求把父原 assess 也限制在输出预留之前，已按上述 SIGALRM 作用域修正，并添加实际慢 mock 验证被中断、原证书保留、原 handler 恢复且无悬挂 timer 的检查。仅重跑 `test_controller.py GuardControls`，14 项合成/加载/清理检查全部通过（0.311 秒）；没有重复真实数学回归或按结果调参。最终新增截止保护交由独立审查及根固定补充评估复核。
