# C3 针对性迭代：冻结对照评估

固定比较冻结 C2 与新 C3，每题三次新 Python 进程，按原题编号在重复间交替先后顺序，全部串行。保留 C2 的 24 条原题对象不变（22 数学任务、2 方法边界控制），另有独立评估代理在求解前封存的 6 条已知数学族邻域题。共 30×2×3=180 次；数学任务总分母28，旧22与新增6分开，边界控制2另列。新增仅为参数迁移检验，不代表新数学族；改变axis属于旋转控制。

完整题单文件 SHA256：`1f286268d9576453b742e4f028991db3f943779f4a30986794c9dac06b61c275`。新增6题规范JSON承诺：`5ba3aca2bfd5d651f22f3b1f40a2aaf2b20b62a55f5741a8a16da6a1a8ba395e`。未对这6题预跑、换题或筛选；源码冻结前未向开发者提供具体参数。共享目录只有纪律隔离，不是权限隔离。正式运行在验证冻结身份之后生成 revealed_neighbors.json，届时参数公开。题单和协议承诺改变均拒绝运行。

每题10秒由独立父进程执行绝对硬截止，包含导入、求解、完整JSON、原始文件刷新和完整重放；同时杀死worker进程组。两版本共用原始冻结 baseline.assess，绑定原kind/function/space/tolerance，再由评估器用精确有理数核对端点宽度或L2误差。不用solver的target_met作裁决。未验证候选保存在原响应中，但不当证书。

记录完整response、certificate、attempts、progress JSONL、失败阶段/原因、真实路线、CPU、wall、峰值RSS、响应与证书UTF-8 JSON bytes。路线区分native_c3、reused_c2、reused_c1、baseline_fallback；复用不自动算错误回退。report.regressions列出稳定成功丢失、逐次成功退步、复用/回退以及共同稳定成功题上的实测成本增加。数学错误/精度不足、明确不支持、其他失败和超时分别保留；不把预期拒绝计作数学成功。

CPU总计来自单个worker的os.wait4，已包含其回收的候选后代，不能重复加children。阶段拆分另给SELF与已回收CHILDREN。RSS是外层与可观测已回收后代的最大单进程高水位，不是同时占用之和；外层硬杀时未回收后代的CPU/RSS可能缺失，记录对应下界/覆盖限制。配对速度只比较同题三次都成功的交集，同时展示真实误差/容差比例。三次重复不提供统计显著性结论。

两个系统显式设置OPENBLAS_NUM_THREADS、OMP_NUM_THREADS、MKL_NUM_THREADS、VECLIB_MAXIMUM_THREADS、NUMEXPR_NUM_THREADS、BLIS_NUM_THREADS为1。Python/NumPy版本、NumPy配置、请求线程限制与可用threadpool信息仅由外部元数据进程在freeze及run前后采样，不向C2 worker额外导入NumPy。当前NumPy为2.3.5；threadpoolctl不可用，因此记录请求限制，不宣称已实测内部线程数。冷启动指新Python进程，系统文件缓存由共享机器管理。

冻结精确保全上一轮619文件并集，再加入本轮实现Python与cases/protocol。evaluation、notes、review、results及缓存目录不进入本轮源码集合，因此自己产生的结果或合成测试文件不会改变源码身份。运行结束重新核验；后置环境探测或manifest异常也保留完整结果并标无效，不丢180行证据。

命令（正式运行仅由获授权的评估代理独占执行）：

```sh
/PYTHON/bin/python3 -B outputs/c_targeted_iteration_20260922/evaluate.py freeze
/PYTHON/bin/python3 -B outputs/c_targeted_iteration_20260922/evaluate.py run
```

每条完成刷新stdout，completed_index.jsonl轻量追加，checkpoint只存进度。完整原始行在各运行目录result.json；末尾一次写report.json（rows、summary、regressions、freeze、身份和环境）。模型API调用为0，实际模型token和缓存token未知，不能以JSON bytes代替。开发自检68项通过，包括硬截止、进程组、CPU/RSS、原题绑定、源码/快照变更拒绝、输出目录排除与重复交替顺序；没有预跑封存题。
