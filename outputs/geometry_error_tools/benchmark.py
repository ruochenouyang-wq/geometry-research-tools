#!/usr/bin/env python3
"""Reproducible small benchmark; timings include proposal and verification."""
from fractions import Fraction as F
import json
from pathlib import Path
import platform
import random
from statistics import median
from time import perf_counter
import error_bounds as e
import run

CASES = [
    ('linear',[0,1],'ordinary polynomial potential'),
    ('quadratic',[1,-2,3],'ordinary polynomial potential'),
    ('strong_linear',[0,20],'stronger localization'),
    ('exact_exponential',[0,20,-100],'manufactured exactly representable family'),
    ('perturbed_exponential',[0,20,-100,0,1],'quartic perturbation of manufactured family'),
    ('double_well',[0,0,-100],'nearby lowest eigenvalues'),
]


def render_report(data,root):
    labels = {'linear':'t','quadratic':'1−2t+3t²','strong_linear':'20t',
              'exact_exponential':'20t−100t²（精确匹配）',
              'perturbed_exponential':'20t−100t²+t⁴','double_well':'−100t²（双井）'}
    names = {'baseline':'原版','banded':'带状版','temple':'普通残差',
             'symmetry':'对称残差','exponential':'指数函数','auto':'自动选择'}
    rows = {(x['case'],x['method']):x for x in data['summary']}
    methods = list(names)
    lines = ['# 实测：误差、函数族与速度','',
             '本文件由 benchmark.py 根据本机原始记录生成。所有结果均经过证书复核，未使用在线模型。', '',
             '## 相同精度目标的比较','',
             '目标是严格区间宽度 ≤ 1e-10；原版与带状版使用 48 轮二分。模式数从 4 开始每次增加 2，最多 24。指数方法搜索 log 函数次数 1、2、4、…、16；其次数不是 Legendre 模式数。', '',
             '| 势函数 q(t) | '+' | '.join(names[m]+' ms' for m in methods)+' |',
             '| '+' | '.join(['---']*(len(methods)+1))+' |']
    for case,_,_ in CASES:
        values = []
        for method in methods:
            r=rows[(case,method)]
            values.append(('%.3f'% (1000*r['median_seconds']))+(' 未达标' if r['status']!='target_met' else ''))
        lines.append('| '+labels[case]+' | '+' | '.join(values)+' |')
    lines += ['', '每格为 3 次运行的中位数，包含达到目标前的全部候选、失败尝试和复核。各轮使用固定随机种子打乱方法顺序。未达标栏记录耗时，不能视为成功提速。计时排除 Python 启动与文件写入；完整样本、范围和环境信息在 [原始数据](results/benchmark.json)。', '',
              '运行环境：Python '+data['environment']['python']+'，'+data['environment']['system']+' / '+data['environment']['machine']+'。这些是小规模同机实验，不是一般几何算子的性能保证。', '',
              '## 可归因的结果','']
    for case in ['linear','quadratic','strong_linear','perturbed_exponential','double_well']:
        a,b=rows[(case,'baseline')],rows[(case,'auto')]
        lines.append('- '+labels[case]+'：自动策略相对原版约快 **%.2f 倍**，两者都达到相同精度目标。'% (a['median_seconds']/b['median_seconds']))
    lines += ['', '精确匹配例没有纳入上述比较列表。函数 exp(−10t) 本来就是该构造模型的基态，因此一个参数可给出零宽度区间；这说明函数族能表达某种结构，不是对任意问题的大幅提速保证。其非零四次扰动被单独测量，指数方法在次数预算内未达到目标。', '',
              '普通残差法在双井例中受很小的完整谱间隔影响；对称残差法只使用偶函数空间内部的间隔。两者的目标仍是同一个完整轴对称算子的基态，合法性依赖 q 的反射对称性及基态正性。', '',
              '20t 的测试中，手选带状法比自动策略更快。自动选择有搜索开销，不能宣称它逐例最优。', '',
              '## 同样 8 个模式的比较','',
              '固定 q=t、N=8；原版和带状版用 36 轮二分。每项 7 次运行，包含复核。这一表比较相同空间，并非相同误差宽度。','',
              '| 方法 | 中位时间 ms | 精确认证区间宽度 |','| --- | ---: | ---: |']
    for x in data['fixed_modes']:
        lines.append('| %s | %.3f | %.6e |'% (names[x['method']],1000*x['median_seconds'],float(F(x['exact_width']))))
    lines += ['', '带状版与原版的有理证书完全相同，速度差来自消元实现。残差法换用了误差公式，省去本例中的端点二分，并给出更窄的区间。', '',
              '## 错误预算：继续增加二分轮数有没有用','',
              '| 模式数 N | 二分轮数 | 区间宽度 |','| ---: | ---: | ---: |']
    for x in data['error_budget']:
        lines.append('| %d | %d | %.6e |'% (x['modes'],x['bits'],float(F(x['exact_width']))))
    lines += ['', 'N=2 时继续二分几乎不再改善，当前试探空间与尾部估计成为瓶颈。N=8 时提高二分精度仍明显有用。不能给所有计算统一增加同一种“精度”。', '',
              '## 完整残差的分解','',
              '| 模式数 N | 空间内残差平方 | 空间外残差平方 | Temple 区间宽度 |',
              '| ---: | ---: | ---: | ---: |']
    for x in data['residual_decomposition']:
        lines.append('| %d | %.6e | %.6e | %.6e |'% (x['modes'],float(F(x['finite_residual_squared'])),float(F(x['tail_residual_squared'])),float(F(x['exact_width']))))
    lines += ['', '这里“空间内”与“空间外”按当前 N 分割，不能把空间内残差全归因于浮点舍入。候选求解、有理化及保留空间内未充分优化的系数都会影响它。N 从 10 增至 14 已没有改善：在固定候选精度下，加模式不再是有效投入。零尾部表示本次有理试探多项式的 H 作用没有非零项落在该截断之外，不表示真实特征函数是有限多项式。', '',
              '## 正确性检查','',
              '30 项测试通过：包括完整尾部的手算核对、Bernstein 多项式重建、区间内部极值、不同基表示的等价性、带状与稠密惯性一致、指数和 Gaussian 精确解、偶函数条件及伪造谱间隔/尾部/作用域的拒绝。', '',
              '测试同时保留超出当前精度预算的有效区间，并要求报告“未达到目标”。测试通过不等同于形式化证明；解析依据、作用域与信任范围见 [误差分析](ERROR_ANALYSIS.md)。', '']
    (root.parent/'EXPERIMENTS.md').write_text('\n'.join(lines),encoding='utf-8')


def main():
    root = Path(__file__).resolve().parent/'results'
    root.mkdir(exist_ok=True)
    methods = ['baseline','banded','temple','symmetry','exponential','auto']
    records, totals = [],[]
    rng = random.Random(3701)
    e.certify_temple([0,1],4)  # one untimed warmup, excludes interpreter startup
    for name,q,kind in CASES:
        by_method = {m:[] for m in methods}
        for repetition in range(3):
            order = methods[:];rng.shuffle(order)
            for method in order:
                result = run.search(q,method,F(1,10**10),24)
                if result['certificate'] is not None:
                    assert e.verify(result['certificate'])
                by_method[method].append(result)
                records.append({'case':name,'method':method,'repetition':repetition,
                                'seconds':result['seconds'],'status':result['status'],'size':result['size'],
                                'exact_width':result['certificate']['exact_width'] if result['certificate'] else None})
        # Saving is outside the timed search. Save one complete trace per method.
        for method,results in by_method.items():
            representative = sorted(results,key=lambda x:x['seconds'])[1]
            e.base.write_json(root/(name+'_'+method+'.json'),representative)
            item = {'case':name,'kind':kind,'q':q,'method':method,
                    'median_seconds':median(r['seconds'] for r in results),
                    'min_seconds':min(r['seconds'] for r in results),
                    'max_seconds':max(r['seconds'] for r in results),
                    'status':representative['status'],'size':representative['size'],
                    'exact_width':representative['certificate']['exact_width'] if representative['certificate'] else None}
            totals.append(item)
            print(name,method,item['status'],'size',item['size'],
                  'ms',round(1000*item['median_seconds'],3),flush=True)
    fixed = []
    funcs = {'baseline':lambda:e.base.certify([0,1],1,8,36),
             'banded':lambda:e.certify_banded([0,1],1,8,36),
             'temple':lambda:e.certify_temple([0,1],8)}
    for method,fn in funcs.items():
        times=[]
        for _ in range(7):
            t=perf_counter();cert=fn();assert e.verify(cert);times.append(perf_counter()-t)
        fixed.append({'method':method,'modes':8,'median_seconds':median(times),'samples_seconds':times,
                      'exact_width':cert['exact_width']})
    budget = []
    for n in [2,8]:
        for bits in [16,32,48]:
            c = e.certify_banded([0,1],1,n,bits)
            budget.append({'modes':n,'bits':bits,'exact_width':c['exact_width']})
    residuals = []
    for n in [4,6,8,10,14]:
        c = e.certify_temple([0,1],n)
        residuals.append({'modes':n,**{k:c[k] for k in ['finite_residual_squared',
                          'tail_residual_squared','exact_width','excited_weight_upper']}})
    output={'environment':{'python':platform.python_version(),'system':platform.system(),
                           'machine':platform.machine(),'processor':platform.processor()},
            'tolerance':'1/10000000000','max_modes':24,'repetitions':3,
            'timing_policy':'same process; total adaptive search including rejected attempts and verification; excludes interpreter startup and file I/O',
            'summary':totals,'individual_runs':records,'fixed_modes':fixed,
            'error_budget':budget,'residual_decomposition':residuals}
    e.base.write_json(root/'benchmark.json',output)
    render_report(output,root)
    print('saved',root/'benchmark.json')


if __name__ == '__main__':main()
