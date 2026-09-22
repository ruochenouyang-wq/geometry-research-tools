"""Deterministic, public benchmark inputs; never adapts tasks to solver results."""
from fractions import Fraction as F
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def function(profile='abs_power', amplitude='-1/2', exponent='-1/4', offset='0'):
    result = {'kind': 'axis_profile', 'profile': profile, 'axis': 2,
              'amplitude': amplitude, 'offset': offset}
    if profile == 'abs_power':
        result['exponent'] = exponent
    return result


def task(q, tolerance='1/1000000', mean_zero=False, kind='spectrum'):
    result = {'kind': kind, 'function': q, 'tolerance': tolerance,
              'budget': {'wall_seconds': 10}}
    if kind == 'spectrum':
        result['mean_zero'] = mean_zero
    return result


def family(q, tolerance, endpoints=None):
    a = F(q['amplitude'])
    endpoints = endpoints or [str(a-abs(a)/100), str(a+abs(a)/100)]
    return {'profile': q['profile'], 'exponent': q['exponent'], 'axis': q['axis'],
            'amplitude_interval': endpoints, 'offset': q['offset'],
            'mean_zero': False, 'tolerance': tolerance}


def inputs():
    cases = []
    def add(dimension, label, q, tolerance='1/1000000', mean_zero=False,
            kind='spectrum'):
        t = task(q, tolerance, mean_zero, kind)
        row = {'id': 'S' + str(len(cases)+1).zfill(2), 'dimension': dimension,
               'label': label, 'task': t}
        if (kind == 'spectrum' and not mean_zero and q['profile'] == 'abs_power'
                and -F(1,2) < F(q['exponent']) < 0 and F(q['amplitude']) < 0):
            row['family'] = family(q, tolerance)
        cases.append(row)

    for tol in ('1/10000', '1/1000000', '1/100000000'):
        add('precision', '固定负幂，目标宽度 ' + tol, function(), tol)
    for alpha in ('-1/6', '-2/5', '-13/27'):
        add('singularity', '改变奇异指数 ' + alpha, function(exponent=alpha))
    for amplitude in ('-1/10', '-3/2'):
        add('amplitude', '改变势的幅度 ' + amplitude, function(amplitude=amplitude))
    add('stress', '保留已知强势失败题', function(amplitude='-3', exponent='-2/5'),
        '1/100000000')
    add('space', '相同负幂改用平均零空间', function(), mean_zero=True)
    add('space', '平均零空间加负偏移', function(offset='-4'), mean_zero=True)
    add('shape', '阶跃势全空间', function('step', '-1'))
    add('shape', '阶跃势平均零空间', function('step', '-1'), mean_zero=True)
    add('zero_crossing', '常势负谱与零模校准', function('step', '0', offset='-6'))
    add('zero_crossing', '平均零最低谱在零上方',
        function('step', '0', offset='-1999999/1000000'), mean_zero=True)
    add('zero_crossing', '平均零最低谱在零下方',
        function('step', '0', offset='-2000001/1000000'), mean_zero=True)
    add('objective_control', '一般负幂函数全域L2逼近',
        function(exponent='-1/3'), kind='approximation')
    add('objective_control', '更奇异的函数L2逼近及更紧精度',
        function(amplitude='-2/3', exponent='-5/12'),
        '1/100000000', kind='approximation')

    # The nested dyadic order is fixed before any solver run.
    ratios = []
    denominator = 2
    while len(ratios) < 16:
        ratios.extend(F(n, denominator) for n in range(1, denominator, 2))
        denominator *= 2
    ratios = ratios[:16]
    families = []
    for index, endpoints in enumerate((['-101/200','-99/200'], ['-3/5','-2/5']), 1):
        q = function()
        spec = family(q, '1/10000', endpoints)
        lo, hi = map(F, endpoints)
        rows = []
        for i, ratio in enumerate(ratios, 1):
            query = function(amplitude=str(lo + ratio*(hi-lo)))
            rows.append({'id': 'F%dQ%02d' % (index, i),
                         'task': task(query, '1/10000')})
        families.append({'id': 'F' + str(index), 'dimension': 'batch_interval_span',
                         'label': '窄幅度区间' if index == 1 else '宽幅度区间',
                         'spec': spec, 'tasks': rows,
                         'query_ratios': [str(r) for r in ratios]})
    protocol = {
        'version': 'parallel_cross_dimension_public_v1',
        'classification': 'public_development_comparison_not_holdout',
        'systems': ['baseline', 'A', 'B', 'C'], 'repetitions': 5,
        'cold_seconds': 10, 'prepare_seconds': 30, 'query_seconds': 10,
        'session_seconds': 200, 'replay_seconds': 10,
        'prefixes': [1,4,16], 'single_case_count': 18, 'batch_family_count': 2,
        'singlepoint_B_interval_rule': 'a plus/minus abs(a)/100, fixed before solving',
        'route_policy': {'A': 'native_negative_spectrum_projected_to_ground_interval',
                         'B': 'native_prepared_family_no_baseline_fallback',
                         'C': 'frozen_adaptive_solver_including_its_explicit_fallback',
                         'baseline': 'frozen_baseline'},
        'unsupported_policy': 'retain and exclude from speed wins; never silently fallback A/B',
        'success_policy': 'external task binding plus complete certificate verification and exact original tolerance, within hard deadline',
        'comparison_policy': 'report fixed denominators; compare costs on identical successful tasks only; show success and errors first',
        'speed_threshold_cost_reduction': 0.20,
        'resource_metrics': ['wall_seconds','cpu_user_seconds','cpu_system_seconds',
                             'peak_rss_bytes','response_json_utf8_bytes',
                             'certificate_json_utf8_bytes','preparation_bytes'],
        'token_policy': 'no model API calls in measured solvers; actual model tokens and research-agent tokens unmeasured; bytes are not tokens',
        'geometry_scope': 'unit S2; dimension means evaluation axes, not a claim of arbitrary manifold dimension',
        'known_case_disclosure': 'S09 is an explicitly retained known development failure; no unseen-test claim',
        'mutation_policy': 'solver sources unchanged; source and public inputs freeze before final serial runs',
    }
    return {'cases.json': cases, 'families.json': families, 'protocol.json': protocol}


if __name__ == '__main__':
    for name, data in inputs().items():
        path = ROOT / name
        content = json.dumps(data, ensure_ascii=False, sort_keys=True, indent=2) + '\n'
        if path.exists() and path.read_text() != content:
            raise RuntimeError('Refusing to silently replace fixed inputs: ' + name)
        path.write_text(content)
    print('Fixed 18 cases and 2 batch families before performance testing.')
