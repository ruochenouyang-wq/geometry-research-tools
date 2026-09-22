"""Public extension cases fixed before any comparative performance run."""
from pathlib import Path
import json
from runtime import baseline

ROOT = Path(__file__).resolve().parent
cases = []


def add(profile, amplitude, mean_zero=False, exponent=None, digits=6, axis=2,
        offset='0', dimension='', classification='new_public', label='', kind='spectrum'):
    q = {'kind': 'axis_profile', 'profile': profile, 'axis': axis,
         'amplitude': amplitude, 'offset': offset}
    if exponent is not None:
        q['exponent'] = exponent
    task = {'kind': kind, 'function': q, 'tolerance': f'1/{10**digits}',
            'budget': {'wall_seconds': 10}}
    if kind == 'spectrum':
        task['mean_zero'] = mean_zero
    cases.append({'id': f'E{len(cases)+1:02}', 'dimension': dimension,
                  'classification': classification, 'label': label,
                  'task': baseline.normalize_task(task)})


for amplitude, digits, axis, offset, name in [
        ('-1/3', 6, 0, '0', '弱阶跃'),
        ('-5/4', 8, 1, '1/2', '偏移阶跃高精度'),
        ('-7/2', 6, 2, '0', '较强阶跃'),
        ('2/3', 6, 0, '-1', '正幅度阶跃')]:
    for mean in (False, True):
        add('step', amplitude, mean, digits=digits, axis=axis, offset=offset,
            dimension='step_strength_space_sign_precision',
            label=name+('，平均零空间' if mean else '，完整空间'))

power_cases = [
    ('-2/7','-2/3',False,6,0,'0','新指数完整空间'),
    ('-2/7','-2/3',False,10,0,'0','同题严格精度'),
    ('-2/7','-2/3',True,8,1,'0','新指数平均零空间'),
    ('-3/8','-3/4',False,8,1,'0','中等奇异性完整空间'),
    ('-3/8','-3/4',True,8,1,'3/2','中等奇异性，正偏移'),
    ('-7/15','-2/5',False,6,0,'0','接近残差可积边界，完整空间'),
    ('-7/15','-2/5',True,6,0,'0','接近残差可积边界，平均零空间'),
    ('-3/8','-2',False,6,1,'0','较强吸引势'),
    ('-2/7','-3',False,6,1,'0','大幅度吸引势'),
    ('-2/7','2/3',True,6,0,'0','正幅度平均零空间')]
for exponent, amplitude, mean, digits, axis, offset, label in power_cases:
    add('abs_power', amplitude, mean, exponent, digits, axis, offset,
        'singularity_strength_space_precision_sign', label=label)

add('abs_power','-3',False,'-2/5',8,2, dimension='old_failure',
    classification='old_failure_regression', label='旧 S09：四路线均失败')
add('abs_power','-1/2',False,'-13/27',8,2, dimension='old_fallback',
    classification='old_failure_regression', label='旧 S06 结构：精度由1e-6提高到1e-8')
add('abs_power','-2/3',exponent='-2/7',axis=1,kind='approximation',
    dimension='approximation_control',label='新逼近题，预计沿用基线')
add('abs_power','-3/4',exponent='-3/8',digits=8,axis=0,kind='approximation',
    dimension='approximation_control',label='新逼近题，更严格误差')
add('abs_power','-1/2',False,'-1/2',6,2,dimension='method_boundary',
    classification='expected_unsupported_control',label='强残差平方不可积边界')
add('abs_power','2/3',False,'-2/7',6,2,dimension='method_boundary',
    classification='expected_unsupported_control',label='完整空间正幅度：当前路线未支持')

assert len(cases) == 24
(ROOT/'cases.json').write_text(json.dumps(cases,ensure_ascii=False,indent=2)+'\n')
protocol = {
    'version': 'c2_public_extension_v1', 'case_count': 24, 'repetitions': 3,
    'cold_seconds': 10, 'systems': ['baseline','C1','C2'],
    'execution': 'serial_cold_process_rotated_order',
    'case_classification': 'public_development_not_holdout',
    'success': 'original task-bound verifier valid AND exact metric <= original tolerance AND process within deadline',
    'stable_success': 'success in all three repetitions',
    'expected_unsupported_controls': 'Report separately; never count expected rejection as a solved math problem.',
    'scope': 'Axis-profile potentials on unit S2; axis changes are rotation controls, not new ambient dimensions.',
    'model_token_measurement': 'No live model API calls; do not infer tokens from certificate bytes.',
    'cases_locked_before_comparative_timing': True,
    'failure_policy': 'Retain every run and trace; no performance tuning or case replacement after freeze.'}
(ROOT/'protocol.json').write_text(json.dumps(protocol,ensure_ascii=False,indent=2)+'\n')
print('24 public cases written: 20 new, 2 old regressions, 2 unsupported controls')
