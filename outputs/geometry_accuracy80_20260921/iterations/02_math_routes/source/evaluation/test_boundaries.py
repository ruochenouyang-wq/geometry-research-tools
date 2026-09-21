"""Evaluator boundary tests; synthetic fixtures never enter accuracy denominators."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import time
import traceback

import evaluate as evaluation


def synthetic_case(kind='spectrum'):
    task = {'kind': kind, 'tolerance': '1/1000000', 'function': {
        'kind': 'axis_profile', 'axis': 2, 'profile': 'step', 'amplitude': '0', 'offset': '0'}}
    if kind == 'spectrum':
        task['mean_zero'] = False
    return {'id': 'BOUNDARY_ONLY', 'task': task, 'budget_seconds': 0.4, 'family': 'selftest_only',
            'classification': 'synthetic_boundary_test_not_an_evaluation_instance'}


def certificate(task):
    result = {'format': 'synthetic_boundary_exact_v1', 'function': deepcopy(task['function']),
              'kind': task['kind'], 'proof': {'identity': 'exact_zero_fixture'}}
    if task['kind'] == 'spectrum':
        result.update(mean_zero=task['mean_zero'], lower='0', upper='0')
    else:
        result.update(error_upper='0', error_squared='0')
    return result


def synthetic_verify(cert, **bindings):
    # An exact fixture for testing evaluator plumbing, not a mathematical verifier.
    expected = {'kind': bindings['expected_kind'], 'function': bindings['expected_function'],
                'tolerance': bindings['expected_tolerance']}
    if expected['kind'] == 'spectrum':
        expected['mean_zero'] = bindings['expected_mean_zero']
    return cert == certificate(expected)


def synthetic_adapter(method):
    def solve(task):
        if method == '_selftest_catches_exception':
            try:
                time.sleep(2)
            except Exception:
                return certificate(task)
        if method == '_selftest_catches_baseexception':
            while True:
                try:
                    time.sleep(2)
                except BaseException:
                    pass
        cert = certificate(task)
        if method == '_selftest_wrong_function':
            cert['function']['offset'] = '1'
        elif method == '_selftest_wrong_space':
            cert['mean_zero'] = not cert['mean_zero']
        elif method == '_selftest_status_forgery':
            cert.pop('proof')
            cert.update(verified=True, certificate_valid=True, target_met=True)
        elif method == '_selftest_unserializable':
            return {'certificate': cert, 'not_json': {1, 2, 3}}
        elif method == '_selftest_nan':
            return {'certificate': cert, 'not_finite': float('nan')}
        return {'certificate': cert}
    def verify(cert, **bindings):
        if method == '_selftest_slow_verify':
            time.sleep(2)
        return synthetic_verify(cert, **bindings)
    return solve, verify, {'solver': 'synthetic test fixture; never an eligible method'}


def run_tests():
    tests = []
    def check(name, function):
        try:
            function()
            tests.append({'name': name, 'passed': True})
        except BaseException as error:
            tests.append({'name': name, 'passed': False, 'error': repr(error),
                          'traceback': traceback.format_exc()})
    def assert_true(value, message='assertion failed'):
        if not value:
            raise AssertionError(message)

    item = synthetic_case()
    def assess_mutation(name, mutate, reason=None):
        cert = certificate(item['task'])
        mutate(cert)
        result = evaluation.assess(cert, item, synthetic_verify)
        assert_true(not result['certificate_valid'], name)
        if reason:
            assert_true(result['reason'] == reason, result)
    check('exact_fixture_accepted', lambda: assert_true(
        evaluation.assess(certificate(item['task']), item, synthetic_verify)['target_met']))
    check('status_only_forgery_rejected', lambda: assess_mutation('forgery', lambda c: c.pop('proof')))
    check('wrong_original_function_rejected', lambda: assess_mutation('function',
        lambda c: c['function'].update(offset='1'), 'original_function_mismatch'))
    check('wrong_original_space_rejected', lambda: assess_mutation('space',
        lambda c: c.update(mean_zero=True), 'original_space_mismatch'))
    check('wrong_conclusion_kind_rejected', lambda: assess_mutation('kind',
        lambda c: c.update(kind='approximation'), 'original_kind_mismatch'))
    check('wrong_embedded_tolerance_rejected', lambda: assess_mutation('target',
        lambda c: c.update(tolerance='1/100'), 'embedded_target_mismatch'))

    def exact_bound_test():
        cert = certificate(item['task'])
        cert['upper'] = item['task']['tolerance']
        assert_true(evaluation.assess(cert, item, lambda *a,**kw: True)['target_met'])
        cert['upper'] = '1000001/1000000000000'
        result = evaluation.assess(cert, item, lambda *a,**kw: True)
        assert_true(result['certificate_valid'] and not result['target_met'])
        cert.update(lower='1', upper='0')
        assert_true(not evaluation.assess(cert, item, lambda *a,**kw: True)['certificate_valid'])
    check('inclusive_exact_tolerance_and_loose_or_reversed_bounds', exact_bound_test)

    def l2_test():
        approx = synthetic_case('approximation')
        cert = certificate(approx['task'])
        cert.update(error_upper='1/1000000', error_squared='1/1000000')
        assert_true(not evaluation.assess(cert, approx, lambda *a,**kw: True)['certificate_valid'])
        cert['error_squared'] = '1/1000000000000'
        assert_true(evaluation.assess(cert, approx, lambda *a,**kw: True)['target_met'])
    check('L2_norm_not_squared_error_and_enclosure_check', l2_test)

    def duplicate_test():
        first = deepcopy(item)
        second = deepcopy(item)
        first['id'], second['id'] = 'A', 'B'
        second['task']['tolerance'] = '1/10000000000'
        second['task']['function']['amplitude'] = '0/2'
        classified = evaluation.classify_cases([first, second], [])
        assert_true(classified[1]['duplicate_of_current_draw'] == 'A')
        known = [{'source': 'test_catalog', 'case_id': 'PREVIOUS', 'task': first['task']}]
        assert_true(evaluation.classify_cases([second], known)[0]['prior_public_overlap'])
        second['task']['mean_zero'] = True
        assert_true(evaluation.instance_key(first['task']) != evaluation.instance_key(second['task']))
    check('duplicate_identity_ignores_tolerance_but_binds_space', duplicate_test)

    def strong_test():
        q = {'kind': 'axis_profile', 'axis': 2, 'profile': 'abs_power', 'amplitude': '-3/5',
             'offset': '0', 'exponent': '-3/7'}
        assert_true(evaluation.eta_squared(q) == evaluation.Fraction(567, 400))
    check('eta_strength_exact_formula', strong_test)
    def pool_test():
        rules = evaluation.parameter_rules_proposal()
        for name, strong in (('power_weak', False), ('power_strong', True)):
            for exponent in rules[name]['exponents']:
                for amplitude in rules[name]['amplitudes']:
                    q = {'kind': 'axis_profile', 'axis': 2, 'profile': 'abs_power',
                         'amplitude': amplitude, 'offset': '0', 'exponent': exponent}
                    assert_true((evaluation.eta_squared(q) >= 1) is strong, q)
        assert_true('-1/4' not in rules['power_approximation']['exponents'])
    check('all_public_parameter_pool_strengths_exactly_match_strata', pool_test)

    def no_seed_test():
        try:
            evaluation.read_json('/never-created/private_seed.json')
        except RuntimeError as error:
            assert_true('never reads' in str(error))
        else:
            raise AssertionError('seed access must be rejected before opening')
    check('seed_reads_forbidden_before_file_open', no_seed_test)

    def run_fixture(method, expected_success, expected_watchdog=None):
        with tempfile.TemporaryDirectory(prefix='accuracy80-boundary-') as directory:
            row = evaluation.run_worker(method, item, 0, Path(directory), selftest=True)
            assert_true(row['success'] is expected_success, row)
            assert_true(row['process']['cpu_seconds'] >= 0 and row['process']['wall_seconds'] > 0)
            if expected_watchdog:
                assert_true(row['process']['watchdog_stage'] == expected_watchdog, row)
            if expected_success:
                assert_true(row['response']['certificate'] == certificate(item['task']))
                assert_true(row['online_wall_seconds'] <= item['budget_seconds'])
            if method in ('_selftest_slow_verify', '_selftest_catches_exception', '_selftest_catches_baseexception'):
                assert_true(row['process']['wall_seconds'] < 1.5, row['process'])
    for method, success, watchdog in (
        ('_selftest_valid', True, None),
        ('_selftest_wrong_function', False, None),
        ('_selftest_wrong_space', False, None),
        ('_selftest_status_forgery', False, None),
        ('_selftest_unserializable', False, None),
        ('_selftest_nan', False, None),
        ('_selftest_slow_verify', False, None),
        ('_selftest_catches_exception', False, None),
        ('_selftest_catches_baseexception', False, 'online_task'),
    ):
        check('cold_worker_'+method, lambda m=method,s=success,w=watchdog: run_fixture(m,s,w))

    def missing_repeat_test():
        rows = []
        for repetition in (0, 1):
            rows.append({'method': 'current', 'case_id': item['id'], 'repetition': repetition,
                         'success': True, 'process': {'wall_seconds': 0.01, 'cpu_seconds': 0.005},
                         'online_wall_seconds': 0.005, 'online_cpu_seconds': 0.003})
        result = evaluation.summarize(rows, [item], 3)
        assert_true(result['current']['robust_successes'] == 0)
        assert_true(not result['current']['formal_acceptance'])
    check('missing_repeat_is_not_stable_success', missing_repeat_test)
    def single_repeat_test():
        row = {'method': 'current', 'case_id': item['id'], 'repetition': 0,
               'success': True, 'process': {'wall_seconds': 0.01, 'cpu_seconds': 0.005},
               'online_wall_seconds': 0.005, 'online_cpu_seconds': 0.003}
        result = evaluation.summarize([row], [item], 1)['current']
        assert_true(result['robust_successes'] == 0 and result['all_observed_runs_successful_cases'] == 1)
        assert_true(not result['stability_evaluated'] and not result['observed_rate_at_least_80_percent'])
    check('single_development_run_never_counts_as_stable_success', single_repeat_test)
    return {'passed': all(t['passed'] for t in tests), 'count': len(tests), 'tests': tests,
            'production_solver_used': False, 'formal_holdout_used': False,
            'synthetic_fixtures_excluded_from_accuracy_denominators': True}


if __name__ == '__main__':
    result = run_tests()
    print(json.dumps(result, ensure_ascii=False, indent=2))
    raise SystemExit(0 if result['passed'] else 1)
