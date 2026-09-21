"""Cross-route original-task binding and public API regression."""
from copy import deepcopy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import solver


def task(kind, profile='abs_power', mean_zero=True):
    q = {'kind': 'axis_profile', 'profile': profile, 'axis': 2,
         'amplitude': '-3/5', 'offset': '0'}
    if profile == 'abs_power':
        q['exponent'] = '-2/11'
    result = {'kind': kind, 'function': q, 'tolerance': '1/1000000'}
    if kind == 'spectrum':
        result['mean_zero'] = mean_zero
    return result


class Integration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.examples = []
        for request in (task('approximation'), task('spectrum'),
                        task('spectrum', mean_zero=False), task('spectrum', 'step')):
            cls.examples.append((request, solver.solve(request)))

    def test_routes_produce_bound_certificates(self):
        for request, result in self.examples:
            with self.subTest(kind=request['kind'], route=result.get('route')):
                self.assertTrue(result['target_met'], result.get('error', result['status']))
                self.assertTrue(solver.assess(result['certificate'], request)['target_met'])

    def test_each_certificate_rejects_a_different_original_function(self):
        for request, result in self.examples:
            wrong = deepcopy(request)
            wrong['function']['amplitude'] = '-2/5'
            self.assertFalse(solver.assess(result['certificate'], wrong)['certificate_valid'])

    def test_kind_space_and_precision_are_bound(self):
        for request, result in self.examples:
            certificate = result['certificate']
            wrong = deepcopy(request)
            wrong['tolerance'] = '1/100000000'
            self.assertFalse(solver.assess(certificate, wrong)['certificate_valid'])
            wrong = deepcopy(request)
            wrong['kind'] = 'spectrum' if request['kind'] == 'approximation' else 'approximation'
            self.assertFalse(solver.assess(certificate, wrong)['certificate_valid'])
            if request['kind'] == 'spectrum':
                wrong = deepcopy(request)
                wrong['mean_zero'] = not request['mean_zero']
                self.assertFalse(solver.assess(certificate, wrong)['certificate_valid'])

    def test_status_only_is_never_evidence(self):
        request = task('spectrum')
        self.assertFalse(solver.assess({'status': 'target_met', 'target_met': True,
                                       'function': request['function']}, request)['target_met'])

    def test_input_and_zero_budget(self):
        request = task('spectrum')
        request['mean_zero'] = 1
        with self.assertRaises(ValueError):
            solver.solve(request)
        request = task('approximation')
        request['budget'] = {'wall_seconds': 0}
        result = solver.solve(request)
        self.assertIsNone(result['certificate'])
        self.assertFalse(result['target_met'])


if __name__ == '__main__':
    unittest.main()
