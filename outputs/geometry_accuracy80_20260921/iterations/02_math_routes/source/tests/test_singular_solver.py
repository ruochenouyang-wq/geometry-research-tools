import copy
from pathlib import Path
import sys
import unittest

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import singular_solver as s


def task(amplitude='-2/5', exponent='-4/9', tolerance='1/100000000'):
    return {'kind': 'spectrum', 'mean_zero': True, 'tolerance': tolerance,
            'function': {'kind': 'axis_profile', 'axis': 2, 'profile': 'abs_power',
                         'amplitude': amplitude, 'offset': '0', 'exponent': exponent}}


class SingularSolverTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.task = task()
        cls.result = s.solve(cls.task)
        cls.cert = cls.result['certificate']

    def test_h07_target_and_independent_binding(self):
        self.assertEqual(self.result['status'], 'target_met')
        self.assertTrue(s.verify(self.cert, self.task['function'], True, self.task['tolerance']))
        self.assertLessEqual(s.F(self.cert['exact_width']), s.F(self.task['tolerance']))
        self.assertTrue(any(a['status'] == 'certified_open' for a in self.result['attempts']))
        self.assertEqual(self.cert['angular_coverage']['first_omitted_m'], 2)

    def test_h09_stricter_target_and_h05_control(self):
        for t in (task('-1/3', tolerance='1/10000000000'), task('-3/5', '-2/11', '1/1000000')):
            with self.subTest(t=t):
                r = s.solve(t)
                self.assertEqual(r['status'], 'target_met')
                self.assertTrue(s.verify(r['certificate'], t['function'], True, t['tolerance']))

    def test_wrong_function_space_tolerance_and_falsebool_rejected(self):
        self.assertFalse(s.verify(self.cert, task('-1/3')['function'], True, self.task['tolerance']))
        for mean in (False, 1, 0):
            self.assertFalse(s.verify(self.cert, self.task['function'], mean, self.task['tolerance']))
        self.assertFalse(s.verify(self.cert, self.task['function'], True, '1/10'))
        with self.assertRaises(ValueError):
            s.solve(dict(self.task, mean_zero=1))

    def test_complete_proof_tampering_rejected(self):
        edits = (
            lambda c: c.update(lower=c['upper'], exact_width='0'),
            lambda c: c.update(mean_zero=1),
            lambda c: c.update(full_residual_not_projected_residual=1),
            lambda c: c['statistics'].update(residual_squared='0'),
            lambda c: c['gap'].update(lower='6', beta='6'),
            lambda c: c['angular_coverage'].update(lower='999'),
            lambda c: c['m0_source'].update(lower='999'),
            lambda c: c['m1_source'].update(lower='999'),
            lambda c: c['powers'].__setitem__(1, '1'),
            lambda c: c['coefficients'].__setitem__(1, '0'),
        )
        for edit in edits:
            c = copy.deepcopy(self.cert)
            edit(c)
            with self.subTest(edit=edit):
                self.assertFalse(s.verify(c, self.task['function'], True, self.task['tolerance']))

    def test_offset_and_axis_are_bound_not_case_id(self):
        t = task('-1/3')
        t['function'] = dict(t['function'], offset='-3', axis=0)
        r = s.solve(t)
        self.assertEqual(r['status'], 'target_met')
        self.assertTrue(s.verify(r['certificate'], t['function'], True, t['tolerance']))
        self.assertLess(s.F(r['certificate']['upper']), 0)


if __name__ == '__main__':
    unittest.main()
