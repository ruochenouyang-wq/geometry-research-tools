import copy
import importlib.util
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

    def test_strong_custom_form_and_every_attempt_retained(self):
        for exponent in ('-5/12', '-13/27'):
            t = task('-22/23', exponent, '1/10000000000')
            r = s.solve(t)
            self.assertEqual(r['status'], 'target_met')
            self.assertTrue(s.verify(r['certificate'], t['function'], True, t['tolerance']))
            self.assertEqual([a['terms'] for a in r['attempts']], [3, 6, 10, 16])
            self.assertTrue(all(a['status'] == 'certified_open' for a in r['attempts'][:-1]))
            self.assertEqual(r['certificate']['m0_source']['format'], s.GROUND_FORMAT)
            for key in ('m0_source', 'm1_source', 'gap'):
                bad = copy.deepcopy(r['certificate'])
                bad[key]['form_certificate']['form']['mass_offset'] = '100'
                self.assertFalse(s.verify(bad, t['function'], True, t['tolerance']))

    def test_private_trial_does_not_relax_shared_frozen_policy(self):
        import adaptive_singular as frozen
        before = frozen.normalize
        q = task('-22/23', '-13/27')['function']
        with self.assertRaises(ValueError):
            frozen.normalize(q)
        engine = s._engine(q)
        self.assertIsNot(engine, frozen)
        self.assertEqual(engine.normalize(q), s.normalize(q))
        self.assertIs(frozen.normalize, before)
        with self.assertRaises(ValueError):
            frozen.normalize(q)
        for exponent in ('-1/2', '-2/3', '0'):
            with self.assertRaises(ValueError):
                engine.normalize(dict(q, exponent=exponent))
        with self.assertRaises(ValueError):
            engine.statistics(q, ['0', '1'], ['1', '1'])

    def test_weak_certificate_is_identical_to_frozen_r02(self):
        root = Path(__file__).resolve().parents[1]
        spec = importlib.util.spec_from_file_location('_test_frozen_r02_singular',
                    root/'iterations'/'02_math_routes'/'source'/'singular_solver.py')
        previous = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(previous)
        expected = previous.solve(self.task)['certificate']
        self.assertEqual(s.canonical(expected), s.canonical(self.cert))


if __name__ == '__main__':
    unittest.main()
