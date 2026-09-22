"""Bounded correctness checks; no performance comparison or acceptance timings."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

PATH = Path(__file__).with_name('b_adapter.py')
SPEC = importlib.util.spec_from_file_location('_cross_b_adapter_tests', PATH)
b = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(b)
FAMILY = {'profile': 'abs_power', 'exponent': '-1/4', 'axis': 2,
          'amplitude_interval': ['-1/2', '-49/100'], 'offset': '0',
          'mean_zero': False, 'tolerance': '1/1000000'}


class AdapterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.engine = b.Engine(FAMILY)
        cls.prepared = cls.engine.prepare(15)
        if not cls.prepared['ok']:
            raise AssertionError(cls.prepared)
        cls.task = cls.engine.module.endpoint_task(FAMILY, '-99/200')
        cls.cert = cls.engine.solve(cls.task)['certificate']

    def test_native_response_and_first_replay_retained(self):
        p = self.prepared
        self.assertTrue(p['verification']['certificate_valid'])
        self.assertIsNotNone(p['preparation']['bank'])
        self.assertTrue(p['preparation']['attempts'])
        self.assertGreaterEqual(p['preparation_wall_seconds'],
                                p['preparation']['preparation_wall_seconds'])
        for a in FAMILY['amplitude_interval']:
            task = self.engine.module.endpoint_task(FAMILY, a)
            out = self.engine.solve(task)
            self.assertTrue(out['ok'])
            self.assertTrue(self.engine.verify(out['certificate'], task, full=False)['ok'])

    def test_absent_family_and_bad_task_never_fallback(self):
        native = self.engine.module
        with patch.object(native, 'solve', side_effect=AssertionError('baseline fallback')):
            absent = b.Engine()
            self.assertEqual(absent.prepare(1)['status'], 'unsupported')
            self.assertEqual(absent.solve(self.task)['status'], 'unsupported')
            self.assertFalse(absent.verify({'format': 'other'}, self.task)['ok'])
            for field, value in [('mean_zero', True), ('mean_zero', 0), ('kind', 'other')]:
                task = dict(self.task, **{field: value})
                self.assertEqual(self.engine.solve(task)['status'], 'unsupported')
            for field, value in [('amplitude', '-3/4'), ('axis', 1), ('offset', '1'),
                                 ('exponent', '-1/3'), ('profile', 'step')]:
                task = copy.deepcopy(self.task)
                task['function'][field] = value
                self.assertEqual(self.engine.solve(task)['status'], 'unsupported')

    def test_family_snapshot_and_returned_bank_isolation(self):
        supplied = copy.deepcopy(FAMILY)
        engine = b.Engine(supplied)
        supplied['exponent'] = '-1/3'
        response = copy.deepcopy(self.prepared['preparation'])
        with patch.object(engine.module, 'prepare_family', return_value=response):
            out = engine.prepare(5)
        self.assertTrue(out['ok'])
        out['preparation']['bank']['cells'][0]['lower_affine']['intercept'] = '999'
        actual = engine.solve(self.task)
        self.assertTrue(actual['ok'])
        self.assertEqual(actual['certificate'], self.cert)

    def test_invalid_bank_and_failed_reprepare_disable_queries(self):
        engine = b.Engine(FAMILY)
        response = copy.deepcopy(self.prepared['preparation'])
        with patch.object(engine.module, 'prepare_family', return_value=response):
            self.assertTrue(engine.prepare(5)['ok'])
        response['bank']['uniform_width_upper'] = '0'
        with patch.object(engine.module, 'prepare_family', return_value=response):
            self.assertFalse(engine.prepare(5)['ok'])
        self.assertEqual(engine.solve(self.task)['status'], 'not_prepared')

    def test_budget_remaining_and_context_charge(self):
        engine = b.Engine(FAMILY)
        response = copy.deepcopy(self.prepared['preparation'])
        clock = [10.0]
        def preparation(spec, preparation_seconds):
            self.assertEqual(spec, FAMILY)
            self.assertEqual(preparation_seconds, 1.0)
            clock[0] += 0.8
            return response
        class Context:
            verification = {'certificate_valid': True, 'target_met': True}
        def context(bank):
            clock[0] += 0.3
            return Context()
        with patch.object(b, 'perf_counter', side_effect=lambda: clock[0]), \
             patch.object(engine.module, 'prepare_family', side_effect=preparation), \
             patch.object(engine.module, 'prepare_context', side_effect=context):
            out = engine.prepare(1)
        self.assertFalse(out['ok'])
        self.assertEqual(out['status'], 'preparation_budget_exceeded')
        self.assertTrue(out['verification']['certificate_valid'])
        self.assertIsNotNone(out['preparation'])
        self.assertEqual(engine.solve(self.task)['status'], 'not_prepared')
        self.assertFalse(engine.prepare(0)['ok'])
        self.assertEqual(engine.prepare(True)['status'], 'invalid_budget')

    def test_target_binding_tampering_and_self_contained_full_replay(self):
        fresh = b.Engine()
        self.assertTrue(fresh.verify(self.cert, self.task)['ok'])
        self.assertEqual(fresh.verify(self.cert, self.task, full=False)['status'], 'not_prepared')
        bad = copy.deepcopy(self.cert)
        bad['lower'] = bad['upper']
        for full in (False, True):
            self.assertFalse(self.engine.verify(bad, self.task, full=full)['ok'])
        tight = dict(self.task, tolerance='1/1000000000000')
        out = self.engine.solve(tight)
        self.assertTrue(out['certificate_valid'])
        self.assertFalse(out['target_met'])
        self.assertEqual(out['certificate']['tolerance'], tight['tolerance'])
        self.assertFalse(self.engine.verify(out['certificate'], self.task)['ok'])
        other = b.Engine(dict(FAMILY, amplitude_interval=['-51/100', '-49/100']))
        self.assertFalse(other.verify(self.cert, self.task)['ok'])
        malformed = b.Engine(dict(FAMILY, mean_zero=True))
        self.assertEqual(malformed.prepare(1)['status'], 'unsupported')
        self.assertFalse(malformed.verify(self.cert, self.task)['ok'])

    def test_independent_process_complete_replay(self):
        script = '''import importlib.util,json,sys
s=importlib.util.spec_from_file_location("_b_cold",sys.argv[1])
m=importlib.util.module_from_spec(s);s.loader.exec_module(m)
p=json.load(sys.stdin)
r=m.Engine().verify(p["cert"],p["task"],full=True)
print(json.dumps(r))
'''
        run = subprocess.run([sys.executable, '-B', '-c', script, str(PATH)],
                             input=json.dumps({'cert': self.cert, 'task': self.task}),
                             text=True, capture_output=True, check=True, timeout=15)
        self.assertTrue(json.loads(run.stdout)['certificate_valid'])


if __name__ == '__main__':
    unittest.main()
