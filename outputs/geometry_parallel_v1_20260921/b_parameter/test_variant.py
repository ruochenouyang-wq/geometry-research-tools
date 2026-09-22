import copy
import importlib.util
import json
from pathlib import Path
import unittest
from unittest.mock import patch

path = Path(__file__).with_name('variant.py')
spec = importlib.util.spec_from_file_location('_b_parameter_tests_variant', path)
b = importlib.util.module_from_spec(spec)
spec.loader.exec_module(b)

FAMILY = {'profile': 'abs_power', 'exponent': '-1/4', 'axis': 2,
          'amplitude_interval': ['-1/2', '-49/100'], 'offset': '0',
          'mean_zero': False, 'tolerance': '1/1000000'}


class ParameterTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.prepared = b.prepare_family(FAMILY)
        cls.bank = cls.prepared['bank']
        if cls.bank is None:
            raise AssertionError(cls.prepared)
        cls.ctx = b.prepare_context(cls.bank)
        cls.task = b.endpoint_task(FAMILY, '-99/200')
        cls.result = cls.ctx.query(cls.task)

    def test_whole_closed_interval_uniform_exact_bound(self):
        self.assertEqual(self.prepared['status'], 'target_met')
        self.assertTrue(b.verify_family(self.bank, FAMILY)['certificate_valid'])
        self.assertLessEqual(b.F(self.bank['uniform_width_upper']), b.F(FAMILY['tolerance']))
        self.assertEqual(self.bank['cells'][0]['amplitude_interval'][0], '-1/2')
        self.assertEqual(self.bank['cells'][-1]['amplitude_interval'][1], '-49/100')
        self.assertEqual(len(self.bank['samples']), 2)
        self.assertGreater(len(self.bank['cells'][0]['critical_points']), 2)

    def test_endpoints_and_16_unprepared_interior_parameters(self):
        amplitudes = [b.F(-1, 2), b.F(-49, 100)] + [b.F(-1, 2)+b.F(j, 1700) for j in range(1, 17)]
        for a in amplitudes:
            task = b.endpoint_task(FAMILY, a)
            result = self.ctx.query(task)
            self.assertIs(result['certificate_valid'], True)
            self.assertIs(result['target_met'], True)
            self.assertTrue(self.ctx.verify(result['certificate'], task)['certificate_valid'])
        self.assertEqual(len(self.bank['samples']), 2)

    def test_cold_query_and_self_contained_json_replay(self):
        response = b.query(self.bank, self.task)
        cert = json.loads(json.dumps(response['certificate']))
        self.assertTrue(b.verify(cert, self.task)['certificate_valid'])
        self.assertIn('family_verification', response)
        self.assertGreater(self.ctx.setup_wall_seconds, 0)

    def test_outside_parameter_box_function_axis_offset_space_rejected(self):
        changes = ({'amplitude': '-51/100'}, {'amplitude': '-48/100'},
                   {'exponent': '-1/3'}, {'axis': 1}, {'offset': '1'})
        for change in changes:
            task = copy.deepcopy(self.task)
            task['function'].update(change)
            with self.assertRaises(ValueError):
                self.ctx.query(task)
            self.assertFalse(b.verify(self.result['certificate'], task)['certificate_valid'])
        for mean in (True, 0, 1):
            task = dict(self.task, mean_zero=mean)
            self.assertFalse(b.verify(self.result['certificate'], task)['certificate_valid'])

    def test_query_preserves_stricter_tolerance_and_honest_open_status(self):
        stricter = dict(self.task, tolerance='1/1000000000000')
        response = self.ctx.query(stricter)
        self.assertTrue(response['certificate_valid'])
        self.assertFalse(response['target_met'])
        self.assertEqual(response['certificate']['tolerance'], stricter['tolerance'])
        self.assertEqual(response['status'], 'certified_open')
        self.assertFalse(self.ctx.verify(response['certificate'], self.task)['certificate_valid'])

    def test_bank_function_bounds_coverage_and_false_bools_rejected(self):
        edits = (
            lambda x: x['family'].update(exponent='-1/3'),
            lambda x: x['family'].update(mean_zero=0),
            lambda x: x.update(uniform_width_upper='0'),
            lambda x: x['cells'][0].update(uniform_width_upper='0'),
            lambda x: x['cells'][0]['critical_points'].pop(1),
            lambda x: x['cells'][0]['lower_affine'].update(intercept='100'),
            lambda x: x['samples'][0]['certificate'].update(lower='100'),
            lambda x: x['samples'][0]['rayleigh_affine'].update(slope='0'),
            lambda x: x['samples'].pop(),
            lambda x: x['proof'].update(covers_every_real_amplitude_in_interval=1),
        )
        for edit in edits:
            bad = copy.deepcopy(self.bank)
            edit(bad)
            self.assertFalse(b.verify_family(bad)['certificate_valid'])
            with self.assertRaises(ValueError):
                b.prepare_context(bad)

    def test_point_certificate_tampering_is_rejected_cold_and_prepared(self):
        edits = (
            lambda x: x.update(lower=x['upper'], exact_width='0'),
            lambda x: x.update(cell_index=False),
            lambda x: x.update(mean_zero=0),
            lambda x: x.update(family_digest='f'*64),
            lambda x: x.update(tolerance='1'),
            lambda x: x['family_certificate'].update(uniform_width_upper='0'),
        )
        for edit in edits:
            bad = copy.deepcopy(self.result['certificate'])
            edit(bad)
            self.assertFalse(self.ctx.verify(bad, self.task)['certificate_valid'])
            self.assertFalse(b.verify(bad, self.task)['certificate_valid'])

    def test_context_owns_snapshot_before_and_after_replay(self):
        caller = copy.deepcopy(self.bank)
        original = b.verify_family
        def mutate_caller_after_check(owned, *args, **kwargs):
            result = original(owned, *args, **kwargs)
            caller['cells'][0]['lower_affine']['intercept'] = '999'
            return result
        with patch.object(b, 'verify_family', side_effect=mutate_caller_after_check):
            context = b.prepare_context(caller)
        response = context.query(self.task)
        self.assertTrue(response['certificate_valid'])
        self.assertEqual(response['certificate']['lower'], self.result['certificate']['lower'])
        self.assertNotEqual(caller, response['certificate']['family_certificate'])

    def test_zero_budget_does_not_issue_a_bank_or_query_proof(self):
        result = b.prepare_family(FAMILY, preparation_seconds=0)
        self.assertIsNone(result['bank'])
        self.assertEqual(result['status'], 'preparation_budget_exceeded')
        query = self.ctx.query(dict(self.task, budget={'wall_seconds': 0}))
        self.assertFalse(query['target_met'])
        self.assertIsNone(query['certificate'])

    def test_failed_uniform_check_never_issues_a_bank(self):
        tighter = dict(FAMILY, tolerance='1/100000000')
        with patch.object(b, 'MAX_SAMPLES', 2):
            result = b.prepare_family(tighter)
        self.assertIsNone(result['bank'])
        self.assertEqual(result['status'], 'preparation_sample_budget_exhausted')
        self.assertTrue(any(a.get('phase') == 'uniform_interval_check' and a['target_met'] is False
                            for a in result['attempts']))

    def test_common_interface_explicit_baseline_fallback(self):
        t = b.endpoint_task(FAMILY, '-1/2', '1/1000000')
        response = b.solve(t)
        self.assertIs(response['fallback'], True)
        self.assertTrue(b.verify(response['certificate'], t)['certificate_valid'])


if __name__ == '__main__':
    unittest.main()
