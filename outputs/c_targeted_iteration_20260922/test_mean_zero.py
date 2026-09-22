"""Small public mathematical controls and deterministic routing regressions."""
from copy import deepcopy
from fractions import Fraction as F
from time import perf_counter
import unittest
from unittest.mock import patch

import mean_zero as route


class Execution:
    def __init__(self, seconds=10):
        self.started, self.seconds = perf_counter(), seconds
        self.attempts, self.counters = [], {}
        self.released = False

    def check(self):
        if perf_counter()-self.started >= self.seconds:
            raise route.Deadline('Functional check deadline')

    def run(self, action, callback, **details):
        self.check()
        return callback()

    def note(self, reason_code, **details):
        self.attempts.append({'reason_code': reason_code, **details})

    def release_reserve(self):
        self.released = True


def task(amplitude='2/3', exponent='-2/7', axis=0, tolerance='1/1000000'):
    return {'kind': 'spectrum', 'mean_zero': True, 'function': {
        'kind': 'axis_profile', 'profile': 'abs_power', 'axis': axis,
        'amplitude': amplitude, 'exponent': exponent, 'offset': '0'},
        'tolerance': tolerance, 'budget': {'wall_seconds': 10}}


class MeanZeroTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.execution = Execution()
        cls.certificate = route.search(task(), cls.execution)

    def test_public_e18_and_frozen_task_bound_replay(self):
        c = self.certificate
        self.assertIsNotNone(c)
        self.assertEqual(c['format'], route.direct.FULL)
        self.assertEqual(c['sectors'][0]['projection'], 'remove_l0')
        self.assertTrue(c['mean_zero'])
        self.assertTrue(route.baseline.assess(c, route.baseline.normalize_task(task()))['target_met'])
        rows = [r for r in self.execution.attempts if r['reason_code']=='mean_zero_refine_active_sector']
        self.assertTrue(rows)
        self.assertTrue(all(r['active_m']==0 for r in rows))
        self.assertEqual(c['sectors'][1]['modes'], 8)

    def test_binding_and_coverage_tamper(self):
        c = self.certificate
        for change in ('space', 'axis', 'tolerance', 'parameter'):
            t = task()
            if change=='space': t['mean_zero']=False
            elif change=='axis': t['function']['axis']=1
            elif change=='parameter': t['function']['amplitude']='3/4'
            else: t['tolerance']='1/100000'
            self.assertFalse(route.baseline.assess(c, route.baseline.normalize_task(t))['certificate_valid'])
        for change in ('lower', 'tail', 'projection', 'remove_sector'):
            bad = deepcopy(c)
            if change=='lower': bad['lower']=bad['upper'];bad['exact_width']='0'
            elif change=='tail': bad['angular_tail_lower']='100000'
            elif change=='projection': bad['sectors'][0]['projection']='none'
            else: bad['sectors'].pop(0)
            self.assertFalse(route.baseline.assess(bad, route.baseline.normalize_task(task()))['certificate_valid'])

    def test_second_positive_control(self):
        t = task('1/3', '-1/4', 2, '1/100000')
        cert = route.search(t, Execution())
        self.assertIsNotNone(cert)
        self.assertTrue(route.baseline.assess(cert, route.baseline.normalize_task(t))['target_met'])

    def test_negative_amplitude_delegates_without_kernel(self):
        with patch.object(route, '_sector', side_effect=AssertionError('Must delegate')):
            self.assertIsNone(route.search(task('-2/3'), Execution()))

    def test_deadline_after_initial_proof_preserves_open(self):
        class Interrupt(Execution):
            def check(self):
                super().check()
                if any(r['reason_code']=='mean_zero_global_interval' for r in self.attempts):
                    raise route.Deadline('Deterministic after initial complete proof')
        ex = Interrupt()
        c = route.search(task(), ex)
        self.assertIsNotNone(c)
        verdict = route.baseline.assess(c, route.baseline.normalize_task(task()))
        self.assertTrue(verdict['certificate_valid'])
        self.assertFalse(verdict['target_met'])
        self.assertTrue(ex.released)

    def test_tight_precision_preserves_valid_open_and_reason(self):
        t = task(tolerance='1/1000000000000000000000000000000')
        with patch.object(route, 'SCHEDULE', ((8, 12),)):
            ex = Execution();c = route.search(t, ex)
        self.assertIsNotNone(c)
        self.assertTrue(route.baseline.assess(c, route.baseline.normalize_task(t))['certificate_valid'])
        self.assertFalse(route.baseline.assess(c, route.baseline.normalize_task(t))['target_met'])
        self.assertTrue(any(r['reason_code']=='mean_zero_active_sector_schedule_exhausted'
                            for r in ex.attempts))

    def test_reserved_assembly_does_not_restart_search(self):
        class InterruptAssembly(Execution):
            def run(self, action, callback, **details):
                if action=='mean_zero_assemble_original_full' and not self.released:
                    raise route.Deadline('Deterministic assembly boundary')
                return super().run(action, callback, **details)
        ex = InterruptAssembly();c = route.search(task(), ex)
        self.assertIsNotNone(c)
        self.assertTrue(ex.released)
        self.assertNotIn('mean_zero_sector_refinements', ex.counters)
        self.assertTrue(any(r['reason_code']=='mean_zero_search_stopped_after_reserved_assembly'
                            for r in ex.attempts))
        self.assertTrue(route.baseline.assess(c, route.baseline.normalize_task(task()))['certificate_valid'])

    def test_dominance_is_not_assumed_from_positive_sign(self):
        c = deepcopy(self.certificate)
        c['sectors'][0]['upper']='10'
        summary = route._summary(c, F(task()['tolerance']))
        self.assertNotIn(0, summary['proved_dominant_sectors'])


if __name__=='__main__':
    unittest.main(verbosity=2)
