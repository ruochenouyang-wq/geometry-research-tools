"""Exact equivalence, counterexamples and per-step measurement helpers."""
from copy import deepcopy
import cProfile
import json
from pathlib import Path
import pstats
import statistics
import time
import unittest

import fast_math as fast
from backend import F, canonical, pieces, direct, enriched, FROZEN

Q = {'kind': 'axis_profile', 'axis': 2, 'profile': 'abs_power',
     'amplitude': '-1', 'offset': '0', 'exponent': '-1/4'}
STEP = {'kind': 'axis_profile', 'axis': 2, 'profile': 'step', 'amplitude': '-3/2', 'offset': '2/7'}


def calls(function, filename, name):
    profiler = cProfile.Profile()
    result = profiler.runcall(function)
    count = sum(nc for (path, _, fun), (_, nc, _, _, _) in pstats.Stats(profiler).stats.items()
                if path.endswith(filename) and fun == name)
    return result, count


def timed(function, clear=None, repeats=3):
    values = []
    for _ in range(repeats):
        if clear:
            clear()
        start = time.perf_counter()
        function()
        values.append(time.perf_counter()-start)
    return {'serial_seconds': values, 'median_seconds': statistics.median(values),
            'not_a_model_speed_claim': True}


def record_step(version, hypothesis, change, measure):
    """Called explicitly after each real source edit, never by test discovery."""
    from campaign import record
    suite = unittest.defaultTestLoader.loadTestsFromModule(__import__(__name__))
    result = unittest.TextTestRunner(verbosity=1).run(suite)
    if not result.wasSuccessful():
        raise AssertionError('Step validation failed; no snapshot recorded')
    comparison = measure()
    entry = record(version, hypothesis=hypothesis, change=change,
                   files=['fast_math.py', 'test_fast_math.py'],
                   validation={'tests_run': result.testsRun, 'failures': 0,
                               'frozen_verifiers_retained': True},
                   comparison=comparison, outcome='implemented_exact_checks_passed',
                   limitations=['Small serial probes; no claim of universal speedup.',
                                'Caches are process-local; cold and warm costs are identified.'])
    path = Path(__file__).resolve().parent/'metrics'/'fast_math.json'
    path.parent.mkdir(exist_ok=True)
    previous = json.loads(path.read_text()) if path.exists() else {'steps': []}
    previous['steps'].append({'version': version, 'comparison': comparison,
                              'source_sha256': entry['source_sha256']})
    path.write_text(json.dumps(previous, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'recorded': version, 'tests': result.testsRun, 'comparison': comparison}, ensure_ascii=False))


class FastMathTests(unittest.TestCase):
    def setUp(self):
        fast.reset_counters(clear_caches=True)

    def test_local_coordinate_cell_exact(self):
        for left,right,degree,amplitude,offset in (
            (F(0),F(1,64),0,F(-1),F(0)),
            (F(1,8),F(1,4),6,F(2,3),F(-1,7)),
            (F(3,4),F(1),8,F(-1),F(2,7)),
            (F(1,2),F(1),4,F(0),F(3))):
            self.assertEqual(fast._cell(left,right,degree,amplitude,offset),
                             pieces._cell(left,right,degree,amplitude,offset))

    def test_piecewise_exact_wire_and_frozen_replay(self):
        for q in (Q, STEP, {**Q, 'amplitude': '2/3', 'offset': '-1/7', 'axis': 0}):
            new = fast.piecewise_model(q, levels=3, degree=4, root_ratio='3/4')
            old = pieces.piecewise_model(q, levels=3, degree=4, root_ratio='3/4')
            self.assertEqual(new, old)
            self.assertTrue(pieces.verify_piecewise(new, expected_function=q))

    def test_shape_cache_is_immutable_and_bound_to_geometry(self):
        a = fast.piecewise_model(Q, levels=2, degree=3)
        a['shape_proof']['error_squared_dt'] = '-1'
        b = fast.piecewise_model({**Q, 'amplitude': '2'}, levels=4, degree=3)
        self.assertTrue(pieces.verify_piecewise(b))
        self.assertEqual(fast.counters()['shape_integrations'], 1)
        fast.piecewise_model(Q, levels=2, degree=3, root_ratio='3/4')
        self.assertEqual(fast.counters()['shape_integrations'], 2)


if __name__ == '__main__':
    unittest.main()
