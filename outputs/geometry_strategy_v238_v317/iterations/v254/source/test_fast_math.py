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


def source_fixture():
    return json.loads((FROZEN/'certificates'/'singular_n4.json').read_text())


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

    def test_batch_evaluation_seams_signs_and_rejection(self):
        for q in (Q,STEP):
            cert = fast.piecewise_model(q,levels=3,degree=3)
            points = [F(-1),F(-1,16),F(0),F(1,4096),F(1,256),F(1,16),F(1)]
            self.assertEqual(fast.evaluate_many(cert,points),
                             [pieces.evaluate(cert,t) for t in points])
            self.assertEqual(fast.evaluate_many(cert,[]),[])
            with self.assertRaises(ValueError): fast.evaluate_many(cert,[F(2)])
            with self.assertRaises(ValueError): fast.evaluate_many(cert,[0.5])
            damaged = deepcopy(cert); damaged['error_squared'] = '77'
            with self.assertRaises(ValueError): fast.evaluate_many(damaged,[F(0)])

    def test_raw_moments_and_full_assembly_exact(self):
        for q in (Q,STEP,{**Q,'amplitude':'2/3','offset':'-1/7'}):
            table = fast.MomentTable(q)
            for power in range(3):
                for degree in range(13):
                    self.assertEqual(table.moment(degree,power),direct.moment(q,degree,power))
            for m,mean_zero,n,near in ((0,True,4,3),(0,False,3,0),(1,True,5,2),(2,False,3,2)):
                self.assertEqual(fast.matrix_assembly(q,m,mean_zero,n,near),
                                 direct.matrix_assembly(q,m,mean_zero,n,near))
        for bad in (-1,257,True,1.0):
            with self.assertRaises(ValueError): fast.moment(Q,bad)
        with self.assertRaises(ValueError): fast.matrix_assembly(Q,mean_zero=1)

    def test_parity_rule_keeps_step_odd_couplings(self):
        step = {**STEP,'amplitude':'2','offset':'-1'}
        for q in (step,{**Q,'amplitude':'0','offset':'2'},Q):
            table = fast.MomentTable(q)
            for m in range(3):
                for l in range(m,m+4):
                    for k in range(m,m+4):
                        for power in range(3):
                            self.assertEqual(table.integral(m,l,k,power),
                                             direct.integral_product(q,m,l,k,power))
            self.assertEqual(fast.matrix_assembly(q,0,True,6,4),
                             direct.matrix_assembly(q,0,True,6,4))
        table = fast.MomentTable(step)
        self.assertNotEqual(table.integral(1,1,2),0)
        self.assertEqual(table.integral(1,1,2,2),0)

    def test_kernel_schur_and_frozen_certificates(self):
        for q in (Q,STEP):
            new = fast.Kernel(q,1,True,4,40,3)
            old = direct.Kernel(q,1,True,4,40,3)
            self.assertEqual(new.evidence(),old.evidence())
            for x in (F(-3),F(0),new.beta-F(1,32)):
                for kind in ('lower','upper'):
                    self.assertEqual(new.matrix(x,kind),old.matrix(x,kind))
                    self.assertEqual(new.count(x,kind),old.count(x,kind))
            with self.assertRaises(ValueError): new.matrix(new.beta)
        for m in (0,1):
            cert = fast.certify_sector(Q,m,True,4,8,40,2)
            self.assertEqual(cert,direct.certify_sector(Q,m,True,4,8,40,2))
            self.assertTrue(direct.verify_sector(cert,Q,m,True))
        cert = fast.full_ground(Q,modes=3,bits=8,max_m=2,near_tail=1)
        self.assertEqual(cert,direct.full_ground(Q,modes=3,bits=8,max_m=2,near_tail=1))
        self.assertTrue(direct.verify_full(cert,Q,True))

    def test_immutable_source_context_and_trial_wire(self):
        from dataclasses import FrozenInstanceError
        source = source_fixture(); original = deepcopy(source)
        context = fast.SourceContext(source)
        source['lower'] = '99'
        copy = context.source(); copy['upper'] = '-99'
        self.assertEqual(context.source(),original)
        with self.assertRaises(FrozenInstanceError): context._text = '{}'
        with self.assertRaises(ValueError): fast.SourceContext(source)
        proposal = fast.propose_trial(enriched.DEFAULT_EXPONENTS[:3],32,3)
        new = fast.trial_certificate(proposal['powers'],proposal['coefficients'],context)
        old = enriched.trial_certificate(proposal['powers'],proposal['coefficients'],original)
        self.assertEqual(new,old)
        self.assertTrue(enriched.verify(new,original))
        run = fast.enrich(context,max_terms=3,precision_bits=32,iterations=3)
        self.assertTrue(enriched.verify(run['certificate'],original))

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
