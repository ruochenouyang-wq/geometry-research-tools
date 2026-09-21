"""Independent adversarial review of request-local replay reuse.

Uses only the named public D06 and original -|z|^-1/4 development problems.
Run directly with Python -B; failures are not suppressed or expected away.
"""
from concurrent.futures import ThreadPoolExecutor
from contextvars import copy_context
from copy import deepcopy
from fractions import Fraction
import importlib
import json
from pathlib import Path
import sys
import tempfile
import threading
import unittest

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
from release_backend import backend, verify_any as frozen_verify
from request_runtime import RequestRuntime, ReplaySession
from release_service import ResearchService


def altered(value, path, replacement):
    out = deepcopy(value)
    node = out
    for key in path[:-1]:
        node = node[key]
    node[path[-1]] = replacement
    return out


def entries_for(session, checker):
    return [json.loads(key) for key in session.entries
            if json.loads(key)['checker'] == checker]


class RuntimeIndependentReview(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.q = {'kind': 'axis_profile', 'profile': 'abs_power', 'exponent': '-1/5',
                 'amplitude': '-1/4', 'offset': '0', 'axis': 2}
        cls.direct = backend.direct.full_ground(cls.q, mean_zero=True, modes=2,
                                               bits=20, max_m=2, tolerance='1/1000')
        cls.original = {'kind': 'spectrum', 'function': {
            'kind': 'axis_profile', 'profile': 'abs_power', 'exponent': '-1/4',
            'amplitude': '-1', 'offset': '0', 'axis': 2}, 'tolerance': '1/100000000'}
        with tempfile.TemporaryDirectory(prefix='geometry-independent-review-') as folder:
            service = ResearchService(evidence_root=folder, cache=False)
            cls.singular = service.evaluate_task(cls.original)['certificate']
        assert frozen_verify(cls.direct) is True
        assert frozen_verify(cls.singular) is True

    def setUp(self):
        self.runtime = RequestRuntime('independent-review-fixed-version')

    def bindings(self, certificate):
        return {'expected_function': certificate['function'],
                'expected_mean_zero': True,
                'expected_tolerance': certificate['tolerance'],
                'expected_kind': 'spectrum'}

    def test_exact_repeat_hits_only_inside_one_request(self):
        with self.runtime.session() as first:
            self.assertTrue(self.runtime.verify(self.direct))
            self.assertTrue(self.runtime.verify(deepcopy(self.direct)))
            self.assertEqual(first.counts['hits'], 1)
        self.assertIsNone(self.runtime.current.get())
        with self.runtime.session() as second:
            self.assertIsNot(first, second)
            self.assertTrue(self.runtime.verify(self.direct))
            self.assertEqual(second.counts['hits'], 0)
            self.assertGreaterEqual(second.counts['misses'], 1)
            self.assertEqual(len(entries_for(second, 'registered_bound_replay')), 1)

    def test_wrong_binding_and_boolean_impostors_do_not_hit(self):
        base = self.bindings(self.direct)
        variations = [
            {'expected_function': dict(self.q, amplitude='-1/3')},
            {'expected_function': dict(self.q, axis=1)},
            {'expected_function': dict(self.q, axis=True)},
            {'expected_mean_zero': False}, {'expected_mean_zero': 1},
            {'expected_mean_zero': 'true'}, {'expected_mean_zero': 1.0},
            {'expected_tolerance': '1/1000000'}, {'expected_tolerance': True},
            {'expected_kind': 'approximation'}, {'_depth': 49},
        ]
        with self.runtime.session() as session:
            self.assertTrue(self.runtime.verify(self.direct, **base))
            for change in variations:
                with self.subTest(change=change):
                    binding = base | change
                    self.assertFalse(frozen_verify(self.direct, **binding))
                    hits = session.counts['hits']
                    self.assertFalse(self.runtime.verify(self.direct, **binding))
                    self.assertEqual(session.counts['hits'], hits)

    def test_complete_direct_mathematical_evidence_tampering_rejected(self):
        changes = [
            (('lower',), '1000'), (('upper',), '-1000'),
            (('exact_width',), '0'), (('angular_tail_lower',), '1000000'),
            (('omitted_azimuth_m_start',), 999),
            (('full_infinite_space_covered',), False),
            (('mean_zero',), 1), (('mean_zero',), False),
            (('sectors', 0, 'lower_inertia'), [0, 0, 0]),
            (('sectors', 0, 'radial_tail_lower'), '1000000'),
            (('sectors', 0, 'near_tail'), True),
            (('sectors', 0, 'modes'), True),
            (('sectors', 0, 'mean_zero'), 1),
            (('sectors', 0, 'function', 'amplitude'), '100'),
            (('sectors',), []),
        ]
        with self.runtime.session() as session:
            self.assertTrue(self.runtime.verify(self.direct))
            for path, replacement in changes:
                with self.subTest(path=path, replacement=replacement):
                    corrupted = altered(self.direct, path, replacement)
                    self.assertFalse(frozen_verify(corrupted))
                    self.assertFalse(self.runtime.verify(corrupted))
                    self.assertFalse(any(backend.canonical(row['certificate']) == backend.canonical(corrupted)
                        for row in entries_for(session, 'registered_bound_replay')))

    def test_singular_source_residual_and_coefficients_tampering_rejected(self):
        self.assertEqual(self.singular['format'], 'adaptive_singular_temple_v1')
        changes = [
            (('lower',), '1000'), (('upper',), '-1000'),
            (('source', 'upper'), '-1000'),
            (('source', 'sectors', 0, 'radial_tail_lower'), '1000000'),
            (('source_digest',), '0' * 64),
            (('statistics', 'residual_squared'), '0'),
            (('statistics', 'rayleigh'), '1000'),
            (('statistics', 'mass'), '0'),
            (('coefficients', 0), '1000'),
            (('matrix_digest',), '0' * 64),
            (('full_residual_not_projected_residual',), False),
            (('mean_zero',), 1),
        ]
        with self.runtime.session():
            self.assertTrue(self.runtime.verify(self.singular, **self.bindings(self.singular)))
            for path, replacement in changes:
                with self.subTest(path=path):
                    corrupted = altered(self.singular, path, replacement)
                    self.assertFalse(frozen_verify(corrupted))
                    self.assertFalse(self.runtime.verify(corrupted, **self.bindings(self.singular)))

    def test_mutating_original_object_cannot_reuse_old_verdict(self):
        certificate = deepcopy(self.direct)
        with self.runtime.session() as session:
            self.assertTrue(self.runtime.verify(certificate))
            certificate['sectors'][0]['upper'] = '-10000'
            self.assertFalse(self.runtime.verify(certificate))
            self.assertEqual(session.counts['hits'], 0)
            self.assertTrue(self.runtime.verify(deepcopy(self.direct)))
            self.assertEqual(session.counts['hits'], 1)

    def test_non_json_tuple_is_not_silently_accepted_on_miss(self):
        certificate = deepcopy(self.direct)
        certificate['sectors'] = tuple(certificate['sectors'])
        self.assertFalse(frozen_verify(certificate))
        with self.runtime.session():
            self.assertFalse(self.runtime.verify(certificate))

    def test_non_json_tuple_cannot_alias_a_list_cache_hit(self):
        certificate = deepcopy(self.direct)
        certificate['sectors'] = tuple(certificate['sectors'])
        with self.runtime.session() as session:
            self.assertTrue(self.runtime.verify(self.direct))
            self.assertFalse(self.runtime.verify(certificate))
            self.assertEqual(session.counts['hits'], 0)

    def test_replay_direct_facade_preserves_boolean_bindings(self):
        facade = self.runtime.adaptive.direct
        with self.runtime.session() as session:
            self.assertTrue(facade.verify_full(self.direct, self.q, True, '1/1000'))
            for bad in (1, 0, 1.0, 'true', False):
                with self.subTest(binding=bad):
                    hits = session.counts['hits']
                    self.assertFalse(facade.verify_full(self.direct, self.q, bad, '1/1000'))
                    self.assertEqual(session.counts['hits'], hits)

    def test_only_literal_true_checker_results_are_cached(self):
        for returned in (False, None, 0, 1, [], [1], 'true'):
            with self.subTest(returned=returned):
                session = ReplaySession('test')
                count = []
                def checker(certificate):
                    count.append(1)
                    return returned
                self.assertFalse(session.call('checker', checker, {}, {}))
                self.assertFalse(session.call('checker', checker, {}, {}))
                self.assertEqual(len(count), 2)
                self.assertEqual(session.entries, {})

    def test_checker_receives_owned_snapshot(self):
        certificate = {'nested': [1]}
        bindings = {'source': {'nested': [2]}}
        def checker(owned, source):
            owned['nested'][0] = 99
            source['nested'][0] = 99
            return True
        session = ReplaySession('test')
        self.assertTrue(session.call('checker', checker, certificate, bindings))
        self.assertEqual(certificate, {'nested': [1]})
        self.assertEqual(bindings, {'source': {'nested': [2]}})

    def test_uncacheable_fraction_binding_runs_original_checker(self):
        binding = self.bindings(self.direct) | {'expected_tolerance': Fraction(1, 1000)}
        with self.runtime.session() as session:
            self.assertTrue(self.runtime.verify(self.direct, **binding))
            self.assertTrue(self.runtime.verify(self.direct, **binding))
            self.assertEqual(entries_for(session, 'registered_bound_replay'), [])
            self.assertEqual(entries_for(session, 'frozen_full'), [])
            self.assertGreaterEqual(session.counts['uncacheable'], 2)

    def test_exception_never_installs_verdict_or_leaks_session(self):
        captured = None
        def raising(certificate):
            raise ArithmeticError('independent deliberate failure')
        with self.assertRaises(ArithmeticError):
            with self.runtime.session() as captured:
                self.runtime.call('deliberate', raising, {}, {})
        self.assertEqual(captured.entries, {})
        self.assertIsNone(self.runtime.current.get())
        with self.runtime.session() as next_session:
            self.assertIsNot(next_session, captured)
            self.assertTrue(self.runtime.verify(self.direct))
            self.assertEqual(next_session.counts['hits'], 0)

    def test_capacity_limits_reduce_reuse_without_changing_validity(self):
        for settings in ({'max_entries': 0}, {'max_bytes': 1}):
            with self.subTest(settings=settings):
                session = ReplaySession('test', **settings)
                for _ in range(2):
                    self.assertTrue(session.call('full', backend.direct.verify_full, self.direct, {}))
                self.assertEqual(session.entries, {})
                self.assertEqual(session.counts['storage_skips'], 2)
                self.assertEqual(session.counts['hits'], 0)

    def test_nested_scope_reuses_outer_and_exception_restores_it(self):
        with self.runtime.session() as outer:
            self.assertTrue(self.runtime.verify(self.direct))
            with self.assertRaises(ValueError):
                with self.runtime.session() as inner:
                    self.assertIs(outer, inner)
                    raise ValueError('nested failure')
            self.assertIs(self.runtime.current.get(), outer)
            self.assertTrue(self.runtime.verify(self.direct))
            self.assertEqual(outer.counts['hits'], 1)
        self.assertIsNone(self.runtime.current.get())

    def test_threads_receive_distinct_sessions(self):
        barrier = threading.Barrier(4)
        def work(index):
            with self.runtime.session() as session:
                barrier.wait(timeout=5)
                self.assertTrue(self.runtime.verify(self.direct))
                self.assertTrue(self.runtime.verify(self.direct))
                snapshot = session.snapshot()
            self.assertIsNone(self.runtime.current.get())
            return session, snapshot
        with ThreadPoolExecutor(max_workers=4) as pool:
            outcomes = list(pool.map(work, range(4)))
        self.assertEqual(len({id(session) for session, _ in outcomes}), 4)
        self.assertTrue(all(row['hits'] == 1 and row['misses'] >= 1 for _, row in outcomes))
        self.assertIsNone(self.runtime.current.get())

    def test_separate_runtime_instances_never_share_sessions(self):
        other = RequestRuntime('independent-review-fixed-version')
        with self.runtime.session() as first, other.session() as second:
            self.assertIsNot(first, second)
            self.assertTrue(self.runtime.verify(self.direct))
            self.assertTrue(other.verify(self.direct))
            self.assertEqual(first.counts['hits'], 0)
            self.assertEqual(second.counts['hits'], 0)

    def test_context_copy_retains_active_request_object_explicitly(self):
        # ContextVar propagates object identity, but a closed session must
        # never provide reusable verdicts in the copied context.
        calls = []
        def checker(certificate):
            calls.append(1)
            return True
        with self.runtime.session() as active:
            self.assertTrue(self.runtime.call('probe', checker, {}, {}))
            context = copy_context()
            self.assertTrue(context.run(self.runtime.call, 'probe', checker, {}, {}))
            self.assertEqual(len(calls), 1)
        self.assertIsNone(self.runtime.current.get())
        self.assertIs(context.run(self.runtime.current.get), active)
        self.assertTrue(active.closed)
        self.assertEqual(active.entries, {})
        self.assertEqual(active.bytes, 0)
        self.assertTrue(context.run(self.runtime.call, 'probe', checker, {}, {}))
        self.assertEqual(len(calls), 2)
        self.assertEqual(active.entries, {})
        def new_request():
            with self.runtime.session() as fresh:
                self.assertIsNot(fresh, active)
                self.assertFalse(fresh.closed)
                self.assertTrue(self.runtime.call('probe', checker, {}, {}))
                self.assertEqual(fresh.counts['hits'], 0)
            return fresh
        fresh = context.run(new_request)
        self.assertTrue(fresh.closed)
        self.assertEqual(len(calls), 3)

    def test_inflight_copied_context_cannot_repopulate_closed_session(self):
        entered, resume = threading.Event(), threading.Event()
        def checker(certificate):
            entered.set()
            if not resume.wait(timeout=5):
                raise TimeoutError('Independent test did not release its checker')
            return True
        with ThreadPoolExecutor(max_workers=1) as pool:
            with self.runtime.session() as active:
                context = copy_context()
                future = pool.submit(context.run, self.runtime.call, 'probe', checker, {}, {})
                self.assertTrue(entered.wait(timeout=5))
            self.assertTrue(active.closed)
            resume.set()
            self.assertTrue(future.result(timeout=5))
        self.assertEqual(active.entries, {})
        self.assertEqual(active.bytes, 0)

    def test_cycle_is_rejected_and_leaves_no_success_entry(self):
        certificate = deepcopy(self.direct)
        certificate['extra'] = certificate
        self.assertFalse(frozen_verify(certificate))
        with self.runtime.session() as session:
            self.assertFalse(self.runtime.verify(certificate))
            self.assertEqual(entries_for(session, 'registered_bound_replay'), [])
            self.assertEqual(entries_for(session, 'frozen_full'), [])

    def test_excessive_nesting_preserves_frozen_false_verdict(self):
        certificate = deepcopy(self.direct)
        certificate['extra'] = {}
        node = certificate['extra']
        for _ in range(1200):
            node['next'] = {}
            node = node['next']
        self.assertFalse(frozen_verify(certificate))
        for reuse in (False, True):
            runtime = RequestRuntime('deep-review', reuse=reuse)
            with runtime.session():
                try:
                    verdict = runtime.verify(certificate)
                except RecursionError:
                    verdict = 'raised RecursionError'
                self.assertIs(verdict, False, 'Excessive nesting must preserve the frozen False verdict')

    def test_sector_binding_m_includes_type_and_value(self):
        sector = self.direct['sectors'][0]
        bindings = {'expected_function': self.q, 'expected_m': 0, 'expected_mean_zero': True}
        with self.runtime.session() as session:
            self.assertTrue(self.runtime.verify_sector(sector, **bindings))
            self.assertTrue(self.runtime.verify_sector(deepcopy(sector), **bindings))
            for change in ({'expected_m': True}, {'expected_m': False},
                           {'expected_m': 0.0}, {'expected_m': 1},
                           {'expected_m': '0'}, {'expected_mean_zero': 1},
                           {'expected_function': dict(self.q, amplitude='-1/2')}):
                with self.subTest(change=change):
                    altered_binding = bindings | change
                    self.assertFalse(backend.direct.verify_sector(sector, **altered_binding))
                    hits = session.counts['hits']
                    self.assertFalse(self.runtime.verify_sector(sector, **altered_binding))
                    self.assertEqual(session.counts['hits'], hits)

    def test_first_sector_replay_builds_frozen_kernel_then_exact_hit_reuses(self):
        sector = self.direct['sectors'][0]
        original = self.runtime.direct.Kernel
        self.assertIsNot(original, self.runtime.fast.Kernel)
        self.assertEqual(original.__init__.__code__.co_code,
                         backend.direct.Kernel.__init__.__code__.co_code)
        calls = []
        def observed(*args, **kwargs):
            calls.append(1)
            return original(*args, **kwargs)
        self.runtime.direct.Kernel = observed
        try:
            with self.runtime.session():
                self.assertTrue(self.runtime.verify_sector(sector))
                self.assertTrue(self.runtime.verify_sector(deepcopy(sector)))
            self.assertEqual(len(calls), 1)
        finally:
            self.runtime.direct.Kernel = original

    def test_sector_reuse_cannot_accept_false_aggregate_or_reordering(self):
        with self.runtime.session():
            for m, sector in enumerate(self.direct['sectors']):
                self.assertTrue(self.runtime.verify_sector(sector, self.q, m, True))
            for certificate in (
                    altered(self.direct, ('sectors',), list(reversed(self.direct['sectors']))),
                    altered(self.direct, ('lower',), '10000'),
                    altered(self.direct, ('angular_tail_lower',), '10000')):
                self.assertFalse(frozen_verify(certificate))
                self.assertFalse(self.runtime.verify_full(certificate))

    def test_private_fast_generation_preserves_frozen_certificate_bytes(self):
        for function in (self.q, dict(self.q, exponent='-1/4', amplitude='-1')):
            with self.subTest(function=function):
                args = {'function': function, 'mean_zero': True, 'modes': 2,
                        'bits': 20, 'max_m': 2, 'tolerance': '1/1000'}
                baseline = backend.direct.full_ground(**args)
                with self.runtime.session():
                    candidate = self.runtime.fast.full_ground(**args)
                self.assertEqual(backend.canonical(baseline), backend.canonical(candidate))
                self.assertTrue(frozen_verify(candidate))

    def test_nonfinite_mathematical_field_cannot_be_cached(self):
        for value in (float('nan'), float('inf'), float('-inf')):
            with self.subTest(value=value):
                certificate = altered(self.direct, ('lower',), value)
                self.assertFalse(frozen_verify(certificate))
                with self.runtime.session() as session:
                    self.assertFalse(self.runtime.verify(certificate))
                    self.assertEqual(entries_for(session, 'registered_bound_replay'), [])

    def test_frozen_global_dependencies_are_not_replaced(self):
        original_adaptive = importlib.import_module('adaptive_singular')
        original_service = importlib.import_module('research_service')
        original_verification = importlib.import_module('verification')
        before = (backend.direct.verify_full, backend.direct.verify_sector,
                  backend.direct.full_ground, backend.direct.certify_sector, backend.direct.Kernel,
                  original_adaptive.direct, original_adaptive.verify,
                  original_service.verify_any, original_service.assess,
                  original_verification.verify_any, importlib.import_module)
        runtime = RequestRuntime('isolation-test')
        with runtime.session():
            self.assertTrue(runtime.verify(self.singular))
        after = (backend.direct.verify_full, backend.direct.verify_sector,
                 backend.direct.full_ground, backend.direct.certify_sector, backend.direct.Kernel,
                 original_adaptive.direct, original_adaptive.verify,
                 original_service.verify_any, original_service.assess,
                 original_verification.verify_any, importlib.import_module)
        self.assertTrue(all(a is b for a, b in zip(before, after)))
        self.assertIsNot(runtime.adaptive, original_adaptive)
        self.assertIsNot(runtime.service, original_service)
        self.assertIsNot(runtime.verification, original_verification)
        self.assertNotIn('_owned_adaptive_singular', sys.modules)
        self.assertNotIn('_owned_direct', sys.modules)
        self.assertNotIn('_owned_fast_math', sys.modules)

    def test_concurrent_full_service_calls_keep_original_problem_bindings(self):
        direct_task = {'kind': 'spectrum', 'function': self.q, 'mean_zero': True,
                       'tolerance': '1/1000'}
        barrier = threading.Barrier(2)
        with tempfile.TemporaryDirectory(prefix='geometry-independent-threads-') as folder:
            service = ResearchService(evidence_root=folder, cache=False)
            def evaluate(task):
                barrier.wait(timeout=5)
                return service.evaluate_task(task)
            with ThreadPoolExecutor(max_workers=2) as pool:
                results = list(pool.map(evaluate, [direct_task, self.original]))
            self.assertIsNone(service.runtime.current.get())
        for task, result in zip([direct_task, self.original], results):
            self.assertTrue(result['certificate_valid'])
            self.assertTrue(result['target_met'])
            self.assertEqual(result['certificate']['function'], task['function'])
            self.assertTrue(frozen_verify(result['certificate'], expected_function=task['function'],
                                         expected_mean_zero=True, expected_tolerance=task['tolerance']))

    def test_original_problem_reuse_toggle_preserves_exact_certificate(self):
        results = []
        for reuse in (False, True):
            with tempfile.TemporaryDirectory(prefix='geometry-independent-toggle-') as folder:
                service = ResearchService(evidence_root=folder, cache=False, replay_reuse=reuse)
                result = service.evaluate_task(self.original)
                self.assertTrue(result['certificate_valid'])
                self.assertTrue(result['target_met'])
                self.assertTrue(frozen_verify(result['certificate'],
                    expected_function=self.original['function'], expected_mean_zero=True,
                    expected_tolerance=self.original['tolerance'], expected_kind='spectrum'))
                results.append(result)
        self.assertEqual(backend.canonical(results[0]['certificate']),
                         backend.canonical(results[1]['certificate']))


if __name__ == '__main__':
    unittest.main(verbosity=2)
