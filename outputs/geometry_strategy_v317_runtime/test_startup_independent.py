"""Independent functional tests for lazy startup; no performance comparisons.

Only public D05 and original -|z|^-1/4 development problems are used.
Filesystem effects are confined to temporary directories.
"""
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
from pathlib import Path
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import release_service
import request_runtime
from release_service import ResearchService
from request_runtime import RequestRuntime
from release_backend import backend, previous_service, verify_any as frozen_verify
import evidence_store
import interaction


STEP = {'kind': 'axis_profile', 'profile': 'step', 'amplitude': '1/4', 'offset': '0', 'axis': 2}
ORIGINAL = {'kind': 'axis_profile', 'profile': 'abs_power', 'exponent': '-1/4',
            'amplitude': '-1', 'offset': '0', 'axis': 2}
CONSTANT = dict(STEP, amplitude='0', offset='-5')
SIMPLE_TASK = {'kind': 'spectrum', 'function': STEP, 'mean_zero': True, 'tolerance': '1/1000'}
SINGULAR_TASK = {'kind': 'spectrum', 'function': ORIGINAL, 'tolerance': '1/100000000'}
APPROX_TASK = {'kind': 'approximation', 'function': ORIGINAL, 'tolerance': '1/1000'}
CONSTANT_REQUEST = {'op': 'research', 'target': 'spectrum', 'function': CONSTANT,
                    'tolerance': '1/1000', 'view': 'full'}


class StartupIndependentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with tempfile.TemporaryDirectory(prefix='geometry-startup-baseline-') as folder:
            old = previous_service.ResearchService(folder, cache=False)
            cls.baseline_simple = old.evaluate_task(SIMPLE_TASK)['certificate']
            cls.baseline_singular = old.evaluate_task(SINGULAR_TASK)['certificate']

    def setUp(self):
        self.folder = tempfile.TemporaryDirectory(prefix='geometry-startup-test-')
        self.addCleanup(self.folder.cleanup)

    def checked(self, result, task):
        self.assertTrue(result['certificate_valid'])
        self.assertTrue(result['target_met'])
        self.assertTrue(frozen_verify(result['certificate'], expected_function=task['function'],
            expected_mean_zero=True if task['kind'] == 'spectrum' else None,
            expected_tolerance=task['tolerance'], expected_kind=task['kind']))

    def test_constructor_defers_math_io_and_shape_objects(self):
        service = ResearchService(self.folder.name, cache=False)
        self.assertFalse(service.runtime._math_ready)
        self.assertFalse(service.service._interaction_ready)
        self.assertIsNone(service.service._shape_cache)
        self.assertNotIn('store', vars(service.service))
        self.assertNotIn('runner', vars(service.service))
        self.assertTrue(Path(self.folder.name).is_dir())

    def test_pure_evaluation_skips_io_but_performs_frozen_math_checks(self):
        with (patch.object(evidence_store.EvidenceStore, '__init__', side_effect=AssertionError('Unused store')),
              patch.object(interaction.RequestRunner, '__init__', side_effect=AssertionError('Unused runner'))):
            service = ResearchService(self.folder.name, cache=False)
            with patch.object(service.runtime, 'uncached_verify', wraps=service.runtime.uncached_verify) as replay:
                result = service.evaluate_task(SIMPLE_TASK)
                self.assertGreater(replay.call_count, 0)
            self.checked(result, SIMPLE_TASK)
            self.assertFalse(service.runtime._math_ready)
            self.assertNotIn('store', vars(service.service))
            self.assertNotIn('runner', vars(service.service))
            broken = deepcopy(result['certificate'])
            broken['upper'] = '-1000'
            with service.runtime.session():
                self.assertFalse(service.runtime.verify(broken))

    def test_approximation_shape_cache_survives_first_handle(self):
        service = ResearchService(self.folder.name, cache=False)
        result = service.evaluate_task(APPROX_TASK)
        self.checked(result, APPROX_TASK)
        shape = service.service._shape_cache
        self.assertIsNotNone(shape)
        self.assertFalse(service.runtime._math_ready)
        self.assertTrue(service.handle(CONSTANT_REQUEST)['certificate_valid'])
        self.assertIs(service.service._shape_cache, shape)

    def test_version_is_computed_once_and_reused_for_io(self):
        with patch.object(ResearchService, 'code_version', return_value='startup-review-version') as version:
            service = ResearchService(self.folder.name, cache=False)
            service.evaluate_task(SIMPLE_TASK)
            service.handle(CONSTANT_REQUEST)
            self.assertEqual(version.call_count, 1)
        self.assertEqual(service.runtime.version, 'startup-review-version')
        self.assertEqual(service.service.store.verifier_version, service.runtime.version)
        self.assertEqual(service.service.runner.namespace, service.runtime.version)

    def test_first_handle_initialization_is_charged_to_wall_budget(self):
        service = ResearchService(self.folder.name, cache=False)
        clock = [10.0]
        original_init = evidence_store.EvidenceStore.__init__
        def initialization_cost(store, *args, **kwargs):
            clock[0] += 0.25
            return original_init(store, *args, **kwargs)
        request = CONSTANT_REQUEST | {'budget': {'wall_seconds': 0.1}}
        with (patch.object(release_service, 'perf_counter', side_effect=lambda: clock[0]),
              patch.object(service.runtime.service, 'perf_counter', side_effect=lambda: clock[0]),
              patch.object(evidence_store.EvidenceStore, '__init__', initialization_cost)):
            first = service.handle(request)
            second = service.handle(request)
        self.assertEqual(first['total_elapsed_seconds'], 0.25)
        self.assertFalse(first['within_budget'])
        self.assertEqual(first['status'], 'budget_exceeded')
        self.assertTrue(first['certificate_valid'])
        self.assertTrue(first['target_met'])
        self.assertEqual(second['total_elapsed_seconds'], 0)
        self.assertTrue(second['within_budget'])

    def test_initialization_failure_retains_shape_cache_and_is_retryable(self):
        service = ResearchService(self.folder.name, cache=False)
        service.evaluate_task(APPROX_TASK)
        shape = service.service._shape_cache
        with patch.object(interaction.RequestRunner, '__init__', side_effect=RuntimeError('Deliberate IO failure')):
            with self.assertRaisesRegex(RuntimeError, 'Deliberate IO failure'):
                service.handle(CONSTANT_REQUEST)
        self.assertFalse(service.service._interaction_ready)
        self.assertIs(service.service._shape_cache, shape)
        self.assertIsNone(service.runtime.current.get())
        result = service.handle(CONSTANT_REQUEST)
        self.assertTrue(result['certificate_valid'])
        self.assertTrue(service.service._interaction_ready)
        self.assertIs(service.service._shape_cache, shape)

    def test_math_initialization_failure_publishes_nothing_and_can_retry(self):
        runtime = RequestRuntime('math-init-failure-test')
        original = request_runtime.private_module
        failed = [False]
        def fail_once(name, filename):
            if name == '_owned_fast_math' and not failed[0]:
                failed[0] = True
                raise RuntimeError('Deliberate math initialization failure')
            return original(name, filename)
        with patch.object(request_runtime, 'private_module', fail_once):
            with self.assertRaisesRegex(RuntimeError, 'Deliberate math initialization failure'):
                runtime._ensure_math()
            self.assertFalse(runtime._math_ready)
            for key in ('direct', 'fast', 'adaptive', 'uncached_full', 'uncached_sector'):
                self.assertNotIn(key, vars(runtime))
            runtime._ensure_math()
        self.assertTrue(runtime._math_ready)
        with runtime.session():
            self.assertTrue(runtime.verify(self.baseline_singular))

    def test_first_singular_certificate_verification_initializes_private_math(self):
        runtime = RequestRuntime('first-certificate-test')
        self.assertFalse(runtime._math_ready)
        with runtime.session():
            self.assertTrue(runtime.verify(self.baseline_singular,
                expected_function=ORIGINAL, expected_mean_zero=True,
                expected_tolerance=SINGULAR_TASK['tolerance'], expected_kind='spectrum'))
        self.assertTrue(runtime._math_ready)

    def test_constructor_rejects_symlink_and_regular_file_like_old_store(self):
        root = Path(self.folder.name)
        target = root / 'target'
        target.mkdir()
        link = root / 'link'
        link.symlink_to(target, target_is_directory=True)
        with self.assertRaises(ValueError):
            ResearchService(link, cache=False)
        regular = root / 'file'
        regular.write_text('not a directory')
        with self.assertRaises(FileExistsError):
            ResearchService(regular, cache=False)

    def test_replaced_storage_root_is_rejected_then_restored_retry_succeeds(self):
        root = Path(self.folder.name)
        storage, saved, elsewhere = root / 'store', root / 'saved', root / 'elsewhere'
        elsewhere.mkdir()
        service = ResearchService(storage, cache=False)
        storage.rename(saved)
        storage.symlink_to(elsewhere, target_is_directory=True)
        with self.assertRaises(ValueError):
            service.handle(CONSTANT_REQUEST)
        self.assertFalse(service.service._interaction_ready)
        self.assertEqual(list(elsewhere.iterdir()), [])
        storage.unlink()
        saved.rename(storage)
        self.assertTrue(service.handle(CONSTANT_REQUEST)['certificate_valid'])

    def test_replaced_parent_is_rejected_without_redirected_storage(self):
        root = Path(self.folder.name)
        parent, saved, elsewhere = root / 'parent', root / 'saved', root / 'elsewhere'
        elsewhere.mkdir()
        service = ResearchService(parent / 'store', cache=False)
        parent.rename(saved)
        parent.symlink_to(elsewhere, target_is_directory=True)
        with self.assertRaises(ValueError):
            service.handle(CONSTANT_REQUEST)
        self.assertFalse(service.service._interaction_ready)
        self.assertFalse((elsewhere / 'store').exists())
        parent.unlink()
        saved.rename(parent)
        self.assertTrue(service.handle(CONSTANT_REQUEST)['certificate_valid'])

    def test_concurrent_first_math_initialization_builds_one_private_set(self):
        runtime = RequestRuntime('concurrent-init-test')
        original = request_runtime.private_module
        labels = []
        barrier = threading.Barrier(4)
        def observed(name, filename):
            labels.append(name)
            return original(name, filename)
        def initialize(index):
            barrier.wait(timeout=5)
            runtime._ensure_math()
            return id(runtime.direct), id(runtime.fast), id(runtime.adaptive)
        with patch.object(request_runtime, 'private_module', observed):
            with ThreadPoolExecutor(max_workers=4) as pool:
                identities = list(pool.map(initialize, range(4)))
        self.assertEqual(len(set(identities)), 1)
        self.assertCountEqual(labels, ['_owned_direct', '_owned_fast_math', '_owned_adaptive_singular'])

    def test_concurrent_first_handle_constructs_io_once(self):
        service = ResearchService(self.folder.name, cache=False)
        original = interaction.RequestRunner.__init__
        count = []
        barrier = threading.Barrier(4)
        def observed(runner, *args, **kwargs):
            count.append(1)
            return original(runner, *args, **kwargs)
        def invoke(index):
            barrier.wait(timeout=5)
            return service.handle(CONSTANT_REQUEST)
        with patch.object(interaction.RequestRunner, '__init__', observed):
            with ThreadPoolExecutor(max_workers=4) as pool:
                results = list(pool.map(invoke, range(4)))
        self.assertEqual(len(count), 1)
        self.assertTrue(all(row['certificate_valid'] and row['target_met'] for row in results))
        self.assertEqual(len({backend.canonical(row['certificate']) for row in results}), 1)
        self.assertIsNone(service.runtime.current.get())

    def test_simple_singular_simple_and_reverse_preserve_exact_certificates(self):
        for order in ((SIMPLE_TASK, SINGULAR_TASK, SIMPLE_TASK),
                      (SINGULAR_TASK, SIMPLE_TASK, SINGULAR_TASK)):
            with self.subTest(first=order[0]['function']['profile']):
                service = ResearchService(self.folder.name, cache=False)
                for task in order:
                    result = service.evaluate_task(task)
                    self.checked(result, task)
                    expected = self.baseline_simple if task is SIMPLE_TASK else self.baseline_singular
                    self.assertEqual(backend.canonical(result['certificate']), backend.canonical(expected))
                self.assertTrue(service.runtime._math_ready)
                self.assertFalse(service.service._interaction_ready)
                self.assertNotIn('store', vars(service.service))

    def test_late_summary_fetch_and_bound_verification_retain_original_behavior(self):
        service = ResearchService(self.folder.name, cache=False)
        service.evaluate_task(SIMPLE_TASK)
        request = {k: v for k, v in CONSTANT_REQUEST.items() if k != 'view'}
        summary = service.handle(request)
        self.assertTrue(summary['certificate_valid'])
        fetched = service.handle({'op': 'fetch', 'reference': summary['evidence_ref']['sha256']})
        self.assertTrue(fetched['certificate_valid'])
        mismatch = service.handle({'op': 'verify', 'certificate': fetched['certificate'],
            'function': CONSTANT, 'mean_zero': False, 'expected_kind': 'spectrum', 'tolerance': '1/1000'})
        self.assertTrue(mismatch['certificate_valid'])
        self.assertFalse(mismatch['verified'])


if __name__ == '__main__':
    unittest.main(verbosity=2)
