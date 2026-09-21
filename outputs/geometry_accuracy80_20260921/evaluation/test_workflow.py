"""Lifecycle tests in temporary directories with fixed synthetic bytes/pools.

Never draw from the real formal parameter pools and never open a real seed.
"""
from copy import deepcopy
import json
from pathlib import Path
import stat
import tempfile
import traceback
from unittest.mock import patch
import zipfile

import evaluate as ev
import workflow as wf

ORIGINAL_DRAW = wf.draw_cases
FIXTURE_BYTES = bytes(range(32))


def synthetic_protocol():
    protocol = deepcopy(ev.draft_protocol())
    rules = protocol['provisional_layout']['parameter_lists_and_sampler']
    # Denominator 101 is disjoint from every real formal amplitude (denominator 23).
    for name, numbers in (('power_weak', range(11, 19)), ('power_strong', range(150, 159)),
                          ('power_approximation', range(21, 38)), ('step_spectrum', range(40, 57))):
        rules[name]['amplitudes'] = [str(ev.Fraction(-k, 101)) for k in numbers]
    return protocol


class Fixture:
    def __init__(self, directory, successes=16):
        self.root = (Path(directory)/'synthetic_geometry_test').resolve()
        self.data = self.root/'evaluation'
        self.data.mkdir(parents=True)
        self.code = self.root/'fixture_solver.py'
        self.code.write_text('FIXTURE = 1\n')
        (self.root/'EVALUATION.md').write_text('Synthetic lifecycle test only.\n')
        ev.write_json(self.data/'protocol.json', ev.draft_protocol())
        old = ev.regression_cases()
        strong = ev.read_json(ev.DATA/'public_strong_stress.json')
        ev.write_json(self.data/'public_strong_stress.json', strong)
        self.entries = [{'source': 'synthetic_public_catalog', 'case_id': c['id'], 'task': c['task']} for c in old+strong]
        self.successes = successes
        self.calls = 0
        self.history = {'runtime_checked_manifest_files': 18, 'runtime_manifest_sha256': 'synthetic',
                        'runtime_preservation': {'checked': 914},
                        'v317_frozen_math': {'checked_files': 1139, 'passed': True}, 'synthetic_only': True}
        self.workflow = wf.Workflow(self.root, catalog_provider=lambda: self.entries,
            source_provider=self.sources, worker_runner=self.worker, history_provider=lambda: deepcopy(self.history))
        reports = []
        for label, cases, phase in (('old', old, 'regression'), ('strong', strong, 'development')):
            directory = self.data/('acceptance_'+label)
            directory.mkdir()
            identities = {str(Path(p).relative_to(self.root.parent)): h for p,h in self.sources().items()}
            for filename in ('source_identity_before.json', 'source_identity_after.json'):
                ev.write_json(directory/filename, identities)
            rows = [{'method': method, 'case_id': c['id'], 'repetition': repetition}
                    for repetition in range(5) for c in cases for method in ev.METHODS]
            report = {'phase': phase, 'repetitions': 5, 'complete': True, 'source_unchanged_during_run': True,
                      'cases': cases, 'rows': rows, 'synthetic_only': True}
            ev.write_json(directory/'report.json', report)
            reports.append({'report_path': str((directory/'report.json').relative_to(self.root)),
                            'sha256': ev.file_digest(directory/'report.json')})
        ev.write_json(self.data/'acceptance_regression.json', {'reports': reports})

    def sources(self):
        return {str(self.code.resolve()): ev.file_digest(self.code)}

    def worker(self, method, item, repetition, directory):
        self.calls += 1
        success = int(item['id'][1:]) <= (self.successes if method == 'current' else 1)
        return {'method': method, 'case_id': item['id'], 'repetition': repetition, 'task': item['task'],
                'success': success, 'process': {'wall_seconds': 0.01, 'cpu_seconds': 0.005},
                'online_wall_seconds': 0.006, 'online_cpu_seconds': 0.003,
                'assessment': {'certificate_valid': success, 'target_met': success},
                'failure_reason': None if success else 'synthetic_retained_failure'}

    def commit(self, batch='B001'):
        with patch.object(wf.os, 'urandom', return_value=FIXTURE_BYTES) as random_source:
            result = self.workflow.commit(batch)
            assert random_source.call_args.args == (32,)
            return result

    def freeze(self, batch='B001'):
        return self.workflow.freeze(batch)

    def authorize(self, batch='B001', authorized=True):
        value = self.workflow.authorization_template(batch)
        value['authorized'] = authorized
        ev.write_json(self.workflow.batch(batch)/'authorization.json', value)

    def reveal(self, batch='B001'):
        synthetic = synthetic_protocol()
        with patch.object(wf, 'draw_cases', side_effect=lambda seed,protocol,catalog: ORIGINAL_DRAW(seed,synthetic,catalog)):
            return self.workflow.reveal(batch)

    def ready(self):
        self.commit()
        self.freeze()
        self.authorize()
        self.reveal()


def run_tests():
    tests = []
    def check(name, operation):
        try:
            with tempfile.TemporaryDirectory(prefix='geometry-lifecycle-synthetic-') as directory:
                operation(Fixture(directory))
            tests.append({'name': name, 'passed': True})
        except BaseException as error:
            tests.append({'name': name, 'passed': False, 'error': repr(error), 'traceback': traceback.format_exc()})
    def expect_failure(operation):
        try:
            operation()
        except (OSError, RuntimeError, ValueError):
            return
        raise AssertionError('Operation unexpectedly accepted')

    def commit_test(f):
        original_open = wf.os.open
        seed_reads = []
        def guarded_open(path, flags, *args):
            if 'private_seed' in str(path) and not flags & wf.os.O_WRONLY:
                seed_reads.append(str(path))
                raise AssertionError('Commit must not read seed')
            return original_open(path, flags, *args)
        with patch.object(wf.os, 'open', side_effect=guarded_open):
            result = f.commit()
        assert not seed_reads
        assert stat.S_IMODE((f.workflow.batch('B001')/'private_seed.bin').stat().st_mode) == 0o600
        assert 'seed_hex' not in ev.canonical(result)
        assert result['commitment']['sha256'] == wf.hashlib.sha256(FIXTURE_BYTES).hexdigest()
        expect_failure(lambda: f.commit('B002'))
    check('commit_writes_0600_random_bytes_without_read_or_disclosure', commit_test)

    def authorization_test(f):
        f.commit(); f.freeze()
        expect_failure(f.reveal)
        f.authorize(authorized=False)
        expect_failure(f.reveal)
        assert not (f.workflow.batch('B001')/'reveal_started.json').exists()
        f.authorize()
        authpath = f.workflow.batch('B001')/'authorization.json'
        auth = ev.read_json(authpath); auth['freeze_sha256'] = '0'*64; ev.write_json(authpath, auth)
        expect_failure(f.reveal)
    check('missing_false_or_wrong_hash_authorization_prevents_seed_read', authorization_test)

    def new_public_test(f):
        f.commit()
        task = deepcopy(f.entries[0]['task']); task['function']['offset'] = '17/101'
        f.entries.append({'source': 'synthetic_new_public', 'case_id': 'NEW', 'task': task})
        expect_failure(f.freeze)
    check('new_public_identity_after_commit_requires_new_batch', new_public_test)

    def old_tolerance_test(f):
        f.commit()
        entry = deepcopy(f.entries[0]); entry['task']['tolerance'] = '1/10000000000'
        f.entries.append(entry)
        f.freeze()
    check('same_public_identity_at_new_tolerance_does_not_change_novelty_set', old_tolerance_test)

    def frozen_source_test(f):
        f.commit(); f.freeze(); f.authorize()
        f.code.write_text('FIXTURE = 2\n')
        expect_failure(f.reveal)
        with zipfile.ZipFile(f.workflow.batch('B001')/'source_snapshot.zip') as archive:
            index = json.loads(archive.read('index.json'))
            assert archive.read(index[str(f.code.resolve())]['member']) == b'FIXTURE = 1\n'
    check('source_mutation_blocks_reveal_and_archived_original_bytes_survive', frozen_source_test)

    def repeated_test(f):
        f.ready()
        expect_failure(f.freeze)
        expect_failure(f.reveal)
        assert len(ev.read_json(f.workflow.batch('B001')/'holdout_cases.json')) == 20
    check('batch_cannot_refreeze_or_reveal_twice', repeated_test)

    def acceptance_test(f):
        pointer = ev.read_json(f.data/'acceptance_regression.json')
        path = f.root/pointer['reports'][0]['report_path']
        report = ev.read_json(path); report['repetitions'] = 1
        ev.write_json(path, report)
        pointer['reports'][0]['sha256'] = ev.file_digest(path)
        ev.write_json(f.data/'acceptance_regression.json', pointer)
        f.commit(); expect_failure(f.freeze)
    check('regression_pointer_requires_five_complete_repeats', acceptance_test)

    def changed_regression_source_test(f):
        f.code.write_text('FIXTURE = 2\n')
        f.commit(); expect_failure(f.freeze)
    check('regression_source_identity_must_match_final_frozen_source', changed_regression_source_test)

    def sampler_test(f):
        protocol = synthetic_protocol()
        first, _ = ORIGINAL_DRAW(FIXTURE_BYTES, protocol, [])
        catalog = wf.identity_catalog([{'source':'synthetic_public', 'case_id':'SEEN', 'task':first[0]['task']}])
        cases, audit = ORIGINAL_DRAW(FIXTURE_BYTES, protocol, catalog)
        assert len(cases) == len({ev.instance_key(c['task']) for c in cases}) == 20
        assert ev.instance_key(first[0]['task']) not in {ev.instance_key(c['task']) for c in cases}
        assert sum(c.get('strength_stratum') == 'power_strong' for c in cases) == 6
        assert all(c['strong_eta_ge_one'] for c in cases if c.get('strength_stratum') == 'power_strong')
        assert sum(c['task']['kind'] == 'approximation' for c in cases) == 6
        assert sum(c['family'] == 'step_spectrum' for c in cases) == 2
        assert all(ev.rational(c['task']['function']['amplitude']).denominator == 101 for c in cases)
        assert audit['solver_called'] is False
    check('synthetic_pool_sampler_keeps_twenty_stratified_new_identities', sampler_test)

    def gate_test(f):
        f.ready()
        result = f.workflow.run('B001')
        assert f.calls == 200 and result['formal_acceptance'] and result['robust_successes'] == 16
        expect_failure(lambda: f.workflow.run('B001'))
        report = ev.read_json(f.workflow.batch('B001')/'formal_report.json')
        assert len(report['rows']) == 200 and sum(not r['success'] for r in report['rows']) == 115
    check('formal_gate_uses_five_repeats_and_sixteen_of_twenty_retaining_failures', gate_test)

    def fail_next_batch_test(f):
        f.successes = 15
        f.ready()
        assert not f.workflow.run('B001')['formal_acceptance']
        revealed = {ev.instance_key(c['task']) for c in ev.read_json(f.workflow.batch('B001')/'holdout_cases.json')}
        f.commit('B002')
        known = wf.identity_set(ev.read_json(f.workflow.batch('B002')/'known_catalog.json'))
        assert revealed <= known
        assert (f.workflow.batch('B001')/'formal_report.json').exists()
    check('failed_batch_is_retained_and_its_revealed_cases_are_public_next_batch', fail_next_batch_test)

    def incomplete_test(f):
        f.ready()
        def crash(*args):
            raise RuntimeError('synthetic worker crash')
        f.workflow.worker_runner = crash
        result = f.workflow.run('B001')
        report = ev.read_json(f.workflow.batch('B001')/'formal_report.json')
        assert result['cases'] == 20 and not result['formal_acceptance'] and not report['complete']
        expect_failure(lambda: f.workflow.run('B001'))
    check('incomplete_formal_run_cannot_pass_or_selectively_restart', incomplete_test)

    def audit_timeout_test(f):
        f.ready()
        original = f.workflow.history_provider
        calls = 0
        def history():
            nonlocal calls
            calls += 1
            if calls >= 2:
                raise wf.subprocess.TimeoutExpired('synthetic historical audit', 60)
            return original()
        f.workflow.history_provider = history
        result = f.workflow.run('B001')
        report = ev.read_json(f.workflow.batch('B001')/'formal_report.json')
        assert len(report['rows']) == 200 and not report['source_unchanged'] and not result['formal_acceptance']
        assert report['error']['type'] == 'TimeoutExpired'
    check('post_run_preservation_timeout_keeps_all_two_hundred_rows_and_fails_gate', audit_timeout_test)

    def protocol_tamper_test(f):
        value = ev.read_json(f.data/'protocol.json'); value['task_seconds'] = 11
        ev.write_json(f.data/'protocol.json', value)
        expect_failure(f.commit)
    check('locked_protocol_bytes_cannot_be_silently_relaxed', protocol_tamper_test)

    return {'passed': all(t['passed'] for t in tests), 'count': len(tests), 'tests': tests,
            'real_commit_freeze_reveal_or_run_invoked': False,
            'real_seed_accessed': False, 'real_formal_parameter_pool_drawn': False,
            'fixture_random_bytes': 'fixed synthetic bytes only, temporary directories',
            'fixture_parameter_pool': 'disjoint denominator 101 amplitudes; never actual heldout denominator 23 pool'}
