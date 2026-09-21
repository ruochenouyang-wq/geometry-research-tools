"""One-way commitment, freeze, authorized reveal, and accuracy evaluation.

This module does nothing on import. Real lifecycle commands require explicit
invocation. Only commit creates private random bytes; only authorized reveal
opens them. The seed is never printed or included in any public artifact.
"""
from copy import deepcopy
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import time
import zipfile

import evaluate as ev

LOCKED_PROTOCOL_SHA256 = 'b8c17b004e02e795051c4d868bd49d3ac41ce8fa8b91af9e157fffb22885d7fc'


def now():
    return datetime.now(timezone.utc).isoformat()


def once(path, value, mode=0o644):
    """O_EXCL is also the transition lock: a transition cannot be repeated."""
    path = ev.forbidden(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, mode)
    with os.fdopen(fd, 'w') as stream:
        stream.write(ev.canonical(value)+'\n')
        stream.flush()
        os.fsync(stream.fileno())


def identity_catalog(entries):
    by_key = {}
    for entry in entries:
        key = ev.instance_key(entry['task'])
        by_key.setdefault(key, {'identity': json.loads(key), 'references': []})['references'].append(
            {'source': entry['source'], 'case_id': entry['case_id']})
    return [by_key[key] for key in sorted(by_key)]


def identity_set(catalog):
    return {ev.canonical(entry['identity']) for entry in catalog}


def slots(protocol):
    rules = protocol['provisional_layout']['parameter_lists_and_sampler']
    for tolerance in protocol['provisional_layout']['tolerances']:
        for mean_zero in (False, True):
            for strength in ('weak', 'strong'):
                yield 'power_'+strength, 'spectrum', mean_zero, tolerance
    for tolerance in protocol['provisional_layout']['tolerances']:
        for _ in range(2):
            yield 'power_approximation', 'approximation', None, tolerance
    for mean_zero in (False, True):
        yield 'step_spectrum', 'spectrum', mean_zero, protocol['provisional_layout']['step_spectrum']['tolerance']


def candidate_pool(protocol, pool_name, kind, mean_zero, tolerance):
    rules = protocol['provisional_layout']['parameter_lists_and_sampler']
    pool = rules[pool_name]
    result = []
    for exponent in pool.get('exponents', [None]):
        for amplitude in pool['amplitudes']:
            function = {'kind': 'axis_profile', 'axis': rules['axis'], 'offset': rules['offset'],
                        'amplitude': amplitude, 'profile': 'step' if exponent is None else 'abs_power'}
            if exponent is not None:
                function['exponent'] = exponent
            task = {'kind': kind, 'function': function, 'tolerance': tolerance}
            if mean_zero is not None:
                task['mean_zero'] = mean_zero
            result.append(task)
    # Exact canonical identity lexicographic order defines "sorted Cartesian".
    return sorted(result, key=ev.instance_key)


def draw_cases(seed, protocol, catalog):
    """Pure deterministic sampler; test fixtures use fixed, non-private bytes."""
    if type(seed) is not bytes or len(seed) != 32:
        raise ValueError('Sampler requires exactly 32 bytes')
    known, selected, cases, audits = identity_set(catalog), set(), [], []
    counter = 0
    for pool_name, kind, scope, tolerance in slots(protocol):
        full = candidate_pool(protocol, pool_name, kind, scope, tolerance)
        excluded = []
        eligible = []
        for task in full:
            key = ev.instance_key(task)
            if key in known or key in selected:
                excluded.append({'identity': json.loads(key),
                                 'reason': 'prior_public_identity' if key in known else 'already_drawn_identity'})
            else:
                eligible.append(task)
        if not eligible:
            raise RuntimeError('A fixed pool is exhausted; no replacement or relaxed rule is allowed')
        size = len(eligible)
        limit = (1 << 256) - ((1 << 256) % size)
        start_counter = counter
        while True:
            word = int.from_bytes(hashlib.sha256(seed + counter.to_bytes(8, 'big')).digest(), 'big')
            counter += 1
            if word < limit:
                index = word % size
                break
        task = eligible[index]
        key = ev.instance_key(task)
        selected.add(key)
        item = {'id': 'N'+str(len(cases)+1).zfill(2), 'task': task,
                'family': 'negative_power_spectrum' if pool_name in ('power_weak', 'power_strong') else pool_name,
                'strength_stratum': pool_name if kind == 'spectrum' and pool_name.startswith('power_') else None,
                'budget_seconds': 10.0, 'classification': 'new_parameters_in_known_mathematical_family',
                'independent_new_instance': True, 'prior_public_overlap': [], 'duplicate_of_current_draw': None}
        if task['function']['profile'] == 'abs_power':
            item['centered_eta_squared'] = str(ev.eta_squared(task['function']))
            item['strong_eta_ge_one'] = ev.eta_squared(task['function']) >= 1
        cases.append(item)
        audits.append({'case_id': item['id'], 'pool': pool_name, 'total_pool_size': len(full),
                       'eligible_pool_size': size, 'excluded_identities': excluded,
                       'counter_start': start_counter, 'counter_end_exclusive': counter,
                       'chosen_index': index})
    if len(cases) != 20 or len(selected) != 20 or selected & known:
        raise RuntimeError('Exactly twenty new unique instances are mandatory')
    return cases, {'ordering': 'canonical JSON identity, Unicode code-point lexicographic order',
                   'counter_initial': 0, 'counter_end_exclusive': counter,
                   'word_rule': 'unsigned big-endian SHA256; reject word >= 2^256-(2^256 mod pool_size)',
                   'draws': audits, 'solver_called': False}


class Workflow:
    def __init__(self, root=None, catalog_provider=None, source_provider=None, worker_runner=None, history_provider=None):
        self.root = Path(root or ev.ROOT).resolve()
        self.data = self.root/'evaluation'
        self.batches = self.data/'batches'
        self.catalog_provider = catalog_provider or ev.public_catalog
        self.source_provider = source_provider or self.default_sources
        self.worker_runner = worker_runner or ev.run_worker
        self.history_provider = history_provider or self.audit_history

    def batch(self, batch_id):
        if not re.fullmatch(r'B[0-9]{3,}', batch_id):
            raise ValueError('Batch ID must be B followed by at least three digits')
        return self.batches/batch_id

    def protocol(self):
        path = self.data/'protocol.json'
        if ev.file_digest(path) != LOCKED_PROTOCOL_SHA256:
            raise RuntimeError('Locked protocol bytes changed; cannot create or assess a formal batch')
        return ev.read_json(path)

    def catalog(self):
        entries = deepcopy(self.catalog_provider())
        # Always remember revealed batches, even after a crash before registry update.
        for path in sorted(self.batches.glob('B*/holdout_cases.json')):
            for item in ev.read_json(path):
                entries.append({'source': str(path), 'case_id': item['id'], 'task': item['task']})
        return identity_catalog(entries)

    def default_sources(self):
        # All current source/tests/reviews plus preserved geometry dependency code.
        paths = set(self.root.rglob('*.py'))
        for directory in self.root.parent.glob('geometry_*'):
            if directory.is_dir():
                paths.update(directory.rglob('*.py'))
        return {str(path.resolve()): ev.file_digest(path) for path in sorted(paths)}

    def audit_history(self):
        """Execute only known read-only audits in a pristine interpreter."""
        code = '''import sys,json,hashlib
from pathlib import Path
runtime,old=sys.argv[1:]
sys.path[:0]=[runtime,old]
import release_evaluation as audit
manifest=audit.checked_manifest()
preserved=audit.verify_old_preservation()
import audit_campaign
legacy=audit_campaign.audit_frozen()
if not legacy['passed']: raise RuntimeError('Historical mathematical baseline changed')
result={'runtime_checked_manifest_files':len(manifest['files']),
        'runtime_manifest_sha256':hashlib.sha256((Path(runtime)/'evaluation/freeze_manifest.json').read_bytes()).hexdigest(),
        'runtime_preservation':preserved,'v317_frozen_math':legacy,
        'old_private_seed_read':False}
print(json.dumps(result,sort_keys=True))
'''
        before = time.monotonic()
        result = subprocess.run([sys.executable, '-B', '-c', code, str(ev.RUNTIME), str(ev.OLD)],
            capture_output=True, text=True, timeout=60,
            env=dict(os.environ, PYTHONDONTWRITEBYTECODE='1'))
        if result.returncode != 0:
            raise RuntimeError('Read-only historical preservation audit failed: '+result.stderr[-2000:])
        report = json.loads(result.stdout)
        report['audit_wall_seconds'] = time.monotonic()-before
        return report

    def commit(self, batch_id):
        protocol = self.protocol()
        self.batches.mkdir(parents=True, exist_ok=True)
        for previous in self.batches.glob('B*'):
            if not (previous/'formal_report.json').exists() and not (previous/'abandoned.json').exists():
                raise RuntimeError('Earlier batch is unfinished; retain it and explicitly abandon or finish it first')
        directory = self.batch(batch_id)
        directory.mkdir(exist_ok=False)
        catalog = self.catalog()
        once(directory/'known_catalog.json', catalog)
        # Seed creation is intentionally below immutable catalog creation.
        secret = os.urandom(32)
        private = directory/'private_seed.bin'
        fd = os.open(private, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
        try:
            os.fchmod(fd, 0o600)
            with os.fdopen(fd, 'wb') as stream:
                stream.write(secret)
                stream.flush()
                os.fsync(stream.fileno())
            commitment = {'batch_id': batch_id, 'created_utc': now(),
                'sha256': hashlib.sha256(secret).hexdigest(), 'seed_bytes': 32,
                'seed_source': 'os.urandom(32)', 'seed_read_during_commit': False,
                'catalog_sha256': ev.file_digest(directory/'known_catalog.json'),
                'protocol_sha256': LOCKED_PROTOCOL_SHA256, 'public_identity_count': len(catalog)}
        finally:
            del secret
        once(directory/'commitment.json', commitment)
        return {'batch_id': batch_id, 'commitment': commitment, 'holdout_generated': False}

    def assert_not_closed(self, directory):
        if (directory/'abandoned.json').exists():
            raise RuntimeError('Batch was abandoned and cannot be resumed')

    def acceptance_regression(self, sources):
        pointer_path = self.data/'acceptance_regression.json'
        pointer = ev.read_json(pointer_path)
        if type(pointer) is not dict or set(pointer) != {'reports'} or len(pointer['reports']) != 2:
            raise RuntimeError('Acceptance pointer must bind two reports: old ten and public four strong tasks')
        expected_old = {c['id']: c['task'] for c in ev.regression_cases()}
        expected_strong = {c['id']: c['task'] for c in ev.read_json(self.data/'public_strong_stress.json')}
        if len(expected_old) != 10 or len(expected_strong) != 4:
            raise RuntimeError('Ten old and four predeclared public strong tasks must be retained')
        artifacts, found = [pointer_path, self.data/'public_strong_stress.json'], set()
        for entry in pointer['reports']:
            if type(entry) is not dict or set(entry) != {'report_path', 'sha256'}:
                raise RuntimeError('Each report pointer needs exactly report_path and sha256')
            relative = Path(entry['report_path'])
            path = (self.root/relative).resolve()
            if relative.is_absolute() or '..' in relative.parts or not path.is_relative_to(self.data) or path.name != 'report.json':
                raise RuntimeError('Report path must stay inside the ROOT-relative evaluation directory')
            if ev.file_digest(path) != entry['sha256']:
                raise RuntimeError('Acceptance regression report digest mismatch')
            report = ev.read_json(path)
            if report.get('phase') not in ('regression', 'development') or report.get('repetitions') != 5 or not report.get('complete') or report.get('source_unchanged_during_run') is not True:
                raise RuntimeError('Acceptance reports require five complete cold repeats on unchanged source')
            actual = {c['id']: c['task'] for c in report['cases']}
            label = 'old_ten' if actual == expected_old else 'public_strong_four' if actual == expected_strong else None
            if label is None or label in found or len(actual) != len(report['cases']):
                raise RuntimeError('Acceptance cases must match the original ten and public strong four exactly')
            found.add(label)
            expected_rows = {(m, cid, r) for m in ev.METHODS for cid in actual for r in range(5)}
            row_keys = [(row['method'], row['case_id'], row['repetition']) for row in report['rows']]
            if set(row_keys) != expected_rows or len(row_keys) != len(expected_rows):
                raise RuntimeError('Acceptance reports must retain every method/case/repeat including failures')
            for name in ('source_identity_before.json', 'source_identity_after.json'):
                identity_path = path.parent/name
                saved = ev.read_json(identity_path)
                absolute = {str((self.root.parent/key).resolve()): value for key, value in saved.items()}
                if absolute != sources:
                    raise RuntimeError('Acceptance regression was run on a different source set or bytes')
                artifacts.append(identity_path)
            artifacts.append(path)
        return pointer, artifacts

    def archive_sources(self, directory, sources, artifacts):
        archive_path = directory/'source_snapshot.zip'
        entries = {**sources, **{str(p.resolve()): ev.file_digest(p) for p in artifacts}}
        index = {}
        fd = os.open(archive_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
        with os.fdopen(fd, 'wb') as stream, zipfile.ZipFile(stream, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
            for number, (path, expected) in enumerate(sorted(entries.items())):
                data = ev.forbidden(path).read_bytes()
                if hashlib.sha256(data).hexdigest() != expected:
                    raise RuntimeError('Source/artifact changed while archiving')
                member = 'files/'+str(number).zfill(6)+'.blob'
                archive.writestr(member, data)
                index[path] = {'member': member, 'sha256': expected, 'bytes': len(data)}
            archive.writestr('index.json', ev.canonical(index)+'\n')
        return archive_path

    def freeze(self, batch_id):
        directory = self.batch(batch_id)
        self.assert_not_closed(directory)
        if (directory/'freeze_started.json').exists() or (directory/'freeze.json').exists():
            raise RuntimeError('This batch already attempted a freeze; never refreeze it')
        self.protocol()
        commitment = ev.read_json(directory/'commitment.json')
        if commitment['batch_id'] != batch_id or commitment['protocol_sha256'] != LOCKED_PROTOCOL_SHA256:
            raise RuntimeError('Commitment identity mismatch')
        catalog_path = directory/'known_catalog.json'
        if ev.file_digest(catalog_path) != commitment['catalog_sha256']:
            raise RuntimeError('Committed novelty catalog changed')
        if identity_set(self.catalog()) != identity_set(ev.read_json(catalog_path)):
            raise RuntimeError('New public identities appeared after commitment; retain/abandon batch and commit a new one')
        if any((directory/name).exists() for name in ('reveal_started.json', 'holdout_cases.json', 'authorization.json')):
            raise RuntimeError('Cannot freeze after an authorization or reveal transition')
        # Inspect metadata only. Do not open/read/hash private seed here.
        metadata = (directory/'private_seed.bin').lstat()
        if not stat.S_ISREG(metadata.st_mode) or stat.S_IMODE(metadata.st_mode) != 0o600 or metadata.st_size != 32:
            raise RuntimeError('Private seed file must be a regular 0600 file of 32 bytes')
        sources = self.source_provider()
        acceptance, regression_artifacts = self.acceptance_regression(sources)
        once(directory/'freeze_started.json', {'batch_id': batch_id, 'utc': now(), 'retry_permitted': False})
        history = self.history_provider()
        once(directory/'historical_preservation.json', history)
        artifacts = [self.data/'protocol.json', directory/'known_catalog.json', directory/'commitment.json',
                     directory/'historical_preservation.json', self.root/'EVALUATION.md'] + regression_artifacts
        if (self.data/'public_test_cases.json').exists():
            artifacts.append(self.data/'public_test_cases.json')
        snapshot = self.archive_sources(directory, sources, artifacts)
        artifacts.append(snapshot)
        manifest = {'batch_id': batch_id, 'frozen_utc': now(), 'sources': sources,
                    'artifacts': {str(path.resolve()): ev.file_digest(path) for path in artifacts},
                    'protocol_sha256': LOCKED_PROTOCOL_SHA256,
                    'acceptance_regression': acceptance,
                    'source_snapshot_archive': str(snapshot.resolve()),
                    'private_seed_read': False, 'reveal_authorization_required': True}
        once(directory/'freeze.json', manifest)
        return {'batch_id': batch_id, 'freeze_sha256': ev.file_digest(directory/'freeze.json'),
                'source_files': len(sources), 'seed_read': False}

    def validate_freeze(self, batch_id):
        directory = self.batch(batch_id)
        self.assert_not_closed(directory)
        self.protocol()
        frozen = ev.read_json(directory/'freeze.json')
        if frozen['batch_id'] != batch_id or frozen['sources'] != self.source_provider():
            raise RuntimeError('Frozen source set or content changed')
        for path, sha256 in frozen['artifacts'].items():
            if ev.file_digest(path) != sha256:
                raise RuntimeError('Frozen artifact changed: '+path)
        history = self.history_provider()
        saved_history = ev.read_json(directory/'historical_preservation.json')
        for key in ('runtime_checked_manifest_files', 'runtime_manifest_sha256', 'runtime_preservation', 'v317_frozen_math'):
            if history[key] != saved_history[key]:
                raise RuntimeError('Historical preservation result changed')
        return frozen

    def authorization_template(self, batch_id):
        directory = self.batch(batch_id)
        self.validate_freeze(batch_id)
        return {'authorized': False, 'authorized_by': 'root', 'action': 'reveal_once',
                'batch_id': batch_id, 'freeze_sha256': ev.file_digest(directory/'freeze.json'),
                'commitment_sha256': ev.file_digest(directory/'commitment.json'),
                'protocol_sha256': LOCKED_PROTOCOL_SHA256,
                'statement': 'Root must explicitly authorize exactly one reveal after reviewing this frozen batch.'}

    def reveal(self, batch_id):
        directory = self.batch(batch_id)
        self.validate_freeze(batch_id)
        expected = self.authorization_template(batch_id)
        authorization = ev.read_json(directory/'authorization.json')
        for key in ('authorized_by', 'action', 'batch_id', 'freeze_sha256', 'commitment_sha256', 'protocol_sha256'):
            if authorization.get(key) != expected[key]:
                raise RuntimeError('Authorization does not bind this exact frozen batch')
        if authorization.get('authorized') is not True:
            raise RuntimeError('Explicit root authorization is required before seed read')
        catalog = ev.read_json(directory/'known_catalog.json')
        if identity_set(self.catalog()) != identity_set(catalog):
            raise RuntimeError('Public task identities changed since commitment')
        once(directory/'reveal_started.json', {'batch_id': batch_id, 'started_utc': now(),
            'authorization_sha256': ev.file_digest(directory/'authorization.json'),
            'freeze_sha256': ev.file_digest(directory/'freeze.json'),
            'retry_permitted': False})
        try:
            fd = os.open(directory/'private_seed.bin', os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0))
            with os.fdopen(fd, 'rb') as stream:
                metadata = os.fstat(stream.fileno())
                if not stat.S_ISREG(metadata.st_mode) or stat.S_IMODE(metadata.st_mode) != 0o600:
                    raise RuntimeError('Private seed metadata changed')
                secret = stream.read(33)
            try:
                commitment = ev.read_json(directory/'commitment.json')
                if len(secret) != 32 or hashlib.sha256(secret).hexdigest() != commitment['sha256']:
                    raise RuntimeError('Seed commitment mismatch')
                cases, audit = draw_cases(secret, self.protocol(), catalog)
            finally:
                del secret
            once(directory/'holdout_cases.json', cases)
            once(directory/'sampling_audit.json', audit)
            once(directory/'revealed.json', {'batch_id': batch_id, 'revealed_utc': now(),
                'holdout_sha256': ev.file_digest(directory/'holdout_cases.json'),
                'sampling_audit_sha256': ev.file_digest(directory/'sampling_audit.json'),
                'case_count': len(cases), 'scope': 'new parameters in known mathematical families only'})
            return {'batch_id': batch_id, 'case_count': len(cases),
                    'holdout_path': str(directory/'holdout_cases.json'), 'seed_disclosed': False}
        except BaseException as error:
            once(directory/'reveal_failed.json', {'error_type': type(error).__name__, 'reason': str(error),
                                                'retry_permitted': False, 'utc': now()})
            raise

    def run(self, batch_id):
        directory = self.batch(batch_id)
        self.validate_freeze(batch_id)
        revealed = ev.read_json(directory/'revealed.json')
        if ev.file_digest(directory/'holdout_cases.json') != revealed['holdout_sha256']:
            raise RuntimeError('Revealed cases changed')
        if ev.file_digest(directory/'sampling_audit.json') != revealed['sampling_audit_sha256']:
            raise RuntimeError('Sampling audit changed')
        cases = ev.read_json(directory/'holdout_cases.json')
        if len(cases) != 20 or len({ev.instance_key(c['task']) for c in cases}) != 20:
            raise RuntimeError('Expected exactly twenty distinct instances')
        if {ev.instance_key(c['task']) for c in cases} & identity_set(ev.read_json(directory/'known_catalog.json')):
            raise RuntimeError('Heldout cases overlap precommitted public identities')
        once(directory/'run_started.json', {'batch_id': batch_id, 'started_utc': now(),
                                          'repetitions': 5, 'methods': list(ev.METHODS), 'retry_permitted': False})
        started, rows = time.monotonic(), []
        error = None
        try:
            for repetition in range(5):
                ordered = cases[repetition:]+cases[:repetition]
                methods = ev.METHODS[repetition % 2:]+ev.METHODS[:repetition % 2]
                for item in ordered:
                    for method in methods:
                        rows.append(self.worker_runner(method, item, repetition, directory/'raw'))
                        if (rows[-1]['method'] != method or rows[-1]['case_id'] != item['id'] or
                                rows[-1]['repetition'] != repetition or rows[-1]['task'] != item['task']):
                            raise RuntimeError('Worker result does not match its assigned original case')
                        ev.write_json(directory/'run_checkpoint.json', {'rows': rows, 'complete': False})
        except BaseException as caught:
            error = {'type': type(caught).__name__, 'reason': str(caught)}
        source_valid = True
        try:
            self.validate_freeze(batch_id)
        except Exception as caught:
            source_valid = False
            error = error or {'type': type(caught).__name__, 'reason': str(caught)}
        summary = ev.summarize(rows, cases, 5)
        complete = len(rows) == 200 and error is None
        for method in ev.METHODS:
            if method in summary:
                summary[method]['formal_acceptance'] = bool(complete and source_valid and summary[method]['robust_successes'] >= 16)
                summary[method]['reason'] = 'locked twenty-case independent parameter batch; five cold repeats, all failures retained'
        passed = bool(complete and source_valid and summary.get('current', {}).get('robust_successes', 0) >= 16)
        report = {'batch_id': batch_id, 'phase': 'independent_holdout', 'repetitions': 5,
            'cases': cases, 'rows': rows, 'summary': summary, 'complete': complete,
            'formal_acceptance': passed, 'source_unchanged': source_valid, 'error': error,
            'whole_evaluation_wall_seconds': time.monotonic()-started,
            'freeze_sha256': ev.file_digest(directory/'freeze.json'),
            'holdout_sha256': revealed['holdout_sha256'], 'prior_failed_batches_retained': True,
            'scope': 'same known mathematical families; no new mathematical structure claim',
            'actual_model_tokens': None, 'token_proxy_used': False,
            'statistics': 'five timing repeats are not independent mathematical cases; descriptive shared-host timings only'}
        once(directory/'formal_report.json', report)
        return {'report_path': str(directory/'formal_report.json'), 'formal_acceptance': passed,
                'robust_successes': summary.get('current', {}).get('robust_successes', 0), 'cases': 20}

    def abandon(self, batch_id, reason):
        directory = self.batch(batch_id)
        if not reason.strip():
            raise ValueError('Abandonment must have an audit reason')
        once(directory/'abandoned.json', {'batch_id': batch_id, 'utc': now(), 'reason': reason,
                                         'files_deleted': False, 'revealed_cases_remain_public': True})
        return {'batch_id': batch_id, 'abandoned': True, 'files_deleted': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=('commit', 'freeze', 'authorization-template', 'reveal', 'run', 'abandon', 'selftest'))
    parser.add_argument('--batch')
    parser.add_argument('--reason')
    args = parser.parse_args()
    if args.command == 'selftest':
        from test_workflow import run_tests
        result = run_tests()
        ev.write_json(ev.DATA/'workflow_selftest.json', result)
        print(ev.canonical(result))
        raise SystemExit(0 if result['passed'] else 1)
    if not args.batch:
        parser.error('--batch is required')
    workflow = Workflow()
    if args.command == 'authorization-template':
        result = workflow.authorization_template(args.batch)
    elif args.command == 'abandon':
        result = workflow.abandon(args.batch, args.reason or '')
    else:
        result = getattr(workflow, args.command)(args.batch)
    print(ev.canonical(result))


if __name__ == '__main__':
    main()
