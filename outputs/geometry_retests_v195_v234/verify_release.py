"""Read-only file integrity, original bindings and certificate replay.

Does not rerun optimization, alter frozen releases or refresh timing records.
The 214-test implementation/independent review run is recorded separately.
"""
import json
from pathlib import Path

from support import ROOT, PREVIOUS, source_hash, original_case
from acceptance import previous_integrity, replay_previous, stages
from research_tool import Store, verify


def check():
    manifest = json.loads((ROOT/'MANIFEST.json').read_text())
    for row in manifest['entries']:
        path = ROOT/row['path']
        assert not path.is_symlink()
        assert path.stat().st_size == row['bytes'], row['path']
        assert source_hash(path) == row['sha256'], row['path']
    dependencies = []
    for dependency in manifest['dependencies']:
        folder = ROOT.parent/dependency['directory']
        assert source_hash(folder/'MANIFEST.json') == dependency['manifest_sha256']
        entries = json.loads((folder/'MANIFEST.json').read_text())['entries']
        for row in entries:
            assert source_hash(folder/row['path']) == row['sha256'], row['path']
        dependencies.append({'release': dependency['directory'], 'files_checked': len(entries)})
    snapshot = previous_integrity()
    accepted = json.loads((ROOT/'ROOT_ACCEPTANCE.json').read_text())
    assert accepted['original_problem_bank_sha256'] == source_hash(PREVIOUS/'challenges.json')
    store = Store()
    original_calls = 0
    for case in accepted['cases']:
        assert case['original'] == original_case(case['case'])
        assert case['original_target_met'] is True
        for row in case['service_calls']:
            record = store.get(row['result']['certificate_id'])
            request = row['request']
            assert record['operation'] == request['op']
            assert record['arguments'] == {k:v for k,v in request.items() if k != 'op'}
            assert verify(record)
            original_calls += 1
    count = 0
    for path in store.directory.glob('*.json'):
        store.get(path.stem)
        count += 1
    prior = replay_previous()
    return {'all_passed': True, 'indexed_release_files': len(manifest['entries']),
            'dependency_manifests': dependencies, 'previous_snapshot': snapshot,
            'original_cases_bound': len(accepted['cases']), 'original_service_calls': original_calls,
            'stored_service_proofs_replayed': count,
            'previous_proofs_replayed': prior['certificates_verified'],
            'executable_increments': len(stages()),
            'optimization_search_rerun': False, 'formal_proof_assistant_checked': False}


if __name__ == '__main__':
    print(json.dumps(check(), ensure_ascii=False, indent=2, sort_keys=True))
