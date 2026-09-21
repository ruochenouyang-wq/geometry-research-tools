"""Read-only release verification: versions, file integrity and challenge proofs.

Run tests separately with the commands in TOOL_GUIDE.md. This script never
regenerates evidence, edits completion records, or contacts a language model.
"""
import hashlib
import json
from pathlib import Path
from common import ROOT
import portal


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_manifest(directory):
    data = json.loads((directory/'MANIFEST.json').read_text())
    for row in data['entries']:
        path = directory/row['path']
        if not path.is_file() or file_hash(path) != row['sha256']:
            raise ValueError('Missing or changed release file: '+str(path))
    return len(data['entries'])


def run():
    progress = json.loads((ROOT/'PROGRESS.json').read_text())
    if progress['status'] != 'complete' or progress['completed_cycles'] != 10:
        raise ValueError('Ten cycles have not all completed')
    versions = []
    for i in range(1,11):
        folder = ROOT/f'round{i:02d}'
        completed = json.loads((folder/'COMPLETED.json').read_text())
        if not completed['validation']['all_passed']:
            raise ValueError('A cycle has no successful acceptance')
        for stage in completed['stages']:
            versions.append(stage['version'])
            for relative in stage['evidence']:
                if not (ROOT/relative).is_file():
                    raise ValueError('Missing stage evidence: '+relative)
        for relative, previous_hash in completed['file_sha256'].items():
            if file_hash(ROOT/relative) != previous_hash:
                recheck = json.loads((folder/'FINAL_REAUDIT.json').read_text())
                if (not recheck['all_passed'] or recheck['current_file_sha256'].get(relative)
                        != file_hash(ROOT/relative)):
                    raise ValueError('Unreviewed post-acceptance change: '+relative)
    if versions != list(range(95,195)):
        raise ValueError('The one hundred executable increments are not consecutive')
    report = {'versions':len(versions),'cycles':10,'release_files':check_manifest(ROOT),
              'frozen_dependency_files':sum(check_manifest(ROOT.parent/name) for name in
                  ('geometry_v36_v90','geometry_v91_v94'))}
    book_path = ROOT/'challenges.json'
    book = json.loads(book_path.read_text())
    results = json.loads((ROOT/'cross_challenge/result.json').read_text())
    if file_hash(book_path) != results['challenge_sha256']:
        raise ValueError('Frozen challenge book changed')
    count = 0
    if [c['id'] for c in results['results']] != [c['id'] for c in book['cases']]:
        raise ValueError('Missing or reordered challenge case')
    for original, case in zip(book['cases'],results['results']):
        if original['input'] != case['original_input'] or original['acceptance'] != case['original_acceptance']:
            raise ValueError('Challenge weakened: '+original['id'])
        for row in case['certificates']:
            path = ROOT/'cross_challenge'/row['path']
            if file_hash(path) != row['sha256'] or not portal.verify(json.loads(path.read_text())):
                raise ValueError('Challenge proof failed: '+str(path))
            count += 1
    report.update(challenges_checked=20,challenge_proofs_replayed=count,all_passed=True)
    return report


if __name__ == '__main__':
    print(json.dumps(run(),ensure_ascii=False,indent=2))
