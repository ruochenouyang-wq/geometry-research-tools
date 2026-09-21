"""Independently inspect iteration snapshots and the immutable baseline."""
from pathlib import Path
import argparse
import hashlib
import json
from backend import ROOT, PROJECT, source_hashes
from campaign import list_records


def audit_records(require_complete=False):
    records = list_records()
    failures = []
    versions = [r.get('version') for r in records]
    if len(versions) != len(set(versions)):
        failures.append('duplicate_version')
    if any(type(v) is not int or not 238 <= v <= 317 for v in versions):
        failures.append('invalid_version')
    if require_complete and versions != list(range(238, 318)):
        failures.append('missing_iterations')
    previous = {}
    for r in records:
        v = r['version']
        changed_production = False
        if not all(r.get(k) for k in ('hypothesis','change','validation','comparison')):
            failures.append(f'v{v}:incomplete_evidence')
        for name, expected in r.get('source_sha256', {}).items():
            snapshot = ROOT/'iterations'/f'v{v}'/'source'/name
            if not snapshot.is_file():
                failures.append(f'v{v}:missing_snapshot:{name}')
                continue
            actual = hashlib.sha256(snapshot.read_bytes()).hexdigest()
            if actual != expected:
                failures.append(f'v{v}:snapshot_hash_mismatch:{name}')
            production = Path(name).suffix == '.py' and not Path(name).name.startswith('test_')
            if production and previous.get(name) != actual:
                changed_production = True
            previous[name] = actual
        if not changed_production:
            failures.append(f'v{v}:no_distinct_production_code_change')
    return {'recorded_iterations':len(records),'versions':versions,
            'complete':versions==list(range(238,318)),
            'failures':failures,'passed':not failures,
            'scope':'Checks provenance and actual source changes; mathematical and performance validation is separate.'}


def audit_frozen():
    saved=json.loads((ROOT/'FROZEN_BASELINE.json').read_text())['files']
    current=source_hashes()
    changes=[name for name,value in saved.items() if current.get(name)!=value]
    additions=sorted(set(current)-set(saved))
    return {'checked_files':len(saved),'changed_or_missing':changes,'new_files_in_frozen_roots':additions,
            'passed':not changes and not additions}


def code_fingerprint():
    return {p.name:hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(ROOT.glob('*.py'))}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--complete',action='store_true')
    args=parser.parse_args()
    result={'iterations':audit_records(args.complete),'frozen':audit_frozen()}
    result['passed']=result['iterations']['passed'] and result['frozen']['passed']
    print(json.dumps(result,ensure_ascii=False,indent=2))
    raise SystemExit(0 if result['passed'] else 1)


if __name__=='__main__':main()
