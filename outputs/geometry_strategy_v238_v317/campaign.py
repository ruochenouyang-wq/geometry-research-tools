"""Small, auditable iteration records; numbered steps are actual code changes."""
from pathlib import Path
import hashlib
import json
import shutil
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent


def record(version, *, hypothesis, change, files, validation, comparison,
           outcome='provisional', limitations=None):
    if type(version) is not int or not 238 <= version <= 317:
        raise ValueError('Version must be in 238..317')
    if not all((hypothesis, change, files, validation, comparison)):
        raise ValueError('Every step requires code, hypothesis, validation and comparison')
    folder = ROOT/'iterations'/f'v{version}'
    if (folder/'record.json').exists():
        raise ValueError('Iteration already recorded; append review evidence separately')
    snapshots = folder/'source'
    snapshots.mkdir(parents=True, exist_ok=True)
    hashes = {}
    for name in files:
        source = (ROOT/name).resolve()
        if not source.is_relative_to(ROOT) or source.suffix != '.py':
            raise ValueError('Snapshot files must be local Python source')
        target = snapshots/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        hashes[name] = hashlib.sha256(source.read_bytes()).hexdigest()
    entry = {'version':version,'recorded_at':datetime.now(timezone.utc).isoformat(),
             'hypothesis':hypothesis,'change':change,'source_sha256':hashes,
             'validation':validation,'comparison':comparison,'outcome':outcome,
             'limitations':limitations or [],
             'note':'A small component iteration; final system claims require integrated evaluation.'}
    (folder/'record.json').write_text(json.dumps(entry, ensure_ascii=False, indent=2, allow_nan=False)+'\n')
    return entry


def list_records():
    return [json.loads(p.read_text()) for p in sorted((ROOT/'iterations').glob('v*/record.json'))]
