"""The completion ledger is driven by saved evidence, never version labels alone."""
from pathlib import Path
from datetime import datetime, timezone
import json
import hashlib
from common import ROOT, save


def stamp():
    return datetime.now(timezone.utc).isoformat()


def progress():
    return json.loads((ROOT/'PROGRESS.json').read_text())


def begin(cycle, problem, weakness, baseline):
    state = progress()
    if cycle != state['completed_cycles']+1 or not 1 <= cycle <= 10:
        raise ValueError('Each cycle starts after the preceding cycle has completed')
    folder = ROOT/f'round{cycle:02d}'
    folder.mkdir(exist_ok=True)
    save(f'round{cycle:02d}/task.json', {'cycle':cycle, 'first_version':95+10*(cycle-1),
         'last_version':104+10*(cycle-1), 'problem':problem, 'observed_weakness':weakness,
         'baseline_evidence':baseline, 'began_utc':stamp()})
    state['active_cycle'] = cycle
    state['status'] = 'running'
    save('PROGRESS.json', state)


def finish(cycle, stages, validation, files):
    """Called only after the root has reviewed the actual mathematical changes."""
    state = progress()
    if cycle != state['completed_cycles']+1:
        raise ValueError('No skipped or repeated completion')
    first = 95+10*(cycle-1)
    if len(stages) != 10 or [s['version'] for s in stages] != list(range(first,first+10)):
        raise ValueError('Exactly ten consecutive reviewed increments required')
    for entry in stages:
        if not entry.get('capability') or not entry.get('callable') or not entry.get('evidence'):
            raise ValueError('An increment needs an executable capability and saved evidence')
        for path in entry['evidence']:
            if not (ROOT/path).is_file():
                raise ValueError('Missing stage evidence: '+path)
    if not validation.get('all_passed'):
        raise ValueError('A cycle cannot complete before its acceptance checks pass')
    hashes = {path:hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in files}
    record = {'cycle':cycle, 'completed_utc':stamp(), 'first_version':first,
              'last_version':first+9, 'stages':stages, 'validation':validation, 'file_sha256':hashes}
    save(f'round{cycle:02d}/COMPLETED.json', record)
    state['cycles'].append({'cycle':cycle,'record':f'round{cycle:02d}/COMPLETED.json'})
    state['completed_cycles'] = cycle
    state['completed_versions'] = 10*cycle
    state['active_cycle'] = None if cycle == 10 else cycle+1
    state['status'] = 'complete' if cycle == 10 else 'running'
    save('PROGRESS.json', state)
    return state
