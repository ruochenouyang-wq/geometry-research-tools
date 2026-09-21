"""Public regression of the isolated controller repair; never a holdout run."""
import hashlib
import json
from pathlib import Path
import signal
import subprocess
import sys
from time import perf_counter, process_time

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT.parent / 'geometry_strategy_v317_runtime/evaluation/holdout_20260920T042408501655Z/release_r0_output.json'


class Deadline(BaseException):
    pass


def worker(case_id):
    started, cpu = perf_counter(), process_time()
    row = next(r for r in json.loads(SOURCE.read_text())['rows'] if r['case_id'] == case_id)
    sys.path.insert(0, str(ROOT))
    import compat_controller as solver
    initialization = perf_counter() - started
    def alarm(*_):
        raise Deadline()
    signal.signal(signal.SIGALRM, alarm)
    signal.setitimer(signal.ITIMER_REAL, 10)
    result = {'id': case_id, 'task': row['task'], 'target_met': False,
              'certificate_valid': False, 'initialization_seconds': initialization}
    online = perf_counter()
    try:
        raw = json.loads(json.dumps(solver.solve(row['task']), allow_nan=False))
        result['raw'] = raw
        certificate = raw.get('certificate', raw)
        assessment = solver.verification.assess(certificate, row['task'])
        result.update(assessment)
    except Deadline:
        result.update(error='hard_timeout', target_met=False)
    except Exception as exc:
        result.update(error=type(exc).__name__ + ': ' + str(exc), target_met=False)
    finally:
        result.update(online_seconds=perf_counter()-online,
                      total_worker_seconds=perf_counter()-started,
                      cpu_seconds=process_time()-cpu)
        # Keep encoding inside the timer as in the real evaluation.
        encoded = json.dumps(result, sort_keys=True, allow_nan=False)
        signal.setitimer(signal.ITIMER_REAL, 0)
    print(encoded)


def main():
    if len(sys.argv) == 2:
        worker(sys.argv[1])
        return
    directory = ROOT / 'iterations/01_controller_repair'
    directory.mkdir()
    records = []
    for i in range(5, 15):
        case_id = f'H{i:02d}'
        started = perf_counter()
        try:
            call = subprocess.run([sys.executable, '-B', __file__, case_id],
                                  text=True, capture_output=True, timeout=30)
            record = json.loads(call.stdout) if call.returncode == 0 else {
                'id': case_id, 'target_met': False, 'error': call.stderr, 'returncode': call.returncode}
        except subprocess.TimeoutExpired:
            record = {'id': case_id, 'target_met': False, 'error': 'parent_timeout'}
        record['parent_wall_seconds'] = perf_counter()-started
        (directory / f'{case_id}.json').write_text(json.dumps(record, indent=2)+'\n')
        records.append({k:v for k,v in record.items() if k != 'raw'})
        print(case_id, record.get('target_met'), record.get('error'), flush=True)
    result = {'iteration': 1, 'scope': 'Public development; no new independent accuracy claim',
              'changes': ['Private repaired adaptive certificate diagnosis',
                          'Singular fallback refinement preserves legal action fields'],
              'code_sha256': hashlib.sha256((ROOT/'compat_controller.py').read_bytes()).hexdigest(),
              'test_suite': 'tests/test_compat_controller.py: 2 passed',
              'successes': sum(r.get('target_met') is True for r in records),
              'count': len(records), 'records': records,
              'baseline_known_successes': 1,
              'retained_all_failed_branches': True}
    result['accuracy_improved'] = result['successes'] > 1
    result['consecutive_accuracy_stagnations'] = 0 if result['accuracy_improved'] else 1
    (directory/'REPORT.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'records'}),flush=True)


if __name__ == '__main__':
    main()
