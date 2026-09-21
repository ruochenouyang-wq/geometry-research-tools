"""Reproducible profiling on named, already public development examples only.

Writes metrics/profile.json. No private seed or unrevealed corpus is consulted.
Each service/case/repetition starts a fresh Python process. First-call costs,
including lazy imports and cache population, remain charged and visible.
"""
from time import perf_counter, process_time
BOOT_WALL, BOOT_CPU = perf_counter(), process_time()
import argparse
import cProfile
import hashlib
import importlib
import json
from pathlib import Path
import platform
import pstats
import resource
import statistics
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
PREVIOUS = ROOT.parent / 'geometry_strategy_v238_v317'
sys.dont_write_bytecode = True


def cases():
    answer = {}
    for pointer, selected in [('latest_development.json', {'D05', 'D06'}),
                              ('latest_holdout.json', {'H04', 'H05', 'H06'})]:
        entry = json.loads((PREVIOUS / 'evaluation' / pointer).read_text())
        report = Path(entry['report_path'])
        # Resolve by the published filename, allowing checkout relocation.
        report = PREVIOUS / 'evaluation' / report.name
        for item in json.loads(report.read_text())['cases']:
            if item['id'] in selected:
                answer[item['id']] = item['task']
    for item in json.loads((PREVIOUS / 'metrics' / 'integrated_challenges.json').read_text()):
        if item['name'] in ('original_function_precision', 'original_spectral_precision'):
            answer[item['name']] = item['task']
    return answer


def clocked(function):
    wall, cpu = perf_counter(), process_time()
    result = function()
    return result, {'cpu_seconds': process_time() - cpu,
                    'wall_seconds': perf_counter() - wall}


def summary(result):
    certificate = result.get('certificate')
    content = json.dumps(certificate, sort_keys=True, separators=(',', ':')).encode()
    return {key: result.get(key) for key in
            ('ok', 'certificate_valid', 'target_met', 'within_budget', 'status',
             'value', 'route', 'cost', 'total_elapsed_seconds')} | {
        'certificate_sha256': hashlib.sha256(content).hexdigest(),
        'certificate_bytes': len(content),
        'certificate_format': certificate.get('format') if certificate else None,
        'attempt_count': len(result.get('attempts', []))}


def worker(args):
    bootstrap = {'cpu_seconds': process_time() - BOOT_CPU,
                 'wall_seconds': perf_counter() - BOOT_WALL}
    sys.path.insert(0, str(PREVIOUS))
    sys.path.insert(0, str(ROOT))
    module, import_cost = clocked(lambda: importlib.import_module(args.service))
    task = cases()[args.case]
    with tempfile.TemporaryDirectory(prefix='geometry-profile-') as folder:
        service, setup_cost = clocked(lambda: module.ResearchService(evidence_root=folder, cache=False))
        profile = cProfile.Profile(timer=process_time) if args.profile else None
        rounds = 1 if profile else 2 if args.audit_replays else 1 + args.warm_repeats
        rows = []
        replay_observations = []
        for index in range(rounds):
            direct = importlib.import_module('backend').direct
            original_verify_full = direct.verify_full
            if args.audit_replays and index == 1:
                def observed_verify_full(certificate, *bindings, **named_bindings):
                    # This observer never skips or substitutes verification.
                    text = json.dumps([certificate, bindings, named_bindings],
                                      sort_keys=True, separators=(',', ':'), allow_nan=False)
                    begin = process_time()
                    valid = original_verify_full(certificate, *bindings, **named_bindings)
                    elapsed = process_time() - begin
                    replay_observations.append({'full_call_sha256': hashlib.sha256(text.encode()).hexdigest(),
                                                'frozen_replay_cpu_seconds': elapsed,
                                                'valid': valid})
                    return valid
                direct.verify_full = observed_verify_full
            if profile:
                profile.enable()
            try:
                result, elapsed = clocked(lambda: service.evaluate_task(task))
            finally:
                direct.verify_full = original_verify_full
            if profile:
                profile.disable()
            rows.append({'phase': 'first_online' if index == 0 else 'repeated_online',
                         **elapsed, 'result': summary(result)})
        # Independent external replay is reported and charged separately.
        verification = importlib.import_module('verification')
        judgment, replay_cost = clocked(lambda: verification.assess(result['certificate'], task))
        assert judgment['certificate_valid'] is True, judgment
        stats = []
        if profile:
            for (filename, line, name), (primitive, calls, own, cumulative, callers) in pstats.Stats(profile).stats.items():
                try:
                    filename = str(Path(filename).resolve().relative_to(ROOT.parent.parent))
                except ValueError:
                    pass
                stats.append({'file': filename, 'line': line, 'function': name,
                              'primitive_calls': primitive, 'calls': calls,
                              'self_cpu_seconds': own, 'cumulative_cpu_seconds': cumulative})
            stats.sort(key=lambda row: -row['cumulative_cpu_seconds'])
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        seen, redundant_cpu = set(), 0.0
        for observation in replay_observations:
            key = observation['full_call_sha256']
            if key in seen:
                redundant_cpu += observation['frozen_replay_cpu_seconds']
            seen.add(key)
        return {'service': args.service, 'case_id': args.case, 'task': task,
                'profiled': args.profile, 'stdlib_bootstrap': bootstrap,
                'service_import': import_cost, 'service_construction': setup_cost,
                'online': rows, 'external_final_replay': replay_cost,
                'external_final_judgment': judgment,
                'process_cpu_seconds_at_end': process_time(),
                'process_peak_rss_bytes': peak if sys.platform == 'darwin' else peak * 1024,
                'replay_audit': {'enabled': args.audit_replays,
                                 'scope': 'Second call only; all original checks still execute; observation overhead included online',
                                 'full_replay_calls': len(replay_observations),
                                 'unique_full_call_keys': len(seen),
                                 'redundant_frozen_replay_cpu_seconds': redundant_cpu,
                                 'observations': replay_observations},
                'profile_all_functions': stats}


def aggregate(rows):
    out = []
    for service, case_id in sorted({(r['service'], r['case_id']) for r in rows}):
        selected = [r for r in rows if r['service'] == service and r['case_id'] == case_id]
        result = {'service': service, 'case_id': case_id, 'fresh_processes': len(selected)}
        for phase in ('stdlib_bootstrap', 'service_import', 'service_construction', 'external_final_replay'):
            for metric in ('cpu_seconds', 'wall_seconds'):
                result[phase + '_' + metric + '_median'] = statistics.median(r[phase][metric] for r in selected)
        for phase in ('first_online', 'repeated_online'):
            observations = [x for r in selected for x in r['online'] if x['phase'] == phase]
            if observations:
                for metric in ('cpu_seconds', 'wall_seconds'):
                    values = [x[metric] for x in observations]
                    result[phase + '_' + metric + '_median'] = statistics.median(values)
                    result[phase + '_' + metric + '_min'] = min(values)
                    result[phase + '_' + metric + '_max'] = max(values)
        result['certificate_hashes'] = sorted({x['result']['certificate_sha256'] for r in selected for x in r['online']})
        result['all_certificates_valid'] = all(x['result']['certificate_valid'] for r in selected for x in r['online'])
        result['all_targets_met'] = all(x['result']['target_met'] for r in selected for x in r['online'])
        out.append(result)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--services', nargs='+', default=['research_service'])
    parser.add_argument('--cases', nargs='+')
    parser.add_argument('--fresh-repeats', type=int, default=3)
    parser.add_argument('--warm-repeats', type=int, default=4)
    parser.add_argument('--worker', action='store_true')
    parser.add_argument('--service')
    parser.add_argument('--case')
    parser.add_argument('--profile', action='store_true')
    parser.add_argument('--audit-replays', action='store_true')
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(args), ensure_ascii=False))
        return
    selected = args.cases or list(cases())
    unknown = set(selected) - set(cases())
    if unknown:
        parser.error('Only fixed public development cases may be profiled: ' + str(unknown))
    runs, profiles, audits = [], [], []
    wall, cpu = perf_counter(), process_time()
    for repeat in range(args.fresh_repeats + 2):
        for case_id in selected:
            # Rotate the service order to limit order bias without an RNG or seed.
            order = args.services[repeat % len(args.services):] + args.services[:repeat % len(args.services)]
            for service in order:
                command = [sys.executable, '-B', str(Path(__file__).resolve()), '--worker',
                           '--service', service, '--case', case_id,
                           '--warm-repeats', str(args.warm_repeats)]
                profiled = repeat == args.fresh_repeats
                audited = repeat == args.fresh_repeats + 1
                if profiled:
                    command.append('--profile')
                if audited:
                    command.append('--audit-replays')
                process_wall = perf_counter()
                completed = subprocess.run(command, text=True, capture_output=True, check=True)
                row = json.loads(completed.stdout)
                row['subprocess_wall_seconds'] = perf_counter() - process_wall
                row['fresh_repeat_index'] = repeat
                (audits if audited else profiles if profiled else runs).append(row)
                print(f'{service} {case_id} {"audit" if audited else "profile" if profiled else "timing"} {repeat + 1}', flush=True)
    report = {
        'format': 'geometry_release_public_development_profile_v1',
        'classification': 'public development only, including formerly held-out H04-H06',
        'private_seed_or_unrevealed_cases_read': False,
        'environment': {'python': sys.version, 'executable': sys.executable,
                        'platform': platform.platform(), 'machine': platform.machine()},
        'protocol': {'case_order': selected, 'service_order': args.services,
                     'fresh_repeats': args.fresh_repeats, 'warm_repeats': args.warm_repeats,
                     'online_adapter': 'ResearchService.evaluate_task; no handle/storage claims',
                     'preparation': 'stdlib bootstrap, module import, constructor separately measured; lazy imports and warmup charged in first_online',
                     'replay': 'external final replay separately charged once per worker; internal replay included in online',
                     'profile': 'Separate fresh-process first-call cProfile using process_time; overlapping cumulative costs are not additive; profiled timings are not speed claims',
                     'replay_audit': 'Separate process first call charged as warmup; second call observes frozen verify_full and exact whole-certificate plus binding call keys, without skipping any check. Repeated-call CPU is an optimization upper bound, not measured savings.',
                     'unprofiled': 'Process CPU and wall; no profiler; first and repeated online shown separately',
                     'machine_isolation': 'Not isolated; other agents/processes may run; descriptive measurements only'},
        'aggregates': aggregate(runs), 'timing_runs': runs, 'profile_runs': profiles,
        'replay_audit_runs': audits,
        'parent_cpu_seconds': process_time() - cpu,
        'experiment_wall_seconds': perf_counter() - wall,
        'total_worker_cpu_seconds': sum(r['process_cpu_seconds_at_end'] for r in runs + profiles + audits)}
    destination = ROOT / 'metrics' / 'profile.json'
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2, ensure_ascii=False) + '\n')
    print(destination)


if __name__ == '__main__':
    main()
