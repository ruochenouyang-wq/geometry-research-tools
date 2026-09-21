#!/usr/bin/env python3
"""Replay v308..317 byte-identical snapshots in their original relative layout.

Usage: python3 -B snapshot_logs/replay_portable.py --python /path/to/python3.9
Only audit output and temporary staging directories are written. The initial
SNAPSHOT_REPLAY.json and its logs are never overwritten.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

CHILD = r'''
import hashlib,importlib,importlib.util,json,sys,unittest
from pathlib import Path
stage=Path(sys.argv[1]).resolve();shared=Path(sys.argv[2]).resolve();version=int(sys.argv[3]);expected=sys.argv[4]
sys.dont_write_bytecode=True
sys.path[:0]=[str(stage),str(shared)]
subject=importlib.import_module('interaction')
assert Path(subject.__file__).resolve()==stage/'interaction.py', 'Wrong production module origin'
source_hash=hashlib.sha256(Path(subject.__file__).read_bytes()).hexdigest()
assert source_hash==expected, 'Staged source is not the recorded snapshot'
sys.path[:]=[str(stage),str(shared)]+[p for p in sys.path if p not in (str(stage),str(shared))]
spec=importlib.util.spec_from_file_location('portable_snapshot_test_v'+str(version),stage/'test_interaction.py')
test=importlib.util.module_from_spec(spec);sys.modules[spec.name]=test;spec.loader.exec_module(test)
result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(test))
summary={'version':version,'python_version':sys.version,'production_module':str(Path(subject.__file__).resolve()),'production_sha256':source_hash,'test_module':str(Path(test.__file__).resolve()),'copied_snapshot_module_loaded':True,'tests_run':result.testsRun,'failures':len(result.failures),'errors':len(result.errors),'skipped':len(result.skipped),'skip_reasons':[{'test':str(t),'reason':r} for t,r in result.skipped],'successful':result.wasSuccessful()}
print('PORTABLE_RESULT_JSON='+json.dumps(summary,ensure_ascii=False),flush=True)
sys.exit(0 if result.wasSuccessful() else 1)
'''


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--python', default='/Library/Developer/CommandLineTools/usr/bin/python3',
                        help='Python 3.9 interpreter compatible with the historical token dependency')
    args = parser.parse_args()
    root = Path(__file__).resolve().parent.parent
    project = root.parents[1]
    logs = root/'snapshot_logs'
    logs.mkdir(exist_ok=True)
    runtime = subprocess.check_output([args.python, '-B', '-c',
        'import platform,sys;print(platform.python_version());print(sys.executable)'], text=True).splitlines()
    if not runtime[0].startswith('3.9.'):
        raise SystemExit('Historical portable replay requires Python 3.9; received '+runtime[0])
    initial = root/'SNAPSHOT_REPLAY.json'
    initial_hash = sha(initial)
    frozen_names = ('geometry_precision_followup','geometry_v36_v90','geometry_v91_v94',
                    'geometry_cycles_v95_v194','geometry_retests_v195_v234',
                    'geometry_generalization_notes','geometry_general_v235_v237')
    links = {name: root.parent/name for name in frozen_names}
    token_cache = root/'token_cache'
    token_deps = project/'work'/'token-deps'
    for path in list(links.values())+[token_cache, token_deps]:
        if not path.is_dir():
            raise SystemExit('Missing existing read dependency: '+str(path))
    started = time.perf_counter()
    report = {'format':'geometry_snapshot_portable_replay_v1',
              'started_at_utc':datetime.now(timezone.utc).isoformat(),
              'historical_python_executable':runtime[1], 'historical_python_version':runtime[0],
              'per_snapshot_timeout_seconds':20, 'sequential_child_processes':True,
              'initial_replay_report':'SNAPSHOT_REPLAY.json', 'prior_harness_failure_report':'snapshot_logs/portable_environment_attempt1.json', 'initial_report_sha256':initial_hash,
              'method':'Copy recorded production and test bytes into a temporary project/outputs layout; link existing historical fixture directories, token_cache and work/token-deps as read dependencies. Each version runs in a fresh process.',
              'limitations':['This historical compatibility run uses Python 3.9 and is distinct from final-framework bundled Python 3.12 tests.',
                             'Only the versioned interaction production module and tests are historical snapshots; shared backend and other current modules are not reconstructed historically.',
                             'Filesystem links are used for reading existing dependencies; this is not an operating-system read-only sandbox.',
                             'Test counts repeat cumulative cases across versions and must not be added as independent test coverage.',
                             'Elapsed times describe audit execution, not algorithm speed or model usage.'],
              'read_dependency_links':{**{k:str(v) for k,v in links.items()},
                                       'token_cache':str(token_cache),'work/token-deps':str(token_deps)},
              'results':[]}
    for version in range(308,318):
        source = root/'iterations'/('v'+str(version))/'source'
        record = json.loads((source.parent/'record.json').read_text())
        before = {name:sha(source/name) for name in record['source_sha256']}
        if before != record['source_sha256']:
            raise AssertionError('Original snapshot mismatch: v'+str(version))
        row = {'version':version,'source_sha256':before,'source_matches_original_record':True}
        with tempfile.TemporaryDirectory(prefix='geometry-portable-v'+str(version)+'-',dir='/tmp') as temp:
            stage_project = Path(temp)
            outputs = stage_project/'outputs';outputs.mkdir()
            stage = outputs/root.name;stage.mkdir()
            for name in before:
                target = stage/name;target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(source/name,target)
            copied = {name:sha(stage/name) for name in before}
            row['copied_sha256'] = copied
            row['copied_bytes_match_snapshot'] = copied == before
            assert copied == before
            for name,target in links.items():
                (outputs/name).symlink_to(target,target_is_directory=True)
            (stage/'token_cache').symlink_to(token_cache,target_is_directory=True)
            (stage_project/'work').mkdir()
            (stage_project/'work'/'token-deps').symlink_to(token_deps,target_is_directory=True)
            env = dict(os.environ,PYTHONDONTWRITEBYTECODE='1',
                       PYTHONPATH=str(stage)+os.pathsep+str(root),TIKTOKEN_CACHE_DIR=str(token_cache))
            command=[args.python,'-B','-c',CHILD,str(stage),str(root),str(version),before['interaction.py']]
            t0=time.perf_counter()
            try:
                process=subprocess.run(command,cwd=stage,env=env,capture_output=True,text=True,timeout=20)
                output=process.stdout+'\n'+process.stderr
                row['exit_code']=process.returncode
                marker=[line[len('PORTABLE_RESULT_JSON='):] for line in process.stdout.splitlines()
                        if line.startswith('PORTABLE_RESULT_JSON=')]
                if marker:
                    row.update(json.loads(marker[-1]))
                row['status']='passed' if process.returncode==0 and row.get('copied_snapshot_module_loaded') else 'failed'
            except subprocess.TimeoutExpired as error:
                def decode(value):
                    return value.decode(errors='replace') if isinstance(value,bytes) else (value or '')
                output=decode(error.stdout)+'\n'+decode(error.stderr)+'\nTIMEOUT: 20 seconds.\n'
                row.update(status='timeout',exit_code=None)
            row['wall_seconds']=time.perf_counter()-t0
            row['staged_source_unchanged']={name:sha(stage/name) for name in before}==copied
        row['original_snapshot_unchanged']={name:sha(source/name) for name in before}==before
        log=logs/('portable_v'+str(version)+'.log');log.write_text(output,encoding='utf-8')
        row['log']=str(log.relative_to(root))
        report['results'].append(row)
        print(version,row['status'],'tests',row.get('tests_run'),'errors',row.get('errors'),
              'skipped',row.get('skipped'),flush=True)
        report['counts']={status:sum(r['status']==status for r in report['results'])
                          for status in ('passed','failed','timeout')}
        report['complete']=len(report['results'])==10
        report['initial_report_unchanged']=sha(initial)==initial_hash
        report['elapsed_wall_seconds']=time.perf_counter()-started
        report['summary']={'test_invocations_including_cumulative_repeats':sum(r.get('tests_run',0) for r in report['results']),
                           'assertion_failures':sum(r.get('failures',0) for r in report['results']),
                           'test_errors':sum(r.get('errors',0) for r in report['results']),
                           'skipped_tests':sum(r.get('skipped',0) for r in report['results']),
                           'all_original_and_copied_sources_unchanged':all(r['original_snapshot_unchanged'] and r['staged_source_unchanged'] for r in report['results'])}
        (root/'SNAPSHOT_REPLAY_PORTABLE.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    if not report['initial_report_unchanged']:
        raise AssertionError('Initial failure report changed during portable replay')
    print('FINAL',report['counts'],report['summary'],flush=True)
    return 0 if report['counts']['passed']==10 and report['summary']['skipped_tests']==0 else 1


if __name__=='__main__':
    sys.exit(main())
