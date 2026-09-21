"""Check preserved work, test logs and the predeclared release outcome."""
import argparse
import json
import hashlib
from pathlib import Path
import subprocess
import sys

from release_backend import ROOT,PREVIOUS
from release_evaluation import verify_old_preservation,checked_manifest,digest,DATA
from audit_campaign import audit_records,audit_frozen


def outcome(report):
    results=report['results']
    methods=results['methods']
    release=methods['release']
    comparisons=results['release_against_every_baseline']
    reliable=all(v['false_accepts']==0 and v['probe_execution_failures']==0
                 and v['worker_failures']==0 for v in methods.values())
    best=max(v['robust_successes'] for name,v in methods.items() if name!='release')
    strongest=[name for name,v in methods.items() if name!='release' and v['robust_successes']==best]
    accuracy=(release['robust_successes']>best and
              all(not comparisons[name]['robust_success_losses'] for name in strongest))
    matched=[name for name,v in methods.items() if name not in ('release','first_release') and
             v['robust_successes']==release['robust_successes']]
    runtime=bool(matched) and all(
        not comparisons[name]['robust_success_losses'] and
        comparisons[name]['median_full_session_wall_ratio_release_over_baseline']<=0.8 and
        comparisons[name]['median_full_session_cpu_ratio_release_over_baseline']<=0.8
        for name in matched)
    no_loss=all(not value['robust_success_losses'] for value in comparisons.values())
    first=nonregression(comparisons['first_release'])
    return {'stable_successes':release['robust_successes'],'cases':release['cases'],
            'reliable_on_declared_probes':reliable,'no_baseline_success_lost':no_loss,
            'accuracy_gate':accuracy,'runtime_gate':runtime,
            'same_success_baselines':matched,'first_release_nonregression':first,
            'passed':reliable and no_loss and first and runtime,
            'scope':'Predeclared finite suite only; no model superiority or universal dominance claim.'}


def nonregression(comparison):
    return not comparison['robust_success_losses'] and all(comparison[key]<=1.0 for key in (
        'median_full_session_wall_ratio_release_over_baseline',
        'median_full_session_cpu_ratio_release_over_baseline',
        'total_experiment_wall_ratio_including_preparation',
        'total_experiment_cpu_ratio_including_preparation'))


def cold_regression():
    pointer=json.loads((DATA/'acceptance_regression.json').read_text())
    path=(ROOT/pointer['report_path']).resolve()
    if not path.is_relative_to(ROOT):raise ValueError('Regression path escapes release')
    if hashlib.sha256(path.read_bytes()).hexdigest()!=pointer['sha256']:
        raise ValueError('Regression evidence changed')
    report=json.loads(path.read_text())
    if report['phase']!='regression':raise ValueError('Wrong regression phase')
    results=report['results'];new=results['methods']['release']
    reliable=not new['false_accepts'] and not new['probe_execution_failures'] and not new['worker_failures']
    reference=results['release_against_every_baseline']
    return {'report':pointer['report_path'],'cases':new['cases'],
            'successes':new['robust_successes'],'reliable':reliable,
            'vs_v317':nonregression(reference['frozen_v317']),
            'vs_first_release':nonregression(reference['first_release']),
            'passed':new['cases']==6 and new['robust_successes']==6 and reliable and
                      nonregression(reference['frozen_v317']) and nonregression(reference['first_release'])}


def run_tests():
    results=[]
    for name,source in [('release',ROOT),('frozen_v317',PREVIOUS)]:
        log=ROOT/('TESTS_'+name+'.log')
        with log.open('w') as handle:
            result=subprocess.run([sys.executable,'-B','-m','unittest','discover','-s',str(source),'-v'],
                                  cwd=source,stdout=handle,stderr=subprocess.STDOUT)
        results.append({'suite':name,'returncode':result.returncode,'log':str(log.relative_to(ROOT))})
    return results


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--tests',action='store_true')
    parser.add_argument('--final',action='store_true')
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    result={'old_preservation':verify_old_preservation(),'iterations':audit_records(True),
            'historical_math_files':audit_frozen()}
    if args.tests:result['tests']=run_tests()
    if args.final:
        manifest=checked_manifest()
        if args.report is None:parser.error('--final requires --report from the released holdout')
        report=json.loads(args.report.read_text())
        if report.get('phase')!='holdout':raise ValueError('Final acceptance requires heldout evaluation')
        if report['protocol_digest']!=manifest['protocol_digest']:
            raise ValueError('Final report protocol mismatch')
        released=json.loads((DATA/'holdout_release.json').read_text())
        if released['freeze_digest']!=digest(manifest) or released['suite_digest']!=digest(report['cases']):
            raise ValueError('Report does not match released frozen suite')
        result['freeze_files']=len(manifest['files'])
        result['heldout']=outcome(report)
        result['cold_regression']=cold_regression()
    result['passed']=(result['iterations']['passed'] and result['historical_math_files']['passed']
       and all(x['returncode']==0 for x in result.get('tests',[]))
       and result.get('heldout',{}).get('passed',True)
       and result.get('cold_regression',{}).get('passed',True))
    target=ROOT/('FINAL_AUDIT.json' if args.final else 'PRE_FREEZE_AUDIT.json')
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'report':str(target),'passed':result['passed'],'heldout':result.get('heldout')},indent=2))
    raise SystemExit(0 if result['passed'] else 1)


if __name__=='__main__':main()
