"""Check preserved work, test logs and the predeclared release outcome."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

from release_backend import ROOT,PREVIOUS
from release_evaluation import verify_old_preservation,checked_manifest
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
    matched=[name for name,v in methods.items() if name!='release' and
             v['robust_successes']==release['robust_successes']]
    runtime=bool(matched) and all(
        not comparisons[name]['robust_success_losses'] and
        comparisons[name]['median_full_session_wall_ratio_release_over_baseline']<=0.8 and
        comparisons[name]['median_full_session_cpu_ratio_release_over_baseline']<=0.8
        for name in matched)
    no_loss=all(not value['robust_success_losses'] for value in comparisons.values())
    return {'stable_successes':release['robust_successes'],'cases':release['cases'],
            'reliable_on_declared_probes':reliable,'no_baseline_success_lost':no_loss,
            'accuracy_gate':accuracy,'runtime_gate':runtime,
            'same_success_baselines':matched,'passed':reliable and no_loss and (accuracy or runtime),
            'scope':'Predeclared finite suite only; no model superiority or universal dominance claim.'}


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
        result['freeze_files']=len(manifest['files'])
        result['heldout']=outcome(report)
    result['passed']=(result['iterations']['passed'] and result['historical_math_files']['passed']
       and all(x['returncode']==0 for x in result.get('tests',[]))
       and result.get('heldout',{}).get('passed',True))
    target=ROOT/('FINAL_AUDIT.json' if args.final else 'PRE_FREEZE_AUDIT.json')
    target.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'report':str(target),'passed':result['passed'],'heldout':result.get('heldout')},indent=2))
    raise SystemExit(0 if result['passed'] else 1)


if __name__=='__main__':main()
