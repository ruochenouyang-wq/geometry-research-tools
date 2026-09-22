"""Synthetic evaluator tests and one optional unchanged-baseline smoke.
No A/B/C mathematical search, input changes, or final comparison is run here.
"""
import argparse
from copy import deepcopy
from pathlib import Path
import json
import benchmark as b


def case(name='synthetic'):
    return {'id':name, 'dimension':'synthetic_only', 'label':'synthetic only',
            'task':{'kind':'spectrum','function':{'kind':'axis_profile','profile':'step',
             'amplitude':'-1','offset':'0','axis':2},'mean_zero':False,
             'tolerance':'1/1000000','budget':{'wall_seconds':10}}}


def run(include_baseline=False):
    destination=b.ROOT/'evaluation'/('selftest_'+b.stamp())
    destination.mkdir(parents=True)
    tests={};records={}
    for behavior in ('echo','raise','unsupported','status_only','sleep','swallow','verify_sleep','memory_sleep'):
        outcome=b.run_process({'mode':'cold','system':'baseline','case':case(),
                               'synthetic':behavior,'limit':.25 if behavior in ('sleep','swallow','verify_sleep','memory_sleep') else 2},
                              destination/behavior)
        records[behavior]=outcome
        tests[behavior]=outcome['success'] if behavior=='echo' else not outcome['success']
        tests[behavior+'_wait4_rss']=outcome['process']['peak_rss_bytes']>0
        tests[behavior+'_wait4_cpu']=outcome['process']['cpu_seconds']>=0
        if behavior in ('sleep','swallow','verify_sleep','memory_sleep'):
            tests[behavior+'_parent_killed']=outcome['process']['watchdog']=='global_deadline' and outcome['process']['returncode']==-9
    tests['partial_certificate_retained_on_verify_timeout']=records['verify_sleep'].get('certificate_json_bytes',0)>0
    tests['timeout_rss_contains_allocated_memory']=records['memory_sleep']['process']['peak_rss_bytes']>=32*1024*1024
    tests['unsupported_did_not_solve']='response' not in records['unsupported']
    original=case();proof=deepcopy(records['echo']['response']['certificate'])
    for name,key,value in [('function','function',{**proof['function'],'offset':'1'}),
                           ('space','mean_zero',True),('tolerance','tolerance','1/2'),
                           ('precision','upper','1'),('negative_width','upper','-1')]:
        changed=deepcopy(proof);changed[key]=value
        tests['reject_'+name]=not b.assess(changed,original['task'],{'certificate_valid':True})['target_met']
    query_cases=[case('Q'+str(i).zfill(2)) for i in range(16)]
    batch=b.run_process({'mode':'batch','system':'baseline','synthetic':'echo','cases':query_cases,
                         'family_metadata':{'id':'synthetic'},'family':{},'repetition':0,'limit':3},destination/'batch_echo')
    tests['batch_fixed_16_success']=len(batch['rows'])==16 and batch['session_complete'] and all(r['success'] for r in batch['rows'])
    failed_batch=b.run_process({'mode':'batch','system':'baseline','synthetic':'swallow','cases':query_cases,
                                'family_metadata':{'id':'synthetic'},'family':{},'repetition':0,'limit':.25},destination/'batch_timeout')
    tests['batch_timeout_keeps_16']=len(failed_batch['rows'])==16 and not any(r['success'] for r in failed_batch['rows'])
    tests['batch_global_deadline']=failed_batch['process']['watchdog']=='global_deadline'
    unsupported_batch=b.run_process({'mode':'batch','system':'baseline','synthetic':'unsupported','cases':query_cases,
                                    'family_metadata':{'id':'synthetic'},'family':{},'repetition':0,'limit':2},destination/'batch_unsupported')
    tests['batch_prepare_failure_keeps_16']=len(unsupported_batch['rows'])==16 and not any(r['success'] for r in unsupported_batch['rows'])
    replay=b.run_process({'mode':'replay','system':'baseline','synthetic':'echo','case':original,
                          'certificate':proof,'limit':2},destination/'replay')
    tests['replay_full_without_prepare']=replay['success'] and replay['preparation'] is None and replay['verification_full'] is True
    # Different capability subsets must not erase a valid baseline/C comparison.
    fake_cases=[case('X'),case('Y')]
    fake_rows=[]
    for system in b.SYSTEMS:
        for item in fake_cases:
            for repeat in range(b.REPEATS):
                row=deepcopy(records['echo']);row.update(system=system,case=b.metadata(item),repetition=repeat)
                row['success']=not(system=='B' and item['id']=='Y')
                fake_rows.append(row)
    summary=b.common_summary(fake_rows,fake_cases)
    tests['all_system_stable_intersection_fixed']=summary['C']['all_four_systems_5_of_5_common_subset']==['X']
    tests['pairwise_preserves_baseline_C_subset']=summary['C']['baseline_pairwise_5_of_5_subset']==['X','Y']
    tests['failure_retained_in_fixed_denominator']=summary['B']['fixed_run_denominator']==10 and summary['B']['successful_runs']==5
    tests['exact_precision_metric_reported']=summary['baseline']['per_case'][0]['actual_metric_exact_by_repetition'][0]['metric']=='1/1000000000'
    records.update(batch=batch,failed_batch=failed_batch,unsupported_batch=unsupported_batch,replay=replay)
    if include_baseline:
        smoke=b.run_process({'mode':'cold','system':'baseline','case':case('baseline_step_smoke')},destination/'baseline_smoke')
        tests['unchanged_baseline_light_smoke']=smoke['success']
        records['baseline_smoke']=smoke
    result={'passed':all(tests.values()),'test_count':len(tests),'tests':tests,
            'synthetic_only':not include_baseline,'only_math_test':'one unchanged baseline step' if include_baseline else None,
            'no_final_comparison':True,'records':records}
    b.write(destination/'selftest_report.json',result)
    print(json.dumps({'passed':result['passed'],'test_count':len(tests),'failed':[n for n,v in tests.items() if not v],
                      'report':str(destination/'selftest_report.json')},ensure_ascii=False,indent=2))
    if not result['passed']:
        raise SystemExit(1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline-smoke',action='store_true')
    run(parser.parse_args().baseline_smoke)
