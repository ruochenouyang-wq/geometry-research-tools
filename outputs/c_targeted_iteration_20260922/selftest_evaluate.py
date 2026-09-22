"""Synthetic evaluator checks only: no C2/C3 mathematical search or tuning."""
from copy import deepcopy
import json
from pathlib import Path
import time
import zipfile
import evaluate as e


def case(name='T01',classification='new_public'):
    return {'id':name,'dimension':'synthetic','label':'synthetic only','classification':classification,
            'task':{'kind':'spectrum','function':{'kind':'axis_profile','profile':'step','amplitude':'-1','offset':'0','axis':2},
                    'mean_zero':False,'tolerance':'1/1000000','budget':{'wall_seconds':10}}}


def main():
    folder=e.ROOT/'evaluation'/('selftest_'+e.stamp())
    folder.mkdir(parents=True)
    tests={};records={}
    for behavior in ('echo','raise','unsupported','status_only','swallow','memory_timeout','verify_timeout','child_busy','child_group_timeout'):
        cutoff=.3 if behavior in ('swallow','memory_timeout','verify_timeout','child_group_timeout') else 3
        marker=folder/'descendant_marker.txt'
        record=e.run_process({'system':'C3','case':case(),'repetition':0,'synthetic':behavior,
                              'synthetic_marker':str(marker),'limit':cutoff},folder/behavior)
        records[behavior]=record
        tests[behavior+'_correct_outcome']=record['success'] if behavior in ('echo','child_busy') else not record['success']
        tests[behavior+'_wait4_rss_positive']=record['process']['peak_rss_bytes']>0
        tests[behavior+'_wait4_cpu_nonnegative']=record['process']['cpu_seconds']>=0
        if cutoff<1:
            tests[behavior+'_parent_kill']=record['process']['watchdog']=='entire_process_deadline' and record['process']['returncode']==-9
    tests['hardkill_progress_survives']=bool(records['swallow']['progress']) and records['swallow']['last_internal_action']=='synthetic_search'
    tests['hardkill_outer_stage_retained']=records['swallow']['failure_stage']=='solve'
    tests['verify_timeout_retains_full_certificate']=records['verify_timeout']['certificate_json_bytes']>0
    tests['verify_timeout_stage_retained']=records['verify_timeout']['failure_stage']=='full_verification'
    tests['timeout_memory_not_dropped']=records['memory_timeout']['process']['peak_rss_bytes']>=32*1024*1024
    tests['unsupported_method_domain_classified']=records['unsupported']['failure_class']=='unsupported'
    tests['timeout_cpu_coverage_flagged']=records['child_group_timeout']['process']['forced_termination_cpu_may_be_lower_bound']
    time.sleep(.4)
    tests['process_group_killed_descendant']=not (folder/'descendant_marker.txt').exists()
    child=records['child_busy']
    resource_record=child['worker_usage_before_final_write']
    outer_cpu=child['process']['cpu_seconds']
    self_cpu=resource_record['self_cpu_seconds'];children_cpu=resource_record['reaped_children_cpu_seconds']
    calibration={'platform':child['process']['platform'],'outer_wait4_cpu':outer_cpu,
        'worker_self_cpu_before_final_write':self_cpu,'worker_reaped_children_cpu':children_cpu,
        'worker_self_plus_reaped_children_cpu':self_cpu+children_cpu,
        'outer_minus_instrumented_sum':outer_cpu-self_cpu-children_cpu,
        'child_reported_cpu':child['response']['synthetic_child_metrics']['cpu'],
        'outer_wait4_peak_rss_bytes':child['process']['outer_wait4_peak_rss_bytes'],
        'worker_self_peak_rss_bytes':resource_record['self_peak_rss_high_water_bytes'],
        'observed_reaped_child_peak_rss_bytes':child['process']['observed_reaped_children_peak_rss_bytes'],
        'reported_largest_single_process_peak_rss_bytes':child['process']['peak_rss_bytes'],
        'interpretation':'Returned/reaped descendant CPU is already included by outer wait4; never add it again.'}
    tests['child_cpu_present_in_worker_children']=children_cpu>.1
    tests['outer_wait4_includes_reaped_child_cpu']=outer_cpu>=self_cpu+children_cpu-.01
    tests['outer_wait4_does_not_require_adding_children_again']=abs(outer_cpu-(self_cpu+children_cpu))<.08
    tests['solve_stage_includes_child_cpu']=child['timings']['solve']['reaped_children_cpu_seconds']>.1
    tests['stage_cpu_self_plus_children']=abs(child['timings']['solve']['cpu_seconds']-
        child['timings']['solve']['self_cpu_seconds']-child['timings']['solve']['reaped_children_cpu_seconds'])<1e-9
    tests['child_rss_exceeds_worker_for_calibration']=child['process']['observed_reaped_children_peak_rss_bytes']>resource_record['self_peak_rss_high_water_bytes']
    tests['reported_rss_covers_observed_child']=child['process']['peak_rss_bytes']>=child['process']['observed_reaped_children_peak_rss_bytes']>=64*1024*1024
    item=case();proof=records['echo']['response']['certificate']
    for name,key,value in [('function','function',{**proof['function'],'offset':'1'}),('space','mean_zero',True),
                           ('kind','kind','approximation'),('tolerance','tolerance','1/2'),
                           ('precision','upper','1'),('negative_width','upper','-1')]:
        changed=deepcopy(proof);changed[key]=value
        tests['reject_'+name]=not e.assess(changed,item['task'],{'certificate_valid':True})['target_met']
    tests['status_is_not_proof']=not e.assess(proof,item['task'],{'status':'target_met','target_met':True})['target_met']
    tests['reused_C2_is_not_unexplained_fallback']=e.provenance('C3',{'route':'reused_c2','native_kind':'existing_demand_threshold','fallback_reason':None})['explanation_complete']
    tests['true_baseline_fallback_identified']=e.provenance('C3',{'route':'baseline_fallback','fallback_reason':'synthetic native failure','counters':{'baseline_fallback_calls':1}})['route']=='baseline_fallback'
    source=folder/'snapshot_source.txt';source.write_text('original source')
    manifest={str(source):e.sha(source)}
    e.assert_files(manifest)
    source.write_text('modified source')
    try:e.assert_files(manifest);tests['changed_source_rejected']=False
    except RuntimeError:tests['changed_source_rejected']=True
    archive=folder/'synthetic_snapshot.zip'
    with zipfile.ZipFile(archive,'w') as z:z.writestr('source','original source')
    archive_hash=e.sha(archive);e.assert_snapshot(archive,archive_hash)
    with archive.open('ab') as stream:stream.write(b'changed')
    try:e.assert_snapshot(archive,archive_hash);tests['changed_snapshot_archive_rejected']=False
    except RuntimeError:tests['changed_snapshot_archive_rejected']=True
    cases=[case('X'),case('Y','old_failure_regression'),case('Z','expected_unsupported_control')]
    fake=[]
    for system in e.SYSTEMS:
        for item in cases:
            for repetition in range(e.REPEATS):
                row=deepcopy(records['unsupported'] if item['id']=='Z' else records['echo'])
                row.update(system=system,case=e.meta(item),repetition=repetition)
                if system=='C2' and item['id']=='Y':row['success']=False
                fake.append(row)
    summary=e.summarize(fake,cases)
    tests['fixed_denominator_retained']=summary['C2']['run_denominator']==9 and summary['C2']['returned_rows']==9
    tests['control_rejection_not_math_success']=summary['C3']['expected_unsupported_control_outcomes']['reported_unsupported_within_deadline']==3 and summary['C3']['expected_unsupported_control_outcomes']['math_successes']==0
    tests['new_and_regression_separate']=summary['C3']['by_classification']['new_public']['case_denominator']==1 and summary['C3']['by_classification']['old_failure_regression']['case_denominator']==1
    tests['pairwise_same_success_only']=summary['C3']['pairwise_comparisons']['C2']['stable_same_success_subset']==['X']
    tests['math_denominator_excludes_controls']=summary['C3']['mathematical_tasks_excluding_method_boundary_controls']['fixed_case_denominator']==2
    # Generated evaluation artifacts must never enter the implementation source set.
    fixture=folder/'source_scope_fixture';(fixture/'evaluation').mkdir(parents=True);(fixture/'notes').mkdir()
    (fixture/'implementation.py').write_text('x = 1\n')
    (fixture/'evaluation'/'generated.py').write_text('x = 2\n')
    (fixture/'notes'/'analysis.py').write_text('x = 3\n')
    tests['outputs_excluded_from_source_freeze']=[x.name for x in e.source_paths(fixture)]==['implementation.py']
    tests['alternating_order_each_repetition']=all(e.system_order(0,i)[0]!=e.system_order(1,i)[0] and e.system_order(0,i)[0]==e.system_order(2,i)[0] for i in range(30))
    # Commitment checks may inspect sealed tasks, but reveal no parameters in results.
    committed_cases,_=e.inputs()
    tests['fixed_30_and_original_24_preserved']=len(committed_cases)==30 and committed_cases[:24]==e.read(e.PREVIOUS/'cases.json')
    tests['six_new_math_identities']=len({e.mathematical_identity(c) for c in committed_cases[24:]})==6
    changes=e.regressions(fake,cases)
    tests['same_success_regression_summary']=len(changes['stable_success_gained'])==1 and not changes['stable_success_lost']
    tests['identical_single_thread_controls']=set(e.THREAD_ENV.values())=={'1'} and len(e.THREAD_ENV)==6
    report={'passed':all(tests.values()),'test_count':len(tests),'tests':tests,'synthetic_only':True,
            'no_mathematical_search':True,'cpu_calibration':calibration,'records':records}
    e.write(folder/'selftest_report.json',report)
    print(json.dumps({'passed':report['passed'],'test_count':len(tests),'failed':[key for key,value in tests.items() if not value],
                      'cpu_calibration':calibration,'report':str(folder/'selftest_report.json')},ensure_ascii=False,indent=2))
    if not report['passed']:raise SystemExit(1)


if __name__=='__main__':main()
