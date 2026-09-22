"""Read-only aggregation of the first fixed results; never reruns candidates."""
from collections import Counter
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from statistics import median

ROOT=Path(__file__).resolve().parent
d=json.loads((ROOT/'RESULTS.json').read_text());rows=d['records']
ARMS=('c2_rational112','numpy_binary64')

def costs(subset):
    out={}
    for phase in ('representation_prepare','candidate','complete_certificate','task_bound_verify',
                  'residual_refinement_subset_of_candidate'):
        values=[r['phases'][phase] for r in subset if phase in r['phases']]
        if values:out[phase]={'count':len(values),**{
            f'{kind}_{metric}':op(x[kind] for x in values)
            for kind in ('wall_seconds','cpu_seconds') for metric,op in (('median',median),('sum',sum))}}
    pipelines=[r['process_pipeline'] for r in subset if r.get('process_pipeline')]
    if pipelines:
        out['process_pipeline']={'count':len(pipelines),**{
            f'{kind}_{metric}':op(x[kind] for x in pipelines)
            for kind in ('wall_seconds','parent_cpu_seconds','child_user_system_cpu_seconds')
            for metric,op in (('median',median),('sum',sum))}}
    return out

summary={'source_results_sha256':hashlib.sha256((ROOT/'RESULTS.json').read_bytes()).hexdigest(),
    'scope':'Fixed warm candidate substitution; differing 53/112bit precision and Ritz/full-residual objectives; no end-to-end/SOTA/token claim.',
    'fixed_denominator':{'configurations_per_arm':12,'repeats':3,'rows_per_arm':36,'total_rows':72},
    'arms':{},'configurations':[],'paired_intersections':{},
    'common_context_wall_seconds':sum(x['phases']['shared_context']['wall_seconds'] for x in d['common_contexts']),
    'common_context_cpu_seconds':sum(x['phases']['shared_context']['cpu_seconds'] for x in d['common_contexts']),
    'common_matrices_wall_seconds':sum(x['phases']['shared_exact_matrices']['wall_seconds'] for x in d['common_matrices']),
    'common_matrices_cpu_seconds':sum(x['phases']['shared_exact_matrices']['cpu_seconds'] for x in d['common_matrices']),
    'common_serialized_matrix_bytes':sum(x['serialized_bytes'] for x in d['common_matrices']),
    'whole_experiment_wall_seconds':d['elapsed_whole_experiment_seconds'],
    'original_c2_equivalence_control':d['single_original_c2_equivalence_check']}
for arm in ARMS:
    rr=[r for r in rows if r['arm']==arm]
    candidates=[r['candidate'] for r in rr if r.get('candidate')]
    summary['arms'][arm]={'row_outcomes':dict(Counter(r['outcome'] for r in rr)),
        'costs':costs(rr),'refinement_attempts':sum(x['refinement_attempts'] for x in candidates),
        'refinement_accepted':sum(x['refinement_accepted'] for x in candidates),
        'refinement_inner_wall_seconds':sum(z['wall_seconds'] for x in candidates for z in x.get('residual_refinement',[])),
        'refinement_inner_cpu_seconds':sum(z['cpu_seconds'] for x in candidates for z in x.get('residual_refinement',[]))}
for ident in ('E09','E14','E16','E19'):
    for count in (6,10,16):
        entry={'case_id':ident,'powers_count':count,'arms':{}}
        for arm in ARMS:
            rr=[r for r in rows if (r['case_id'],r['powers_count'],r['arm'])==(ident,count,arm)]
            widths=[Fraction(r['certificate']['exact_width']) for r in rr if r.get('certificate')]
            entry['arms'][arm]={'outcomes':dict(Counter(r['outcome'] for r in rr)),
                'costs':costs(rr),'median_certified_width':str(median(widths)) if widths else None,
                'median_certified_width_float':float(median(widths)) if widths else None,
                'widths_identical_across_repeats':len(set(widths))<=1 if widths else None,
                'refinement_attempts':[r['candidate']['refinement_attempts'] if r.get('candidate') else None for r in rr],
                'errors':list({r.get('reason') for r in rr if r.get('reason')})}
        summary['configurations'].append(entry)
for arm in ARMS:
    configs=[x['arms'][arm] for x in summary['configurations']]
    summary['arms'][arm]['configuration_outcomes']=dict(Counter(next(iter(x['outcomes']))
        if len(x['outcomes'])==1 else 'mixed_repeats' for x in configs))
for label,accepted in (('both_target_met',{'target_met'}),('both_certificate_valid',{'target_met','valid_open'})):
    selected=[(x['case_id'],x['powers_count']) for x in summary['configurations'] if all(
        set(x['arms'][arm]['outcomes'])<=accepted for arm in ARMS)]
    summary['paired_intersections'][label]={'configurations':selected,'configurations_count':len(selected),
        'arms':{arm:costs([r for r in rows if r['arm']==arm and (r['case_id'],r['powers_count']) in selected]) for arm in ARMS}}
summary['consistency']={'all_rows_unique':len({(r['case_id'],r['powers_count'],r['repeat'],r['arm']) for r in rows})==72,
    'all_expected_rows':len(rows)==72,'all_returned_certificates_taskbound_valid':all(
        r['verification']['certificate_valid'] for r in rows if r.get('certificate')),
    'certificate_count':sum(bool(r.get('certificate')) for r in rows),
    'all_repeated_widths_identical':all(x['arms'][arm]['widths_identical_across_repeats'] is not False
        for x in summary['configurations'] for arm in ARMS),
    'all_source_hashes_unchanged':all(hashlib.sha256(Path(path).read_bytes()).hexdigest()==value for path,value in d['source_hashes'].items()),
    'all_processes_within_deadline':all(r['process_pipeline']['within_deadline'] for r in rows if r.get('process_pipeline'))}
(ROOT/'SUMMARY.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ('configurations','arms','paired_intersections')},indent=2))
for x in summary['configurations']:
    print(x['case_id'],x['powers_count'],[(arm,v['outcomes'],v['costs']['candidate']['wall_seconds_median'],
        v['costs']['process_pipeline']['wall_seconds_median'],v['median_certified_width_float']) for arm,v in x['arms'].items()])
for name,x in summary['paired_intersections'].items():print(name,json.dumps(x))
