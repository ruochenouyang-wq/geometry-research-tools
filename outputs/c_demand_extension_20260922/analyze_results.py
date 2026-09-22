"""Derive compact, fixed-denominator comparisons from the frozen run ledger."""
import csv
from fractions import Fraction
import json
from pathlib import Path
import statistics
import sys

ROOT=Path(__file__).resolve().parent


def main(path):
    report=json.loads(Path(path).read_text())
    assert report['valid_frozen_comparison'] and report['complete_fixed_denominator']
    rows=report['rows'];summary=report['summary']
    cases=json.loads((ROOT/'cases.json').read_text())
    bysystem={s:{r['id']:r for r in summary[s]['per_case']} for s in summary}
    def runs(s,ids):return [r for r in rows if r['system']==s and r['case']['id'] in ids]
    def stable(s,c):return bysystem[s][c]['stable_3_of_3']
    comparison={}
    for reference in ('baseline','C1'):
        ids=[c['id'] for c in cases if stable(reference,c['id']) and stable('C2',c['id'])]
        data={'common_stable_cases':ids,
              'gained_stable_cases':[c['id'] for c in cases if not stable(reference,c['id']) and stable('C2',c['id'])],
              'lost_stable_cases':[c['id'] for c in cases if stable(reference,c['id']) and not stable('C2',c['id'])]}
        for metric in ('wall_seconds','cpu_seconds'):
            ours=sum(r['process'][metric] for r in runs('C2',ids))
            prior=sum(r['process'][metric] for r in runs(reference,ids))
            data[metric]={'C2':ours,reference:prior,'reduction_fraction':1-ours/prior if prior else None}
        comparison[reference]=data
    compact={'source_report':str(Path(path).resolve()),'valid_frozen_comparison':True,
             'classification':'public_development_not_holdout','comparisons':comparison,
             'systems':summary,'row_count':len(rows)}
    destination=ROOT/'results';destination.mkdir(exist_ok=True)
    (destination/'summary.json').write_text(json.dumps(compact,ensure_ascii=False,indent=2)+'\n')
    fields=['id','classification','label','kind','tolerance']
    for system in ('baseline','C1','C2'):
        fields += [system+'_'+x for x in ('successes','wall_median','cpu_median','peak_rss_median','certificate_bytes_median','metric','route')]
    table=[]
    for c in cases:
        row={k:c[k] for k in ('id','classification','label')}
        row.update(kind=c['task']['kind'],tolerance=c['task']['tolerance'])
        for s in ('baseline','C1','C2'):
            rr=runs(s,[c['id']]);first=rr[0]
            row[s+'_successes']=sum(r['success'] for r in rr)
            for label,metric in [('wall','wall_seconds'),('cpu','cpu_seconds'),('peak_rss','peak_rss_bytes')]:
                row[s+'_'+label+'_median']=statistics.median(r['process'][metric] for r in rr)
            row[s+'_certificate_bytes_median']=statistics.median(r['certificate_json_bytes'] for r in rr)
            metric=first.get('assessment',{}).get('metric')
            row[s+'_metric']=float(Fraction(metric)) if metric is not None else ''
            row[s+'_route']=';'.join(sorted({r['provenance']['route'] for r in rr}))
        table.append(row)
    with (destination/'case_comparison.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fields);writer.writeheader();writer.writerows(table)
    traces=[]
    for row in rows:
        if row['system']!='C2':continue
        response=row.get('response',{})
        interesting=[a for a in response.get('attempts',[]) if a.get('status') in ('failed','budget_stop') or a.get('reason_code') or a.get('action') in ('candidate_assessment','stop')]
        traces.append({'case':row['case']['id'],'repetition':row['repetition'],
            'success':row['success'],'status':row['status'],'failure_stage':row.get('failure_stage'),
            'failure_reason':row.get('failure_reason'),'provenance':row['provenance'],
            'last_progress':row.get('last_progress'),'attempts':interesting})
    (destination/'c2_decisions_and_failures.json').write_text(json.dumps(traces,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'row_count':len(rows),'comparisons':comparison,
        'counts':{s:{'stable':summary[s]['stable_success_cases'],'strata':summary[s]['by_classification']} for s in summary}},ensure_ascii=False,indent=2))


if __name__=='__main__':main(sys.argv[1])
