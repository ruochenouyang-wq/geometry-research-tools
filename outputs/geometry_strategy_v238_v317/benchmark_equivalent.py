"""Same-certificate, warmed serial probes; no external model or speed guarantee."""
from copy import deepcopy
import json
import platform
from statistics import median
import sys
from time import perf_counter,process_time

import backend
import fast_math


def paired(name,before,after,verify,repetitions=21):
    old=before();new=after()
    if old!=new or verify(new) is not True:
        raise AssertionError('Same-result prerequisite failed: '+name)
    samples={'before':[],'after':[]}
    # Explicit warm-up, then alternate ordering to reduce directional drift.
    for index in range(repetitions):
        for label,call in ([('before',before),('after',after)] if index%2==0 else
                           [('after',after),('before',before)]):
            wall,cpu=perf_counter(),process_time()
            value=call()
            samples[label].append({'wall_seconds':perf_counter()-wall,'cpu_seconds':process_time()-cpu})
            if value!=old:raise AssertionError('Output changed in repeat')
    medians={key:{metric:median(row[metric] for row in rows)
                  for metric in ('wall_seconds','cpu_seconds')} for key,rows in samples.items()}
    return {'name':name,'exact_same_certificate':True,'frozen_replay_passed':True,
            'repetitions_each':repetitions,'samples':samples,'medians':medians,
            'wall_ratio_before_over_after':medians['before']['wall_seconds']/medians['after']['wall_seconds']}


def main():
    q=deepcopy(backend.enriched.FUNCTION)
    source=json.loads((backend.FROZEN/'certificates'/'singular_n4.json').read_text())
    kwargs={'mean_zero':True,'modes':6,'bits':28,'max_m':2,'near_tail':6}
    results=[
        paired('full_spectrum_n6_near6',lambda:backend.direct.full_ground(q,**kwargs),
               lambda:fast_math.full_ground(q,**kwargs),backend.direct.verify_full),
        paired('piecewise_16_levels_degree6',lambda:backend.pieces.piecewise_model(q,levels=16,degree=6),
               lambda:fast_math.piecewise_model(q,levels=16,degree=6),backend.pieces.verify_piecewise),
        paired('enriched_six_terms',lambda:backend.enriched.enrich(source,max_terms=6)['certificate'],
               lambda:fast_math.enrich(source,max_terms=6)['certificate'],backend.enriched.verify)]
    result={'python':sys.version,'platform':platform.platform(),'scope':'fixed inputs; warmed serial process; imports excluded',
            'limitations':['No isolated-machine guarantee; no statistical generalization.',
                           'Shared mathematical caches are warm; no cold-start claim.',
                           'Same-result component probes do not imply automatic research is faster on every task.'],
            'results':results,'actual_model_usage':None}
    target=backend.ROOT/'metrics'/'equivalent_benchmark.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'report':str(target),'ratios':{r['name']:r['wall_ratio_before_over_after'] for r in results}},indent=2))


if __name__=='__main__':main()
