"""Small public development probe, never an independent or cold benchmark."""
from pathlib import Path
import cProfile
import json
import pstats
import statistics
import sys
from time import perf_counter
sys.dont_write_bytecode = True
sys.path.insert(0,str(Path(__file__).resolve().parent))
import variant


def tasks():
    return [dict(kind='spectrum',function=dict(kind='axis_profile',profile='abs_power',
                 amplitude=a,exponent=alpha,offset='0',axis=2),mean_zero=mz,
                 tolerance='1/100000000',budget={'wall_seconds':10})
            for a,alpha in (('-1','-1/4'),('-1/3','-2/5'),('-2','-1/3'))
            for mz in (False,True)]


def run():
    rows=[]
    for task in tasks():
        records={'baseline':[],'C':[]}
        for repeat in range(3):
            for name in (('baseline','C') if repeat%2==0 else ('C','baseline')):
                call=variant.baseline.solve if name=='baseline' else variant.solve
                start=perf_counter();result=call(task)
                encoded=json.dumps(result,allow_nan=False)
                check=variant.verify(json.loads(encoded)['certificate'],task) if result.get('certificate') else {}
                elapsed=perf_counter()-start
                records[name].append({'seconds_with_json_and_external_replay':elapsed,
                    'target_met':result['target_met'],'independent_target_met':check.get('target_met',False),
                    'status':result['status'],'metric':result.get('metric'),
                    'search_counters':result.get('counters')})
        profiles={}
        for name,call in (('baseline',variant.baseline.solve),('C',variant.solve)):
            profiler=cProfile.Profile();profiler.runcall(call,task)
            stats=pstats.Stats(profiler)
            profiles[name]=[{'file':Path(path).name,'function':fun,'calls':nc}
                for (path,line,fun),(cc,nc,tt,ct,callers) in stats.stats.items()
                if fun in ('inertia','banded_inertia')]
        rows.append({'task':task,'runs':records,'profiles_separate_from_timing':profiles,
            'median_seconds':{name:statistics.median(r['seconds_with_json_and_external_replay'] for r in rs)
                              for name,rs in records.items()}})
    report={'scope':'Six explicitly public development inputs; same process, shared baseline caches, alternating order; not a cold or held-out benchmark.',
            'model_calls':0,'model_token_savings':None,'rows':rows}
    path=Path(__file__).resolve().parent/'PUBLIC_PROBE.json'
    path.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({'report':str(path),'rows':[{'function':r['task']['function'],
          'mean_zero':r['task']['mean_zero'],'medians':r['median_seconds'],
          'success':{n:all(x['target_met'] and x['independent_target_met'] for x in xs)
                     for n,xs in r['runs'].items()},'profiles':r['profiles_separate_from_timing']}
          for r in rows]},ensure_ascii=False))


if __name__=='__main__':run()
