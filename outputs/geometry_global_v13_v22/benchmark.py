"""Reproducible sequential benchmark, including failures and slower strategies."""
import argparse
import json
import platform
import statistics
from pathlib import Path
from fractions import Fraction as F
from time import perf_counter
import global_core as core
import global_search as baseline
import point_oracle
import bernstein_v17 as bernstein
import equality_v19 as equality
import robust_v20 as robust
import constrained_v21 as constrained

HERE=Path(__file__).resolve().parent


def cases():
    result=[]
    p=core.problem([0],[[0,1]],[[-6,6]],[F(1,50)],[[F(1,5)]])
    result.append(('v11_fixed12_tilt',lambda:baseline.search(p,modes=12)))
    for v in (13,14,15,22):result.append((f'v{v}_tilt',lambda v=v:core.search(p,version=v)))
    result.append(('v13_fixed4_floor',lambda:core.search(p,F(1,10**8),version=13,max_modes=4,max_leaves=16)))
    m=core.all_real_problem([0,8],[[0,1],[0,0,1]],[0,0],[[F(1,5),0],[0,1]])
    for cuts in (False,True):
        result.append(('v16_'+('cuts' if cuts else 'independent_moments'),lambda cuts=cuts:core.search(m,F(1,10000),version=16,domain='all_real',moment_cuts=cuts)))
    for d in (3,4):
        high=core.problem([0],[[0]*i+[1] for i in range(1,d+1)],[[-1,1]]*d,[0]*d,[[2*int(i==j) for j in range(d)] for i in range(d)])
        for v in (17,18,22):result.append((f'v{v}_{d}d',lambda high=high,v=v:core.search(high,F(1,10),version=v,max_leaves=64)))
    eq=core.problem([0],[[0,1],[0,0,1],[0,0,0,1]],[[-1,1]]*3,[0]*3,[[int(i==j) for j in range(3)] for i in range(3)])
    result.append(('v19_three_to_two',lambda:equality.search(eq,[[1,0,1]],[0])))
    rp=core.problem([0],[[0,1]],[[-1,1]],[0],[[1]])
    result.append(('v20_symmetric',lambda:robust.search(rp,[[0,1]],[[-1,1]])))
    result.append(('v20_asymmetric',lambda:robust.search(rp,[[0,1]],[[-1,2]])))
    result.append(('v20_budget_open',lambda:robust.search(rp,[[0,1]],[[-1,2]],F(1,10**6),max_steps=1)))
    rp2=core.problem([0],[[0,1],[0,0,1]],[[-1,1]]*2,[0]*2,[[2,0],[0,2]])
    result.append(('v20_two_design_two_uncertainty',lambda:robust.search(rp2,[[0,1],[0,0,1]],[[F(-1,2),F(1,2)]]*2,F(1,100))))
    for target in (F(1,2),F(9,10),F(2)):
        result.append(('v21_target_'+str(target).replace('/','_'),lambda target=target:constrained.search([0],[0,1],target)))
    result.append(('v21_feasibility_not_found',lambda:constrained.search([0],[0,1],F(1,2),max_steps=1)))
    extra={
        'weak':core.problem([0],[[0,F(1,10)]],[[-1,1]],[0],[[1]]),
        'asymmetric_quadratic':core.problem([2,-3,F(3,2)],[[0,1]],[[-4,4]],[F(3,100)],[[F(3,10)]]),
        'even_quartic_2d':core.problem([0,0,-8,0,2],[[0,1],[0,0,1]],[[-3,3],[-1,1]],[0,0],[[F(1,2),0],[0,1]]),
        'spectral_action':core.problem([0],[[0,1]],[[-4,4]],[0],[[F(1,5)]])}
    for name,problem in extra.items():
        for v in (15,22):result.append((f'v{v}_extra_{name}',lambda problem=problem,v=v:core.search(problem,F(1,1000),version=v,max_leaves=128)))
    return result


def microbenchmarks():
    p=core.problem([0],[[0,1]],[[-20,20]])
    points={}
    for name,maximum in [('fixed4',4),('adaptive',24)]:
        start=perf_counter();o=point_oracle.Oracle(p,F(1,10**5),max_modes=maximum);r=o.records[o.point([20])]
        points[name]={'seconds':perf_counter()-start,'width':str(point_oracle.width(r)),**o.statistics(),'record':r}
    reuse={}
    p=core.problem([0,2],[[1]],[[-2,2]])
    for enabled in (False,True):
        start=perf_counter();o=point_oracle.Oracle(p,F(1,100000),orbits=enabled)
        for i in range(-16,17):o.point([F(i,8)])
        reuse[str(enabled)]={'seconds':perf_counter()-start,**o.statistics()}
    poly={(0,):F(1,4),(1,):F(-1),(2,):F(1)}
    low,tree=bernstein.refine(poly,1,2)
    return {'v13_spectral_floor':points,'v14_shift_orbits':reuse,
            'v18_same_polynomial':{'polynomial':'(x-1/2)^2 on [0,1]','before':str(min(bernstein.coefficients(poly,1).values())),
                                  'after':str(low),'new_spectral_evaluations':0,'tree':tree}}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repeats',type=int,default=3);parser.add_argument('--output',type=Path,default=HERE/'results');args=parser.parse_args()
    if not 1<=args.repeats<=10:raise ValueError('Repeat budget 1..10')
    args.output.mkdir(parents=True,exist_ok=True);certdir=args.output/'certificates';certdir.mkdir(exist_ok=True)
    rows=[]
    for name,fn in cases():
        times=[];runs=[]
        for repetition in range(args.repeats):
            result=fn();times.append(result['seconds']);c=result['certificate']
            runs.append({key:c.get(key) for key in ('global_lower','candidate_upper','global_gap','status')})
        (certdir/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
        if any(row!=runs[0] for row in runs):raise AssertionError('Proof result changed between repetitions')
        row={'case':name,'version':result['algorithm_version'],'seconds_runs':times,'median_seconds':statistics.median(times),
             'statistics':result.get('statistics',{'points':result.get('point_evaluations')}),'actions':result.get('actions'),
             'leaf_count':c.get('leaf_count'),'moment_excluded_leaves':c.get('moment_excluded_leaves'),
             'tolerance':c.get('tolerance'),**runs[-1],'certificate':'certificates/'+name+'.json'}
        rows.append(row);print(name,c['status'],round(row['median_seconds'],4),'seconds',flush=True)
        (args.output/'benchmark.json').write_text(json.dumps({'python':platform.python_version(),'platform':platform.platform(),
            'repeats':args.repeats,'timing':'sequential; includes proof construction and independent dense verification; excludes JSON serialization',
            'cases':rows},indent=2)+'\n')
    micros=[microbenchmarks() for _ in range(args.repeats)]
    (args.output/'microbenchmarks.json').write_text(json.dumps(micros,indent=2)+'\n')

if __name__=='__main__':main()
