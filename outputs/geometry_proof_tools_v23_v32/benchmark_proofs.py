"""Sequential component/structured-goal ablations. No model performance claims."""
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
import argparse,json,statistics,platform
import compiler_v23 as compiler
import reductions
import decision_v28 as decision
import family_v29 as family
import matrix_v30
import transport_v31
import planner_v32 as planner
import global_core

HERE=Path(__file__).resolve().parent


def goal(potential='a*t',penalty='a*a/4',parameters=None,box=None,threshold='0',scope='full_sphere'):
    return {'geometry':'unit_sphere','function_space':scope,'eigenvalue_index':1,'parameters':['a'] if parameters is None else parameters,
            'potential':potential,'penalty':penalty,'domain':{'kind':'all_real'} if box is None else {'kind':'box','box':box},'threshold':threshold}


def routed(raw,rules=(),optimize=False,max_leaves=64):
    start=perf_counter();compiled=compiler.compile_goal(raw);p=compiled['problem'];path=[]
    if p['function_space']=='full_sphere':
        edge=reductions.sphere_ground(p);path.append(edge);p=edge['reduced']
    for fn in rules:
        edge=fn(p);path.append(edge);p=edge['reduced']
    if optimize:
        cp=decision.as_core(p);run=global_core.search(cp,F(1,10**6),version=22,domain=p['domain']['kind'],max_leaves=max_leaves)
        proof=decision.result(p,compiled['threshold'],run['certificate'],'partition');stats=run['statistics']
    else:
        run=decision.search(p,compiled['threshold'],max_leaves=max_leaves);proof=run['certificate'];stats=run['statistics']
    cert=planner.certificate(compiled,path,{'kind':'decision','proof':proof})
    return {'certificate':cert,'statistics':stats,'seconds':perf_counter()-start}


def measured(fn):
    start=perf_counter();c=fn();return {'certificate':c,'statistics':{},'seconds':perf_counter()-start}


def cases():
    out=[]
    out.append(('v23_expression_compilation',23,lambda:measured(lambda:compiler.compile_goal(goal('(a+1)*t-a*t+t*t','(a+1)**2-a*a')))))
    out.append(('v24_full_sphere_ground',24,lambda:routed(goal('2*t*t','0',[],[],threshold='0'))))
    passive=goal('a*t','a*a/4+b*b+c*c+d*d',['a','b','c','d'],[[-2,2]]*4,threshold='-1/10')
    out.append(('v25_without_elimination',25,lambda:routed(passive)))
    out.append(('v25_with_elimination',25,lambda:routed(passive,[reductions.absorb_and_eliminate])))
    reflection=goal('a*t','a*a/10',box=[[-6,6]],threshold='-1/5')
    out.append(('v26_full_box',26,lambda:routed(reflection)))
    out.append(('v26_half_box',26,lambda:routed(reflection,[reductions.reflect])))
    rank=goal('(a+2*b-c+3*d)*t','(a*a+b*b+c*c+d*d)/2',['a','b','c','d'],threshold='-4')
    out.append(('v27_four_parameters',27,lambda:routed(rank)))
    out.append(('v27_one_effective_parameter',27,lambda:routed(rank,[reductions.compress])))
    threshold=goal('a*t','a*a/10+a/50',box=[[-6,6]],threshold='-1/10')
    out.append(('v28_optimize_then_decide',28,lambda:routed(threshold,optimize=True,max_leaves=128)))
    out.append(('v28_targeted_decision',28,lambda:routed(threshold)))
    out.append(('v29_scalar_linear_lemma',29,lambda:measured(lambda:family.synthesize([0,1]))))
    out.append(('v29_scalar_quadratic_lemma',29,lambda:measured(lambda:family.synthesize([F(-1,3),0,1]))))
    out.append(('v29_scalar_degree_six_lemma',29,lambda:measured(lambda:family.synthesize([0,0,0,0,0,0,1]))))
    out.append(('v30_matrix_lemma',30,lambda:measured(lambda:matrix_v30.synthesize([[0,1],[F(-1,3),0,1]]))))
    out.append(('v30_nonidentity_weight',30,lambda:measured(lambda:matrix_v30.synthesize([[0,1],[0,0,1]],[[2,F(1,2)],[F(1,2),1]]))))
    out.append(('v30_scale_budget_open',30,lambda:measured(lambda:matrix_v30.synthesize([[0,0,0,0,0,0,1]],tolerance=F(1,10**12),max_leaves=1))))
    library=[family.synthesize([0,1]),matrix_v30.synthesize([[0,1],[F(-1,3),0,1]])]
    transfer=goal('(a+2*b)*t+(b-a)*(t*t-1/3)+t**6','((a+2*b)**2+(b-a)**2)/4',['a','b'])
    out.append(('v31_affine_dominance_transfer',31,lambda:planner.solve(transfer,library=library,synthesize=False)))
    for name,raw in [
        ('uniform_linear',goal()),
        ('uniform_two_directions',goal('a*t+b*(t*t-1/3)','(a*a+b*b)/4',['a','b'])),
        ('counterexample',goal('a*t','a*a/10',box=[[-4,4]],threshold='-1/10')),
        ('too_weak_lemma',goal('a*t','a*a/5',box=[[-2,2]])),
        ('four_parameter_counterexample',goal('(a+b+c+d)*t','(a*a+b*b+c*c+d*d)/2',['a','b','c','d'])),
        ('generic_asymmetric',goal('2-3*t+3*t*t/2+a*t','3*a*a/20+3*a/100',box=[[-4,4]],threshold='-1'))]:
        out.append(('v32_'+name,32,lambda raw=raw:planner.solve(raw,max_leaves=32)))
        if name in ('uniform_linear','uniform_two_directions','generic_asymmetric'):
            out.append(('v32_no_lemmas_'+name,32,lambda raw=raw:planner.solve(raw,synthesize=False,max_leaves=32)))
    return out


def transfer_batch(output,repeats):
    timings=[];synthesis=[];counts=[]
    for repeat in range(repeats):
        start=perf_counter();library=[family.synthesize([0,1])];synthesis.append(perf_counter()-start);start=perf_counter();points=0
        for index in range(16):
            scale=index%4+1;offset=index%5-2;positive=index%3+1
            raw=goal(f'({scale}*a+({offset}))*t+{positive}*t**6',f'({scale}*a+({offset}))**2/4')
            result=planner.solve(raw,library=library,synthesize=False);assert result['certificate']['status']=='proved';points+=result['statistics']['points']
            if repeat==repeats-1:(output/f'v31_transfer_{index:02d}.json').write_text(json.dumps(result,indent=2)+'\n')
        timings.append(perf_counter()-start);counts.append(points)
    return {'tasks':16,'precompute_seconds':synthesis,'online_batch_seconds':timings,'spectral_points_per_batch':counts,
            'protocol':'Constructed transfer probes, not a held-out model benchmark. Online time includes compilation, route selection, applicability checks and proof replay.'}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repeats',type=int,default=3);parser.add_argument('--output',type=Path,default=HERE/'results');args=parser.parse_args()
    if not 1<=args.repeats<=10:raise ValueError('Repeat budget')
    args.output.mkdir(parents=True,exist_ok=True);proofdir=args.output/'proofs';proofdir.mkdir(exist_ok=True);rows=[]
    for name,version,fn in cases():
        times=[];certs=[]
        for repeat in range(args.repeats):
            result=fn();times.append(result['seconds']);certs.append(result['certificate'])
        if any(c!=certs[0] for c in certs):raise AssertionError('Mathematical output changed between repetitions')
        (proofdir/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n');c=result['certificate']
        status=c.get('status',c.get('range_proof',{}).get('status','verified_component'))
        row={'case':name,'version':version,'median_seconds':statistics.median(times),'seconds_runs':times,'status':status,
             'lower_bound':c.get('lower_bound'),'upper_bound':c.get('upper_bound'),'statistics':result['statistics'],'proof':'proofs/'+name+'.json'}
        rows.append(row);print(name,status,result['statistics'].get('points','—'),round(row['median_seconds'],6),flush=True)
        (args.output/'benchmark.json').write_text(json.dumps({'python':platform.python_version(),'platform':platform.platform(),'repeats':args.repeats,
            'model_calls':0,'protocol':'Sequential structured-goal/component ablations. Timing includes compile, reductions and proof reconstruction where applicable. This does not measure model reasoning time.',
            'cases':rows},indent=2)+'\n')
    (args.output/'transfer_batch.json').write_text(json.dumps(transfer_batch(proofdir,args.repeats),indent=2)+'\n')

if __name__=='__main__':main()
