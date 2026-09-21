"""Sequential development probes: cost/gap improvement and goal-aware ablation.
No model calls; no hidden or held-out problem claim. All proof replay is timed.
"""
from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
import argparse,json,platform,statistics
import matrix_v30
import dual_v33 as dual
import optimize_v34 as optimizer
import planner_v35 as new
import planner_v32 as old
import run

HERE=Path(__file__).resolve().parent
DIRS=[[0,1],[F(-1,3),0,1]]


def goal(penalty='a*a/4+b*b/9'):
    return {'geometry':'unit_sphere','function_space':'full_sphere','eigenvalue_index':1,
            'parameters':['a','b'],'potential':'a*t+b*(t*t-1/3)','penalty':penalty,
            'domain':{'kind':'all_real'},'threshold':'0'}


def cases():
    out=[]
    for name,dirs,cost in [
        ('linear',[[0,1]],None),('odd_even',DIRS,None),
        ('nonorthogonal',[[0,1,1],[0,2,-1]],None),
        ('quartic',[[0,1],[0,0,0,0,1]],None),
        ('degree_six',[[0,1],[0,0,0,0,0,0,1]],None),
        ('three_directions',[[0,1],[0,0,1],[0,0,0,1]],None),
        ('four_directions',[[0,1],[0,0,1],[0,0,0,1],[0,0,0,0,1]],None),
        ('weighted_cost',DIRS,[[2,F(1,2)],[F(1,2),1]])]:
        _,_,c=dual.setup(dirs,cost)
        out.append(('v30_identity_'+name,lambda dirs=dirs:{'certificate':matrix_v30.synthesize(dirs),'statistics':{'model_calls':0}},c))
        out.append(('v34_optimized_'+name,lambda dirs=dirs,cost=cost:optimizer.solve(dirs,cost),c))
    out.append(('v34_one_round',lambda:optimizer.solve(DIRS,max_rounds=1),None))
    out.append(('v34_range_budget',lambda:optimizer.solve([[0,1],[0,0,0,0,0,0,1]],max_rounds=1,max_range_leaves=1),None))
    raw=goal();mapped=goal('(a+2*b)**2/4+(b-a)**2/9');mapped['potential']='(a+2*b)*t+(b-a)*(t*t-1/3)'
    for name,task in [('anisotropic_goal',raw),('mapped_goal',mapped)]:
        out.append(('v32_'+name,lambda task=task:old.solve(task,max_leaves=32),None))
        out.append(('v35_'+name,lambda task=task:new.solve(task,max_leaves=32),None))
    weak=goal();weak.update({'parameters':['a'],'potential':'a*t','penalty':'a*a/5','domain':{'kind':'box','box':[[-2,2]]}})
    counter=goal();counter.update({'parameters':['a'],'potential':'a*t','penalty':'a*a/10','threshold':'-1/10','domain':{'kind':'box','box':[[-4,4]]}})
    generic=goal();generic.update({'parameters':['a'],'potential':'2-3*t+3*t*t/2+a*t','penalty':'3*a*a/20+3*a/100','threshold':'-1','domain':{'kind':'box','box':[[-4,4]]}})
    out.append(('v35_weak_poisson_ansatz',lambda:new.solve(weak,max_leaves=16),None))
    out.append(('v35_counterexample',lambda:new.solve(counter,max_leaves=16),None))
    out.append(('v32_generic',lambda:old.solve(generic,max_leaves=32),None))
    out.append(('v35_generic',lambda:new.solve(generic,max_leaves=32),None))
    return out


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--repeats',type=int,default=3);args=parser.parse_args()
    if not 1<=args.repeats<=10:raise ValueError('Repeat budget')
    output=HERE/'results';proofdir=output/'proofs';proofdir.mkdir(parents=True,exist_ok=True);rows=[]
    for name,fn,cost in cases():
        timings=[];certs=[]
        for _ in range(args.repeats):
            start=perf_counter();result=fn();cert=result['certificate']
            assert run.verify(cert)
            assert all(new.verify_budget(d) for d in result.get('diagnostics',[]))
            timings.append(perf_counter()-start);certs.append(cert)
        if any(c!=certs[0] for c in certs):raise AssertionError('Exact certificate changed between repeated runs')
        (proofdir/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
        upper=cert.get('upper_bound');lower=cert.get('lower_bound');gap=cert.get('gap')
        if cert['method']=='poisson_matrix_v30':
            upper=str(dual.inner(cost,[[F(z) for z in row] for row in cert['quadratic_bound']]))
        row={'case':name,'status':cert.get('status','valid_envelope'),'objective_cost':None if cost is None else dual.encode(cost),
             'lower':lower,'upper':upper,'gap':gap,'median_seconds':statistics.median(timings),'seconds_runs':timings,
             'statistics':result.get('statistics',{}),'diagnostic_statuses':[d['status'] for d in result.get('diagnostics',[])],
             'proof':'proofs/'+name+'.json'}
        rows.append(row)
        (output/'benchmark.json').write_text(json.dumps({'python':platform.python_version(),'platform':platform.platform(),'repeats':args.repeats,
            'model_calls':0,'protocol':'Sequential development probes. Time includes construction/search AND explicit replay of the final certificate and diagnostics. Result file I/O is excluded. Not model or held-out evaluation.',
            'cases':rows},indent=2)+'\n')
        print(name,row['status'],float(F(upper)) if upper is not None else None,round(row['median_seconds'],6),flush=True)

if __name__=='__main__':main()
