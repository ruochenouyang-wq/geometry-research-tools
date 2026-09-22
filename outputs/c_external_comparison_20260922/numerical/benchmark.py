"""One fixed, two-arm candidate substitution experiment. No result tuning."""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[_key]='1'
from contextlib import redirect_stdout
from fractions import Fraction as F
import hashlib
import importlib.util
import io
import json
import math
import multiprocessing as mp
from pathlib import Path
import platform
import resource
import sys
from time import perf_counter, process_time

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
OUTPUTS=ROOT.parents[1]
C2=OUTPUTS/'c_demand_extension_20260922'
IDS=('E09','E14','E16','E19')
COUNTS=(6,10,16)
ARMS=('c2_rational112','numpy_binary64')
PROTOCOL_HASH='f66ba99664b42f5e9780faf58d4a3906db4ce03d19a4f63d5ade060f4af97997'
_extension=None
_full=None


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load():
    global _extension,_full
    if _extension is None:
        spec=importlib.util.spec_from_file_location('_external_comparison_frozen_c2',C2/'fullspace_extension.py')
        _extension=importlib.util.module_from_spec(spec);spec.loader.exec_module(_extension)
        _full=_extension._modules()[1]
    return _extension,_full


def timed(phases,name,operation):
    w,c=perf_counter(),process_time()
    try:return operation()
    finally:phases[name]={'wall_seconds':perf_counter()-w,'cpu_seconds':process_time()-c}


def immutable(matrix):return tuple(tuple(x for x in row) for row in matrix)


def prepare(case,count,context):
    _,full=load();q=context['q'];ss=full.generated_powers(q,count)
    M,H,R=tuple(immutable(x) for x in full.trial_matrices(q,ss))
    C=immutable(full.cancellation_transform(q,ss))
    mr,hr,rr=tuple(immutable(full.congruence(x,C)) for x in (M,H,R))
    return {'q':q,'ss':tuple(ss),'M':M,'H':H,'R':R,'C':C,'mr':mr,'hr':hr,'rr':rr,
            'gap':context['gap'],'ground':context['ground'],'tolerance':case['task']['tolerance']}


def rational_candidate(p,phases):
    """Exact C2 loop; only common validation/matrix preparation moved out."""
    _,full=load();q,ss=p['q'],p['ss'];M,H,R=p['M'],p['H'],p['R']
    C,mr,hr,rr=p['C'],p['mr'],p['hr'],p['rr']
    beta=F(p['gap']['beta']);shift=F(p['ground']['lower'])-F(1,1024)
    epsilon=F(p['tolerance'])*F(15,16)
    solve=full.backend.enriched._ldl_solver([
        [h-shift*mr[i][j] for j,h in enumerate(row)] for i,row in enumerate(hr)])
    denominator=1<<112;v=[F(i==0) for i in range(len(mr))]
    def rounded(x):
        scale=max(map(abs,x))
        if not scale:raise ArithmeticError('Zero inverse iteration candidate')
        return [F(round(z/scale*denominator),denominator) for z in x]
    def stats(w):return full.statistics(q,ss,full.matvec(C,w),(M,H,R))
    used=0;trace=[]
    for used in range(1,25):
        v=rounded(solve(full.matvec(mr,v)))
        if used%4==0:
            current=stats(v);mu,variance=F(current['rayleigh']),F(current['residual_squared'])
            if mu<beta and variance/(beta-mu)<=epsilon:break
    rw,rc=perf_counter(),process_time()
    for _ in range(2):
        before=stats(v);mu,variance=F(before['rayleigh']),F(before['residual_squared'])
        if mu<beta and variance/(beta-mu)<=epsilon:break
        w,c=perf_counter(),process_time()
        kappa=max(F(1,1<<24),variance)
        residual=[[rr[i][j]-2*mu*hr[i][j]+(mu*mu+kappa)*mr[i][j]
                   for j in range(len(mr))] for i in range(len(mr))]
        candidate=rounded(full.backend.enriched._ldl_solver(residual)(full.matvec(mr,v)))
        after=stats(candidate);amu,avar=F(after['rayleigh']),F(after['residual_squared'])
        accepted=amu<beta and (mu>=beta or avar/(beta-amu)<variance/(beta-mu))
        trace.append({'before_residual_squared':str(variance),'after_residual_squared':str(avar),
                      'accepted':accepted,'wall_seconds':perf_counter()-w,'cpu_seconds':process_time()-c})
        if not accepted:break
        v=candidate
    phases['residual_refinement_subset_of_candidate']={'wall_seconds':perf_counter()-rw,'cpu_seconds':process_time()-rc}
    return {'powers':list(map(str,ss)),'coefficients':list(map(str,full.matvec(C,v))),
            'statistics':stats(v),'iterations':used,'precision_bits':112,
            'inverse_shift':str(shift),'verified_search_beta':str(beta),
            'residual_refinement':trace,'refinement_attempts':len(trace),
            'refinement_accepted':sum(x['accepted'] for x in trace)}


def numpy_prepare(p):
    import numpy as np
    return np.array(p['mr'],dtype=np.float64),np.array(p['hr'],dtype=np.float64)


def numpy_candidate(p,arrays,state):
    import numpy as np
    mass,energy=arrays
    state['stage']='cholesky'
    L=np.linalg.cholesky(mass)
    state['stage']='whitening'
    left=np.linalg.solve(L,energy)
    whitened=np.linalg.solve(L,left.T).T
    whitened=(whitened+whitened.T)/2
    state['stage']='eigh'
    values,vectors=np.linalg.eigh(whitened)
    state['stage']='coefficient_recovery'
    vector=np.linalg.solve(L.T,vectors[:,0])
    scale=np.max(np.abs(vector))
    if not np.isfinite(scale) or scale==0:raise ArithmeticError('Nonfinite or zero floating candidate')
    vector=vector/scale
    if not np.all(np.isfinite(vector)):raise ArithmeticError('Nonfinite floating coefficients')
    rational=tuple(F.from_float(float(x)) for x in vector)
    coefficients=_full.matvec(p['C'],rational)
    return {'powers':list(map(str,p['ss'])),'coefficients':list(map(str,coefficients)),
            'free_coefficients_binary64_as_exact_rationals':list(map(str,rational)),
            'approximate_ritz_value_not_a_bound':float(values[0]),
            'precision_bits':53,'refinement_attempts':0,'refinement_accepted':0}


def _worker(output_path,arm,p,task):
    wall,cpu=perf_counter(),process_time();phases={};state={'stage':'prepare'}
    result={'arm':arm,'phases':phases,'candidate':None,'certificate':None,'verification':None}
    try:
        if arm==ARMS[0]:
            phases['representation_prepare']={'wall_seconds':0.0,'cpu_seconds':0.0,
                'note':'References common immutable exact matrices; preparation counted once outside both arms.'}
            result['candidate']=timed(phases,'candidate',lambda:rational_candidate(p,phases))
        else:
            arrays=timed(phases,'representation_prepare',lambda:numpy_prepare(p))
            result['candidate']=timed(phases,'candidate',lambda:numpy_candidate(p,arrays,state))
        state['stage']='complete_certificate'
        candidate=result['candidate']
        result['certificate']=timed(phases,'complete_certificate',lambda:_full.certificate(
            p['q'],candidate['powers'],candidate['coefficients'],p['gap'],p['tolerance']))
        state['stage']='task_bound_verify'
        result['verification']=timed(phases,'task_bound_verify',lambda:_extension.verify(result['certificate'],task))
        check=result['verification']
        result['outcome']=('target_met' if check['target_met'] else
                           'valid_open' if check['certificate_valid'] else 'verification_failed')
    except Exception as error:
        message=str(error)
        if state['stage']=='cholesky':outcome='cholesky_failed'
        elif 'Rayleigh quotient >= second-eigenvalue bound' in message:outcome='mu_not_below_proved_beta'
        elif state['stage'] in ('complete_certificate','task_bound_verify'):outcome='residual_or_certificate_failed'
        else:outcome='candidate_failed'
        result.update(outcome=outcome,error=type(error).__name__,reason=message,failed_stage=state['stage'])
    result['child_total']={'wall_seconds':perf_counter()-wall,'cpu_seconds':process_time()-cpu,
                          'scope':'Excludes final JSON serialization/write; included in process_pipeline.'}
    encoded=json.dumps(result,allow_nan=False)
    temporary=output_path.with_suffix('.partial')
    temporary.write_text(encoded+'\n');temporary.replace(output_path)


def safe_candidate(arm,p,task,record_id):
    ctx=mp.get_context('fork')
    output_path=ROOT/'worker_responses'/(record_id+'.json')
    before=resource.getrusage(resource.RUSAGE_CHILDREN)
    wall,cpu=perf_counter(),process_time()
    deadline=wall+10
    process=ctx.Process(target=_worker,args=(output_path,arm,p,task));process.start()
    process.join(max(0.0,deadline-perf_counter()))
    result=None
    expired=process.is_alive() or perf_counter()>deadline
    if process.is_alive():
        process.terminate();process.join(0.25)
        if process.is_alive():process.kill();process.join(0.25)
    if not expired and process.exitcode==0 and output_path.exists():
        try:result=json.loads(output_path.read_text())
        except Exception as error:
            result={'arm':arm,'outcome':'worker_response_failed','candidate':None,'certificate':None,
                    'verification':None,'phases':{},'reason':type(error).__name__+': '+str(error)}
    expired=expired or perf_counter()>deadline
    if result is None or expired:
        result={'arm':arm,'outcome':'safety_timeout' if expired else 'worker_exit_failed',
                'candidate':None,'certificate':None,'verification':None,'phases':{},
                'reason':'Ten-second absolute process deadline exceeded.' if expired else 'Worker exited without a complete successful response.',
                'raw_response_path':str(output_path) if output_path.exists() else None}
    after=resource.getrusage(resource.RUSAGE_CHILDREN)
    result['process_pipeline']={'wall_seconds':perf_counter()-wall,'parent_cpu_seconds':process_time()-cpu,
        'child_user_system_cpu_seconds':after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
        'exit_code':process.exitcode,'safety_limit_seconds':10,'within_deadline':not expired,
        'cleanup_note':'Deadline enforcement may include up to 0.5 seconds of termination/reaping overhead; such rows cannot count as target_met.'}
    return result


def context_for(case):
    extension,full=load();_,_,sg,forms=extension._modules()
    q=full.normalize(case['task']['function']);run=extension.Run(120)
    form=forms.select_form(q,16)
    kernel=sg._kernel(q,0,False,16,16,40,form,full.verify_form)
    gap,ground=extension._coarse_gap_and_shift(kernel,form,run)
    return {'q':q,'gap':gap,'ground':ground,'form':form,'search_counters':run.counts}


def diagnostics(p):
    import numpy as np
    m,h=numpy_prepare(p);result={}
    for name,mat in (('mass',m),('energy',h)):
        value=float(np.linalg.cond(mat))
        result[name]={'condition_2':value if math.isfinite(value) else str(value),
                      'numerical_rank':int(np.linalg.matrix_rank(mat)),'dimension':len(mat)}
    try:
        L=np.linalg.cholesky(m);x=np.linalg.solve(L,h);k=np.linalg.solve(L,x.T).T
        result['whitening_asymmetry_max_abs']=float(np.max(np.abs(k-k.T)))
    except Exception as error:result['whitening_diagnostic_error']=type(error).__name__+': '+str(error)
    return result


def environment():
    import numpy as np
    out=io.StringIO()
    with redirect_stdout(out):np.show_config()
    return {'python_executable':sys.executable,'python_version':sys.version,'platform':platform.platform(),
        'numpy_version':np.__version__,'numpy_blas_lapack_configuration':out.getvalue(),
        'numpy_runtime_note':'numpy.linalg routines use the NumPy build shown above; no SciPy or arbitrary-precision backend.',
        'multiprocessing_start_method':'fork','threads':{k:os.environ[k] for k in
            ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS')}}


def run():
    if sha(ROOT/'PROTOCOL.md')!=PROTOCOL_HASH:raise ValueError('Protocol changed after lock')
    if (ROOT/'RESULTS.json').exists() or (ROOT/'records.jsonl').exists():raise ValueError('First complete run artifacts already exist; do not overwrite')
    start=perf_counter();extension,full=load()
    cases_all=json.loads((C2/'cases.json').read_text())
    cases=[next(x for x in cases_all if x['id']==ident) for ident in IDS]
    hashes={str(p):sha(p) for p in (ROOT/'PROTOCOL.md',Path(__file__),C2/'fullspace_extension.py',C2/'cases.json')}
    report={'protocol_sha256':PROTOCOL_HASH,'source_hashes':hashes,'environment':environment(),
            'cases':cases,'common_contexts':[],'common_matrices':[],'records':[],
            'scope':'Same reduced exact trial space; 53bit Ritz vs112bit residual-refined candidate; no end-to-end or SOTA ranking.'}
    (ROOT/'ENVIRONMENT.json').write_text(json.dumps(report['environment'],indent=2)+'\n')
    (ROOT/'worker_responses').mkdir(exist_ok=False)
    smoke_done=False
    diagnostic_inputs=[]
    for case_index,case in enumerate(cases):
        phases={};context=None;preparation_error=None
        try:context=timed(phases,'shared_context',lambda:context_for(case))
        except Exception as error:preparation_error={'outcome':'shared_context_failed','error':type(error).__name__,'reason':str(error)}
        report['common_contexts'].append({'case_id':case['id'],'phases':phases,'context':context,'failure':preparation_error})
        for power_index,count in enumerate(COUNTS):
            phases={};prepared=None;matrix_error=preparation_error;encoded=''
            if context is not None:
                try:prepared=timed(phases,'shared_exact_matrices',lambda:prepare(case,count,context))
                except Exception as error:matrix_error={'outcome':'shared_matrices_failed','error':type(error).__name__,'reason':str(error)}
            if prepared is not None:
                encoded=json.dumps({k:[[str(x) for x in row] for row in prepared[k]] for k in ('M','H','R','C','mr','hr','rr')})
            report['common_matrices'].append({'case_id':case['id'],'powers_count':count,'phases':phases,
                'prepared_matrix_sha256':hashlib.sha256(encoded.encode()).hexdigest() if prepared else None,
                'serialized_bytes':len(encoded.encode()),'failure':matrix_error,
                'powers':list(map(str,prepared['ss'])) if prepared else None,
                'reduced_dimension':len(prepared['mr']) if prepared else None})
            (ROOT/'PREPARATION.json').write_text(json.dumps({'common_contexts':report['common_contexts'],
                'common_matrices':report['common_matrices']},indent=2)+'\n')
            if prepared is not None:diagnostic_inputs.append((case['id'],count,prepared))
            if prepared is not None and not smoke_done:
                smoke_start=perf_counter();local=rational_candidate(prepared,{})
                original=extension._candidate_worker({'function':prepared['q'],'powers':list(map(str,prepared['ss'])),
                    'gap':prepared['gap'],'shift_source':prepared['ground'],'tolerance':prepared['tolerance']})
                check={'powers_equal':local['powers']==original['powers'],
                    'coefficients_equal':local['coefficients']==original['coefficients'],
                    'statistics_equal':local['statistics']==original['statistics'],
                    'iterations_equal':local['iterations']==original['iterations'],
                    'elapsed_seconds':perf_counter()-smoke_start,'case_id':case['id'],'powers_count':count}
                if not all(check[k] for k in ('powers_equal','coefficients_equal','statistics_equal','iterations_equal')):
                    raise ArithmeticError('Adapted C2 candidate differs from original')
                report['single_original_c2_equivalence_check']=check;smoke_done=True
            for repeat in range(3):
                order=ARMS if (case_index+power_index+repeat)%2==0 else tuple(reversed(ARMS))
                for arm in order:
                    record_id=f'{case["id"]}_{count}_{repeat+1}_{arm}'
                    result=(safe_candidate(arm,prepared,case['task'],record_id) if prepared is not None else
                        {'arm':arm,'candidate':None,'certificate':None,'verification':None,'phases':{},
                         'process_pipeline':None,'candidate_not_started':True,**matrix_error})
                    result.update(case_id=case['id'],powers_count=count,repeat=repeat+1,
                        order_in_pair=order.index(arm)+1,shared_matrix_sha256=report['common_matrices'][-1]['prepared_matrix_sha256'])
                    report['records'].append(result)
                    with (ROOT/'records.jsonl').open('a') as stream:stream.write(json.dumps(result,allow_nan=False)+'\n')
                    print(json.dumps({'completed':len(report['records']),'case':case['id'],'powers':count,'arm':arm,'outcome':result['outcome']}),flush=True)
    # Auxiliary numerical ranks/conditions are evaluated only after all
    # candidate results are locked, and cannot influence their selection.
    report['auxiliary_diagnostics']=[]
    for ident,count,p in diagnostic_inputs:
        phases={}
        try:diag=timed(phases,'auxiliary_only',lambda:diagnostics(p))
        except Exception as error:diag={'error':type(error).__name__,'reason':str(error)}
        report['auxiliary_diagnostics'].append({'case_id':ident,'powers_count':count,'phases':phases,'diagnostics':diag})
    assert len(report['records'])==72
    assert all(sha(path)==value for path,value in hashes.items())
    report['elapsed_whole_experiment_seconds']=perf_counter()-start
    (ROOT/'RESULTS.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':run()
