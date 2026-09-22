"""C3 strong search: optional float proposals, rational fallback, deferred assembly.

No frozen verifier is replaced. The radial tail barrier is never supplied as
the second-eigenvalue bound. This component does not run a baseline fallback.
"""
from fractions import Fraction as F
import importlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys
from time import perf_counter, process_time

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
FROZEN = ROOT.parent/'geometry_accuracy80_20260921'
VERSION = 'c3-hybrid-proved-gap-fullspace'
_LOADED = None


def _modules():
    global _LOADED
    if _LOADED is None:
        if str(FROZEN) not in sys.path:
            sys.path.insert(0,str(FROZEN))
        spec = importlib.util.spec_from_file_location('_c3_strong_frozen_solver',FROZEN/'solver.py')
        baseline = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(baseline)
        _LOADED = (baseline,importlib.import_module('fullspace_singular'),
                   importlib.import_module('spectral_gap'),importlib.import_module('power_forms'))
    return _LOADED


class SearchBudget(Exception):
    pass


class Run:
    def __init__(self,seconds,progress=None):
        self.started = perf_counter()
        self.seconds = seconds
        self.reserve = 0.0
        self.search_end = self.started+seconds
        self.progress = progress
        self.attempts = []
        self.prepared_trials = {}
        self.task = None
        self.counts = {'kernel_constructions':0,'search_inertia_checks':0,
                       'coarse_gap_steps':0,'threshold_gap_checks':0,
                       'candidate_processes':0,'candidate_timeouts':0,
                       'candidate_iterations':0,'final_frozen_replays':0,
                       'budget_recovery_assemblies':0,'gap_priority_upgrades':0,
                       'fast_candidate_attempts':0,'fast_candidate_accepted':0,
                       'rational_after_fast':0,'deferred_assemblies':0}

    def remaining(self):
        return max(0.0,self.seconds-(perf_counter()-self.started))

    def available(self):
        return max(0.0,self.search_end-perf_counter())

    def check(self):
        if self.available() <= 0:
            raise SearchBudget('Component allocation exhausted; caller retains its separate fallback reserve')

    def emit(self,event,action,**details):
        if self.progress is not None:
            self.progress({'event':event,'action':action,'component':'fullspace_extension',
                'elapsed_seconds':perf_counter()-self.started,'remaining_seconds':self.remaining(),**details})

    def action(self,name,callback,*,permit_expired=False,**details):
        if not permit_expired:self.check()
        start=perf_counter();before=self.counts['search_inertia_checks']
        row={'action':name,'permit_expired_allocation':permit_expired,**details}
        self.emit('begin',name,**details)
        try:
            result=callback();row['status']='returned';return result
        except Exception as error:
            row.update(status='failed',error=type(error).__name__,reason=str(error));raise
        finally:
            row['elapsed_seconds']=perf_counter()-start
            row['search_inertia_checks']=self.counts['search_inertia_checks']-before
            self.attempts.append(row)
            self.emit('error' if row.get('status')=='failed' else 'end',name,detail=row)


def _count(k,x,run):
    run.check()
    matrix=k.matrix(x,'lower')
    run.counts['search_inertia_checks']+=1
    return _modules()[2].direct.base.inertia(matrix)


def _gap_at(k,lower,form,run):
    if F(lower)>=k.beta:
        raise ValueError('Requested threshold reaches the strict radial tail barrier')
    run.check()
    run.counts['search_inertia_checks']+=1
    return _modules()[2]._certificate(k,lower,form)


def _coarse_gap_and_shift(k,form,run):
    """Eight exact second-index steps plus a separately proved ground shift."""
    _,full,sg,_=_modules()
    analytic=sg.direct.tail_bound(k.form,0)
    distance=F(1)
    for _ in range(32):
        low=analytic-distance
        if low<k.beta and _count(k,low,run)[0]==0:
            break
        distance*=2
    else:
        raise ValueError('No coarse positive comparison bracket')
    ground_low=low
    high=k.beta
    for _ in range(8):
        mid=(low+high)/2
        negative=_count(k,mid,run)[0]
        run.counts['coarse_gap_steps']+=1
        if negative<=1:
            low=mid
            if negative==0:
                ground_low=max(ground_low,mid)
        else:
            high=mid
    gap=_gap_at(k,low,form,run)
    if not sg.verify_gap(gap,k.function,0,False,form_verifier=full.verify_form):
        raise ArithmeticError('Frozen verifier rejected the coarse second-eigenvalue proof')
    # A second search only needs a reasonably close certified ground bound.
    ground_high=min(low,min(k.data['A'][i][i]/k.data['mass'][i] for i in range(k.n)))
    for _ in range(6):
        mid=(ground_low+ground_high)/2
        if _count(k,mid,run)[0]==0:
            ground_low=mid
        else:
            ground_high=mid
    upper=min(k.data['A'][i][i]/k.data['mass'][i] for i in range(k.n))
    # The frozen constructor recounts lower and upper independently.
    ground=sg._ground_certificate(k,ground_low,upper,form)
    if not sg.verify_ground(ground,k.function,0,False,form_verifier=full.verify_form):
        raise ArithmeticError('Frozen verifier rejected the ground shift source')
    return gap,ground


def _ceil_dyadic(value,bits=20):
    value=F(value);d=1<<bits
    return F(-((-value.numerator*d)//value.denominator),d)


def _candidate_worker(payload):
    """Search only; all accepted results are later rebuilt by the old verifier."""
    _,full,_,_=_modules()
    q=full.normalize(payload['function']);ss=full.powers(payload['powers'])
    gap=payload['gap'];form=gap['form_certificate']
    sg=_modules()[2]
    if not sg.verify_gap(gap,q,0,False,form_verifier=full.verify_form):
        raise ValueError('Worker requires a matching proved second-eigenvalue bound')
    ground=payload['shift_source']
    if not sg.verify_ground(ground,q,0,False,form_verifier=full.verify_form):
        raise ValueError('Worker requires a matching proved ground shift source')
    beta=F(gap['beta']);shift=F(ground['lower'])-F(1,1024)
    epsilon=F(payload['tolerance'])*F(15,16)
    M,H,R=full.trial_matrices(q,ss)
    C=full.cancellation_transform(q,ss)
    mr,hr,rr=(full.congruence(mat,C) for mat in (M,H,R))
    solve=full.backend.enriched._ldl_solver([
        [h-shift*mr[i][j] for j,h in enumerate(row)] for i,row in enumerate(hr)])
    denominator=1<<112
    v=[F(i==0) for i in range(len(mr))]
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
            if mu<beta and variance/(beta-mu)<=epsilon:
                break
    for _ in range(2):
        before=stats(v);mu,variance=F(before['rayleigh']),F(before['residual_squared'])
        if mu<beta and variance/(beta-mu)<=epsilon:break
        kappa=max(F(1,1<<24),variance)
        residual=[[rr[i][j]-2*mu*hr[i][j]+(mu*mu+kappa)*mr[i][j]
                   for j in range(len(mr))] for i in range(len(mr))]
        candidate=rounded(full.backend.enriched._ldl_solver(residual)(full.matvec(mr,v)))
        after=stats(candidate);amu,avar=F(after['rayleigh']),F(after['residual_squared'])
        accepted=amu<beta and (mu>=beta or avar/(beta-amu)<variance/(beta-mu))
        trace.append({'before_residual_squared':str(variance),'after_residual_squared':str(avar),
                      'accepted':accepted})
        if not accepted:break
        v=candidate
    return {'powers':list(map(str,ss)),'coefficients':list(map(str,full.matvec(C,v))),
            'statistics':stats(v),'iterations':used,'precision_bits':112,
            'inverse_shift':str(shift),'verified_search_beta':str(beta),
            'residual_refinement':trace}


def _candidate_seconds(available, terminal=False):
    # The caller already reserves 30% for replay/fallback. A final candidate
    # need not reserve half of its own allocation for nonexistent later trials.
    usable = available-0.15 if terminal else available*0.5
    return min(2.0,max(0.0,usable))


def _run_rational_candidate(q,powers,gap,ground,tolerance,run,*,terminal=False):
    timeout=_candidate_seconds(run.available(),terminal)
    if timeout<0.08:raise SearchBudget('Insufficient protected budget to start another candidate')
    run.counts['candidate_processes']+=1
    payload={'function':q,'powers':list(map(str,powers)),'gap':gap,
             'shift_source':ground,'tolerance':str(tolerance)}
    try:
        process=subprocess.run([sys.executable,'-B',str(Path(__file__).resolve()),'--candidate'],
            input=json.dumps(payload),text=True,capture_output=True,timeout=float(timeout),check=False)
    except subprocess.TimeoutExpired as error:
        run.counts['candidate_timeouts']+=1
        raise SearchBudget('Candidate hard deadline expired after '+str(float(timeout))+' seconds; child terminated') from error
    if process.returncode!=0:
        raise ValueError('Candidate child failed: '+process.stderr[-500:])
    result=json.loads(process.stdout)
    if 'error' in result:raise ValueError(result['error'])
    run.counts['candidate_iterations']+=result['iterations']
    result['hard_timeout_seconds']=float(timeout)
    result['candidate_kind']='c2_rational112'
    return result


def _run_candidate(q,powers,gap,ground,tolerance,run,*,terminal=False):
    """Only the final original replay can accept either kind of candidate."""
    if len(powers) in (6,10):
        import fast_candidate
        run.check()
        run.counts['fast_candidate_attempts'] += 1
        result = run.action('try_fast_candidate', lambda: fast_candidate.try_candidate(
            q, powers, gap, tolerance, prepared=run.prepared_trials.get(tuple(powers)),
            full=_modules()[1], deadline=run.search_end), terms=len(powers))
        prepared = result.get('prepared')
        if prepared is not None:
            run.prepared_trials[tuple(powers)] = prepared
        row = {'action':'fast_candidate_decision','status':'search_only',
               'reason':result.get('reason'),'accepted_for_assembly':bool(result.get('accepted')),
               'terms':len(powers),'phases':result.get('phases',{}),
               'candidate':result.get('candidate'),'elapsed_seconds':0.0,
               'is_certificate':False}
        run.attempts.append(row)
        run.emit('end','fast_candidate_decision',detail=row)
        if result.get('accepted') and result.get('candidate') is not None:
            run.counts['fast_candidate_accepted'] += 1
            candidate = dict(result['candidate'])
            candidate.update(candidate_kind='numpy_binary64',iterations=0,
                             hard_timeout_seconds=None,
                             worker_cpu_seconds=result.get('cpu_seconds'))
            return candidate
        run.counts['rational_after_fast'] += 1
        run.attempts.append({'action':'rational_fallback_after_fast','status':'decision',
                            'reason':result.get('reason','fast_candidate_not_accepted'),
                            'terms':len(powers),'elapsed_seconds':0.0})
    return _run_rational_candidate(q,powers,gap,ground,tolerance,run,terminal=terminal)


def verify(certificate,task):
    baseline=_modules()[0]
    return baseline.assess(certificate,baseline.normalize_task(task))


def _assemble_best(best,q,tolerance,run,permit_expired=False):
    import certification
    _,candidate,gap=best
    def assemble():
        run.counts['deferred_assemblies'] += 1
        context = certification.prepare_fullspace(run.task,gap)
        return context.assemble(candidate['powers'],candidate['coefficients'],
            prepared=run.prepared_trials.get(tuple(map(F,candidate['powers']))),
            statistics=candidate['statistics'])
    pending=run.action('recover_best_candidate_assembly' if permit_expired else 'assemble_pending_certificate',
        assemble,
        permit_expired=permit_expired)
    cert=json.loads(json.dumps(pending['certificate'],allow_nan=False))
    met=F(cert['exact_width'])<=F(tolerance)
    return cert,{'certificate_valid':None,'target_met':False,
        'candidate_target_met':met,'verification_pending':True,
        'metric':cert['exact_width'],'tolerance':tolerance}


def solve(raw_task,*,wall_seconds=None,progress=None):
    entry=perf_counter()
    if type(raw_task) is not dict:raise ValueError('Task must be a dictionary')
    seconds=(raw_task.get('budget',{}).get('wall_seconds',10) if wall_seconds is None else wall_seconds)
    if isinstance(seconds,bool) or not math.isfinite(float(seconds)) or not 0<=float(seconds)<=86400:
        raise ValueError('Invalid component wall_seconds')
    if progress is not None and not callable(progress):raise ValueError('progress must be callable')
    run=Run(float(seconds),progress)
    # Include module loading and input normalization in this component budget.
    run.started=entry;run.search_end=entry+float(seconds)
    cert=None;assessment={'certificate_valid':False,'target_met':False}
    reason='stage_budget';best=None
    try:
        run.check()
        baseline,full,sg,power_forms=run.action('load_frozen_modules',_modules)
        task=baseline.normalize_task(raw_task)
        run.task=task
        if task['kind']!='spectrum' or task.get('mean_zero') is not False:
            raise ValueError('Full-space extension requires a spectrum task with mean_zero=false')
        q=full.normalize(task['function']);tol=F(task['tolerance'])
        # This component explicitly extends the formerly skipped Holder regime.
        for modes in (12,16):
            form=run.action('verified_holder_form',lambda:power_forms.select_form(q,modes),modes=modes)
            def construct():
                run.counts['kernel_constructions']+=1
                return sg._kernel(q,0,False,modes,16,40,form,full.verify_form)
            kernel=run.action('construct_complete_kernel',construct,modes=modes)
            gap,ground=run.action('coarse_gap_and_ground_proofs',
                lambda:_coarse_gap_and_shift(kernel,form,run),modes=modes)
            # The external component experiment showed that ten powers can
            # suffice with a stronger proved gap. Try them before sixteen;
            # this is a general schedule choice, never a task-id exception.
            for terms in ((6,10,14,16) if modes==12 else (10,16)):
                candidate=run.action('bounded_candidate',lambda:_run_candidate(
                    q,full.generated_powers(q,terms),gap,ground,tol,run,
                    terminal=(modes==16 and terms==16)),terms=terms,modes=modes,
                    terminal_candidate=(modes==16 and terms==16))
                mu=F(candidate['statistics']['rayleigh']);variance=F(candidate['statistics']['residual_squared'])
                selected_gap=gap;beta=F(gap['beta'])
                width=variance/(beta-mu) if mu<beta else None
                # Commit a reconstructible open trial BEFORE any additional
                # threshold work can exhaust the allocated search budget.
                if width is not None and (best is None or width<best[0]):
                    best=(width,candidate,gap)
                if width is None or width>tol:
                    for portion in (F(15,16),F(1)):
                        required=mu+max(variance/(tol*portion),F(1,1<<24))
                        threshold=_ceil_dyadic(required)
                        try:
                            def check_threshold():
                                run.counts['threshold_gap_checks']+=1
                                return _gap_at(kernel,threshold,form,run)
                            selected_gap=run.action('demand_gap_threshold',check_threshold,
                                required_beta=str(required),threshold=str(threshold),
                                radial_tail_lower=str(kernel.beta),allocated_width=str(tol*portion))
                            beta=F(selected_gap['beta']);width=variance/(beta-mu)
                            break
                        except (ValueError,ArithmeticError):
                            continue
                row={'action':'candidate_assessment','status':'valid_search_statistics',
                     'terms':terms,'modes':modes,'mu':str(mu),'full_residual_squared':str(variance),
                     'verified_beta':str(beta),'predicted_width':str(width) if width is not None else None,
                     'reason':'candidate_above_proved_gap' if width is None else 'width_checked',
                     'elapsed_seconds':0.0,'worker_cpu_seconds':candidate.get('worker_cpu_seconds'),
                     'hard_timeout_seconds':candidate['hard_timeout_seconds'],'iterations':candidate['iterations'],
                     'candidate_kind':candidate.get('candidate_kind')}
                run.attempts.append(row)
                if width is not None and (best is None or width<best[0]):
                    best=(width,candidate,selected_gap)
                if width is not None and width<=tol:
                    reason='target_candidate';break
                if modes==12 and terms==6 and mu>=beta:
                    # A search-order heuristic, not an impossibility proof:
                    # a larger trial might lower mu enough, but first try to
                    # improve the proved comparison that blocked this trial.
                    decision='candidate_above_current_proved_gap_prioritize_stronger_gap'
                    run.counts['gap_priority_upgrades']+=1
                    run.attempts.append({'action':'prioritize_stronger_gap',
                        'status':'heuristic_schedule_change','reason':decision,
                        'mu':str(mu),'verified_beta':str(beta),'from_modes':12,'to_modes':16,
                        'skipped_terms_at_current_modes':[10,14,16],
                        'not_a_proof_that_other_candidates_are_impossible':True,
                        'elapsed_seconds':0.0})
                    run.emit('end','prioritize_stronger_gap',reason=decision,from_modes=12,to_modes=16)
                    break
            if reason=='target_candidate':break
        if best is not None:
            # Assembly reuses exact search data without claiming validity.
            # The common controller serializes and performs the original
            # full task-bound replay; this component cannot declare success.
            cert,assessment=_assemble_best(best,q,task['tolerance'],run)
            reason='candidate_ready' if assessment['candidate_target_met'] else 'certified_candidate_open'
    except SearchBudget as error:
        reason='protected_budget_stop'
        run.attempts.append({'action':'stop','status':'budget_stop','reason':str(error),'elapsed_seconds':0.0})
    except (ValueError,ArithmeticError,TypeError,KeyError) as error:
        reason='extension_failed'
        run.attempts.append({'action':'stop','status':'failed','error':type(error).__name__,
                             'reason':str(error),'elapsed_seconds':0.0})
    search_stop_reason=reason
    if reason=='protected_budget_stop' and best is not None and cert is None:
        # A later candidate timeout must not discard an earlier complete trial.
        # The caller explicitly owns a separate reserve. This one recovery may
        # exceed the component allocation, is timed, and remains subject to its
        # enclosing hard deadline. No search or retry is restarted here.
        run.counts['budget_recovery_assemblies']+=1
        try:
            cert,assessment=_assemble_best(best,q,task['tolerance'],run,permit_expired=True)
            reason='protected_budget_stop_candidate_retained'
        except Exception as error:
            reason='best_candidate_recovery_failed'
            run.attempts.append({'action':'recovery_failed','status':'failed',
                'error':type(error).__name__,'reason':str(error),'elapsed_seconds':0.0})
    recoverable=None
    if best is not None and cert is None:
        _,candidate,selected_gap=best
        recoverable={'function':q,'mean_zero':False,'tolerance':task['tolerance'],
            'powers':candidate['powers'],'coefficients':candidate['coefficients'],
            'gap':selected_gap,'requires_original_certificate_construction_and_verification':True}
    elapsed=perf_counter()-entry
    if elapsed>float(seconds):
        assessment['target_met']=False
        reason='component_budget_exceeded_candidate_retained' if cert is not None else 'budget_exceeded'
    fallback=not assessment.get('candidate_target_met',False)
    fallback_reason=(run.attempts[-1].get('reason',reason) if fallback and run.attempts else
                     reason if fallback else None)
    run.emit('end','component_complete',stop_reason=reason,fallback_recommended=fallback)
    return {'version':VERSION,'certificate':cert,**assessment,'status':reason,'stop_reason':reason,
            'search_stop_reason':search_stop_reason,'best_candidate':recoverable,
            'attempts':run.attempts,'counters':run.counts,'elapsed_seconds':elapsed,
            'remaining_seconds':max(0.0,float(seconds)-elapsed),'reserved_seconds':run.reserve,
            'fallback_recommended':fallback,'fallback_reason':fallback_reason,
            'within_budget':elapsed<=float(seconds) and float(seconds)>0,
            'counter_scope':'Local explicit search comparisons; excludes constructor and frozen replay internals.',
            'budget_scope':'Caller supplies allocated budget and owns fallback reserve; candidate child hard timeout; parent algebra cooperative; outer hard deadline required',
            'actual_model_calls':0,'actual_model_tokens':None}


if __name__=='__main__':
    if sys.argv[1:]==['--candidate']:
        start=process_time()
        try:
            result=_candidate_worker(json.load(sys.stdin));result['worker_cpu_seconds']=process_time()-start
            print(json.dumps(result,allow_nan=False))
        except Exception as error:
            print(json.dumps({'error':type(error).__name__+': '+str(error)}))
    else:
        for line in sys.stdin:
            if line.strip():print(json.dumps(solve(json.loads(line)),allow_nan=False),flush=True)
