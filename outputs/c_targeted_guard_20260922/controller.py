"""C3.1: one protected baseline attempt may improve a verified open result.

No mathematical algorithm or verifier is replaced. The fallback is restricted
to the frozen positive mean-zero abs_power route, which creates no children.
"""
from contextlib import contextmanager
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
import resource
import signal
import subprocess
import sys
import threading
from time import perf_counter, process_time

sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
C3_ROOT=ROOT.parent/'c_targeted_iteration_20260922'
BASELINE_ROOT=ROOT.parent/'geometry_accuracy80_20260921'
C1_ROOT=ROOT.parent/'geometry_parallel_v1_20260921'
VERSION='c3.1-positive-open-budget-guard'
FINAL_RESERVE=1.1
OUTPUT_RESERVE=0.25
MIN_FALLBACK=0.2
MAX_RESPONSE_CHARS=16_000_000
_C3=None
_SCOPED={}
_MISSING=object()
_LOCAL_NAMES=('runtime','certification','fast_candidate','weak_search','strong_search','mean_zero',
              'shared_baseline','solver')


class ParentReplayDeadline(BaseException):
    """Bypass frozen verifiers' ordinary Exception rejection handlers."""
    pass


def _load_absolute(path,name):
    existing=sys.modules.get(name)
    if existing is not None:
        if Path(existing.__file__).resolve()!=path.resolve():raise ImportError('Private frozen module path mismatch')
        return existing
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);sys.modules[name]=module
    spec.loader.exec_module(module)
    return module


@contextmanager
def _frozen_scope():
    """Keep C3's absolute imports local without clobbering caller bindings."""
    prior_path=list(sys.path);prior={name:sys.modules.get(name,_MISSING) for name in _LOCAL_NAMES}
    sys.path[:]=[str(C3_ROOT),str(BASELINE_ROOT),str(C1_ROOT)]+prior_path
    try:
        for name in _LOCAL_NAMES:
            if name in _SCOPED:sys.modules[name]=_SCOPED[name]
            else:sys.modules.pop(name,None)
        runtime=_load_absolute(C3_ROOT/'runtime.py','_c31_frozen_c3_runtime')
        sys.modules['runtime']=runtime
        yield
    finally:
        for name in _LOCAL_NAMES:
            if name in sys.modules:_SCOPED[name]=sys.modules[name]
            if prior[name] is _MISSING:sys.modules.pop(name,None)
            else:sys.modules[name]=prior[name]
        sys.path[:]=prior_path


def _modules():
    global _C3
    if _C3 is None:_C3=_load_absolute(C3_ROOT/'controller.py','_c31_frozen_c3_controller')
    if Path(_C3.baseline.__file__).resolve()!=BASELINE_ROOT/'solver.py':
        raise ImportError('Frozen baseline must be the original absolute solver path')
    return _C3,_C3.baseline


def _stop_direct_child(process):
    """The fixed fallback has no descendants and inherits the outer group."""
    if process.poll() is None:
        process.terminate()
        try:process.wait(timeout=0.1)
        except subprocess.TimeoutExpired:
            process.kill();process.wait(timeout=0.2)


def _baseline_child(task,seconds):
    """Start, JSON transport and return within one absolute child deadline."""
    started=perf_counter();deadline=started+seconds;process=None
    result={'outcome':'worker_failed','response':None,'error':None}
    try:
        payload=json.dumps(task,allow_nan=False)
        process=subprocess.Popen([sys.executable,'-B',str(Path(__file__).resolve()),'--baseline-child'],
            stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,
            start_new_session=False)
        output,stderr=process.communicate(payload,timeout=max(0.0,deadline-perf_counter()))
        if process.returncode!=0:
            result.update(outcome='worker_failed',error='Nonzero baseline worker exit',
                          exit_code=process.returncode,stderr=stderr)
        else:
            if len(output)>MAX_RESPONSE_CHARS:raise ValueError('Baseline response exceeds the finite JSON transport limit')
            response=json.loads(output)
            if type(response) is not dict:raise ValueError('Baseline response must be a JSON object')
            result.update(outcome='returned',response=response,stderr=stderr,exit_code=process.returncode)
        if perf_counter()>deadline:
            result.update(outcome='timeout',error='Baseline response exceeded the protected absolute deadline')
    except subprocess.TimeoutExpired:
        result.update(outcome='timeout',error='Protected baseline worker deadline expired')
    except Exception as error:
        result.update(outcome='worker_failed',error=type(error).__name__+': '+str(error))
    finally:
        if process is not None:
            try:_stop_direct_child(process)
            except Exception as error:result.update(outcome='worker_cleanup_failed',cleanup_error=type(error).__name__+': '+str(error))
            for pipe in (process.stdin,process.stdout,process.stderr):
                if pipe is not None:pipe.close()
            result['exit_code']=process.returncode
        result.update(elapsed_seconds=perf_counter()-started,allocated_seconds=seconds,
            process_scope='Frozen single-process fallback inherits outer process group; local timeout kills/reaps direct child.')
    return result


def _eligible(task,result):
    return (task['kind']=='spectrum' and task.get('mean_zero') is True and
            task['function']['profile']=='abs_power' and F(task['function']['amplitude'])>0 and
            result.get('certificate_valid') is True and result.get('target_met') is False and
            type(result.get('certificate')) is dict)


def _parent_assess(baseline,certificate,task,deadline):
    """Bound this original replay on the Unix main thread; never patch it."""
    if threading.current_thread() is not threading.main_thread():
        raise RuntimeError('Guard replay requires the Unix main thread')
    if any(signal.getitimer(signal.ITIMER_REAL)):
        raise RuntimeError('An existing process timer prevents isolated guard replay')
    seconds=deadline-perf_counter()
    if seconds<=0:raise ParentReplayDeadline('No protected time remains for parent replay')
    previous=signal.getsignal(signal.SIGALRM)
    def expired(signum,frame):
        raise ParentReplayDeadline('Parent frozen replay reached the output-reserved deadline')
    signal.signal(signal.SIGALRM,expired)
    try:
        signal.setitimer(signal.ITIMER_REAL,seconds)
        return baseline.assess(certificate,task)
    finally:
        signal.setitimer(signal.ITIMER_REAL,0)
        signal.signal(signal.SIGALRM,previous)


def solve(raw_task,progress=None):
    started,cpu_started=perf_counter(),process_time()
    child_started=resource.getrusage(resource.RUSAGE_CHILDREN)
    if threading.current_thread() is not threading.main_thread():
        raise RuntimeError('C3.1 guard requires sequential calls on the Unix main thread')
    wrapper_attempts=[];guard={'attempted':False,'replaced':False,'reason':'not_eligible',
                              'reserved_seconds':FINAL_RESERVE,'selected_origin':'frozen_c3'}
    def emit(event,action,**details):
        if progress is not None:progress({'component':'c3_guard','event':event,'action':action,
            'elapsed_seconds':perf_counter()-started,**details})
    def stage(name,operation,**details):
        w,c=perf_counter(),process_time();row={'action':name,**details}
        emit('begin',name,**details)
        try:
            value=operation();row['status']='returned';return value
        except ParentReplayDeadline as error:
            row.update(status='failed',error=type(error).__name__,reason=str(error));raise
        except Exception as error:
            row.update(status='failed',error=type(error).__name__,reason=str(error));raise
        finally:
            row.update(elapsed_seconds=perf_counter()-w,parent_cpu_seconds=process_time()-c)
            wrapper_attempts.append(row);emit('error' if row.get('status')=='failed' else 'end',name,detail=row)
    with _frozen_scope():
        c3,baseline=stage('load_frozen_c3',_modules)
        task=baseline.normalize_task(raw_task);limit=task['budget']['wall_seconds']
        deadline=started+limit
        c3_task={**task,'budget':{'wall_seconds':max(0.0,deadline-perf_counter())}}
        original=stage('frozen_c3_complete_solve',lambda:c3.solve(c3_task,progress=progress))
        result=dict(original);original_attempts=list(original.get('attempts',[]))
        guard.update(eligible=_eligible(task,original),base_status=original.get('status'),
                     base_elapsed_seconds=original.get('total_elapsed_seconds'))
        if guard['eligible']:
            try:
                original_width=F(original['certificate']['exact_width'])
                guard['original_exact_width']=str(original_width)
                available=max(0.0,deadline-perf_counter()-FINAL_RESERVE)
                if available<MIN_FALLBACK:
                    guard['reason']='insufficient_remaining_budget_for_protected_fallback'
                else:
                    guard.update(attempted=True,allocated_seconds=available,reason='protected_baseline_attempt')
                    reduced_task={**task,'budget':{'wall_seconds':available}}
                    child=stage('protected_baseline_remaining_budget',lambda:_baseline_child(reduced_task,available),
                                allocated_seconds=available,reserved_seconds=FINAL_RESERVE)
                    guard['fallback_worker']=child
                    if child['outcome']!='returned':
                        guard['reason']='fallback_'+child['outcome']
                    elif perf_counter()>=deadline:
                        guard['reason']='no_budget_for_parent_frozen_replay'
                    else:
                        response=child['response'];certificate=response.get('certificate')
                        # Child status/validity fields are not acceptance evidence.
                        candidate=json.loads(json.dumps(certificate,allow_nan=False))
                        replay_deadline=deadline-OUTPUT_RESERVE
                        assessment=stage('guard_parent_frozen_replay',lambda:_parent_assess(baseline,candidate,task,replay_deadline),
                                         output_reserve_seconds=OUTPUT_RESERVE)
                        if perf_counter()>=replay_deadline:
                            raise ParentReplayDeadline('Replay returned after the output-reserved deadline')
                        guard['fallback_assessment']=assessment
                        if assessment.get('certificate_valid') is not True:
                            guard['reason']='fallback_original_verifier_rejected'
                        else:
                            width=F(candidate['exact_width']);guard['fallback_exact_width']=str(width)
                            if width<original_width:
                                if perf_counter()>=replay_deadline:
                                    raise ParentReplayDeadline('Verified improvement could not be selected before the output-reserved deadline')
                                result.update(certificate=candidate,**assessment,
                                    route='c3_guard_baseline_improvement',failure_stage=None if assessment.get('target_met') else 'precision_gate')
                                guard.update(replaced=True,reason='strictly_narrower_verified_certificate',selected_origin='baseline')
                            else:guard['reason']='fallback_no_strict_width_improvement'
            except ParentReplayDeadline as error:
                guard.update(reason='fallback_parent_replay_timeout',error=str(error))
            except Exception as error:
                # The selected result was untouched until complete independent
                # acceptance; any failure preserves the original verified best.
                guard.update(reason='fallback_exception',error=type(error).__name__+': '+str(error))
        emit('end','guard_decision',reason=guard['reason'],attempted=guard['attempted'],replaced=guard['replaced'])
        counters=dict(original.get('counters',{}))
        counters.update(guard_baseline_calls=int(guard['attempted']),
            guard_parent_replays=sum(x['action']=='guard_parent_frozen_replay' for x in wrapper_attempts))
        result.update(version=VERSION,base_version=original.get('version'),guard=guard,
            attempts=original_attempts+wrapper_attempts,counters=counters,
            actual_model_calls=0,actual_model_tokens=None)
        # Count a complete final JSON interchange check inside solve. The CLI
        # write itself is also measured by the external whole-process evaluator.
        result=stage('response_json_boundary',lambda:json.loads(json.dumps(result,allow_nan=False)))
        result['attempts']=original_attempts+wrapper_attempts
    elapsed=perf_counter()-started;child_end=resource.getrusage(resource.RUSAGE_CHILDREN)
    parent_cpu=process_time()-cpu_started
    child_cpu=max(0.0,(child_end.ru_utime-child_started.ru_utime)+(child_end.ru_stime-child_started.ru_stime))
    within=limit>0 and elapsed<=limit
    result.update(total_elapsed_seconds=elapsed,total_cpu_seconds=parent_cpu+child_cpu,
        parent_cpu_seconds=parent_cpu,child_cpu_seconds=child_cpu,within_budget=within,
        cpu_scope='Parent process_time plus reaped-child rusage; fixed fallback is single-process.',
        budget_scope='Whole wrapper includes frozen C3, protected fallback, failures, parent replay and JSON check; outer group hard deadline still required.')
    if not within:result.update(target_met=False,status='budget_exceeded')
    return result


if __name__=='__main__':
    if sys.argv[1:]==['--baseline-child']:
        with _frozen_scope():
            _,baseline=_modules()
            print(json.dumps(baseline.solve(json.load(sys.stdin)),allow_nan=False),flush=True)
    else:
        for line in sys.stdin:
            if line.strip():print(json.dumps(solve(json.loads(line)),allow_nan=False),flush=True)
