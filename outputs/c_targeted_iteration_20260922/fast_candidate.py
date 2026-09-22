"""Optional binary64 trial proposals; never a certificate or proof verifier.

The controller owns rational fallback and independent frozen final replay.
Importing this module does not import NumPy or the frozen mathematical stack.
"""
from dataclasses import dataclass
from fractions import Fraction as F
import importlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from time import perf_counter, process_time

VERSION='c3-optional-binary64-candidate'
ALLOWED_COUNTS=(6,10)
_FULL=None


@dataclass(frozen=True)
class PreparedTrial:
    function_key: str
    powers: tuple
    M: tuple
    H: tuple
    R: tuple
    C: tuple
    mr: tuple
    hr: tuple


class CandidateDeadline(Exception):
    pass


class InvalidFloatCandidate(ArithmeticError):
    pass


def _resolve_full(full=None):
    global _FULL
    if full is not None:return full
    if _FULL is None:
        source=Path(__file__).resolve().parent.parent/'c_demand_extension_20260922'/'fullspace_extension.py'
        spec=importlib.util.spec_from_file_location('_c3_fast_candidate_frozen_c2',source)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        _FULL=module._modules()[1]
    return _FULL


def _exact(value):
    if isinstance(value,bool) or not isinstance(value,(F,int,str)):
        raise TypeError('Exact matrix/scalar inputs require Fraction, integer or rational string')
    return F(value)


def _matrix(raw,rows,columns,*,symmetric=False):
    if not isinstance(raw,(list,tuple)) or len(raw)!=rows:
        raise ValueError('Matrix row count does not match the trial')
    if any(not isinstance(row,(list,tuple)) or len(row)!=columns for row in raw):
        raise ValueError('Matrix column count does not match the trial')
    result=tuple(tuple(_exact(x) for x in row) for row in raw)
    if symmetric and any(result[i][j]!=result[j][i] for i in range(rows) for j in range(columns)):
        raise ValueError('Exact trial matrices must be symmetric')
    return result


def prepare_exact(function,powers,*,full=None,matrices=None,transform=None):
    """Prepare request-local search data, not independently verified evidence.

Caller-supplied exact matrices are structurally checked, not authenticated.
The frozen final verifier must rebuild them from the original function.
"""
    full=_resolve_full(full);q=full.normalize(function);ss=full.powers(powers);n=len(ss)
    raw=full.trial_matrices(q,ss) if matrices is None else matrices
    if not isinstance(raw,(list,tuple)) or len(raw)!=3:raise ValueError('Require M, H and complete R')
    M,H,R=(_matrix(x,n,n,symmetric=True) for x in raw)
    expected=full.cancellation_transform(q,ss)
    C=_matrix(expected if transform is None else transform,n,len(expected[0]))
    if C!=_matrix(expected,n,len(expected[0])):
        raise ValueError('Transform differs from the original cancellation space')
    k=len(C[0])
    mr,hr=(_matrix(full.congruence(x,C),k,k,symmetric=True) for x in (M,H))
    return PreparedTrial(full.canonical(q),tuple(ss),M,H,R,C,mr,hr)


def prepare_trial(function,powers,*,full=None):
    """Public exact-only preparation for either old or new candidate paths."""
    return prepare_exact(function,powers,full=full)


def _bind_prepared(prepared,q,ss,full):
    if not isinstance(prepared,PreparedTrial):raise TypeError('Require a PreparedTrial from this component')
    if prepared.function_key!=full.canonical(q) or prepared.powers!=tuple(ss):
        raise ValueError('Prepared matrices belong to a different function or powers')
    return prepared


def _check_deadline(deadline):
    if deadline is not None and perf_counter()>=deadline:
        raise CandidateDeadline('Fast candidate exhausted its allotted time; rational fallback belongs to the controller')


def _numpy():
    return importlib.import_module('numpy')


def _float_candidate(np,prepared,arrays,state,full):
    mass,energy=arrays
    state['substage']='cholesky';L=np.linalg.cholesky(mass)
    state['substage']='whitening'
    left=np.linalg.solve(L,energy);whitened=np.linalg.solve(L,left.T).T
    whitened=(whitened+whitened.T)/2
    state['substage']='eigh';values,vectors=np.linalg.eigh(whitened)
    state['substage']='coefficient_recovery';vector=np.linalg.solve(L.T,vectors[:,0])
    scale=np.max(np.abs(vector))
    if not np.isfinite(scale) or scale==0 or not np.isfinite(values[0]):
        raise InvalidFloatCandidate('Nonfinite eigenvalue or nonfinite/zero floating vector')
    vector=vector/scale
    if not np.all(np.isfinite(vector)):
        raise InvalidFloatCandidate('Nonfinite normalized floating vector')
    reduced=tuple(F.from_float(float(x)) for x in vector)
    coefficients=full.matvec(prepared.C,reduced)
    return {'powers':list(map(str,prepared.powers)),'coefficients':list(map(str,coefficients)),
            'precision_bits':53,'method':'numpy_cholesky_whitening_eigh',
            'approximate_ritz_value_not_a_bound':float(values[0]),
            'iterations':0,'residual_refinement':[]}


def try_candidate(function,powers,gap,tolerance,*,prepared=None,full=None,progress=None,deadline=None):
    """Return a proposal and explicit fallback decision, never certificate_valid.

accepted means exact statistics meet tolerance *conditional on caller beta*.
No gap proof is established here. All failures retain available candidate data.
deadline is an absolute perf_counter value; outer hard process limits remain
necessary because running imports/LAPACK/exact algebra are not interruptible.
"""
    started,cpu_started=perf_counter(),process_time();phases={};state={'phase':'input','substage':None}
    result={'version':VERSION,'accepted':False,'fallback_required':True,'reason':None,
            'candidate':None,'prepared':prepared,'phases':phases,
            'numpy_preloaded_at_entry':'numpy' in sys.modules,
            'screening_only':True,'gap_proof_checked_here':False}
    def emit(event,action,**details):
        if callable(progress):progress({'component':'fast_candidate','event':event,'action':action,
            'elapsed_seconds':perf_counter()-started,**details})
    def stage(name,operation):
        _check_deadline(deadline);state.update(phase=name,substage=None)
        w,c=perf_counter(),process_time();row={'status':'returned'}
        try:
            emit('begin',name);return operation()
        except Exception as error:
            row.update(status='failed',error=type(error).__name__,reason=str(error));raise
        finally:
            row.update(wall_seconds=perf_counter()-w,cpu_seconds=process_time()-c);phases[name]=row
            emit('error' if row['status']=='failed' else 'end',name,detail=row)
    try:
        if progress is not None and not callable(progress):raise TypeError('progress must be callable')
        if deadline is not None:
            if isinstance(deadline,bool) or not isinstance(deadline,(float,int)) or not math.isfinite(deadline):
                raise ValueError('deadline must be a finite absolute perf_counter value')
        _check_deadline(deadline)
        if not isinstance(powers,(list,tuple)):raise ValueError('powers must be a list or tuple')
        if len(powers) not in ALLOWED_COUNTS:
            result['reason']='unsupported_power_count';return result
        full=stage('load_frozen_full',lambda:_resolve_full(full))
        q=full.normalize(function);ss=full.powers(powers);tol=_exact(tolerance)
        if not F(1,10**30)<=tol<=1:raise ValueError('Tolerance must lie in 1e-30..1')
        if (not isinstance(gap,dict) or full.normalize(gap['function'])!=q or
                gap.get('azimuth_m')!=0 or gap.get('mean_zero') is not False):
            raise ValueError('Caller gap context must match the original function and full m=0 space')
        beta=_exact(gap['beta'])
        result['prepared']=stage('exact_prepare',lambda:
            prepare_exact(q,ss,full=full) if prepared is None else _bind_prepared(prepared,q,ss,full))
        p=result['prepared']
        np=stage('numpy_import',_numpy)
        arrays=stage('representation_prepare',lambda:(np.array(p.mr,dtype=np.float64),np.array(p.hr,dtype=np.float64)))
        result['candidate']=stage('floating_candidate',lambda:_float_candidate(np,p,arrays,state,full))
        candidate=result['candidate']
        candidate['statistics']=stage('exact_statistics',lambda:full.statistics(q,ss,candidate['coefficients'],(p.M,p.H,p.R)))
        mu=_exact(candidate['statistics']['rayleigh']);variance=_exact(candidate['statistics']['residual_squared'])
        candidate.update(screening_beta=str(beta),requested_tolerance=str(tol),screening_width=None)
        _check_deadline(deadline)
        if variance<0:raise ArithmeticError('Complete residual variance is negative')
        if mu>=beta:
            result['reason']='mu_not_below_screening_beta';return result
        width=variance/(beta-mu);candidate['screening_width']=str(width)
        if width>tol:
            result['reason']='exact_width_exceeds_tolerance';return result
        result.update(accepted=True,fallback_required=False,reason='exact_statistics_meet_requested_tolerance')
    except Exception as error:
        if isinstance(error,CandidateDeadline):reason='deadline_exceeded'
        elif isinstance(error,ImportError) and state['phase']=='numpy_import':reason='numpy_unavailable'
        elif state.get('substage')=='cholesky':reason='cholesky_failed'
        elif isinstance(error,InvalidFloatCandidate):reason='invalid_float_candidate'
        elif state['phase']=='exact_statistics':reason='exact_statistics_failed'
        else:reason=state['phase']+'_failed'
        result.update(reason=reason,error=type(error).__name__,error_detail=str(error),
                      failed_phase=state['phase'],failed_substage=state.get('substage'))
    finally:
        result.update(elapsed_seconds=perf_counter()-started,cpu_seconds=process_time()-cpu_started,
            budget_scope='Cooperative stage checks; controller retains fallback budget and outer hard deadline.',
            prepared_scope='Request-local exact search matrices; final frozen verifier must independently rebuild.')
        try:emit('end','fast_candidate_complete',accepted=result['accepted'],reason=result['reason'])
        except Exception as error:
            result.update(accepted=False,fallback_required=True,reason='progress_callback_failed',
                          progress_error=type(error).__name__+': '+str(error))
    return result
