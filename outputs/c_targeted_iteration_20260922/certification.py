"""Task-local, unverified assembly of the unchanged full-space certificate.

This is a candidate builder, deliberately not a verifier.  Search matrices and
gaps are snapshotted and bound to the task, but their mathematical correctness
is established ONLY by the controller's serialized, frozen baseline.assess.
There is no accepted/valid shortcut, certificate cache, or global monkeypatch.
"""
from collections.abc import Mapping
from dataclasses import dataclass, field
from fractions import Fraction as F
import importlib
import importlib.util
import json
from pathlib import Path
import sys

sys.dont_write_bytecode=True
FROZEN=Path(__file__).resolve().parent.parent/'geometry_accuracy80_20260921'
_MODULES=None


def _modules():
    """Cache code imports only; no task, gap, matrix or certificate is cached."""
    global _MODULES
    if _MODULES is None:
        if str(FROZEN) not in sys.path:sys.path.insert(0,str(FROZEN))
        spec=importlib.util.spec_from_file_location('_c3_certification_frozen_solver',FROZEN/'solver.py')
        baseline=importlib.util.module_from_spec(spec);spec.loader.exec_module(baseline)
        _MODULES=(baseline,importlib.import_module('fullspace_singular'),
                  importlib.import_module('spectral_gap'))
    return _MODULES


def _canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)


def _r(value):
    if isinstance(value,F):return value
    return _modules()[0].rational(value)


def _normalize_task(task):
    baseline,full,_=_modules();t=baseline.normalize_task(task)
    if t['kind']!='spectrum' or t['mean_zero'] is not False:
        raise ValueError('This builder requires a full-space spectrum task with mean_zero=False')
    t['function']=full.normalize(t['function'])
    return t


def _identity(task):
    return _canonical({key:task[key] for key in ('kind','function','mean_zero','tolerance')})


def _matrix(raw,n,name):
    if type(raw) not in (tuple,list) or len(raw)!=n:
        raise ValueError(name+' must have one row per trial power')
    rows=[]
    for row in raw:
        if type(row) not in (tuple,list) or len(row)!=n:
            raise ValueError(name+' must be square')
        rows.append(tuple(_r(value) for value in row))
    matrix=tuple(rows)
    if any(matrix[i][j]!=matrix[j][i] for i in range(n) for j in range(i)):
        raise ValueError(name+' must be exactly symmetric')
    return matrix


def _field(prepared,name):
    try:return prepared[name] if isinstance(prepared,Mapping) else getattr(prepared,name)
    except (KeyError,AttributeError) as exc:raise ValueError('Prepared trial missing '+name) from exc


@dataclass(frozen=True,eq=False)
class PendingEvaluation:
    """Immutable arithmetic receipt, never evidence of mathematical validity.

    Its owner binds it to a single context.  Even a structurally valid receipt
    must be submitted to the frozen verifier in the assembled certificate.
    """
    _owner:object=field(repr=False)
    powers:tuple
    coefficients:tuple
    matrices:tuple=field(repr=False)
    statistics_items:tuple
    reused_search_matrices:bool

    @property
    def statistics(self):return dict(self.statistics_items)


@dataclass(frozen=True,init=False,eq=False)
class UnverifiedAssemblyContext:
    _task_json:str=field(repr=False)
    _identity_key:str=field(repr=False)
    _gap_json:str=field(repr=False)
    _owner:object=field(repr=False)

    def __init__(self,task,gap):
        t=_normalize_task(task);_,full,sg=_modules()
        if not sg._json_native(gap) or type(gap) is not dict:
            raise ValueError('Gap must be a native immutable-snapshot-compatible certificate candidate')
        # Freeze before inspecting. No caller-owned dict/list survives here.
        encoded=_canonical(gap);snapshot=json.loads(encoded)
        if snapshot.get('format')!=sg.FORMAT or snapshot.get('function')!=t['function']:
            raise ValueError('Gap format or original function/axis mismatch')
        if snapshot.get('mean_zero') is not False or type(snapshot.get('azimuth_m')) is not int or snapshot['azimuth_m']!=0:
            raise ValueError('Full-space trial requires an m=0, mean_zero=False gap candidate')
        if snapshot.get('geometry')!='unit_S2' or snapshot.get('eigenvalue_index')!=2:
            raise ValueError('Expected a unit-S2 second-eigenvalue gap candidate')
        if snapshot.get('scope')!='one_real_radial_fourier_component':
            raise ValueError('Wrong gap scope')
        if _r(snapshot['beta'])!=_r(snapshot['lower']):
            raise ValueError('Inconsistent gap endpoints')
        object.__setattr__(self,'_task_json',_canonical(t))
        object.__setattr__(self,'_identity_key',_identity(t))
        object.__setattr__(self,'_gap_json',encoded)
        object.__setattr__(self,'_owner',object())

    @property
    def task(self):return json.loads(self._task_json)

    @property
    def gap(self):return json.loads(self._gap_json)

    def _bind(self,task=None):
        t=self.task if task is None else _normalize_task(task)
        if _identity(t)!=self._identity_key:
            raise ValueError('Context belongs to a different mathematical task or tolerance')
        return t

    def evaluate(self,powers,coefficients,*,prepared=None,statistics=None,task=None):
        """Snapshot matrices and exactly recompute the three scalar quadratics.

        A prepared object must bind q and powers, but that is only structural
        provenance. Its matrix entries are NOT trusted as proof input.
        """
        t=self._bind(task);_,full,_=_modules();q=t['function']
        ss=full.powers(powers)
        if type(coefficients) not in (list,tuple) or len(coefficients)!=len(ss):
            raise ValueError('One exact coefficient per trial power is required')
        vector=tuple(_r(x) for x in coefficients)
        if not any(vector):raise ValueError('Nonzero trial required')
        if prepared is None:
            raw=full.trial_matrices(q,ss)
        else:
            if _field(prepared,'function_key')!=full.canonical(q):
                raise ValueError('Prepared matrices belong to another function or axis')
            if full.powers(_field(prepared,'powers'))!=ss:
                raise ValueError('Prepared matrices belong to another trial basis')
            raw=tuple(_field(prepared,name) for name in ('M','H','R'))
        matrices=tuple(_matrix(mat,len(ss),name) for mat,name in zip(raw,('M','H','R')))
        stats=full.statistics(q,ss,vector,matrices)
        if statistics is not None:
            if type(statistics) is not dict or _canonical(statistics)!=_canonical(stats):
                raise ValueError('Candidate statistics disagree with exact quadratic recomputation')
        return PendingEvaluation(self._owner,ss,vector,matrices,tuple(stats.items()),prepared is not None)

    def assemble_evaluation(self,evaluation,*,task=None):
        """Assemble from a same-context receipt; no math-validity claim is made."""
        t=self._bind(task);_,full,_=_modules()
        if type(evaluation) is not PendingEvaluation or evaluation._owner is not self._owner:
            raise ValueError('Statistics receipt belongs to another context')
        ss,vector=evaluation.powers,evaluation.coefficients
        stats=evaluation.statistics;gap=self.gap;q=t['function']
        mu,residual,beta=F(stats['rayleigh']),F(stats['residual_squared']),F(gap['beta'])
        if residual<0 or mu>=beta:
            raise ValueError('Candidate does not meet strict Temple scalar conditions')
        lower=mu-residual/(beta-mu);width=mu-lower;tol=F(t['tolerance'])
        encoded=[[[str(x) for x in row] for row in matrix] for matrix in evaluation.matrices]
        # Deliberately byte-for-byte field-compatible with the old constructor.
        # The old status is only a candidate field until baseline.assess passes.
        cert={'format':full.FORMAT,'function':q,'geometry':'unit_S2',
              'measure':'d_sigma/(4*pi)','scope':'all_real_H1_on_unit_S2',
              'mean_zero':False,'eigenvalue_index':1,
              'trial':{'azimuth_m':0,'radial_parity':'even',
                       'formula':'sum(a_s*abs(z)^s)','powers':list(map(str,ss)),
                       'coefficients':list(map(str,vector)),
                       'strong_domain':'s=0_or_s>3/2','matrix_measure':'dz'},
              'matrix_digest':full.digest(encoded),'statistics':stats,
              'gap':gap,'gap_digest':full.digest(gap),
              'full_space_domination':json.loads(_canonical(full.DOMINATION)),
              'lower':str(lower),'upper':str(mu),'exact_width':str(width),
              'tolerance':str(tol),'status':'target_met' if width<=tol else 'certified_open',
              'full_infinite_space_covered':True,'original_function_moments_exact':True,
              'function_approximation_error':'0','formal_proof_assistant_checked':False}
        return {'certificate':cert,'certificate_valid':None,'target_met':False,
                'verification_pending':True,'candidate_target_met':width<=tol,
                'metric':str(width),'tolerance':str(tol),
                'construction_metadata':{
                    'context':'UnverifiedAssemblyContext',
                    'reused_search_matrices':evaluation.reused_search_matrices,
                    'local_gap_verifications':0,'local_final_verifications':0,
                    'external_verified_flags_trusted':False,
                    'required_acceptance':'serialize then unchanged baseline.assess against the original task'}}

    def assemble(self,powers,coefficients,prepared=None,statistics=None,*,task=None):
        evaluation=self.evaluate(powers,coefficients,prepared=prepared,statistics=statistics,task=task)
        return self.assemble_evaluation(evaluation,task=task)


def prepare_fullspace(task,gap):
    """Create a new task-local candidate context; never verify or accept it."""
    return UnverifiedAssemblyContext(task,gap)
