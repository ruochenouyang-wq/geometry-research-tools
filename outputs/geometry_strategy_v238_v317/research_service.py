"""Integrated local research tool: exact evidence, automatic routes, brief results.

One JSON request per line. No network or external model calls are performed.
The integration is reviewed separately from the 80 component iterations.
"""
from copy import deepcopy
from pathlib import Path
from time import perf_counter
import hashlib
import importlib
import json
import math
import os
import sys

import backend
from verification import verify_any, assess


def parse_request(text):
    """Reject ambiguous JSON before routing old and new operations."""
    import request_contract
    value=json.loads(text,object_pairs_hook=request_contract._unique_object,
                     parse_constant=request_contract._reject_constant)
    if type(value) is not dict:
        raise ValueError('Request must be an object')
    request_contract._finite_values(value)
    return value


def validate_spectral_controls(task):
    """Validate controller controls even when an exact shortcut needs no search."""
    import research_controller
    request=research_controller.normalize({k:deepcopy(v) for k,v in task.items() if k!='kind'} |
                                         {'target':'spectrum'})
    budget=request.get('budget',{})
    policy=request.get('policy',{})
    backend.direct.integer(budget.get('max_attempts',8),1,1024,'max_attempts')
    backend.direct.integer(budget.get('work_units',2000000),1,10**12,'work_units')
    wall=budget.get('wall_seconds',60)
    if isinstance(wall,bool) or not math.isfinite(float(wall)) or not 0<=float(wall)<=86400:
        raise ValueError('wall_seconds must be finite and in 0..86400')
    backend.direct.integer(policy.get('stagnation_patience',2),1,64,'stagnation_patience')
    gain=backend.pieces.rational(policy.get('min_relative_gain','1/1000'))
    if not 0<=gain<1:
        raise ValueError('min_relative_gain must be in [0,1)')
    return request


class ResearchService:
    def __init__(self, evidence_root=None, cache=True):
        import evidence_store
        import interaction
        import adaptive_approx
        self.shape_cache = adaptive_approx.ShapeCache()
        self.store = evidence_store.EvidenceStore(
            evidence_root or backend.ROOT/'runtime_evidence', verifier=verify_any,
            verifier_version=self.code_version())
        self.runner = interaction.RequestRunner(self._execute,
            namespace=self.code_version(), cache=cache)

    @staticmethod
    def code_version():
        files={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
               for p in sorted(backend.ROOT.glob('*.py'))}
        return 'geometry-v317:'+backend.digest(files)

    @staticmethod
    def _spectrum(action):
        args={k:v for k,v in action.items() if k!='op'}
        fast=importlib.import_module('fast_math')
        solver=getattr(fast,'full_ground',backend.direct.full_ground)
        return {'ok':True,'certificate':solver(**args)}

    @staticmethod
    def _singular(action):
        import adaptive_singular
        args={k:v for k,v in action.items() if k!='op'}
        # A certified shift is a proposal aid. Final acceptance still replays
        # the original operator and every omitted component.
        args.setdefault('shift_policy','certified')
        return adaptive_singular.solve(**args)

    @staticmethod
    def _singular_supports(request):
        if request.get('mean_zero',True) is not True:
            return False
        try:
            import adaptive_singular
            q=adaptive_singular.normalize(request['function'])
            backend.direct.form_bound(q)
            return True
        except (ValueError,TypeError,KeyError,ArithmeticError):
            return False

    def evaluate_task(self, task):
        """Evaluator adapter. It returns full evidence, never a summary alone."""
        if type(task) is not dict:
            raise ValueError('Task must be an object')
        allowed={'kind','function','tolerance','mean_zero','budget','policy','id'}
        if not set(task)<=allowed:
            raise ValueError('Unknown task field')
        started=perf_counter()
        result=self._research({k:v for k,v in task.items() if k!='id'})
        return self._total_budget(result,task,started)

    @staticmethod
    def _total_budget(result,request,started):
        elapsed=perf_counter()-started
        result['total_elapsed_seconds']=elapsed
        if (request.get('target',request.get('kind'))=='spectrum' and
                request.get('op','research')=='research'):
            limit=float(request.get('budget',{}).get('wall_seconds',60))
            result['within_budget']=result.get('within_budget',True) and elapsed<=limit
            result['budget_scope']='whole_service_call_including_replay_and_storage; cooperative'
            if result['within_budget'] is False:
                result['status']='budget_exceeded'
        return result

    def _research(self, task):
        import adaptive_approx
        import research_controller
        import proof_transport
        kind=task.get('kind')
        if kind not in ('approximation','spectrum'):
            raise ValueError('Choose approximation or spectrum explicitly')
        q=backend.direct.normalize(task['function'])
        tol=backend.pieces.rational(task.get('tolerance','1/100000000'))
        if not backend.F(1,10**30)<=tol<=1:
            raise ValueError('Tolerance must be in 1e-30..1')
        mean_zero=task.get('mean_zero',True)
        if type(mean_zero) is not bool:
            raise ValueError('mean_zero must be boolean')
        budget=task.get('budget',{})
        if type(budget) is not dict:
            raise ValueError('Budget must be an object')
        if kind=='spectrum':
            spectral_request=validate_spectral_controls(task)
        if kind=='approximation':
            if 'mean_zero' in task or 'policy' in task:
                raise ValueError('Function approximation does not have a spectral constraint or controller policy')
            if 'shape_cache' in budget:
                raise ValueError('Planning cache is managed by the service')
            raw=adaptive_approx.solve(q,tolerance=str(tol),shape_cache=self.shape_cache,**budget)
            raw['route']='automatic_function_approximation'
        elif backend.F(q['amplitude'])==0:
            started=perf_counter()
            certificate=proof_transport.constant_for_function(q,mean_zero=mean_zero)
            elapsed=perf_counter()-started
            within_budget=elapsed<=float(budget.get('wall_seconds',60))
            raw={'certificate':certificate,'attempts':[],
                 'route':'exact_constant_spectrum','cost':{'solver_calls':0,'elapsed_seconds':elapsed},
                 'within_budget':within_budget,
                 'status':'target_met' if within_budget else 'budget_exceeded'}
        else:
            # Public development evidence: loose targets often need only the
            # cheaper direct probe. Keep the enriched route available if it fails.
            direct_first=(tol>=backend.F(1,1000) and
                          spectral_request.get('policy',{}).get('compare_routes',True))
            if direct_first:
                spectral_request['route']='spectrum'
            solvers={'spectrum':self._spectrum,'singular':self._singular,
                     'approximate':backend.legacy.handle,
                     'singular_supports':self._singular_supports}
            raw=research_controller.solve(spectral_request,solvers=solvers,verifier=verify_any)
            raw['route']='evidence_driven_spectral_controller'
            raw['initial_route_policy']=('direct_first_for_loose_target' if direct_first else
                                         'structure_selected')
        certificate=raw.get('certificate')
        if certificate is None:
            return {'ok':False,'execution_ok':True,'certificate_valid':False,
                    'target_met':False,**raw}
        judgment=assess(certificate,{'kind':kind,'function':q,'tolerance':str(tol),
                                    **({'mean_zero':mean_zero} if kind=='spectrum' else {})})
        if raw.get('within_budget') is False and judgment['certificate_valid']:
            judgment['status']='budget_exceeded'
        return {**raw,**judgment,'ok':judgment['certificate_valid'],
                'execution_ok':True,'certificate':certificate}

    def _legacy_route(self, request):
        op=request['op']
        if op=='spectrum':
            return self._spectrum(request)
        if op=='approximate':
            import fast_math
            certificate=fast_math.piecewise_model(**{k:v for k,v in request.items() if k!='op'})
            return {'ok':True,'kind':'function_L2_approximation_only','certificate':certificate}
        if op=='precise_singular':
            return {'ok':True,**self._singular(request)}
        return backend.legacy.handle(request)

    def _execute(self, request):
        import request_contract
        if type(request) is not dict:
            raise ValueError('Request must be an object')
        if request.get('op')=='verify':
            from verification import verify_request
            return verify_request(request)
        if request.get('op')=='research':
            allowed={'op','target','function','tolerance','mean_zero','budget','policy','request_id'}
            if not set(request)<=allowed:
                raise ValueError('Unknown research request field')
            for key in ('function','tolerance','mean_zero','budget','policy'):
                if key in request and request[key] is None:
                    raise ValueError('Explicit null research binding: '+key)
            if 'request_id' in request:
                identifier=request['request_id']
                if type(identifier) is not str or not 1<=len(identifier)<=128:
                    raise ValueError('Invalid request_id')
            task={k:deepcopy(v) for k,v in request.items() if k not in ('op','target','request_id')}
            task['kind']=request.get('target')
            result=self._research(task)
            result['request_digest']=backend.digest(request)
            if 'request_id' in request:
                result['request_id']=request['request_id']
            return result
        request_contract.CERTIFICATE_KINDS['adaptive_singular_temple_v1']='spectrum'
        return request_contract.execute(request,executor=self._legacy_route,verifier=verify_any)

    def handle(self, request):
        started=perf_counter()
        result=self._handle(request)
        return self._total_budget(result,request,started)

    def _handle(self, request):
        import interaction
        if type(request) is not dict:
            raise ValueError('Request must be an object')
        import request_contract
        request_contract._finite_values(request)
        if request.get('op')=='fetch':
            ledger=self.runner.ledger
            try:
                ledger.admit_request()
                ledger.start_execution()
                if not set(request)<={'op','reference'}:
                    raise ValueError('Unknown fetch field')
                with ledger.stage('evidence_read'):
                    certificate=self.store.get(request['reference'])
                with ledger.stage('verification'):
                    valid=verify_any(certificate)
                ledger.counts['successes' if valid else 'failures']+=1
                return {'ok':valid,'certificate':certificate,'certificate_valid':valid,
                        'binding':'unbound_evidence_replay','target_met':None,'usage':ledger.snapshot()}
            except Exception:
                ledger.counts['failures']+=1
                raise
        if request.get('op')=='batch':
            if set(request)!={'op','requests'} or type(request['requests']) is not list or len(request['requests'])>64:
                raise ValueError('Batch requires at most 64 requests')
            responses=[]
            for child in request['requests']:
                try:
                    if type(child) is dict and child.get('op')=='batch':
                        raise ValueError('Nested batches are unsupported')
                    responses.append(self.handle(child))
                except Exception as error:
                    responses.append({'ok':False,'error':type(error).__name__,'reason':str(error)})
            return {'ok':True,'responses':responses,'usage':self.runner.ledger.snapshot()}
        payload=deepcopy(request)
        view=payload.pop('view','summary')
        if view not in ('summary','full'):
            raise ValueError('view must be summary or full')
        envelope=self.runner.execute(payload)
        result=envelope['result']
        if view=='full' or 'certificate' not in result or result['certificate'] is None:
            return {**result,'interaction':{k:v for k,v in envelope.items() if k!='result'}}
        effective_request=deepcopy(payload)
        if payload.get('op') in ('research','spectrum','precise_singular'):
            effective_request.setdefault('tolerance','1/100000000')
        def bound_verifier(cert):
            kind=('approximation' if payload.get('op')=='approximate' else
                  payload.get('target') if payload.get('op')=='research' else 'spectrum')
            return verify_any(cert,expected_function=payload.get('function'),
                expected_mean_zero=payload.get('mean_zero',True) if kind=='spectrum' else None,
                expected_tolerance=effective_request.get('tolerance'),expected_kind=kind)
        def save(key,cert):
            ref=self.store.put(cert)
            if ref!=key:
                raise ArithmeticError('Evidence reference mismatch')
            return ref
        summary=interaction.present_result(result,request=effective_request,verifier=bound_verifier,
                                           store=save,ledger=self.runner.ledger)
        for key in ('execution_ok','certificate_valid','target_met','stop_reason','cost','request_id','request_digest'):
            if key in result:
                summary[key]=result[key]
        if result.get('within_budget') is not None:
            summary['within_budget']=result['within_budget']
        if result.get('status')=='budget_exceeded':
            summary['status']='budget_exceeded'
        summary['usage']=self.runner.ledger.snapshot()
        summary['interaction']={k:v for k,v in envelope.items() if k!='result'}
        return summary


def main():
    import request_contract
    service=ResearchService(evidence_root=os.environ.get('GEOMETRY_EVIDENCE_ROOT'))
    for line in sys.stdin:
        if not line.strip():
            continue
        started=perf_counter()
        try:
            result=service.handle(parse_request(line))
        except Exception as error:
            result={'ok':False,'execution_ok':False,'error':type(error).__name__,
                    'reason':str(error),'status':'request_failed'}
        result['total_elapsed_seconds']=perf_counter()-started
        print(json.dumps(result,sort_keys=True,ensure_ascii=False,allow_nan=False),flush=True)


if __name__=='__main__':main()
