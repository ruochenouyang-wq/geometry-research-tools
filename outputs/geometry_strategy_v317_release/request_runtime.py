"""Request-local exact replay reuse without mutating any shared backend module.

Private module instances execute the frozen source unchanged. Only their
dependency references are routed through this owned runtime. The cache keys
contain complete canonical evidence and all supplied bindings, not object ids.
"""
from contextlib import contextmanager
from contextvars import ContextVar
from copy import deepcopy
import importlib
import importlib.util
import json
from time import perf_counter
from types import SimpleNamespace

from release_backend import PREVIOUS,backend


def json_native(value):
    """Cache only exact JSON types; never silently coerce tuples or dict keys."""
    if type(value) in (str,int,float,bool,type(None)):return True
    if type(value) is list:return all(json_native(item) for item in value)
    if type(value) is dict:
        return all(type(key) is str and json_native(item) for key,item in value.items())
    return False


class ReplaySession:
    def __init__(self,version,max_entries=256,max_bytes=16*1024*1024):
        self.version=version
        self.max_entries=max_entries
        self.max_bytes=max_bytes
        self.entries={}
        self.bytes=0
        self.closed=False
        self.counts={'requests':0,'hits':0,'misses':0,'uncacheable':0,'rejected':0,
                     'storage_skips':0,'key_seconds':0.0,'replay_seconds':0.0}

    def call(self,name,checker,certificate,bindings):
        self.counts['requests']+=1
        started=perf_counter()
        try:
            if not json_native(certificate) or not json_native(bindings):
                raise TypeError('Non-JSON values require original-type replay')
            key=backend.canonical({'version':self.version,'checker':name,
                                   'certificate':certificate,'bindings':bindings})
        except (TypeError,ValueError,RecursionError):
            self.counts['uncacheable']+=1
            try:
                owned_certificate,owned_bindings=deepcopy(certificate),deepcopy(bindings)
            except RecursionError:
                self.counts['rejected']+=1
                return False
            return checker(owned_certificate,**owned_bindings) is True
        finally:
            self.counts['key_seconds']+=perf_counter()-started
        if key in self.entries:
            self.counts['hits']+=1
            return True
        self.counts['misses']+=1
        owned=json.loads(key)
        started=perf_counter()
        try:
            valid=checker(owned['certificate'],**owned['bindings']) is True
        finally:
            self.counts['replay_seconds']+=perf_counter()-started
        if not valid:
            self.counts['rejected']+=1
            return False
        size=len(key.encode('utf-8'))
        if not self.closed and len(self.entries)<self.max_entries and self.bytes+size<=self.max_bytes:
            self.entries[key]=True
            self.bytes+=size
        else:
            self.counts['storage_skips']+=1
        return True

    def snapshot(self):
        return {**self.counts,'entries':len(self.entries),'canonical_key_bytes':self.bytes,
                'scope':'current request only; complete content and bindings; successful replays only',
                'memory_limit_scope':'canonical keys, excludes Python object overhead',
                'timing_note':'Nested checker times overlap; use complete service/evaluator elapsed time for totals.',
                'version':self.version}


def private_module(name,filename):
    spec=importlib.util.spec_from_file_location(name,PREVIOUS/filename)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ImportFacade:
    def __init__(self,adaptive):self.adaptive=adaptive
    def import_module(self,name,*args,**kwargs):
        return self.adaptive if name=='adaptive_singular' else importlib.import_module(name,*args,**kwargs)


class RequestRuntime:
    def __init__(self,version,reuse=True):
        self.version=version
        self.reuse=reuse
        self.current=ContextVar('geometry_replay_session',default=None)
        self.direct=private_module('_owned_direct',backend.direct.__file__)
        self.uncached_full=self.direct.verify_full
        self.uncached_sector=self.direct.verify_sector
        self.direct.verify_full=self.verify_full
        self.direct.verify_sector=self.verify_sector
        self.fast=private_module('_owned_fast_math','fast_math.py')
        self.fast.direct=self.direct
        # Accelerated generation, with the unchanged legacy Kernel kept in
        # direct.verify_sector for every first replay of a sector.
        self.direct.full_ground=self.fast.full_ground
        self.direct.certify_sector=self.fast.certify_sector
        self.adaptive=private_module('_owned_adaptive_singular','adaptive_singular.py')
        self.adaptive.direct=self.direct
        self.verification=private_module('_owned_verification','verification.py')
        self.verification.backend=SimpleNamespace(**{**vars(backend),'direct':self.direct})
        self.verification.importlib=ImportFacade(self.adaptive)
        self.uncached_verify=self.verification.verify_any
        self.verification.verify_any=self.verify
        self.service=private_module('_owned_research_service','research_service.py')
        self.service.verify_any=self.verify
        self.service.assess=self.verification.assess

    @contextmanager
    def session(self):
        # Nested operations within one request reuse the same owned session.
        active=self.current.get()
        if active is not None and not active.closed:
            yield active
            return
        active=ReplaySession(self.version)
        token=self.current.set(active)
        try:
            yield active
        finally:
            active.closed=True
            active.entries.clear()
            active.bytes=0
            self.current.reset(token)

    def call(self,name,checker,certificate,bindings):
        session=self.current.get()
        if session is None or session.closed or not self.reuse:
            try:
                owned_certificate,owned_bindings=deepcopy(certificate),deepcopy(bindings)
            except RecursionError:
                return False
            return checker(owned_certificate,**owned_bindings) is True
        return session.call(name,checker,certificate,bindings)

    def verify_full(self,certificate,expected_function=None,expected_mean_zero=None,
                    expected_tolerance=None):
        return self.call('frozen_full',self.uncached_full,certificate,
             {'expected_function':expected_function,'expected_mean_zero':expected_mean_zero,
              'expected_tolerance':expected_tolerance})

    def verify_sector(self,certificate,expected_function=None,expected_m=None,
                      expected_mean_zero=None):
        return self.call('frozen_sector',self.uncached_sector,certificate,
             {'expected_function':expected_function,'expected_m':expected_m,
              'expected_mean_zero':expected_mean_zero})

    def verify(self,certificate,expected_function=None,expected_mean_zero=None,
               expected_tolerance=None,expected_kind=None,_depth=0):
        # Kind checking is cheap and always performed, including on cache hits.
        kind=self.verification.conclusion_kind(certificate)
        if kind is None or (expected_kind is not None and expected_kind!=kind):return False
        return self.call('registered_bound_replay',self.uncached_verify,certificate,
            {'expected_function':expected_function,'expected_mean_zero':expected_mean_zero,
             'expected_tolerance':expected_tolerance,'_depth':_depth})

    def solve_singular(self,action):
        args={k:v for k,v in action.items() if k!='op'}
        args.setdefault('shift_policy','certified')
        return self.adaptive.solve(**args)

    def solve_spectrum(self,action):
        return {'ok':True,'certificate':self.fast.full_ground(**{k:v for k,v in action.items() if k!='op'})}
