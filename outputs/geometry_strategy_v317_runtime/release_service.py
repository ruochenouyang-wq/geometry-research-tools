"""V317 acceptance revision: same mathematical tools, scoped replay reuse."""
import hashlib
import json
import os
import sys
from pathlib import Path
from threading import RLock
from time import perf_counter

from release_backend import ROOT,PREVIOUS,backend,previous_service
from request_runtime import RequestRuntime


class ResearchService:
    def __init__(self,evidence_root=None,cache=True,replay_reuse=True):
        version=self.code_version()
        self.runtime=RequestRuntime(version,reuse=replay_reuse)
        base=self.runtime.service.ResearchService
        class DemandService(base):
            def __init__(owned):
                owned._shape_cache=None
                owned._interaction_ready=False
                owned._initialization_lock=RLock()

            @property
            def shape_cache(owned):
                with owned._initialization_lock:
                    if owned._shape_cache is None:
                        import adaptive_approx
                        owned._shape_cache=adaptive_approx.ShapeCache()
                    return owned._shape_cache

            @shape_cache.setter
            def shape_cache(owned,value):
                owned._shape_cache=value

            @staticmethod
            def code_version():return version

            def handle(owned,request):
                with owned._initialization_lock:
                    if not owned._interaction_ready:
                        if storage.is_symlink() or storage.resolve()!=storage:
                            raise ValueError('Store root location changed')
                        shape_cache=owned._shape_cache
                        try:
                            base.__init__(owned,evidence_root=storage,cache=cache)
                        finally:
                            if shape_cache is not None:owned._shape_cache=shape_cache
                        owned._interaction_ready=True
                return super().handle(request)

        # Preserve constructor path validation and creation; unused evidence and
        # interaction machinery is constructed at its first measured use.
        storage=Path(evidence_root or ROOT/'runtime_evidence')
        if storage.is_symlink():raise ValueError('Store root must not be a symlink')
        storage=storage.resolve()
        storage.mkdir(parents=True,exist_ok=True)
        self.service=DemandService()
        self.service._singular=self.runtime.solve_singular
        self.service._spectrum=self.runtime.solve_spectrum

    @staticmethod
    def code_version():
        files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.glob('*.py'))}
        files['v317_freeze']=hashlib.sha256((PREVIOUS/'evaluation/freeze_manifest.json').read_bytes()).hexdigest()
        return 'geometry-v317-runtime:'+backend.digest(files)

    def _execute(self,method,request,include_replay_statistics=True):
        started=perf_counter()
        with self.runtime.session() as session:
            result=method(request)
            if include_replay_statistics:
                result['replay_reuse']=session.snapshot()
                result['replay_reuse']['enabled']=self.runtime.reuse
        return self.service._total_budget(result,request,started)

    def evaluate_task(self,task):
        return self._execute(self.service.evaluate_task,task)

    def handle(self,request):
        return self._execute(self.service.handle,request,
                             include_replay_statistics=isinstance(request,dict) and request.get('view')=='full')


def main():
    service=ResearchService(evidence_root=os.environ.get('GEOMETRY_EVIDENCE_ROOT'))
    for line in sys.stdin:
        if not line.strip():continue
        started=perf_counter()
        try:
            result=service.handle(previous_service.parse_request(line))
        except Exception as error:
            result={'ok':False,'execution_ok':False,'status':'request_failed',
                    'error':type(error).__name__,'reason':str(error)}
        result['total_elapsed_seconds']=perf_counter()-started
        print(json.dumps(result,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)


if __name__=='__main__':main()
