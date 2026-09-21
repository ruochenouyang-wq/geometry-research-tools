"""V317 acceptance revision: same mathematical tools, scoped replay reuse."""
import hashlib
import json
import os
import sys
from time import perf_counter

from release_backend import ROOT,PREVIOUS,backend,previous_service
from request_runtime import RequestRuntime


class ResearchService:
    def __init__(self,evidence_root=None,cache=True,replay_reuse=True):
        self.runtime=RequestRuntime(self.code_version(),reuse=replay_reuse)
        self.runtime.service.ResearchService.code_version=staticmethod(self.code_version)
        self.service=self.runtime.service.ResearchService(evidence_root or ROOT/'runtime_evidence',cache=cache)
        self.service._singular=self.runtime.solve_singular
        self.service._spectrum=self.runtime.solve_spectrum

    @staticmethod
    def code_version():
        files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(ROOT.glob('*.py'))}
        files['v317_freeze']=hashlib.sha256((PREVIOUS/'evaluation/freeze_manifest.json').read_bytes()).hexdigest()
        return 'geometry-v317-release:'+backend.digest(files)

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
