"""Controlled H14 diagnostic comparisons; frozen algorithms remain unchanged."""
import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter

ROOT=Path(__file__).resolve().parent
RELEASE=ROOT.parent/'geometry_strategy_v317_runtime'
sys.path.insert(0,str(RELEASE))
from release_backend import backend
from release_evaluation import timed_call
import fast_math


def main():
    raw=json.loads((RELEASE/'evaluation/holdout_20260920T042408501655Z/release_r0_output.json').read_text())
    row=next(x for x in raw['rows'] if x['case_id']=='H14')
    task=row['task']
    configs=[('recorded_configuration',8,16,20),('endpoint_precision_only',8,16,48),
             ('near_tail_only',8,32,48),('basis_size_only',16,16,48)]
    records=[]
    for name,modes,near,bits in configs:
        args={'function':task['function'],'mean_zero':True,'modes':modes,'near_tail':near,
              'bits':bits,'max_m':2,'tolerance':task['tolerance']}
        def work():
            cert=fast_math.full_ground(**args)
            assert backend.direct.verify_full(cert,task['function'],True,task['tolerance'])
            return cert
        result=timed_call(work,10)
        cert=result.pop('response',None)
        rec={'name':name,'configuration':args,'execution':result}
        if cert:
            path=ROOT/('step_'+name+'.json')
            path.write_text(json.dumps(cert,indent=2)+'\n')
            rec.update(certificate_file=path.name,certificate_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                       width=cert['exact_width'],width_float=float(backend.F(cert['exact_width'])),
                       lower_float=float(backend.F(cert['lower'])),upper_float=float(backend.F(cert['upper'])),
                       angular_tail_lower_float=float(backend.F(cert['angular_tail_lower'])),
                       sectors=[{'m':x['azimuth_m'],'lower':float(backend.F(x['lower'])),
                                 'upper':float(backend.F(x['upper'])),
                                 'width':float(backend.F(x['exact_width']))} for x in cert['sectors']],
                       independently_replayed=True)
        records.append(rec)
        print(name,rec.get('width_float'),result.get('wall_seconds'),flush=True)
    out={'scope':'Post-holdout development diagnosis, not new independent evaluation or an algorithm version.',
         'task':task,'original_stop_reason':row['attempts'][0]['response']['stop_reason'],
         'single_call_timeout_seconds':10,'records':records}
    (ROOT/'step_diagnostics.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
