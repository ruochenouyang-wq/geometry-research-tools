from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from statistics import median
import v05_range as v
import v04_search as previous

def main():
    root=Path(__file__).resolve().parent/'results'/'v05';root.mkdir(parents=True,exist_ok=True)
    rows=[]
    cases=[('positive_quadratic',[0,0,100]),
           ('shifted_quartic',['100/81','-400/27','200/3','-400/3',100]),
           ('perturbed',[0,20,-100,0,1])]
    for name,q in cases:
        for method in ['coefficient_range','certified_range']:
            samples=[]
            for _ in range(3):
                start=perf_counter()
                r=previous.search(q,max_modes=24) if method=='coefficient_range' else v.search(q,max_modes=24)
                assert previous.v3.verify(r['certificate']);samples.append(perf_counter()-start)
            v.e.base.write_json(root/(name+'_'+method+'.json'),r)
            c=r['certificate'];row={'case':name,'q':q,'method':method,'status':r['status'],'modes':c['modes'],
                                   'median_seconds':median(samples),'samples':samples,'exact_width':c['exact_width']}
            rows.append(row);print(name,method,r['status'],c['modes'],round(1000*median(samples),3),flush=True)
    q=[F(1,9),0,F(-2,3),0,1];proof=v.make_range(q)
    v.e.base.write_json(root/'irrational_extrema.json',{'method':'polynomial_range_document_v5',
                        'q_coefficients':[str(x) for x in q],'range_proof':proof})
    v.e.base.write_json(root/'benchmark.json',{'version':5,'rows':rows,'repetitions':3,
        'range_example':{'polynomial':'(t^2-1/3)^2','leaves':len(proof['intervals']),
                         'lower':proof['lower'],'upper':proof['upper'],'target_met':proof['range_target_met']}})

if __name__=='__main__':main()
