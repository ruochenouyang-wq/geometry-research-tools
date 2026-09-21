from pathlib import Path
from time import perf_counter
from statistics import median
import v06_parity as v

def main():
    root=Path(__file__).resolve().parent/'results'/'v06';root.mkdir(parents=True,exist_ok=True)
    rows=[]
    for name,q,k in [('double_well_first',[0,0,-100],1),('double_well_second',[0,0,-100],2),
                     ('sextic_third',[0,0,100,0,-200,0,100],3)]:
        for method in ['v5_full','v6_parity']:
            samples=[]
            for _ in range(3):
                start=perf_counter()
                r=v.v5.search(q,k,max_modes=24) if method=='v5_full' else v.search(q,k,max_modes=24)
                assert v.v3.verify(r['certificate']);samples.append(perf_counter()-start)
            v.v3.base.write_json(root/(name+'_'+method+'.json'),r)
            rows.append({'case':name,'q':q,'k':k,'method':method,'status':r['status'],
                         'modes':r['certificate']['modes'],'median_seconds':median(samples),'samples':samples,
                         'exact_width':r['certificate']['exact_width']})
            print(name,method,r['status'],round(1000*median(samples),3),flush=True)
    v.v3.base.write_json(root/'benchmark.json',{'version':6,'repetitions':3,'rows':rows})

if __name__=='__main__':main()
