from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from statistics import median
import v03_tail as v

def main():
    root=Path(__file__).resolve().parent/'results'/'v03';root.mkdir(parents=True,exist_ok=True)
    rows=[]
    cases=[('linear_coarse',[0,1],1,2),('quadratic',[1,-2,3],1,6),
           ('quartic',[0,20,-100,0,1],1,12),('broken_symmetry_cluster',[0,'1/1000000',-100],2,20)]
    for name,q,k,n in cases:
        for policy in v.POLICIES:
            samples=[]
            for _ in range(3):
                t=perf_counter();c=v.certify(q,k,n,48,policy);assert v.verify(c);samples.append(perf_counter()-t)
            v.base.write_json(root/(name+'_'+policy+'.json'),c)
            rows.append({'case':name,'q':q,'k':k,'modes':n,'policy':policy,
                         'median_seconds':median(samples),'samples':samples,'exact_width':c['exact_width']})
            print(name,policy,'width',float(F(c['exact_width'])),'ms',round(1000*median(samples),3),flush=True)
    v.base.write_json(root/'benchmark.json',{'version':3,'bits':48,'repetitions':3,
                       'timing':'generation and verification, same process, excludes file I/O','rows':rows})

if __name__=='__main__':main()
