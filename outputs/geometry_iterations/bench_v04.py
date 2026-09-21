from fractions import Fraction as F
from pathlib import Path
from time import perf_counter
from statistics import median
import v04_search as v

def main():
    root=Path(__file__).resolve().parent/'results'/'v04';root.mkdir(parents=True,exist_ok=True)
    rows=[]
    cases=[('linear',[0,1],1,8,F(1,10**10)),('quartic',[0,20,-100,0,1],1,20,F(1,10**10)),
           ('cluster_second',[0,'1/1000000',-100],2,20,F(1,10**10)),
           ('high_precision',[0,1],1,14,F(1,10**30))]
    for name,q,k,n,tolerance in cases:
        for method in ['v3_exact','v4_guided']:
            samples=[];counts=None;path=None
            for _ in range(3):
                start=perf_counter()
                if method=='v3_exact':c=v.v3.certify(q,k,n,112 if tolerance<F(1,10**20) else 48)
                else:
                    r=v.guided(q,k,n,tolerance);c=r['certificate'];counts=r['exact_count_calls'];path=r['path']
                assert v.v3.verify(c)
                samples.append(perf_counter()-start)
            v.v3.base.write_json(root/(name+'_'+method+'.json'),c)
            row={'case':name,'q':q,'k':k,'modes':n,'method':method,'tolerance':str(tolerance),
                 'target_met':F(c['exact_width'])<=tolerance,'median_seconds':median(samples),
                 'samples':samples,'exact_width':c['exact_width'],'search_count_calls':counts,'path':path}
            rows.append(row);print(name,method,'ms',round(1000*median(samples),3),'counts',counts,'met',row['target_met'],flush=True)
    v.v3.base.write_json(root/'benchmark.json',{'version':4,'rows':rows,'repetitions':3,
                        'timing':'generation plus verification; same process; no file I/O'})

if __name__=='__main__':main()
