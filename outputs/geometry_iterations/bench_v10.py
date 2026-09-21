from pathlib import Path
from statistics import median
from time import perf_counter
from fractions import Fraction as F
import v10_adaptive as v

def main():
    root=Path(__file__).resolve().parent/'results'/'v10';root.mkdir(parents=True,exist_ok=True);rows=[]
    cases=[('weak_linear',[0,'1/2'],1,24),('cubic',[2,-7,12,-3],1,28),
           ('asymmetric_quartic',[5,8,-70,2,4],2,28),('even_sextic',[0,0,-40,0,10,0,2],3,28),
           ('near_pair',[0,'1/10000',-150],2,28),('budget_failure',[0,40,-400,0,5],1,12)]
    for name,q,k,budget in cases:
        for method in ['fixed_grid','adaptive']:
            samples=[]
            for _ in range(3):
                start=perf_counter()
                even=not any(F(x) for x in q[1::2])
                if method=='fixed_grid':r=(v.v6.search if even else v.v5.search)(q,k,max_modes=budget)
                else:r=v.adaptive(q,k,max_modes=budget)
                assert v.v3.verify(r['certificate']);samples.append(perf_counter()-start)
            v.v3.base.write_json(root/(name+'_'+method+'.json'),r)
            row={'case':name,'q':q,'k':k,'max_modes':budget,'tolerance':'1/10000000000','method':method,
                 'baseline':'V6' if even else 'V5','status':r['status'],
                 'visited_modes':[x['modes'] for x in r['attempts']],
                 'exact_width':r['certificate']['exact_width'],'median_seconds':median(samples),'samples':samples}
            rows.append(row);print(name,method,row['status'],row['visited_modes'],round(1000*median(samples),3),flush=True)
    for name,q,tolerance in [('precision_route',[0,1],F(1,10**40)),('spatial_route',[0,1],F(1,10**12))]:
        r=v.function(q,tolerance,start_modes=4,max_modes=24)
        assert r['certificate'] is not None and v.v7.verify(r['certificate'])
        v.v3.base.write_json(root/(name+'.json'),r)
    v.v3.base.write_json(root/'benchmark.json',{'version':10,'repetitions':3,'rows':rows,
        'case_selection':'New fixed potentials after the scheduling rule was written; no tuning on these results. Small diagnostic suite, not a broad performance guarantee.'})

if __name__=='__main__':main()
