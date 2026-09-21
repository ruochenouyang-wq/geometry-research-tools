from pathlib import Path
from statistics import median
from fractions import Fraction as F
import v09_family as v

def main():
    root=Path(__file__).resolve().parent/'results'/'v09';root.mkdir(parents=True,exist_ok=True);rows=[]
    for name,threshold,budget in [('true_margin',F(-3,10),64),('false_threshold',F(-1,10),64),
                                  ('insufficient_budget',F(-3,10),1)]:
        for adaptive in [False,True]:
            results=[v.family([0],[0,1],threshold=threshold,max_cells=budget,adaptive=adaptive) for _ in range(3)]
            r=results[-1];assert v.verify(r['certificate'])
            method='adaptive' if adaptive else 'uniform';v.v3.base.write_json(root/(name+'_'+method+'.json'),r)
            row={'case':name,'method':method,'max_cells':budget,'modes':10,'threshold':str(threshold),
                 'status':r['certificate']['status'],'point_evaluations':r['point_evaluations'],
                 'median_seconds':median(x['seconds'] for x in results),'samples':[x['seconds'] for x in results]}
            rows.append(row);print(name,method,row['status'],row['point_evaluations'],round(1000*row['median_seconds'],3),flush=True)
    v.v3.base.write_json(root/'benchmark.json',{'version':9,'repetitions':3,'rows':rows})

if __name__=='__main__':main()
