"""Compare retained V2 against both V10 routes at identical ground-state targets."""
from pathlib import Path
from statistics import median
from time import perf_counter
from fractions import Fraction as F
import previous_run as old
import v10_adaptive as v
from check_document import verify_document


def main():
    root=Path(__file__).resolve().parent/'results'/'final';root.mkdir(parents=True,exist_ok=True);rows=[]
    cases=[('linear',[0,1],F(1,10**10)),('even_well',[0,0,-100],F(1,10**10)),
           ('near_pair',[0,'1/1000000',-100],F(1,10**10)),
           ('quartic',[0,20,-100,0,1],F(1,10**10)),('precision',[0,1],F(1,10**40))]
    for name,q,tolerance in cases:
        for method in ['v2_auto','v10_spectrum','v10_function']:
            samples=[]
            for _ in range(3):
                start=perf_counter()
                if method=='v2_auto':r=old.search(q,tolerance=tolerance,max_modes=24)
                elif method=='v10_spectrum':r=v.adaptive(q,tolerance=tolerance,max_modes=24)
                else:r=v.function(q,tolerance,start_modes=4,max_modes=24)
                assert verify_document(r);samples.append(perf_counter()-start)
            old.e.base.write_json(root/(name+'_'+method+'.json'),r)
            row={'case':name,'q':q,'method':method,'k':1,'max_modes':24,'tolerance':str(tolerance),
                 'status':r['status'],'exact_width':r['certificate']['exact_width'],
                 'median_seconds':median(samples),'samples':samples}
            rows.append(row);print(name,method,row['status'],round(1000*median(samples),3),flush=True)
    old.e.base.write_json(root/'benchmark.json',{'version':'final_comparison','repetitions':3,'rows':rows,
        'comparison':'Identical ground-state eigenvalue width targets and maximum Legendre index budget. V2 auto also tries its original degree-1 exponential structural candidate; charged to cost. Function error is not the target metric here.'})


if __name__=='__main__':main()
