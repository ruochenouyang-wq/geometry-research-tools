from pathlib import Path
from statistics import median
import v07_precision as v

def main():
    root=Path(__file__).resolve().parent/'results'/'v07';root.mkdir(parents=True,exist_ok=True)
    rows=[]
    for name,q,n in [('linear_low_space',[0,1],4),('linear_high_precision',[0,1],14),
                     ('quadratic_high_precision',[1,-2,3],20)]:
        for steps in [0,1,2]:
            results=[v.refine(q,n,112,steps) for _ in range(3)]
            r=results[-1];assert v.verify(r['certificate'])
            v.v3.base.write_json(root/(name+'_steps'+str(steps)+'.json'),r)
            row={'case':name,'q':q,'modes':n,'steps':steps,'coefficient_bits':112,
                 'median_seconds':median(x['seconds'] for x in results),
                 'samples':[x['seconds'] for x in results],
                 'exact_width':r['certificate']['exact_width'],
                 'inside_residual_squared':r['certificate']['inside_residual_squared'],
                 'outside_residual_squared':r['certificate']['outside_residual_squared']}
            rows.append(row)
            print(name,steps,round(row['median_seconds']*1000,3),float(v.F(row['exact_width'])),flush=True)
    v.v3.base.write_json(root/'benchmark.json',{'version':7,'repetitions':3,'rows':rows})

if __name__=='__main__':main()
