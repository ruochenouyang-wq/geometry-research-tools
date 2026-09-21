from pathlib import Path
from statistics import median
import v08_cluster as v

def main():
    root=Path(__file__).resolve().parent/'results'/'v08';root.mkdir(parents=True,exist_ok=True)
    rows=[]
    for name,q in [('separated',[0,1]),('near_doublet',[0,'1/1000000',-100])]:
        for n in [12,16,20]:
            for m in [1,2]:
                try:
                    results=[v.cluster(q,m,n) for _ in range(3)];r=results[-1];c=r['certificate']
                    assert v.verify(c);v.v3.base.write_json(root/(name+'_n'+str(n)+'_m'+str(m)+'.json'),r)
                    row={'case':name,'q':q,'modes':n,'cluster_size':m,'status':'certified',
                         'median_seconds':median(x['seconds'] for x in results),
                         'samples':[x['seconds'] for x in results],
                         'projector_squared_bound':c['projector_hilbert_schmidt_squared_upper'],
                         'external_gap':str(v.F(c['external_lower'])-v.F(c['ritz_upper']))}
                except ValueError as exc:
                    row={'case':name,'q':q,'modes':n,'cluster_size':m,'status':'no_certified_gap','reason':str(exc)}
                rows.append(row);print(name,n,m,row['status'],float(v.F(row.get('projector_squared_bound','0'))),flush=True)
    v.v3.base.write_json(root/'benchmark.json',{'version':8,'repetitions':3,'rows':rows,
        'comparison':'Same potential and mode budgets. m=1 and m=2 concern DIFFERENT projectors; no equal-task speedup claimed.'})

if __name__=='__main__':main()
