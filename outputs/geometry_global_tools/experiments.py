"""Reproducible global-vs-local and global-bound comparisons."""
from fractions import Fraction as F
from pathlib import Path
from statistics import median
from time import perf_counter
import global_search as g
import coercive
import moment_v12


def check(c):
    if c.get('method')==g.METHOD:return g.verify(c)
    if c.get('method')==coercive.METHOD:return coercive.verify(c)
    if c.get('method')==moment_v12.METHOD:return moment_v12.verify(c)
    return False


def local_search(p,start,iterations=30,modes=10):
    """Deliberately local coordinate polling, not a global optimum certificate."""
    t=perf_counter();oracle=g.Oracle(p,modes,F(1,10**8));x=list(map(F,start));box=[g.rationals(pair) for pair in p['box']]
    steps=[(hi-lo)/16 for lo,hi in box];index=oracle.point(x)
    for _ in range(iterations):
        best=index
        for i in range(len(x)):
            for sign in [-1,1]:
                y=x[:];y[i]+=sign*steps[i]
                if g.inside(y,box):
                    idx=oracle.point(y)
                    if oracle.values[idx][1]<oracle.values[best][1]:best=idx
        if best==index:steps=[step/2 for step in steps]
        else:index=best;x=g.rationals(oracle.records[index]['at'])
    return {'candidate_parameters':[str(z) for z in x],'candidate_upper':str(oracle.values[index][1]),
            'point_evaluations':len(oracle.records),'seconds':perf_counter()-t,
            'global_optimality_proved':False,'candidate_record':oracle.records[index]}


def main():
    root=Path(__file__).resolve().parent/'results';root.mkdir(exist_ok=True);rows=[]
    full=g.problem([0],[[0,1]],[[-6,6]],linear=['1/50'],hessian=[['1/5']])
    two=g.problem([0],[[0,1],[0,0,1]],[[-4,4],[-1,1]],linear=['1/50','-1/6'],hessian=[['1/5','1/50'],['1/50',1]])
    cases=[('tilted_double_well',full,F(1,1000),1024),('two_parameters',two,F(1,100),256),
           ('strict_budget',full,F(1,1000),1)]
    for name,p,tolerance,budget in cases:
        for policy in ['lipschitz','concavity']:
            samples=[]
            for _ in range(3):
                r=g.search(p,tolerance,modes=10,max_leaves=budget,policy=policy);samples.append(r['seconds'])
            c=r['certificate'];assert check(c);g.v3.base.write_json(root/(name+'_'+policy+'.json'),r)
            row={'case':name,'policy':policy,'parameter_domain':'declared_box','modes':10,'max_leaves':budget,
                 'tolerance':str(tolerance),'status':c['status'],'evaluations':r['point_evaluations'],
                 'global_lower':c['global_lower'],'candidate_upper':c['candidate_upper'],'global_gap':c['global_gap'],
                 'candidate_parameters':c['candidate_parameters'],'samples_seconds':samples,'median_seconds':median(samples)}
            rows.append(row);print(name,policy,c['status'],r['point_evaluations'],float(F(c['global_gap'])),round(median(samples),4),flush=True)
    # Same global objective, now on ALL real parameters, with two domain proofs.
    for method in ['growth_at_infinity','moment_elimination']:
        samples=[]
        for _ in range(3):
            fn=coercive.search if method=='growth_at_infinity' else moment_v12.search
            r=fn([0],[[0,1]],['1/50'],[['1/5']],modes=10,tolerance=F(1,1000));samples.append(r['seconds'])
        c=r['certificate'];assert check(c);g.v3.base.write_json(root/('all_real_'+method+'.json'),r)
        row={'case':'all_real_double_well','policy':method,'parameter_domain':'all_real','modes':10,'max_leaves':256,
             'tolerance':'1/1000','status':c['status'],'evaluations':r['point_evaluations'],
             'global_lower':c['global_lower'],'candidate_upper':c['candidate_upper'],'global_gap':c['global_gap'],
             'candidate_parameters':c['candidate_parameters'],'search_box':c['boxed_certificate']['problem']['box'],
             'samples_seconds':samples,'median_seconds':median(samples)}
        rows.append(row);print('all_real',method,c['status'],r['point_evaluations'],float(F(c['global_gap'])),round(median(samples),4),flush=True)
    local=local_search(full,[3]);g.v3.base.write_json(root/'local_search.json',{'problem':full,**local})
    samples=[]
    for _ in range(3):
        cost=two['penalty']
        r=moment_v12.search(two['q0'],two['directions'],cost['linear'],cost['hessian'],modes=10,tolerance=F(1,100))
        samples.append(r['seconds'])
    c=r['certificate'];assert check(c);g.v3.base.write_json(root/'all_real_two_parameters.json',r)
    rows.append({'case':'all_real_two_parameters','policy':'moment_elimination','parameter_domain':'all_real',
                 'modes':10,'max_leaves':256,'tolerance':'1/100','status':c['status'],
                 'evaluations':r['point_evaluations'],'global_lower':c['global_lower'],
                 'candidate_upper':c['candidate_upper'],'global_gap':c['global_gap'],
                 'candidate_parameters':c['candidate_parameters'],'search_box':c['boxed_certificate']['problem']['box'],
                 'samples_seconds':samples,'median_seconds':median(samples)})
    print('all_real_two_parameters',c['status'],r['point_evaluations'],float(F(c['global_gap'])),round(median(samples),4),flush=True)
    positive=g.problem([0],[[0,1]],[[0,6]],linear=['1/50'],hessian=[['1/5']])
    r=g.search(positive,F(1,1000),modes=10);assert check(r['certificate'])
    g.v3.base.write_json(root/'positive_basin_global.json',r)
    full_cert=__import__('json').loads((root/'tilted_double_well_concavity.json').read_text())['certificate']
    positive_cert=r['certificate'];endpoints=[]
    for a in [F(0),F(6)]:
        record=next(x for x in positive_cert['points'] if x['at']==[str(a)])
        lower=F(record['spectral_certificate']['lower'])+g.penalty(positive,[a])
        assert lower>F(positive_cert['candidate_upper']);endpoints.append(str(lower))
    difference=F(positive_cert['global_lower'])-F(full_cert['candidate_upper']);assert difference>0
    trap={'local_candidate':local['candidate_parameters'],'local_upper':local['candidate_upper'],
          'positive_basin_lower':positive_cert['global_lower'],'positive_basin_upper':positive_cert['candidate_upper'],
          'positive_basin_boundary_lowers':endpoints,'full_domain_upper':full_cert['candidate_upper'],
          'certified_improvement_over_every_positive_basin_point':str(difference),
          'interpretation':'Positive-side restricted minimum is interior and hence a true local minimum, but the full-domain candidate is strictly better.'}
    g.v3.base.write_json(root/'local_trap_evidence.json',trap)
    g.v3.base.write_json(root/'benchmark.json',{'repetitions':3,'rows':rows,'local_search':{k:v for k,v in local.items() if k!='candidate_record'},
        'timing':'Total search plus final dense certificate reconstruction; excludes interpreter startup and file writes.',
        'case_selection':'Constructed diagnostic cases, including a rigorously witnessed inferior local minimum; not a representative geometry benchmark.'})


if __name__=='__main__':main()
