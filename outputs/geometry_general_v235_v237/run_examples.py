"""Bounded, unchanged-function comparisons for the three increments."""
import json
import time
from dependencies import ROOT, F, save, sha
import function_models as f
import spectral_transfer as s

ANALYTIC={'kind':'analytic_sum','terms':[
    {'function':'exp','argument':{'0,0,1':'1/2'}},
    {'function':'sin','argument':{'0,0,1':'1/3'}},
    {'function':'log1p','argument':{'0,0,1':'1/4'}}]}
MULTI={'kind':'analytic_sum','polynomial':{'2,0,0':'1','0,2,0':'2','0,0,2':'3'},
       'terms':[{'function':'exp','argument':{'1,0,0':'1','0,1,0':'1','0,0,1':'1'},'coefficient':'1/1000000'}]}
STEP={'kind':'axis_profile','axis':2,'profile':'step'}
SINGULAR={'kind':'axis_profile','axis':2,'profile':'abs_power','exponent':'-1/4','amplitude':'-1'}


def run():
    rows=[];proofs=[];start=time.perf_counter()
    def record(name,cert,elapsed):
        assert s.verify(cert)
        path=f'certificates/{name}.json';save(path,cert)
        proofs.append({'path':path,'sha256':sha(ROOT/path),'verified':True,
                       'certificate_format':cert['format']})
        row={'case':name,'function':cert['function'],'method':cert['format'],
             'lower':cert['lower'],'upper':cert['upper'],'width':cert['exact_width'],
             'width_decimal':float(F(cert['exact_width'])),'status':cert['status'],
             'error_budget':cert['error_budget'],'seconds_before_final_saved_replay':elapsed}
        if cert['coercivity']:
            row.update(fixed_reference=cert['coercivity']['reference'],
                       shift=cert['coercivity']['shift'],theta=cert['coercivity']['distance_upper'])
        rows.append(row)
        print(name,cert['status'],row['width_decimal'],flush=True)
    for order in [4,8,12]:
        tick=time.perf_counter()
        c=s.solve(ANALYTIC,order=order,source_options={'modes':8,'max_modes':12,'bits':40})
        record(f'V235_analytic_order{order}',c,time.perf_counter()-tick)
    assert rows[-1]['status']=='target_met'
    tick=time.perf_counter()
    c=s.solve(MULTI,order=2,source_options={'L':3,'bits':32})
    record('V235_multiaxis_original_function',c,time.perf_counter()-tick)
    assert c['source']['backend']=='cartesian_constraints'
    for degree in [0,4,8,12]:
        tick=time.perf_counter()
        c=s.solve(STEP,method='l2',degree=degree,
                  source_options={'modes':6,'max_modes':8,'bits':32,'tolerance':'1/1000000'})
        record(f'V236_step_degree{degree}',c,time.perf_counter()-tick)
    reference=f.l2_model(SINGULAR,degree=4)['polynomial']
    reference_ranges=[]
    for degree in [4,8,12]:
        tick=time.perf_counter()
        model=f.l2_model(SINGULAR,degree=degree)
        source=s.polynomial_source(model['polynomial'],modes=4,max_modes=4,max_m=2,
                                   bits=12,tolerance='1/1000')
        direct=s.transfer_l2(model,source)
        record(f'V236_negative_singularity_degree{degree}',direct,time.perf_counter()-tick)
        tick=time.perf_counter()
        fixed=s.transfer_fixed(model,source,reference)
        record(f'V237_fixed_p4_degree{degree}',fixed,time.perf_counter()-tick)
        reference_ranges.append(fixed['coercivity']['range_evidence'])
    assert all(v==reference_ranges[0] for v in reference_ranges)
    # Norm-only family evidence demonstrates a uniformly stable fixed reference
    # through degree24 even when the frozen source solver rejects that degree.
    norms=[]
    for degree in [0,4,8,12,16,20,24]:
        model=f.l2_model(SINGULAR,degree=degree)
        fixed=s.fixed_reference(model['polynomial'],{'0,0,0':'-4/3'})
        assert F(fixed['distance_squared'])+F(model['error_squared'])==F(2,9)
        assert F(fixed['distance_upper']) < F(1,2)
        norms.append({'degree':degree,'model':model,'coercivity':fixed,
                      'identity':'distance_squared+projection_error_squared=2/9',
                      'spectrum_computed_in_this_record':False})
    save('FIXED_REFERENCE_NORMS.json',norms)
    failures=[]
    def rejected(label,action):
        try:action()
        except (ValueError,TypeError,ArithmeticError) as e:
            failures.append({'case':label,'status':'not_certified','exception':type(e).__name__,'reason':str(e)})
        else:raise AssertionError(label+' unexpectedly accepted')
    rejected('degree24_source_coefficient_budget',lambda:s.polynomial_source(f.l2_model(SINGULAR,degree=24)['polynomial']))
    rejected('non_L2_boundary_alpha_minus_half',lambda:f.l2_model({**SINGULAR,'exponent':'-1/2'}))
    rejected('log_domain_not_certified',lambda:f.analytic_model({'kind':'analytic_sum','terms':[{'function':'log1p','argument':{'0,0,1':'1'}}]}))
    big=f.l2_model({**SINGULAR,'amplitude':'-4'},degree=0)
    big_source=s.polynomial_source(big['polynomial'],modes=2,max_modes=2,bits=12)
    rejected('relative_L2_error_at_least_one',lambda:s.transfer_l2(big,big_source))
    step0=f.l2_model(STEP,degree=0)
    step_source=s.polynomial_source(step0['polynomial'],modes=2,max_modes=2,bits=12)
    rejected('fixed_relative_error_equal_one',lambda:s.transfer_fixed(step0,step_source,{}))
    save('FAILURES.json',failures)
    for proof in proofs:
        assert sha(ROOT/proof['path'])==proof['sha256']
        assert s.verify(json.loads((ROOT/proof['path']).read_text()))
    result={'versions':[235,236,237],'results':rows,'certificates':proofs,
            'saved_certificates_verified':len(proofs),'preserved_failures':len(failures),
            'all_original_function_bindings_verified':True,
            'wall_seconds':time.perf_counter()-start,
            'timing_scope':'local Python observations; not a model or universal speed benchmark',
            'interpretation':'Rough potential intervals remain open; function errors are included, not sampled.'}
    save('RESULTS.json',result)
    return result


if __name__=='__main__':run()
