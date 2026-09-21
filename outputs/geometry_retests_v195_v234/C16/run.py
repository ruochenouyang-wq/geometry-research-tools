"""Fixed original C16 retest, exact proof replay, and a small same-input timing run."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
from statistics import median
from time import perf_counter
import io
import json
import sys
import unittest
from unittest.mock import patch
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import compact_gn as g
import support


def save(name,obj):
    (HERE/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2,sort_keys=True)+'\n')


def size(obj):return len(json.dumps(obj,sort_keys=True,separators=(',',':')).encode())


def run():
    baseline={'original_case':support.original_case('C16'),
        'previous_source':'outputs/geometry_cycles_v95_v194/gn_variation.py',
        'previous_source_sha256':support.source_hash(support.PREVIOUS/'gn_variation.py'),
        'input_scope_unchanged':'u_n=(1+t)^n-2^n/(n+1), n=4,16,64, dmu=dt/2',
        'probes':[]}
    originals={}
    for n in (4,16,64):
        p=g.original_cap(n);originals[n]=p;start=perf_counter()
        try:
            result=support.old_gn.moments(p)
            row={'n':n,'outcome':'certificate_returned','certificate':result,
                 'verified':support.old_gn.verify(result,p),'seconds':perf_counter()-start}
        except Exception as e:
            row={'n':n,'outcome':'rejected_as_outside_previous_scope',
                 'exception':type(e).__name__,'message':str(e),'seconds':perf_counter()-start}
        row['original_t_coefficients']=list(map(str,p));baseline['probes'].append(row)
    baseline['diagnosis']='The previous degree36 t-polynomial adapter rejects the original n64 input. The new moment core changes representation and kernels; no previous adapter is modified.'
    save('baseline.json',baseline)

    cases=[]
    for n,p in originals.items():
        c=g.bridge_original(p,n)
        accepted=g.verify(c,expected_n=n,expected_original_coefficients=p)
        assert accepted
        cases.append({'n':n,'certificate':c,'verified':accepted})
    save('original_cases.json',cases)
    terms=g.centered_function([[64,1]])['terms']
    kernels={'mean':g.mean_certificate([[0,2],[4,-1],[64,3]]),
        'centered':g.centered_function([[4,1],[16,F(-2,3)],[64,F(1,5)]],F(-7,3)),
        'mass_form':g.mass_form([[4,2],[16,-3],[64,F(1,2)]]),
        'energy_form':g.energy_form([[4,2],[16,-3],[64,F(1,2)]]),
        'quartic':g.quartic_moment(terms),'ratio_original64':g.ratio(g.centered_function([[64,1]],2**64)['terms'])}
    assert all(g.verify(c) for c in kernels.values());save('kernel_evidence.json',kernels)
    residuals={str(n):{'residual':g.residual(g.centered_function([[n,1]])['terms']),
                       'direction':g.direction(g.centered_function([[n,1]])['terms'])} for n in (4,16,64,1024)}
    assert all(g.verify(c) for row in residuals.values() for c in row.values());save('residuals.json',residuals)
    planes=[]
    for a,b,original in [(4,16,False),(16,64,False),(4,64,True)]:
        basis=[g.centered_function([[n,1]],2**n if original else 1)['terms'] for n in (a,b)]
        c=g.plane(basis,F(1,10**8));assert g.verify_plane(c,basis)
        planes.append({'n_pair':[a,b],'original_scale':original,'certificate':c,'verified':True})
    save('planes.json',planes)
    family={'law':g.family_law(),'values':[g.family_value(n) for n in (1,4,16,64,1024,4096,10**6,10**12)]}
    assert g.verify(family['law']) and all(g.verify(c) for c in family['values']);save('family.json',family)

    benchmark={'measurement':'Five bounded repeats on preconstructed identical original functions; moment call time includes each implementation own certificate construction, excludes input conversion and verification.',
        'timing_limit':'Small local timing samples only; no universal speed or model-token claim. n64 has no successful old path, hence no comparable speed ratio.',
        'operation_count_limit':'Old convolution coefficient products are instrumented. New convolution products and energy-kernel evaluations are reported separately; these are not total machine instructions.',
        'rows':[]}
    for n,p in originals.items():
        terms=g.centered_function([[n,1]],2**n)['terms'];new_times=[];old_times=[]
        new=g.moments(terms)
        for _ in range(5):
            start=perf_counter();g.moments(terms);new_times.append(perf_counter()-start)
        row={'n':n,'new_seconds':new_times,'new_median_seconds':median(new_times),
            'new_operation_counts':new['operation_counts'],
            'dense_original_input_bytes':size(list(map(str,p))),'compact_original_input_bytes':size(terms),
            'new_serialized_certificate_bytes':size(new)}
        if n<=36:
            old=support.old_gn.moments(p)
            for _ in range(5):
                start=perf_counter();support.old_gn.moments(p);old_times.append(perf_counter()-start)
            counter={'convolution_calls':0,'convolution_coefficient_products':0}
            original_mul=support.old_gn.mul
            def counted(a,b):
                counter['convolution_calls']+=1;counter['convolution_coefficient_products']+=len(a)*len(b)
                return original_mul(a,b)
            with patch.object(support.old_gn,'mul',counted): measured=support.old_gn.moments(p)
            assert all(F(old[k])==F(new[k]) for k in ('mean','mass','energy','quartic'))
            assert support.old_gn.verify(measured,p)
            row.update(old_seconds=old_times,old_median_seconds=median(old_times),
                old_operation_counts=counter,old_serialized_certificate_bytes=size(old),exact_moments_agree=True)
        else:
            row.update(old_outcome='degree36_input_rejection',speed_comparison_available=False)
        benchmark['rows'].append(row)
    save('benchmark.json',benchmark)

    # Retain both the old scope failure and the actual first implementation
    # search failure. Neither was converted into a success proof.
    failures=[{'case':'original_n64_previous_adapter','status':'rejected','message':baseline['probes'][-1]['message']},
        {'case':'first_compact_plane_original_u4_u64','tolerance':'1/1000000',
         'status':'eight_unconditioned_ceiling_attempts_failed','exception':'ArithmeticError',
         'message':'No verified plane ceiling within the search budget',
         'resolution':'Record positive independent basis scales and replay the transformed charts and original basis witness coefficients. Original input functions remain bound.'},
        {'case':'independent_audit_two_33_term_plane_bases',
         'f_weights':'powers 1..32: coefficient1 at power1, otherwise1/1024; centered',
         'g_weights':'powers 33..64: coefficient1 at power33, otherwise1/1024; centered',
         'status':'previous_compact_draft_returned_an_unreplayable_65_term_witness',
         'resolution':'Enforce the complete union support budget at the plane entrance; reject with ValueError before returning evidence. Verify every generated plane certificate before return. Original C16 two-term functions remain supported.'},
        {'case':'continuous_real_n_monotonicity','status':'false_hypothesis_rejected',
         'counterexample':'J(11/10)<J(1); derivative at n=1 is -3/100',
         'retained_claim':'Only integer n>=1 forward differences are proved positive.'}]
    save('failure_history.json',failures)
    import test_compact_gn
    stream=io.StringIO();tests=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(test_compact_gn))
    save('tests.json',{'tests_run':tests.testsRun,'passed':tests.wasSuccessful(),'failures':len(tests.failures),'errors':len(tests.errors),'output':stream.getvalue()})
    assert tests.wasSuccessful()
    stages=[dict(row,outcome='implemented_and_evidence_replayed') for row in g.STAGES]
    save('STAGES.json',stages)
    summary={'original_n_values':[4,16,64],'all_original_inputs_supported':True,
        'original_polynomial_coefficients_verified':[len(originals[n]) for n in (4,16,64)],
        'original_n64_probability_ratio':cases[-1]['certificate']['evaluation']['ratio']['probability_ratio'],
        'original_n64_pi_times_area_K':cases[-1]['certificate']['evaluation']['ratio']['pi_times_area_K'],
        'tests_passed':tests.testsRun,'three_planes_replayed':True,
        'large_n_closed_values_replayed':[1024,4096,10**6,10**12],
        'full_function_space_GN_sharpness_claimed':False,
        'all_evidence_replayed':True,'implementation_sha256':support.source_hash(HERE.parent/'compact_gn.py')}
    save('RESULTS.json',summary);print(json.dumps(summary,indent=2));return summary


if __name__=='__main__':run()
