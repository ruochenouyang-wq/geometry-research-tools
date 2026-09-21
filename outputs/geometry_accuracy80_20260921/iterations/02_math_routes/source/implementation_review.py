"""Independent, bounded mathematical checks on public development inputs.

This script changes no implementation or evidence files and reads no holdout.
The strongest check integrates the merged full differential residual directly,
rather than contracting the implementation's strong-operator Gram matrix.
"""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
from unittest.mock import patch

ROOT = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(ROOT))
import fullspace_singular as fs
import spectral_gap as sg


def add(left, right):
    out = dict(left)
    for p, c in right.items():
        out[p] = out.get(p, F(0))+c
    return {p:c for p,c in out.items() if c}


def multiply(left, right):
    out = {}
    for p,a in left.items():
        for q,b in right.items():
            out[p+q] = out.get(p+q, F(0))+a*b
    return {p:c for p,c in out.items() if c}


def derivative(poly):
    return {p-1:p*c for p,c in poly.items() if p*c}


def scale(poly, factor):
    return {p:c*factor for p,c in poly.items() if c*factor}


def integral_even(poly):
    if any(p <= -1 for p in poly):
        raise AssertionError('Independent merged integral is not integrable')
    return sum((2*c/(p+1) for p,c in poly.items()), F(0))


def independent_statistics(cert):
    q = cert['function']
    u = dict(zip(map(F, cert['trial']['powers']), map(F, cert['trial']['coefficients'])))
    potential = add({F(q['exponent']):F(q['amplitude'])}, {F(0):F(q['offset'])})
    du = derivative(u)
    flux = multiply({F(0):F(1), F(2):F(-1)}, du)
    hu = add(scale(derivative(flux), F(-1)), multiply(potential, u))
    mass = integral_even(multiply(u,u))
    energy = integral_even(multiply(u,hu))
    weak = integral_even(add(multiply({F(0):F(1),F(2):F(-1)},multiply(du,du)),
                             multiply(potential,multiply(u,u))))
    assert energy == weak
    mu = energy/mass
    residual = add(hu,scale(u,-mu))
    variance = integral_even(multiply(residual,residual))/mass
    hnorm = integral_even(multiply(hu,hu))
    assert variance == hnorm/mass-mu*mu
    actual = {'mass':mass,'operator_form':energy,'operator_norm_squared':hnorm,
              'rayleigh':mu,'residual_squared':variance}
    for key,value in actual.items():
        assert value == F(cert['statistics'][key]), key
    beta = F(cert['gap']['beta'])
    assert mu < beta
    lower = mu-variance/(beta-mu)
    assert lower == F(cert['lower']) and mu == F(cert['upper'])
    assert mu-lower == F(cert['exact_width'])
    return {'merged_trial_terms':len(u),'merged_Hu_terms':len(hu),
            'merged_full_residual_terms':len(residual),
            'independent_mass_weak_energy_Hu_norm_and_full_residual_exact':True,
            'temple_endpoints_exact':True,'width_approx':float(mu-lower),
            'mu_approx':float(mu),'beta_approx':float(beta)}


def independent_clip(form):
    q, proof = form['function'], form['proof']
    assert form['method'] == 'clipped_negative_L2'
    alpha, strength = F(q['exponent']), -F(q['amplitude'])
    base = proof['clip_base'];r,s=-alpha.numerator,alpha.denominator
    delta, cap = F(1,base**s), strength*base**r
    assert delta == F(proof['delta']) and cap == F(proof['cap'])
    # Integrate k² z^(2 alpha) - 2 k C z^alpha + C² on [0,delta].
    first = strength**2*F(1,base**(s-2*r))/(1+2*alpha)
    cross = -2*strength*cap*F(1,base**(s-r))/(1+alpha)
    constant = cap**2*delta
    residual = first+cross+constant
    eta=F(proof['eta_upper'])
    assert residual == F(proof['residual_L2_squared'])
    assert residual >= 0 and eta*eta >= residual and 0 <= eta < 1
    assert F(form['form']['energy_factor']) == 1-eta
    assert F(form['form']['mass_offset']) == F(q['offset'])-cap-eta
    return {'delta':str(delta),'cap':str(cap),'independent_residual_squared':str(residual),
            'eta_upper':str(eta),'positive_energy_factor':True}


def hashes():
    return {name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
            for name in ('spectral_gap.py','fullspace_singular.py')}


def main():
    before=hashes();checks=[]
    q={'kind':'axis_profile','axis':2,'profile':'abs_power','exponent':'-3/7',
       'amplitude':'-3/5','offset':'0'}  # Public old H10, already a regression case.
    run=fs.solve(q,tolerance='1/10000000000',mean_zero=False,max_terms=10,wall_seconds=10)
    cert=run['certificate'];assert cert is not None and run['target_met']
    with patch.object(fs,'proposal',side_effect=AssertionError('Verifier called search')):
        assert fs.verify(cert,q,False,'1/10000000000')
    checks.append({'case':'public_H10','exact_width':cert['exact_width'],
                   'stages':[a['terms'] for a in run['attempts']],
                   'verification_without_candidate_search':True,
                   'clip':independent_clip(cert['gap']['form_certificate']),
                   **independent_statistics(cert)})
    rejections=[]
    cases=[]
    bad=deepcopy(cert);bad['gap']['form_certificate']['proof']['cap']='0';cases.append(('clipping_cap',bad))
    bad=deepcopy(cert);bad['gap']['form_certificate']['proof']['residual_L2_squared']='0';cases.append(('clipping_integral',bad))
    bad=deepcopy(cert);bad['gap']['form_certificate']['form']['energy_factor']='1';cases.append(('clipping_energy_factor',bad))
    bad=deepcopy(cert);bad['gap']['comparison_inertia']=[0,0,bad['gap']['modes']];cases.append(('gap_inertia',bad))
    bad=deepcopy(cert);bad['gap']['azimuth_m']=1;cases.append(('gap_wrong_sector',bad))
    bad=deepcopy(cert);bad['gap']['mean_zero']=True;cases.append(('gap_wrong_space',bad))
    bad=deepcopy(cert);bad['mean_zero']=True;cases.append(('full_to_mean_zero',bad))
    bad=deepcopy(cert);bad['statistics']['residual_squared']='0';cases.append(('full_residual',bad))
    bad=deepcopy(cert);bad['full_space_domination']['requires_mean_zero_false']=False;cases.append(('domination_claim',bad))
    for field,value in [('axis',0),('exponent','-1/4'),('offset','1/7'),('amplitude','-1/2')]:
        bad=deepcopy(cert);bad['function'][field]=value;cases.append(('function_'+field,bad))
        assert not fs.verify(cert,expected_function=dict(q,**{field:value}))
    for label,bad in cases:
        assert not fs.verify(bad), label
        rejections.append(label)
    assert not fs.verify(cert,expected_mean_zero=True)
    assert not fs.verify(cert,expected_mean_zero=0)
    for bad_power in ('1','3/2'):
        try:fs.trial_matrices(q,['0',bad_power])
        except ValueError:pass
        else:raise AssertionError('Unsafe strong-domain power accepted')
    checks.append({'case':'H10_adversarial_mutations','rejected':rejections,
                   'unsafe_strong_domain_powers_rejected':['1','3/2']})
    q2={'kind':'axis_profile','axis':0,'profile':'abs_power','exponent':'-1/4',
        'amplitude':'-1/2','offset':'2/7'}
    cert2=fs.solve(q2,tolerance='1/1000000',mean_zero=False,max_terms=6,wall_seconds=10)['certificate']
    assert cert2 and fs.verify(cert2,q2,False,'1/1000000')
    checks.append({'case':'public_axis0_offset_example',**independent_statistics(cert2)})
    zero={'kind':'axis_profile','profile':'step','amplitude':'0','offset':'0'}
    k=sg._kernel(zero,0,False,2,0,40)
    endpoint=sg._certificate(k,F(2))
    assert endpoint['comparison_inertia']==[1,1,0]
    assert sg.verify_gap(endpoint,zero,0,False)
    too_high=deepcopy(endpoint);too_high['lower']=too_high['beta']='5/2'
    assert not sg.verify_gap(too_high,zero,0,False)
    boundary=deepcopy(endpoint);boundary['lower']=boundary['beta']=str(k.beta)
    assert not sg.verify_gap(boundary,zero,0,False)
    checks.append({'case':'analytic_zero_potential_second_value',
                   'index2_endpoint':'2','inertia':[1,1,0],
                   'one_negative_plus_one_zero_accepted':True,
                   'too_high_and_nonstrict_tail_rejected':True})
    after=hashes()
    assert before==after, 'Authors changed inspected implementation during this run; rerun review'
    output={'source_sha256':after,'implementation_unchanged_during_check':True,
            'all_independent_checks_passed':True,'checks':checks,
            'limits':['Public development inputs only; no new holdout used.',
                      'Kernel assembly and its exact inertia routine remain inherited trust dependencies.',
                      'Merged differential residual is independently integrated, without using the implementation Gram contraction.',
                      'No formal proof assistant or broad performance claim.']}
    print(json.dumps(output,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
