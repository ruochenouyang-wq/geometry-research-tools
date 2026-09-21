"""Read-only R03 audit on explicitly public strong-potential development cases.

No implementation or result files are written. This is not a holdout run.
"""
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys
from unittest.mock import patch

ROOT=Path(__file__).resolve().parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT))
import power_forms as pf
import spectral_gap as sg
import singular_solver as ss
import fullspace_singular as fs
from implementation_review import independent_statistics


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def code_hashes():
    return {name:sha(ROOT/name) for name in ('power_forms.py','spectral_gap.py','singular_solver.py','fullspace_singular.py')}


def plus(a,b):
    out=dict(a)
    for p,c in b.items():out[p]=out.get(p,F(0))+c
    return {p:c for p,c in out.items() if c}


def times(a,b):
    out={}
    for p,c in a.items():
        for q,d in b.items():out[p+q]=out.get(p+q,F(0))+c*d
    return {p:c for p,c in out.items() if c}


def scaled(a,k):
    return {p:c*k for p,c in a.items() if c*k}


def diff(a):
    return {p-1:p*c for p,c in a.items() if p*c}


def weighted_integral(a):
    assert all(p>-1 for p in a)
    return sum((2*c/(p+1)-2*c/(p+3) for p,c in a.items()),F(0))


def independent_m1_residual(cert):
    q=cert['function'];f=dict(zip(map(F,cert['powers']),map(F,cert['coefficients'])))
    df=diff(f);flux=times({F(0):F(1),F(2):F(-1)},df)
    potential=plus({F(q['exponent']):F(q['amplitude'])},{F(0):F(q['offset'])})
    # H(sqrt(1-z²) cos(phi) f) divided by that angular/polar factor.
    hf=plus(plus(scaled(diff(flux),-1),times({F(1):F(2)},df)),
            plus(scaled(f,2),times(potential,f)))
    mass=weighted_integral(times(f,f));energy=weighted_integral(times(f,hf))
    norm=weighted_integral(times(hf,hf));mu=energy/mass
    residual=plus(hf,scaled(f,-mu));variance=weighted_integral(times(residual,residual))/mass
    expected={'mass':mass,'operator_form':energy,'operator_norm_squared':norm,
              'rayleigh':mu,'residual_squared':variance}
    assert variance==norm/mass-mu*mu
    for key,value in expected.items():assert value==F(cert['statistics'][key]),key
    beta=F(cert['gap']['beta']);assert mu<beta
    assert F(cert['temple_lower'])==mu-variance/(beta-mu)
    lower=min(F(cert['m0_source']['lower']),max(F(cert['temple_lower']),F(cert['m1_source']['lower'])))
    upper=min(mu,F(cert['m0_source']['upper']),F(cert['m1_source']['upper']))
    assert lower==F(cert['lower']) and upper==F(cert['upper'])
    return {'independent_merged_m1_residual_exact':True,'merged_operator_terms':len(hf),
            'width_approx':float(upper-lower),'status':cert['status']}


def independent_form(cert):
    q=cert['function'];raw=cert.get('selected_form',cert);proof=raw['proof']
    p,A,b=-F(q['exponent']),-F(q['amplitude']),F(q['offset'])
    m,s=p.numerator,p.denominator;base=F(proof['base'])
    delta,cap=base**(-s),A*base**m
    # Substitute t=delta*y^(2s) in the three-term core majorant.
    D=F(2*s,2*s-3*m)-F(3*s,2*s-m)+F(s,2*s+m)
    simplified=p*(2+5*p)/(4*(1-F(3,2)*p)*(1-p*p/4))
    assert D==simplified and D>0
    norm_cube=cap**3*delta**2*D**2
    eta=F(proof['eta_upper']);unit=F(1,2**raw['root_bits'])
    assert F(proof['delta'])==delta and F(proof['cap'])==cap
    assert F(proof['integral_coefficient_D'])==D
    assert F(proof['residual_norm_cube_upper'])==norm_cube
    assert eta**3>=norm_cube and (eta-unit)**3<norm_cube
    assert 0<eta<F(1,2)
    assert F(raw['form']['energy_factor'])==1-2*eta
    assert F(raw['form']['mass_offset'])==b-cap-eta
    if 'selection' in cert:
        scores=[]
        for attempt in cert['selection']['attempts']:
            e,k=F(attempt['energy_factor']),F(attempt['mass_offset'])
            n=cert['selection']['tail_start']
            score=e*n*(n+1)+k
            assert score==F(attempt['initial_tail_lower'])
            if attempt['admissible']:scores.append((score,-F(attempt['base'])))
        best=max(scores)
        assert -best[1]==base
        assert cert['form']==raw['form']
        assert cert['selection']['optimal_over_all_bases_claimed'] is False
    return {'base':str(base),'D':str(D),'eta_cube_enclosure_exact':True,
            'energy_factor':str(raw['form']['energy_factor']),
            'selection_is_only_between_declared_candidates':True}


def main():
    before=code_hashes();results=[]
    # Exact dyadic cube boundaries, including tiny and huge input scales.
    total=0
    for bits in (8,48,160):
        for q in (F(0),F(1),F(8),F(1,7),F(1,2**900),F(2**900)):
            u=pf.cube_root_upper(q,bits);step=F(1,2**bits)
            assert u**3>=q and (not u or (u-step)**3<q)
            total+=1
    for k in (1,7,255,256,300):
        c=F(k,256)**3;eps=F(1,2**40)
        for q,expected in ((c,F(k,256)),(c-eps,F(k,256)),(c+eps,F(k+1,256))):
            assert pf.cube_root_upper(q,8)==expected;total+=1
    for j in range(65):
        x=F(j,64);S=1-x/2
        assert S>=0 and S*S>=1-x
        assert (1-x)*S==1-F(3,2)*x+x*x/2
    results.append({'case':'cube_root_and_majorant','exact_root_cases':total,
                    'rational_majorant_points':65,'all_outward_and_minimal':True})
    public=json.loads((ROOT/'evaluation'/'public_strong_stress.json').read_text())
    assert all(row['classification']=='public_development_stress' for row in public)
    functions=[]
    for row in public:
        q=row['task']['function']
        if q not in functions:functions.append(q)
    forms=[]
    for i,q in enumerate(functions):
        fc=pf.select_form(q,13);forms.append(fc)
        with patch.object(pf,'select_form',side_effect=AssertionError('Verifier called proposal')):
            assert pf.verify_form(fc,q)
        result={'case':'public_strong_form_'+str(i+1),**independent_form(fc)}
        rejected=[]
        edits=[('D',lambda c:c['selected_form']['proof'].update(integral_coefficient_D='0')),
               ('cube',lambda c:c['selected_form']['proof'].update(residual_norm_cube_upper='0')),
               ('eta',lambda c:c['selected_form']['proof'].update(eta_upper='0')),
               ('selected_base',lambda c:c['selection'].update(selected_base='2')),
               ('optimality',lambda c:c['selection'].update(optimal_over_all_bases_claimed=True)),
               ('selected_energy',lambda c:c['form'].update(energy_factor='1')),
               ('axis',lambda c:c['function'].update(axis=0)),
               ('offset',lambda c:c['function'].update(offset='1'))]
        for name,edit in edits:
            bad=deepcopy(fc);edit(bad);assert not pf.verify_form(bad,q),name;rejected.append(name)
        result['rejected_mutations']=rejected;results.append(result)
    # Check custom ground lower/upper inertia and source binding on both strong profiles.
    for i,(q,fc) in enumerate(zip(functions,forms)):
        ground=sg.certified_ground(q,0,True,12,12,16,form_certificate=fc,form_verifier=pf.verify_form)
        assert sg.verify_ground(ground,q,0,True,pf.verify_form)
        assert ground['lower_inertia'][0]==0
        assert ground['upper_inertia'][0]+ground['upper_inertia'][1]>=1
        assert F(ground['lower'])<F(ground['radial_tail_lower'])
        assert not sg.verify_ground(ground,q,0,False,pf.verify_form)
        assert not sg.verify_ground(ground,q,1,True,pf.verify_form)
        assert not sg.verify_ground(ground,q,0,True)
        for path,value in [('lower','100'),('mean_zero',False),('eigenvalue_index',2)]:
            bad=deepcopy(ground);bad[path]=value;assert not sg.verify_ground(bad,q,0,True,pf.verify_form)
        bad=deepcopy(ground);bad['form_certificate']['selected_form']['proof']['eta_upper']='0'
        assert not sg.verify_ground(bad,q,0,True,pf.verify_form)
        results.append({'case':'public_strong_custom_ground_'+str(i+1),
                        'full_radial_lower_and_upper_replayed':True,
                        'space_sector_form_and_endpoint_tampering_rejected':True,
                        'width_approx':float(F(ground['upper'])-F(ground['lower']))})
    # Strong private algebra must not relax the old module's normalization policy.
    frozen=ss.trial
    old_globals=dict(frozen.__dict__)
    old_normalize=frozen.normalize;old_form=sg.direct.form_bound
    old_cache=frozen._matrices.cache_info()
    old_file_hash=sha(frozen.__file__)
    for q in functions:
        try:frozen.normalize(q)
        except ValueError:pass
        else:raise AssertionError('Expected old eta limitation did not apply')
        private=ss._engine(q)
        assert private is not frozen and private.__dict__ is not frozen.__dict__
        assert private._matrices is not frozen._matrices
        assert private._matrices.__wrapped__.__globals__ is private.__dict__
        assert private.normalize is ss.normalize
        assert private.generated_powers(q,6)
    for row in public:
        if row['task']['mean_zero'] is not True:continue
        run=ss.solve(row['task']);cert=run['certificate']
        if cert is None:
            results.append({'case':row['id'],'status':run['status'],'certificate':None,
                            'meaning':'An open method outcome, not a successful mathematical certificate'})
            continue
        assert ss.verify(cert,row['task']['function'],True,row['task']['tolerance'])
        result={'case':row['id'],**independent_m1_residual(cert),
                'attempt_count':len(run['attempts'])}
        bad=deepcopy(cert);bad['m0_source']['mean_zero']=False;assert not ss.verify(bad)
        bad=deepcopy(cert);bad['gap']['mean_zero']=False;assert not ss.verify(bad)
        bad=deepcopy(cert);bad['statistics']['residual_squared']='0';assert not ss.verify(bad)
        result['wrong_source_space_gap_space_and_residual_rejected']=True
        results.append(result)
    assert frozen.normalize is old_normalize and sg.direct.form_bound is old_form
    assert all(frozen.__dict__.get(key) is value for key,value in old_globals.items())
    assert set(frozen.__dict__)==set(old_globals)
    assert frozen._matrices.cache_info()==old_cache
    assert sha(frozen.__file__)==old_file_hash
    assert sys.modules['adaptive_singular'] is frozen
    for q in functions:
        try:frozen.normalize(q)
        except ValueError:pass
        else:raise AssertionError('Private normalization leaked into the old module')
    results.append({'case':'private_trial_isolation','old_normalize_and_direct_form_identity_preserved':True,
                    'all_old_trial_global_bindings_preserved':True,'old_trial_matrix_cache_unchanged':True,
                    'old_source_sha256_unchanged':old_file_hash,'sys_modules_keeps_original_trial':True})
    for row in public:
        if row['task']['mean_zero'] is not False:continue
        task=row['task']
        run=fs.solve(task['function'],task['tolerance'],False,max_terms=10,wall_seconds=30)
        cert=run['certificate'];assert cert is not None and run['target_met']
        with patch.object(fs,'proposal',side_effect=AssertionError('Verifier called proposal')):
            assert fs.verify(cert,task['function'],False,task['tolerance'])
        assert cert['gap']['form_certificate']['format']==pf.SELECTED_FORMAT
        assert not fs.verify(cert,expected_mean_zero=True)
        bad=deepcopy(cert);bad['gap']['mean_zero']=True;assert not fs.verify(bad)
        bad=deepcopy(cert);bad['gap']['form_certificate']['function']['axis']=0;assert not fs.verify(bad)
        bad=deepcopy(cert);bad['statistics']['residual_squared']='0';assert not fs.verify(bad)
        results.append({'case':row['id'],**independent_statistics(cert),
                        'fullspace_source_problem_and_strong_residual_binding_checked':True,
                        'verification_without_proposal':True,
                        'terms':[attempt['terms'] for attempt in run['attempts']],
                        'gap_modes':[attempt['modes'] for attempt in run['gap_attempts']]})
    after=code_hashes();assert before==after,'Inspected code changed during review; rerun against final hashes'
    print(json.dumps({'all_independent_checks_passed':True,'source_sha256':after,
          'source_unchanged_during_review':True,'checks':results,
          'limitations':['Only public strong development parameters were used; no seed or holdout file was read.',
                         'Whole frozen Kernel/inertia implementation remains a shared mathematical dependency.',
                         'The private module reuses pure/shared mathematical helpers, but does not replace old globals.',
                         'No continuous or all-base optimality is certified.',
                         'Full-space reduction only applies with mean_zero=False.']},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
