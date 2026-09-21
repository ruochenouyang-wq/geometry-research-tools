"""Reproduce this round without modifying frozen external-trial files."""
from pathlib import Path
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import sys
from time import perf_counter
sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
import gn_variation as g


def save(name,value):
    (HERE/name).write_text(json.dumps(value,ensure_ascii=False,sort_keys=True,indent=2)+'\n')


def run():
    old=HERE.parent.parent/'geometry_external_trials'
    source=old/'interpolation_results.json'
    result=json.loads(source.read_text())
    spec=importlib.util.spec_from_file_location('frozen_interpolation_trials',old/'interpolation_trials.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    baseline={'source':str(source),'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'scope':result['global_scope'],'previous_cases':[],'limits':[
        'Exactly two configured 2D trial spaces, each basis degree <=6.',
        'No complete Euler-Lagrange residual or automatic variation loop in the frozen module.',
        'The old plane upper bounds never purported to bound all sphere H1 functions.']}
    for row in result['cases']:
        cert=json.loads((old/row['certificate']).read_text())
        baseline['previous_cases'].append({'case':row['case'],
            'verified':module.verify(cert,expected_basis=row['basis']),
            'lower':cert['lower'],'upper':cert['upper']})
    try:
        module.basis_input([g.centered_cap(8),g.centered_cap(12)])
        baseline['degree8_12_probe']={'status':'unexpectedly_accepted'}
    except ValueError as e:
        baseline['degree8_12_probe']={'status':'rejected_as_outside_old_scope','message':str(e)}
    save('baseline.json',baseline)
    U=[F(0),F(2,5),F(0),F(1)]
    simple={'fixed_initial_label':'P1+(2/5)P3','u':g.strings(U),
        'degree36_moments':g.moments(g.centered_cap(36)),
        'normalization':g.normalize(g.add(g.scale(U,-3),[7])),
        'initial_ratio':g.ratio(U),'full_residual':g.residual(U),
        'direction':g.direction(U),'projected_residual':g.project_residual(U,5),
        'projected_zero_full_nonzero':g.project_residual([0,1],1),
        'zero_iteration':g.refine(U,iterations=0,degree_budget=3)}
    simple['all_certificates_verified']=all(g.verify(c) for c in simple.values() if isinstance(c,dict))
    save('results.json',simple)
    start=perf_counter();refined=g.refine(U,iterations=2,degree_budget=12,tolerance=F(1,10**7))
    replay=g.verify(refined,U)
    save('refinement.json',refined)
    summary={'refinement_verified':replay,'refinement_seconds_including_replay':perf_counter()-start,
        'initial_pi_K':simple['initial_ratio']['pi_K_lower'],
        'final_pi_K':refined['final']['pi_K_lower'],'accepted_steps':sum(x['accepted'] for x in refined['steps']),
        'stop':refined['stop'],'full_global_optimum_claimed':False}
    caps=g.screen_caps();save('caps.json',caps)
    hp=g.plane([g.centered_cap(8),g.centered_cap(12)],F(1,10**7))
    save('high_degree_plane.json',hp)
    summary['caps_verified']=g.verify(caps)
    summary['high_degree_plane_verified']=g.verify_plane(hp,[g.centered_cap(8),g.centered_cap(12)])
    diagnostics=[]
    for name,u in [('fixed_initial',U),('cap_degree12',g.centered_cap(12))]:
        start=perf_counter()
        c=g.spectral_diagnostic(u,modes=4,max_modes=4,bits=18,max_m=3,tolerance=F(1,1000))
        diagnostics.append({'case':name,'certificate':c,'verified':g.verify(c,u),
            'seconds_including_replay':perf_counter()-start})
    save('spectral_diagnostics.json',diagnostics)
    summary['diagnostics']=[{'case':r['case'],'verified':r['verified'],
        'status':r['certificate']['status'],
        'spectral_status':r['certificate'].get('spectral',{}).get('status'),
        'exact_width':r['certificate'].get('spectral',{}).get('exact_width')} for r in diagnostics]
    summary['failure_records']=[{'case':'degree8_12_old_adapter','outcome':baseline['degree8_12_probe']},
        {'case':'P1_degree_budget1','outcome':'projected_stationary_with_nonzero_full_residual'},
        {'case':'fixed_seed_zero_iteration','outcome':simple['zero_iteration']['stop']}]
    save('summary.json',summary);save('STAGES.json',g.STAGES)
    print(json.dumps(summary,indent=2))
    return summary


if __name__=='__main__': run()
