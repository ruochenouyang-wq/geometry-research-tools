"""Reproduce the ten continuum-family capability demonstrations."""
from pathlib import Path
import sys
import json
import time
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import F, save
import parameter_family as p
import witness


def run():
    started = time.perf_counter()
    family = p.validate_family([0, 0, 2], [[0, 1]], [[-1, 1]])
    cache, stages, certificates = {}, [], []
    def record(version, capability, callable_name, certs, outcome):
        paths = []
        for name, certificate in certs.items():
            if not p.verify(certificate):
                raise ArithmeticError('Rejected '+name)
            path = save('round04/certificates/'+name+'.json', certificate)
            paths.append(path)
            certificates.append({'path': path, 'verified': True, 'format': certificate['format']})
        stages.append({'version': version, 'capability': capability, 'callable': callable_name,
                       'evidence': paths, 'outcome': outcome})
    cheap = p.cheap_family_enclosure(family)
    record(125, 'Compute a search-free whole-family spectral enclosure from coefficient boxes and a fixed admissible harmonic.',
           'parameter_family.cheap_family_enclosure', {'v125_family': family, 'v125_cheap_enclosure': cheap},
           'q_s(t)=2t²+s t, every real s in [-1,1]: initial entire-range spectral enclosure ['+cheap['lower']+', '+cheap['upper']+'].')
    norms = p.direction_norms(family)
    cancellation = p.validate_family([0], [[0, 1, 0, -1]], [[-1, 1]])
    cancellation_norms = p.direction_norms(cancellation)
    record(126, 'Certify multiplier L-infinity bounds by rational polynomial extrema.',
           'parameter_family.direction_norms', {'v126_norms': norms, 'v126_cancellation_norms': cancellation_norms},
           'For t-t³ the norm upper is '+cancellation_norms['norm_upper'][0]+' instead of the coefficient-sum bound 2.')
    anchor = p.anchor_cell(family, cache=cache, modes=6, bits=36, norms=norms)
    record(127, 'Propagate a verified point spectrum to every point of a whole cell.',
           'parameter_family.anchor_cell', {'v127_anchor': anchor},
           'Certified midpoint Lipschitz interval; uniform lower='+anchor['lower'])
    partition = p.partition_certificate(family, [[[-1, 0]], [[0, 1]]])
    record(128, 'Prove exact rational box coverage with disjoint interiors.',
           'parameter_family.partition_certificate', {'v128_partition': partition},
           'Both closed halves cover the original continuum; missing/overlapping cells rejected in tests.')
    lower = p.cell_lower_bound(family, cache=cache, modes=6, bits=36)
    record(129, 'Use concavity to certify a uniform lower from every vertex.',
           'parameter_family.cell_lower_bound', {'v129_concavity': lower},
           'Parameter minimum in ['+lower['lower']+', '+lower['minimum_upper']+']; vertex maxima are not claimed as upper bounds.')
    center_proposal = witness.mass_normalized_proposal([0, 0, 2], m=1, modes=6)
    center_trial = witness.quantize_proposal(center_proposal)['certificate']
    majorant = p.affine_majorant(family, center_trial)
    plus = p.affine_majorant(family, witness.rayleigh_certificate([0], [1, F(1, 8)], m=1))
    minus = p.affine_majorant(family, witness.rayleigh_certificate([0], [1, -F(1, 8)], m=1))
    record(130, 'Produce exact affine Rayleigh majorants from fixed explicit functions.',
           'parameter_family.affine_majorant', {'v130_center_majorant': majorant,
             'v130_plus_majorant': plus, 'v130_minus_majorant': minus},
           'Mean-zero m=1 trials cover nonaxisymmetric functions; center majorant has zero s slope by parity.')
    envelope = p.envelope_1d(family, [plus, minus])
    best_envelope = p.envelope_1d(family, [plus, minus, majorant])
    record(131, 'Certify the maximum of a one-dimensional minimum affine envelope.',
           'parameter_family.envelope_1d', {'v131_intersection_envelope': envelope,
             'v131_tighter_envelope': best_envelope},
           'All pairwise rational line intersections and endpoints tested; maximum upper='+best_envelope['maximum_upper'])
    family2 = p.validate_family([0, 0, 2], [[0, 1], [1]], [[-1, 1], [F(-1, 4), F(1, 4)]])
    majorant2 = p.affine_majorant(family2, center_trial)
    upper2 = p.envelope_box(family2, [majorant2])
    record(132, 'Provide a safe multi-parameter cell upper from affine majorants.',
           'parameter_family.envelope_box', {'v132_two_parameter_upper': upper2},
           'Two-parameter full-box upper='+upper2['maximum_upper']+'; no false maximum-of-vertex-minima rule.')
    budget = p.adaptive_family(family, F(1, 10**10), max_cells=2, modes=4, bits=28)
    adaptive = p.adaptive_family(family, F(1, 5), max_cells=8, modes=4, bits=28)
    accelerated = p.adaptive_family(family, F(1, 10**7), max_cells=8, majorants=[majorant],
                                    modes=6, bits=36, cache=cache)
    record(133, 'Adaptively subdivide the global-supremum interval using reusable anchor proofs.',
           'parameter_family.adaptive_family', {'v133_budget_open': budget,
             'v133_subdivided': adaptive, 'v133_majorant_accelerated': accelerated},
           'Two-cell tight request remains budget_open; ordinary subdivision meets 1/5 at 8 cells; affine majorant meets 1e-7 at '+str(len(accelerated['cells']))+' cell(s).')
    proved = p.universal_inequality(family, F(23, 10), lower_proof=lower, modes=6)
    refuted = p.universal_inequality(family, F(12, 5), lower_proof=lower, modes=6)
    unresolved = p.universal_inequality(family, F(lower['lower'])+F(1, 10**20), lower_proof=lower, modes=6)
    record(134, 'Resolve a universally quantified weighted inequality or exhibit parameter and function.',
           'parameter_family.universal_inequality', {'v134_proved': proved, 'v134_refuted': refuted,
             'v134_undetermined': unresolved},
           'Threshold 2.3 proved; threshold 2.4 explicitly refuted; a boundary-straddling threshold preserved as undetermined.')
    summary = {'problem': 'Continuum weighted Poincare: q_s(t)=2t²+s t, s∈[-1,1]',
               'stages': 10, 'certificate_count': len(certificates), 'certificates': certificates,
               'midpoint_lipschitz_lower': anchor['lower'], 'concavity_uniform_lower': lower['lower'],
               'parameter_minimum_upper': lower['minimum_upper'], 'parameter_minimum_gap': lower['minimum_gap'],
               'parameter_supremum_lower': accelerated['maximum_lower'],
               'parameter_supremum_upper': accelerated['maximum_upper'],
               'parameter_supremum_gap': accelerated['maximum_gap'],
               'vertex_maximum_upper_is_invalid': F(anchor['anchor_proof']['lower']) > max(F(v['spectrum']['upper']) for v in lower['vertices']),
               'universal_threshold_2_3': proved['status'], 'universal_threshold_2_4': refuted['status'],
               'counterexample_parameter': refuted['counterexample']['parameter'],
               'counterexample_rayleigh': refuted['counterexample']['trial']['rayleigh_quotient'],
               'counterexample_harmonic_terms': refuted['counterexample']['trial']['harmonic_terms'],
               'development_run_seconds': time.perf_counter()-started,
               'not_a_new_GN_theorem': True, 'not_an_LLM_token_or_model_speed_benchmark': True}
    save('round04/STAGES.json', stages)
    save('round04/RESULTS.json', summary)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('certificates','counterexample_harmonic_terms')}, indent=2))
    return summary


if __name__ == '__main__':
    run()
