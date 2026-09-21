"""Rebuild the ten radius-stage evidence records with full certificate replay."""
import json
import sys
from pathlib import Path
from time import perf_counter
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import sphere_scale as s
import ordered_spectrum as ordered
from common import F,save

FAST=dict(modes=4,max_modes=8,max_m=4,bits=28)
proofs=[]
stages=[]

def record(name,cert):
    if not s.verify(cert):raise AssertionError('Replay failed: '+name)
    path=save('round10/certificates/'+name+'.json',cert)
    proofs.append({'path':path,'verified':True,'format':cert['format']})
    return path

def stage(version,capability,callable_,evidence,outcome):
    evidence_path=save('round10/evidence/v'+str(version)+'.json',
                       {'version':version,'actual_outcome':outcome,'mathematical_evidence':evidence})
    stages.append({'version':version,'capability':capability,'callable':callable_,
                   'evidence':evidence+[evidence_path],'actual_outcome':outcome})

start=perf_counter()
zero=s.backend_ground([0],2,**FAST)
p185=record('V185_R2_laplace',zero)
negative=s.backend_ground([-6],2,**FAST)
stage(185,'Exact sphere-dilation transfer of full mean-zero ground bounds',
      'sphere_scale.scale_ground',[p185,record('V185_R2_negative_shift',negative)],
      {'R2_q0_ground':'1/2','R2_q_minus6_ground':'-11/2','unit_only_baseline':'unsupported radius field'})

pullback=s.physical_pullback([1,-2,3],2)
stage(186,'Ambient polynomial coordinate pullback and operator scaling',
      'sphere_scale.physical_pullback',[record('V186_ambient_pullback',pullback)],
      {'pulled_back_q':pullback['pulled_back_potential'],'unit_operator_q':pullback['unit_operator_potential']})

normalized=s.scale_trial([0,0,1],2,[1])
physical=s.scale_trial([0,0,1],2,[1],coordinates='ambient_x3')
wide_trial=s.scale_trial([0]*12+[1],2,[1])
stage(187,'Natural-area explicit trial mass, Dirichlet and wide-potential energy scaling',
      'sphere_scale.scale_trial',[record('V187_normalized_t2',normalized),record('V187_physical_x3square',physical),
                                 record('V187_degree12_trial',wide_trial)],
      {'same_trial_u':'x3/R','normalized_q_rayleigh':normalized['rayleigh_quotient'],
       'physical_q_rayleigh':physical['rayleigh_quotient'],'degree12_rayleigh':wide_trial['rayleigh_quotient']})

unit_ordered=ordered.adaptive_spectrum([-6],k=8,modes=2,max_modes=3,max_radial=3,bits=20)['certificate']
scaled_ordered=s.scale_ordered(['-3/2'],2,unit_ordered)
unit_ordered2=ordered.adaptive_spectrum([-24],k=24,modes=2,max_modes=4,max_radial=4,bits=20)['certificate']
scaled_ordered2=s.scale_ordered([-6],2,unit_ordered2)
stage(188,'Preserve ordered eigenvalue multiplicities and scale the complete negative trace',
      'sphere_scale.scale_ordered',[record('V188_scaled_minus6_unit_spectrum',scaled_ordered),
                                   record('V188_R2_physical_minus6_spectrum',scaled_ordered2)],
      {'q_minus3over2_R2_negative_trace':scaled_ordered['negative_trace_lower'],
       'q_minus3over2_R2_strict_negative_count':scaled_ordered['negative_count_lower'],
       'q_minus6_R2_negative_trace':scaled_ordered2['negative_trace_lower'],
       'q_minus6_R2_strict_negative_count':scaled_ordered2['negative_count_lower'],'zero_levels_excluded':True})

zero_cell=s.radius_cell([0],[1,2],s.backend_ground([0],'3/2',**FAST))
negative_cell=s.radius_cell([-6],[1,2],s.backend_ground([-6],'3/2',**FAST))
negative_variable=s.radius_cell([-6,1],[1,2],s.backend_ground([-6,1],'3/2',**FAST))
stage(189,'Uniform radius perturbation enclosure with sign-aware interval products',
      'sphere_scale.radius_cell',[record('V189_laplace_radius_cell',zero_cell),
                                 record('V189_negative_radius_cell',negative_cell),
                                 record('V189_nonconstant_negative_cell',negative_variable)],
      {'R_in_1_2_laplace_range':[zero_cell['lower'],zero_cell['upper']],
       'negative_constant_range':[negative_cell['lower'],negative_cell['upper']],
       'negative_nonconstant_bound':[negative_variable['lower'],negative_variable['upper']],
       'uses_all_four_signed_products':True})

baseline=s.adaptive_minimum([0,0,1],[1,2],tolerance='1/10000',max_splits=0,**FAST)
refined=s.adaptive_minimum([0,0,1],[1,2],tolerance='1/10000',max_splits=18,**FAST)
physical_open=s.adaptive_minimum([0,0,1],[1,2],coordinates='ambient_x3',tolerance='1/10000',max_splits=8,**FAST)
zero_global=s.adaptive_minimum([0],[1,2],max_splits=0,**FAST)
stage(190,'Complete-coverage adaptive global minimization over positive radii',
      'sphere_scale.adaptive_minimum',[record('V190_initial_nonconstant_global',baseline['certificate']),
                                      record('V190_refined_nonconstant_global',refined['certificate']),
                                      record('V190_physical_nonconstant_open',physical_open['certificate']),
                                      record('V190_exact_laplace_minmax',zero_global['certificate'])],
      {'initial_width':baseline['certificate']['width'],'refined_width':refined['certificate']['width'],
       'width_reduction_factor':str(F(baseline['certificate']['width'])/F(refined['certificate']['width'])),
       'nonconstant_global_bracket':[refined['certificate']['lower'],refined['certificate']['upper']],
       'status':refined['certificate']['status'],'target':'1/10000',
       'physical_coordinate_search_status':physical_open['certificate']['status'],
       'retained_partition_cells':len(refined['certificate']['cells'])})
save('round10/V190_search_history.json',{'normalized':refined['search_history'],
      'normalized_work':refined['work'],'physical':physical_open['search_history'],'physical_work':physical_open['work']})

first=s.adaptive_minimum([0,0,1],[1,2],max_splits=1,tolerance='1/10000',**FAST)
cp=s.checkpoint(first['certificate'])
resumed=s.resume_minimum(cp,max_splits=2,**FAST)
repeated=s.resume_minimum(s.checkpoint(refined['certificate']),max_splits=0,**FAST)
stage(191,'Verified cache continuation retaining old incumbents and monotone parent bounds',
      'sphere_scale.resume_minimum',[record('V191_checkpoint',cp),record('V191_resumed_global',resumed['certificate']),
                                    record('V191_replayed_global',repeated['certificate'])],
      {'resumed_work':resumed['work'],'zero_step_resume_work':repeated['work'],
       'old_lower':first['certificate']['lower'],'new_lower':resumed['certificate']['lower'],
       'old_upper':first['certificate']['upper'],'new_upper':resumed['certificate']['upper'],
       'review_fix':'Retain incomparable old/new proofs at one radius; width alone cannot replace an incumbent.'})

ce=s.transfer_counterexample(s.scale_trial([0,0,'1/2'],2,[1,'-1/36'],m=1,degrees=[1,3]),'3/5')
stage(192,'Explicit scaled function counterexamples with threshold binding',
      'sphere_scale.transfer_counterexample',[record('V192_explicit_R2_counterexample',ce)],
      {'rayleigh':ce['rayleigh_quotient'],'threshold':ce['threshold'],'strict_margin':ce['strict_margin'],
       'status':ce['status'],'function':'[P_1^1(x3/R)-P_3^1(x3/R)/36]*cos(phi)'})

precise=s.backend_ground([0,0,'1/2'],2,tolerance=F(1,10**18),modes=4,max_modes=16,max_m=4)
broad=s.backend_ground([0]*12+[1],1,**FAST)
stage(193,'Choose and bind suitable exact backends by degree and required scaled accuracy',
      'sphere_scale.backend_ground',[record('V193_precision_backend',precise),record('V193_wide_backend',broad)],
      {'precision_backend':precise['unit_backend'],'precision_physical_width':precise['width'],
       'wide_backend':broad['unit_backend'],'wide_degree':12,
       'precision_target_met':F(precise['width'])<=F(1,10**18)})

proved=s.research_job([0,0,1],'3/5',domain=[1,2],saved=s.checkpoint(refined['certificate']),max_splits=0,**FAST)
refuted=s.research_job([0,0,1],'7/10',radius=2,trial={'coefficients':[1,'-1/36'],'m':1,'degrees':[1,3]},**FAST)
undetermined=s.research_job([0,0,1],'3/5',domain=[1,2],max_splits=0,**FAST)
stage(194,'Unified radius research jobs with content-bound proof replay and honest tri-state decisions',
      'sphere_scale.research_job',[record('V194_proved_interval_inequality',proved['certificate']),
                                  record('V194_refuted_explicit_inequality',refuted['certificate']),
                                  record('V194_undetermined_interval_inequality',undetermined['certificate'])],
      {'decisions':[proved['certificate']['status'],refuted['certificate']['status'],undetermined['certificate']['status']],
       'proved_resume_work':proved['work'],'model_token_measurement':None,
       'work_measurement_is_only_local_anchor_computations':True})

for item in proofs:
    cert=json.loads((ROOT/item['path']).read_text())
    if not s.verify(cert):raise AssertionError('Saved proof failed replay '+item['path'])
summary={'cycle':10,'versions':list(range(185,195)),'stage_count':len(stages),
         'proof_records':proofs,'all_saved_proofs_replayed':True,
         'elapsed_seconds_including_generation_and_replay':perf_counter()-start,
         'no_model_api_calls':True,'model_token_measurement':None,
         'mathematical_scope':'Full real mean-zero H1(S2_R), fixed polynomial potential or explicitly quantified radius interval',
         'remaining_gap':'The nonconstant radius global minimum target 1e-4 remains open at 18 splits.'}
save('round10/STAGES.json',stages)
save('round10/RESULTS.json',summary)
print(json.dumps({'stages':len(stages),'saved_proofs':len(proofs),'all_replayed':True,
      'elapsed_seconds':summary['elapsed_seconds_including_generation_and_replay']},indent=2))
