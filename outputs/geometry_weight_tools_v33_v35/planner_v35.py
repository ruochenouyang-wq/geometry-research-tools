"""Goal-aware weights, a scoped obstruction diagnostic, and verified fallback."""
from fractions import Fraction as F
from time import perf_counter
import compiler_v23 as compiler
import reductions
import matrix_v30
import family_v29 as family
import transport_v31 as transport
import planner_v32 as baseline
import dual_v33 as dual


def budget_certificate(p,template):
    p=compiler.normalized(p)
    if not family.verify(template) or template['directions']!=p['directions']:raise ValueError('Budget template must describe the current parameter directions')
    _,_,_,_,h=reductions.data(p);weight=[[z/2 for z in row] for row in h]
    if template['weight']!=dual.encode(weight):raise ValueError('Template weight is not the goal penalty Hessian / 2')
    low=F(template['range_proof']['maximum_lower']);upper=F(template['range_proof']['maximum_upper'])
    status='poisson_envelope_fits' if upper<=1 else ('poisson_envelope_incompatible' if low>1 else 'budget_test_unresolved')
    return {'method':'penalty_budget_v35','problem':p,'template':template,'status':status,
            'scale_lower':str(low),'scale_upper':str(upper),
            'meaning':'Tests A(t) <= penalty Hessian / 2 for this fixed Poisson family; does not decide the original spectral goal.',
            'formal_assistant_checked':False}


def verify_budget(cert):
    try:return cert==budget_certificate(cert['problem'],cert['template'])
    except (ValueError,KeyError,TypeError,IndexError,ZeroDivisionError):return False


def solve(raw,library=None,use_penalty_weights=True,synthesize=True,max_states=16,max_leaves=64,max_modes=24):
    start=perf_counter();compiled=compiler.compile_goal(raw);p=compiled['problem'];threshold=F(compiled['threshold'])
    library=list(library or [])
    if len(library)>32 or any(not family.verify(t) for t in library):raise ValueError('Invalid lemma library')
    diagnostics=[];attempts=0;seen=set()
    if use_penalty_weights and synthesize:
        states,trace=baseline.reduction_search(p,max_states)
        for state,path in states:
            _,dirs,_,_,h=reductions.data(state)
            if not dirs:continue
            weight=[[z/2 for z in row] for row in h];signature=(tuple(map(tuple,dirs)),tuple(map(tuple,weight)))
            if signature in seen:continue
            seen.add(signature)
            try:dual.setup(dirs,weight)
            except ValueError:continue
            attempts+=1;template=matrix_v30.synthesize(dirs,weight)
            diagnostic=budget_certificate(state,template);diagnostics.append(diagnostic)
            try:proof=transport.apply(state,template)
            except ValueError:continue
            if F(proof['lower_bound'])>=threshold:
                cert=baseline.certificate(compiled,path,{'kind':'transport','proof':proof})
                return {'algorithm_version':35,'certificate':cert,'strategy':'penalty_weighted_lemma',
                        'statistics':{'weighted_attempts':attempts,'expanded_states':len(states),'points':0,'spectral_searches':0,'model_calls':0},
                        'diagnostics':diagnostics,'generated_library':[template],'seconds':perf_counter()-start}
    result=baseline.solve(raw,library=library,synthesize=synthesize,max_states=max_states,max_leaves=max_leaves,max_modes=max_modes)
    result['algorithm_version']=35;result['strategy']='v32_fallback';result['statistics']['weighted_attempts']=attempts
    result['diagnostics']=diagnostics;result['seconds']=perf_counter()-start
    return result
