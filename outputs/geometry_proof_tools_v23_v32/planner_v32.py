"""V32: best-first exact reductions, lemma synthesis/transport, then decision search.

All routing is heuristic. The returned proof chain is reconstructed independently
of the search history. This program makes no model/API calls.
"""
from fractions import Fraction as F
from time import perf_counter
import copy,heapq,json
import compiler_v23 as compiler
import reductions
import family_v29 as family
import matrix_v30
import transport_v31 as transport
import decision_v28 as decision
import algebra
import error_bounds as e


def key(p):return json.dumps(p,sort_keys=True,separators=(',',':'))


def score(p):
    d=len(p['directions']);volume=F(1)
    for lo,hi in p['domain'].get('box',[]):volume*=F(hi)-F(lo)
    return d,int(p['function_space']=='full_sphere'),volume,sum(map(len,p['directions']))


def reduction_search(p,max_states=16):
    if type(max_states) is not int or not 1<=max_states<=64:raise ValueError('Reduction search budget')
    queue=[(score(p),0,p,[])];seen={key(p)};states=[];trace=[];serial=0
    while queue and len(states)<max_states:
        _,_,current,path=heapq.heappop(queue);states.append((current,path))
        for name,fn in reductions.RULES.items():
            try:edge=fn(current)
            except ValueError as exc:
                trace.append({'rule':name,'outcome':'precondition_not_met','detail':str(exc)});continue
            nxt=edge['reduced'];k=key(nxt)
            if k in seen:continue
            seen.add(k);serial+=1;heapq.heappush(queue,(score(nxt),serial,nxt,path+[edge]));trace.append({'rule':name,'outcome':'new_state','dimension':len(nxt['directions'])})
    return sorted(states,key=lambda item:score(item[0])),trace


def template_directions(p):
    _,dirs,_,_,_=reductions.data(p);stripped=[]
    for q in dirs:
        q=q[:];q[0]=F(0)
        if any(q):stripped.append(e.trim(q))
    if not stripped:return []
    n=max(map(len,stripped));rows=[q+[F(0)]*(n-len(q)) for q in stripped]
    return [stripped[i] for i in algebra.independent_rows(rows)]


def terminal_bounds(terminal,p,threshold):
    proof=terminal['proof'];kind=terminal['kind']
    if kind=='transport':
        if proof['problem']!=p or not transport.verify(proof):raise ValueError('Invalid transported lemma')
        return F(proof['lower_bound']),None,None,None
    if kind=='decision':
        if proof['problem']!=p or proof['threshold']!=str(threshold) or not decision.verify(proof):raise ValueError('Invalid terminal decision')
        low=F(proof['lower_bound']);upper=None if proof['upper_bound'] is None else F(proof['upper_bound']);x=proof['candidate_parameters'];trial=None
        if upper is not None:
            inner=proof['inner']
            if proof['proof_kind']=='fixed':trial=inner['record']['trial_legendre']
            else:trial=inner['points'][inner['candidate_index']]['trial_legendre']
        return low,upper,None if x is None else list(map(F,x)),trial
    raise ValueError('Unknown terminal rule')


def certificate(compiled,path,terminal):
    if not compiler.verify(compiled) or not isinstance(path,list) or len(path)>12:raise ValueError('Invalid compiled task or proof path')
    original=compiled['problem'];p=original;threshold=F(compiled['threshold'])
    for edge in path:
        if edge['original']!=p or not reductions.verify(edge):raise ValueError('Broken reduction chain')
        p=edge['reduced']
    low,upper,x,trial=terminal_bounds(terminal,p,threshold)
    if x is not None:
        for edge in reversed(path):x=reductions.lift(edge,x)
        q0,dirs,c,l,h=reductions.data(original);q=q0[:]
        if len(x)!=len(dirs):raise ValueError('Witness dimension mismatch')
        if original['domain']['kind']=='box' and any(not F(lo)<=v<=F(hi) for v,(lo,hi) in zip(x,original['domain']['box'])):raise ValueError('Lifted witness outside original domain')
        for z,direction in zip(x,dirs):q=e.add(q,e.scale(direction,z))
        actual=e.residual_statistics(q,list(map(F,trial)))[0]+algebra.value(c,l,h,x)
        if actual!=upper:raise ValueError('Lifted witness does not have the claimed original energy')
    status='proved' if low>=threshold else ('refuted' if upper is not None and upper<threshold else 'unresolved')
    return {'method':'proof_search_v32','compiled_goal':compiled,'reductions':path,'terminal':terminal,
            'status':status,'lower_bound':str(low),'upper_bound':None if upper is None else str(upper),
            'candidate_parameters':None if x is None else list(map(str,x)),
            'original_function_space':original['function_space'],'threshold':str(threshold),'formal_assistant_checked':False}


def verify(cert):
    try:return cert==certificate(cert['compiled_goal'],cert['reductions'],cert['terminal'])
    except (ValueError,KeyError,TypeError,IndexError,ZeroDivisionError,RecursionError):return False


def solve(raw,library=None,synthesize=True,max_states=16,max_leaves=64,max_modes=24):
    start=perf_counter();compiled=compiler.compile_goal(raw);p=compiled['problem'];threshold=F(compiled['threshold']);library=copy.deepcopy(library or [])
    if len(library)>32 or any(not family.verify(t) for t in library):raise ValueError('Invalid lemma library')
    states,trace=reduction_search(p,max_states);stats={'expanded_states':len(states),'template_attempts':0,'template_syntheses':0,'points':0,'spectral_searches':0,'model_calls':0}
    attempted=set();synthesized=set()
    for state,path in states:
        if synthesize:
            directions=template_directions(state);signature=tuple(tuple(q) for q in directions)
            if directions and signature not in synthesized:
                synthesized.add(signature)
                template=family.synthesize(directions[0]) if len(directions)==1 else matrix_v30.synthesize(directions)
                if all(t['template_id']!=template['template_id'] for t in library):library.append(template);stats['template_syntheses']+=1
        for template in library:
            signature=(key(state),template['template_id'])
            if signature in attempted:continue
            attempted.add(signature);stats['template_attempts']+=1
            try:proof=transport.apply(state,template)
            except ValueError as exc:
                trace.append({'route':'lemma','outcome':'not_applicable_or_too_weak','detail':str(exc)});continue
            trace.append({'route':'lemma','outcome':'bound','lower':proof['lower_bound']})
            if F(proof['lower_bound'])>=threshold:
                cert=certificate(compiled,path,{'kind':'transport','proof':proof})
                return {'algorithm_version':32,'certificate':cert,'statistics':stats,'trace':trace,'generated_library':library,'seconds':perf_counter()-start}
    # This prototype executes the smallest admissible numerical representation;
    # it does not claim the ranking is a globally optimal allocation of work.
    failures=[]
    for state,path in states:
        if state['function_space']!='axisymmetric':continue
        try:result=decision.search(state,threshold,max_leaves=max_leaves,max_modes=max_modes)
        except ValueError as exc:failures.append(str(exc));continue
        cert=certificate(compiled,path,{'kind':'decision','proof':result['certificate']})
        stats.update({k:result['statistics'].get(k,0) for k in ('points','spectral_searches')})
        trace.append({'route':'numerical_decision','outcome':cert['status'],'dimension':len(state['directions'])})
        return {'algorithm_version':32,'certificate':cert,'statistics':stats,'trace':trace,'generated_library':library,'seconds':perf_counter()-start}
    raise ValueError('No implemented route supports this goal: '+'; '.join(failures))
