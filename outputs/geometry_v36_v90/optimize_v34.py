"""Adaptive SDP proposals + certified all-interval envelopes + global dual gap."""
from fractions import Fraction as F
from time import perf_counter
import math
import algebra as a
import dual_v33 as dual
import float_sdp
import family_v29 as family
import matrix_v30


def propose(directions,cost,samples,tolerance):
    dirs,v,c=dual.setup(directions,cost);d=len(dirs);aa=[dual.at(v,t) for t in samples]
    mu=F(tolerance)/(16*d*len(samples))
    raw,actual_mu,stats=float_sdp.solve([[[float(z) for z in row] for row in m] for m in aa],[[float(z) for z in row] for row in c],float(mu))
    g=[[F(round(z*10**12),10**12) for z in row] for row in raw]
    # Quantization may cross a finite constraint. Add a bounded exact diagonal repair.
    shift=F(0)
    for _ in range(48):
        candidate=[[g[i][j]+(shift if i==j else 0) for j in range(d)] for i in range(d)]
        if all(a.inertia([[candidate[i][j]-m[i][j] for j in range(d)] for i in range(d)])==[0,0,d] for m in aa):g=candidate;break
        shift=F(1,10**12) if not shift else 2*shift
    else:raise ValueError('Could not certify a positive finite-sample candidate')
    atoms=[];rational_mu=F(str(actual_mu))
    for t,m in zip(samples,aa):
        inv=a.inverse([[g[i][j]-m[i][j] for j in range(d)] for i in range(d)])
        atoms.append({'t':str(t),'matrix':dual.encode([[rational_mu*z for z in row] for row in inv])})
    return g,dual.repair_dual(dirs,c,atoms),stats


def separation_sample(template):
    proof=template['range_proof'];location=F(proof['witness'])
    if proof['coordinate']=='t_squared':
        location=F(round(math.sqrt(float(location))*10**8),10**8)
    return max(F(-1),min(F(1),location))


def solve(directions,cost=None,tolerance=F(1,10**5),max_rounds=8,max_range_leaves=128):
    start=perf_counter();tolerance=F(tolerance)
    if not F(1,10**8)<=tolerance<=F(1,10):raise ValueError('Optimizer accuracy budget: 1e-8 to 1e-1')
    if type(max_rounds) is not int or not 1<=max_rounds<=24:raise ValueError('One to 24 exchange rounds')
    dirs,v,c=dual.setup(directions,cost);d=len(dirs);samples={F(i,4) for i in range(-4,5)};trace=[];best_template=None;best_upper=None;best_dual=None
    # Known-valid fallback even if all floating proposals fail.
    initial=matrix_v30.synthesize(dirs,tolerance=F(1,10**8),max_leaves=max_range_leaves)
    best_template=initial;best_upper=dual.inner(c,[[F(z) for z in row] for row in initial['quadratic_bound']])
    zero={'t':'0','matrix':[['0']*d for _ in range(d)]};best_dual=dual.lower_certificate(dirs,c,[zero]);total_steps=0
    for index in range(max_rounds):
        try:g,lower,stats=propose(dirs,c,sorted(samples),tolerance)
        except (ValueError,OverflowError,ZeroDivisionError) as exc:
            trace.append({'round':index+1,'outcome':'proposal_failed','reason':str(exc)});break
        total_steps+=stats['newton_steps']
        # Relative envelope accuracy adapts to the objective scale.
        scale=max(F(1),dual.inner(c,g));range_tol=max(F(1,10**12),min(F(1,10**6),tolerance/(32*scale)))
        template=matrix_v30.synthesize(dirs,g,tolerance=range_tol,max_leaves=max_range_leaves)
        upper=dual.inner(c,[[F(z) for z in row] for row in template['quadratic_bound']])
        if upper<best_upper:best_upper=upper;best_template=template
        if F(lower['lower_bound'])>F(best_dual['lower_bound']):best_dual=lower
        cert=dual.certificate(best_template,best_dual,tolerance)
        trace.append({'round':index+1,'sample_count':len(samples),'lower':cert['lower_bound'],'upper':cert['upper_bound'],'gap':cert['gap'],
                      'envelope_scale':template['scale'],'range_status':template['range_proof']['status'],**stats})
        if cert['status']=='global_gap_closed':break
        t=separation_sample(template);before=len(samples);samples.add(t)
        # Opposite signs cover parity pairs; rational samples are proposals, not coverage proofs.
        samples.add(-t)
        if len(samples)==before:
            trace[-1]['stop_reason']='no_new_separation_sample';break
    cert=dual.certificate(best_template,best_dual,tolerance)
    return {'algorithm_version':34,'certificate':cert,'statistics':{'exchange_rounds':len(trace),'newton_steps':total_steps,'sample_count':len(samples),'model_calls':0},
            'trace':trace,'seconds':perf_counter()-start}
