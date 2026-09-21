"""V10: budget allocation driven by certified errors; predictions remain untrusted."""
from fractions import Fraction as F
from math import ceil, log, isfinite
from time import perf_counter
import v03_tail as v3
import v04_search as v4
import v05_range as v5
import v06_parity as v6
import v07_precision as v7


def next_size(history,tolerance,max_modes):
    current=history[-1]['modes'];step=2
    if len(history)>=2:
        previous=history[-2]
        w0,w1=F(previous['exact_width']),F(history[-1]['exact_width'])
        if 0<w1<w0:
            slope=(log(float(w1))-log(float(w0)))/(current-previous['modes'])
            proposed=(log(float(tolerance/4))-log(float(w1)))/slope
            if isfinite(proposed):step=min(8,max(2,2*ceil(proposed/2)))
    return min(max_modes,current+step)


def adaptive(q,k=1,tolerance=F(1,10**10),max_modes=32):
    v3.sizes(k,max_modes);tolerance=F(tolerance)
    if not F(1,10**40)<=tolerance<=1:raise ValueError('Tolerance must be in [1e-40,1]')
    n=max(4,k+2)
    if n>max_modes:raise ValueError('Mode budget too small')
    start=perf_counter();q=v3.base.potential([str(x) for x in q]);proof=v5.make_range(q)
    parity=not any(q[i] for i in range(1,len(q),2));factory=v6.ParityKernel if parity else v3.Kernel
    history=[];best=None
    while True:
        r=v4.guided(q,k,n,tolerance,range_proof=proof,kernel_factory=factory);c=r['certificate']
        history.append({'modes':n,'exact_width':c['exact_width'],'path':r['path'],
                        'exact_count_calls':r['exact_count_calls'],'seconds':r['seconds']})
        if best is None or F(c['exact_width'])<F(best['exact_width']):best=c
        if F(best['exact_width'])<=tolerance or n==max_modes:break
        proposed=next_size(history,tolerance,max_modes)
        # Prediction errors affect cost only. Enforce progress and the hard budget.
        n=min(max_modes,max(n+1,proposed)) if type(proposed) is int else min(max_modes,n+2)
    if not v3.verify(best):raise ArithmeticError('Final adaptive certificate rejected')
    return {'algorithm_version':10,'route':'spectral','status':'target_met' if F(best['exact_width'])<=tolerance else 'target_not_met',
            'tolerance':str(tolerance),'max_modes':max_modes,'parity_blocks':parity,
            'certificate':best,'attempts':history,'seconds':perf_counter()-start}


def function(q,tolerance=F(1,10**40),start_modes=8,max_modes=24,bits=160):
    """Choose between Newton coefficient refinement and spatial enrichment.

    The requested target is the eigenvalue enclosure produced by the function;
    the independent excited-weight bound is also retained, not equated to it.
    """
    v3.base.check_sizes(1,max_modes)
    if type(start_modes) is not int or not 4<=start_modes<=max_modes:raise ValueError('Invalid initial mode count')
    if type(bits) is not int or not 48<=bits<=224:raise ValueError('Coefficient bits must be 48..224')
    tolerance=F(tolerance)
    if not F(1,10**80)<=tolerance<=1:raise ValueError('Function-derived eigenvalue tolerance must be [1e-80,1]')
    start=perf_counter();q=v3.base.potential([str(x) for x in q]);best=None;history=[];n=start_modes
    while True:
        try:
            result=v7.refine(q,n,bits,0);c=result['certificate']
            history.append({'modes':n,'action':'initial_trial','exact_width':c['exact_width']})
            if best is None or F(c['exact_width'])<F(best['exact_width']):best=c
            for step in range(2):
                if F(c['exact_width'])<=tolerance:break
                inside=F(c['inside_residual_squared']);outside=F(c['outside_residual_squared'])
                if inside<outside:break
                coefficients=[F(x) for x in c['trial_legendre_coefficients']]
                proposal=v7.newton_step(q,coefficients,bits)
                refined=v7.certificate(q,proposal,c['gap_proof'])
                history.append({'modes':n,'action':'coefficient_refinement','exact_width':refined['exact_width']})
                if F(refined['exact_width'])>=F(c['exact_width']):break
                c=refined
                if F(c['exact_width'])<F(best['exact_width']):best=c
        except ValueError as exc:
            history.append({'modes':n,'action':'unavailable_at_this_size','reason':str(exc)})
        if (best is not None and F(best['exact_width'])<=tolerance) or n==max_modes:break
        n=min(max_modes,n+4)
    if best is not None and not v7.verify(best):raise ArithmeticError('Function certificate rejected')
    met=best is not None and F(best['exact_width'])<=tolerance
    return {'algorithm_version':10,'route':'function','status':'target_met' if met else 'target_not_met',
            'tolerance':str(tolerance),'max_modes':max_modes,'coefficient_bits':bits,
            'certificate':best,'attempts':history,'seconds':perf_counter()-start}
