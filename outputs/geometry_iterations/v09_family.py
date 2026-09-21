"""V9: rigorous coverage of affine potential families, not sampled assertions."""
from fractions import Fraction as F
from collections import deque
from time import perf_counter
import v05_range as v5
import v03_tail as v3
import v04_search as v4


def at_parameter(q0,direction,a):
    q=v3.old.add(q0,v3.old.scale(direction,a))
    return v3.base.potential([str(x) for x in q])


def setup(q0,direction,interval,k,threshold,range_proof):
    q0=v3.base.potential([str(x) for x in q0]);direction=v3.base.potential([str(x) for x in direction])
    if not isinstance(interval,(list,tuple)) or len(interval)!=2:raise ValueError('Two parameter endpoints required')
    a,b=[v3.base.rational(str(x)) for x in interval]
    if a>=b:raise ValueError('Parameter interval must have positive length')
    if type(k) is not int or not 1<=k<=32:raise ValueError('Eigenvalue index must be 1..32')
    threshold=v3.base.rational(str(threshold));lo,hi=v3.range_bounds(direction,range_proof)
    norm=max(abs(lo),abs(hi))
    metadata={'model':v3.base.MODEL,'scope':v3.base.SCOPE,'q0_coefficients':[str(x) for x in q0],
              'direction_coefficients':[str(x) for x in direction],
              'parameter_interval':[str(a),str(b)],'eigenvalue_index':k,'threshold':str(threshold),
              'direction_range_proof':range_proof,'lipschitz_upper':str(norm),'formal_assistant_checked':False}
    return q0,direction,a,b,threshold,norm,metadata


def checked_point(q0,direction,k,anchor,cert):
    if not v3.verify(cert) or cert['eigenvalue_index']!=k or v3.base.potential(cert['q_coefficients'])!=at_parameter(q0,direction,anchor):
        raise ValueError('Point certificate must match parameter, potential, and index')
    return F(cert['lower']),F(cert['upper'])


def certificate(q0,direction,interval,k,threshold,range_proof,leaves):
    q0,direction,a,b,threshold,norm,meta=setup(q0,direction,interval,k,threshold,range_proof)
    if not isinstance(leaves,list) or not 1<=len(leaves)<=256:raise ValueError('At most 256 leaves')
    expected=[];previous=a
    for leaf in leaves:
        left,right=[v3.base.rational(str(x)) for x in leaf['interval']]
        if left!=previous or right<=left or right>b:raise ValueError('Leaves must cover interval without gaps or overlap')
        anchor=(left+right)/2;lo,hi=checked_point(q0,direction,k,anchor,leaf['point_certificate'])
        lower=lo-norm*(right-left)/2
        expected.append({'interval':[str(left),str(right)],'anchor':str(anchor),
                         'point_certificate':leaf['point_certificate'],'uniform_lower':str(lower),
                         'proves_threshold':lower>=threshold})
        previous=right
    if previous!=b:raise ValueError('Incomplete parameter cover')
    return {'method':'uniform_affine_family_v9',**meta,'leaves':expected,
            'status':'proved' if all(x['proves_threshold'] for x in expected) else 'unresolved',
            'uniform_lower':str(min(F(x['uniform_lower']) for x in expected))}


def counterexample(q0,direction,interval,k,threshold,range_proof,anchor,point_certificate):
    q0,direction,a,b,threshold,norm,meta=setup(q0,direction,interval,k,threshold,range_proof)
    anchor=v3.base.rational(str(anchor))
    if not a<=anchor<=b:raise ValueError('Counterexample outside parameter interval')
    lo,hi=checked_point(q0,direction,k,anchor,point_certificate)
    if hi>=threshold:raise ValueError('Point does not disprove the threshold')
    return {'method':'affine_family_counterexample_v9',**meta,'status':'disproved',
            'parameter':str(anchor),'point_certificate':point_certificate}


def verify(cert):
    try:
        args=[cert['q0_coefficients'],cert['direction_coefficients'],cert['parameter_interval'],
              cert['eigenvalue_index'],cert['threshold'],cert['direction_range_proof']]
        if cert['method']=='uniform_affine_family_v9':expected=certificate(*args,cert['leaves'])
        elif cert['method']=='affine_family_counterexample_v9':expected=counterexample(*args,cert['parameter'],cert['point_certificate'])
        else:return False
        return expected==cert
    except (KeyError,TypeError,ValueError,ZeroDivisionError,IndexError):return False


def family(q0,direction,interval=(-1,1),k=1,threshold=F(-3,10),modes=10,max_cells=64,adaptive=True):
    if type(max_cells) is not int or not 1<=max_cells<=256:raise ValueError('Cell budget must be 1..256')
    if type(adaptive) is not bool:raise ValueError('Adaptive flag must be boolean')
    v3.sizes(k,modes);start=perf_counter()
    q0=v3.base.potential([str(x) for x in q0]);direction=v3.base.potential([str(x) for x in direction])
    rp=v5.make_range(direction)
    q0,direction,a,b,threshold,norm,_=setup(q0,direction,interval,k,threshold,rp)
    arguments=[q0,direction,interval,k,threshold,rp]
    pending=deque([(a,b)] if adaptive else [(a+(b-a)*i/max_cells,a+(b-a)*(i+1)/max_cells) for i in range(max_cells)])
    leaves=[];evaluations=0
    while pending:
        left,right=pending.popleft();anchor=(left+right)/2;q=at_parameter(q0,direction,anchor)
        point=v4.guided(q,k,modes,F(1,10**10),range_proof=v5.make_range(q))['certificate'];evaluations+=1
        if F(point['upper'])<threshold:
            result=counterexample(*arguments,anchor,point)
            break
        leaf={'interval':[str(left),str(right)],'point_certificate':point}
        lower=F(point['lower'])-norm*(right-left)/2
        if adaptive and lower<threshold and len(pending)+len(leaves)+2<=max_cells:
            pending.extend([(left,anchor),(anchor,right)])
        else:leaves.append(leaf)
    else:
        leaves.sort(key=lambda x:F(x['interval'][0]));result=certificate(*arguments,leaves)
    if not verify(result):raise ArithmeticError('Parameter-family certificate rejected')
    return {'algorithm_version':9,'certificate':result,'point_evaluations':evaluations,
            'max_cells':max_cells,'adaptive':adaptive,'modes':modes,'seconds':perf_counter()-start}
