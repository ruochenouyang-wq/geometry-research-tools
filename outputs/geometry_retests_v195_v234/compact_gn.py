"""V215--224: compact exact GN trials in powers of y=(1+t)/2.

No dense t expansion is used by the moment, residual or plane kernels.  The
original C16 input is retained and checked separately by bridge_original.
Only trial lower bounds / maxima in a recorded two-function plane are claimed.
"""
from fractions import Fraction as F
from math import comb
import re
import sys
sys.dont_write_bytecode=True
from support import arithmetic as poly, exact_extrema as ex, old_gn, same

MAX_POWER=10**6
MAX_TERMS=64
MAX_PAIRS=262144
COORDINATE='y=(1+t)/2,t=x3_on_unit_S2'
MEASURE='sphere_probability_mu;axisymmetric_dmu=dy'
SCOPE='one_axisymmetric_mean_zero_trial_not_a_full_space_GN_upper_bound'
CONVENTION='probability_ratio=N/(M*E);pi_times_area_K=N/(4*M*E)'


def exact(x):
    if type(x) in (int,F): return F(x)
    if isinstance(x,str) and re.fullmatch(r'-?\d+(?:/[1-9]\d*)?',x): return F(x)
    raise ValueError('Exact integer or rational value required')


def integer(x,low=0,high=MAX_POWER):
    if type(x) is not int or not low<=x<=high: raise ValueError('Integer exponent/budget outside supported range')
    return x


def sparse(raw,max_power=MAX_POWER,max_terms=MAX_TERMS):
    """Canonical sorted [integer exponent, exact coefficient] sparse powers."""
    if not isinstance(raw,(list,tuple)) or not len(raw)<=max_terms: raise ValueError('Sparse term budget exceeded')
    result={}
    for term in raw:
        if not isinstance(term,(list,tuple)) or len(term)!=2: raise ValueError('Require [exponent,coefficient] pairs')
        n,c=integer(term[0],0,max_power),exact(term[1])
        result[n]=result.get(n,F(0))+c
    return tuple((n,c) for n,c in sorted(result.items()) if c)


def wire(p): return [[n,str(c)] for n,c in p]


def _add(*ps):
    out={}
    for p in ps:
        for n,c in p: out[n]=out.get(n,F(0))+c
    return tuple((n,c) for n,c in sorted(out.items()) if c)


def _scale(p,c): return tuple((n,c*a) for n,a in p if c*a)


def _mul(p,q,counter=None):
    if len(p)*len(q)>MAX_PAIRS: raise ValueError('Convolution pair budget exceeded')
    if counter is not None: counter['coefficient_products']=counter.get('coefficient_products',0)+len(p)*len(q)
    result={}
    for i,a in p:
        for j,b in q: result[i+j]=result.get(i+j,F(0))+a*b
    return tuple((n,c) for n,c in sorted(result.items()) if c)


def _mean(p): return sum((c/F(n+1) for n,c in p),F(0))


def _inner(p,q): return _mean(_mul(p,q))


def _lap(p):
    # -d/dy[y(1-y) d/dy] y^n = n(n+1)y^n - n^2 y^(n-1).
    return _add(tuple((n,n*(n+1)*c) for n,c in p if n),
                tuple((n-1,-n*n*c) for n,c in p if n))


def mean_certificate(raw):
    """V215: exact sparse-power integration using integral_0^1 y^n=1/(n+1)."""
    p=sparse(raw)
    return {'format':'compact_mean_v215','coordinate':COORDINATE,'measure':MEASURE,
        'terms':wire(p),'integrated_terms':[str(c/F(n+1)) for n,c in p],
        'mean':str(_mean(p)),'kernel_evaluations':len(p)}


def centered_function(weights,scale=1):
    """V216: sum_a c_a*(y^a-1/(a+1)), with an explicit nonzero scale."""
    weights=sparse(weights);scale=exact(scale)
    if scale==0: raise ValueError('Function representative scale must be nonzero')
    p=_scale(_add(weights,((0,-_mean(weights)),)),scale)
    if len(p)>MAX_TERMS: raise ValueError('Centering exceeds the sparse function term budget')
    return {'format':'compact_centered_v216','coordinate':COORDINATE,'measure':MEASURE,
        'weights':wire(weights),'scale':str(scale),'terms':wire(p),
        'mean':str(_mean(p)),'zero_function':not bool(p)}


def mass_kernel(a,b):
    """V217: covariance of centered powers h_a and h_b, including h_0=0."""
    a,b=integer(a),integer(b)
    return F(1,a+b+1)-F(1,(a+1)*(b+1))


def mass_form(weights):
    weights=sparse(weights)
    gram=[[mass_kernel(a,b) for b,_ in weights] for a,_ in weights]
    value=sum((ca*cb*gram[i][j] for i,(_,ca) in enumerate(weights) for j,(_,cb) in enumerate(weights)),F(0))
    return {'format':'compact_mass_form_v217','weights':wire(weights),
        'gram':[[str(x) for x in row] for row in gram],'value':str(value),
        'function':centered_function(weights)}


def energy_kernel(a,b):
    """V218: ab*integral_0^1 y^(a+b-1)(1-y)dy, zero for constants."""
    a,b=integer(a),integer(b)
    return F(a*b,(a+b)*(a+b+1)) if a and b else F(0)


def energy_form(weights):
    weights=sparse(weights)
    gram=[[energy_kernel(a,b) for b,_ in weights] for a,_ in weights]
    value=sum((ca*cb*gram[i][j] for i,(_,ca) in enumerate(weights) for j,(_,cb) in enumerate(weights)),F(0))
    return {'format':'compact_energy_form_v218','weights':wire(weights),
        'gram':[[str(x) for x in row] for row in gram],'value':str(value),
        'function':centered_function(weights)}


def quartic_moment(raw):
    """V219: exponent-sum convolution, never a dense coefficient expansion."""
    p=sparse(raw);count={};square=_mul(p,p,count);fourth=_mul(square,square,count)
    return {'format':'compact_quartic_v219','coordinate':COORDINATE,'measure':MEASURE,
        'terms':wire(p),'square':wire(square),'fourth_power':wire(fourth),
        'quartic':str(_mean(fourth)),'operation_counts':count,
        'support_sizes':{'input':len(p),'square':len(square),'fourth':len(fourth)},
        'dense_t_polynomial_expanded':False}


def moments(raw):
    p=sparse(raw);four=quartic_moment(p)
    M=_mean(tuple((n,F(c)) for n,c in four['square']))
    energy=sum((ca*cb*energy_kernel(a,b) for a,ca in p for b,cb in p),F(0))
    return {'format':'compact_moments_v220','coordinate':COORDINATE,'measure':MEASURE,
        'terms':wire(p),'mean':str(_mean(p)),'mass':str(M),'energy':str(energy),
        'quartic':four['quartic'],'quartic_certificate':four,
        'operation_counts':{'convolution_coefficient_products':four['operation_counts'].get('coefficient_products',0),
                            'energy_kernel_evaluations':len(p)**2,
                            'mass_mean_terms':len(four['square']),'quartic_mean_terms':len(four['fourth_power'])}}


def ratio(raw):
    """V220: exact probability ratio and pi times the natural-area quotient."""
    p=sparse(raw);m=moments(p);M,E,N=[F(m[k]) for k in ('mass','energy','quartic')]
    if _mean(p)!=0 or M<=0 or E<=0: raise ValueError('Require a nonzero mean-zero positive-energy function')
    return {'format':'compact_ratio_v220','terms':wire(p),'moments':m,
        'probability_ratio':str(N/(M*E)),'pi_times_area_K':str(N/(4*M*E)),
        'convention':CONVENTION,'scope':SCOPE}


def residual(raw):
    """V221: complete EL residual in sparse powers; highest degree is retained."""
    p=sparse(raw);rc=ratio(p);M,E,N=[F(rc['moments'][k]) for k in ('mass','energy','quartic')]
    cubic=_mul(_mul(p,p),p)
    before=_add(_lap(p),_scale(p,E/M),_scale(cubic,-2*E/N))
    r=_add(before,((0,-_mean(before)),))
    return {'format':'compact_residual_v221','terms':wire(p),'ratio':rc,
        'laplacian':wire(_lap(p)),'cubic':wire(cubic),'residual':wire(r),
        'removed_mean':str(_mean(before)),'residual_mean':str(_mean(r)),
        'inner_trial':str(_inner(r,p)),'residual_norm_squared':str(_inner(r,r)),
        'full_degree':max((n for n,_ in r),default=0),'residual_support_size':len(r),
        'full_stationary':not bool(r),'degree_truncation_used':False}


def direction(raw,v=None):
    p=sparse(raw);rc=residual(p);r=sparse(rc['residual'],3*MAX_POWER,MAX_PAIRS)
    v=_scale(r,-1) if v is None else sparse(v,3*MAX_POWER,MAX_PAIRS)
    if _mean(v)!=0: raise ValueError('Direction must preserve the mean-zero constraint')
    M,E,N=[F(rc['ratio']['moments'][k]) for k in ('mass','energy','quartic')]
    J=N/(4*M*E);dm=2*_inner(p,v);de=2*_inner(_lap(p),v);dn=4*_inner(_mul(_mul(p,p),p),v)
    derivative=J*(dn/N-dm/M-de/E);identity=-2*J*_inner(r,v)/E
    if derivative!=identity: raise ArithmeticError('Sparse EL derivative identity failed')
    return {'format':'compact_direction_v221','terms':wire(p),'direction':wire(v),
        'residual_certificate':rc,'derivative_pi_times_K':str(derivative),
        'identity_derivative':str(identity),'strict_local_improvement':derivative>0,
        'finite_step_improvement_claimed':False}


def _basis(raw):
    if not isinstance(raw,(list,tuple)) or len(raw)!=2: raise ValueError('Two sparse functions required')
    f,g=map(sparse,raw)
    if len({n for n,_ in f}|{n for n,_ in g})>MAX_TERMS:
        raise ValueError('The union of plane supports exceeds the sparse witness term budget')
    ratio(f);ratio(g)
    if _inner(f,f)*_inner(g,g)-_inner(f,g)**2<=0: raise ValueError('Plane basis is linearly dependent')
    return f,g


def chart_moments(f,g):
    M=[_inner(f,f),2*_inner(f,g),_inner(g,g)]
    E=[_inner(f,_lap(f)),2*_inner(f,_lap(g)),_inner(g,_lap(g))]
    N=[]
    for j in range(5):
        p=((0,F(1)),)
        for _ in range(4-j): p=_mul(p,f)
        for _ in range(j): p=_mul(p,g)
        N.append(comb(4,j)*_mean(p))
    return {'mass':poly.trim(M),'energy':poly.trim(E),'numerator':poly.trim(N),
            'denominator':poly.scale(poly.mul(M,E),4)}


def _sign(p,strict=False):
    factor=1/max(map(abs,p)) if any(p) else F(1)
    proof=ex.maximize(poly.scale(p,factor),tolerance=F(1,10**20),max_nodes=128)
    b=F(proof['maximum_upper'])
    return {'positive_multiplier':str(factor),'proof':proof}, b<0 if strict else b<=0


def _verify_sign(c,p,strict=False):
    factor=exact(c['positive_multiplier']);proof=c['proof']
    if factor<=0 or proof['polynomial']!=list(map(str,poly.scale(p,factor))) or proof['interval']!=['-1','1'] or not ex.verify(proof): return False
    b=F(proof['maximum_upper'])
    return b<0 if strict else b<=0


def plane(raw_basis,tolerance=F(1,10**8)):
    """V222: global two-function projective optimization, scalar degree <=4.

    Trial powers may be 64 or higher; the Sturm problem is only in the chart
    parameter. Frozen proposal code receives these four-degree moment arrays,
    never an old degree-limited trial-function adapter.
    """
    original_basis=_basis(raw_basis);tolerance=exact(tolerance)
    if not F(1,10**18)<=tolerance<=1: raise ValueError('Tolerance outside [1e-18,1]')
    # The original C16 functions have factors 2^n. Equalizing their sparse
    # coefficient magnitudes is an invertible positive basis transformation,
    # recorded below and replayed exactly; it does not change the plane.
    basis_scales=[1/max(abs(c) for _,c in p) for p in original_basis]
    basis=tuple(_scale(p,c) for p,c in zip(original_basis,basis_scales))
    data=[chart_moments(*basis),chart_moments(*reversed(basis))]
    proposals=[old_gn._propose(m) for m in data]
    which=max(range(2),key=lambda j:proposals[j][1]);parameter,lower=proposals[which]
    upper=lower+tolerance;attempts=[]
    for _ in range(8):
        proofs=[_sign(poly.add(m['numerator'],poly.scale(m['denominator'],-upper))) for m in data]
        attempts.append({'upper':str(upper),'accepted_charts':[x[1] for x in proofs]})
        if all(x[1] for x in proofs): break
        upper=lower+10*(upper-lower)
    else: raise ArithmeticError('No verified plane ceiling within the search budget')
    charts=[]
    for j,m in enumerate(data):
        mp,mok=_sign(poly.scale(m['mass'],-1),True);ep,eok=_sign(poly.scale(m['energy'],-1),True)
        if not mok or not eok: raise ArithmeticError('Plane positivity proof failed')
        charts.append({'order':[j,1-j],'moments':{k:list(map(str,v)) for k,v in m.items()},
            'mass_positive':mp,'energy_positive':ep,'upper_nonpositive':proofs[j][0]})
    witness=_add(basis[which],_scale(basis[1-which],parameter))
    original_coefficients=[F(0),F(0)]
    original_coefficients[which]=basis_scales[which]
    original_coefficients[1-which]=parameter*basis_scales[1-which]
    certificate={'format':'compact_plane_v222','basis':[wire(p) for p in original_basis],
        'positive_basis_scales':list(map(str,basis_scales)),'chart_basis':[wire(p) for p in basis],
        'scope':'global_in_recorded_two_function_plane_only','convention':CONVENTION,
        'lower_pi_times_K':str(lower),'upper_pi_times_K':str(upper),'gap':str(upper-lower),
        'requested_gap':str(tolerance),'target_met':upper-lower<=tolerance,
        'charts':charts,'witness':{'chart':which,'parameter':str(parameter),'terms':wire(witness),
                                  'original_basis_coefficients':list(map(str,original_coefficients))},
        'ceiling_attempts':attempts,'parameter_polynomial_degree_at_most':4,
        'function_degree_limit_does_not_come_from_Sturm_degree_limit':True,
        'universal_GN_upper_bound':False}
    if not verify_plane(certificate,original_basis):
        raise ArithmeticError('Generated compact plane certificate failed independent replay')
    return certificate


def verify_plane(c,expected_basis=None):
    try:
        if c['format']!='compact_plane_v222' or c['scope']!='global_in_recorded_two_function_plane_only' or c['convention']!=CONVENTION or c['universal_GN_upper_bound'] is not False or c['parameter_polynomial_degree_at_most']!=4 or c['function_degree_limit_does_not_come_from_Sturm_degree_limit'] is not True: return False
        original_basis=_basis(c['basis'])
        if expected_basis is not None and original_basis!=_basis(expected_basis): return False
        basis_scales=[1/max(abs(v) for _,v in p) for p in original_basis]
        basis=tuple(_scale(p,v) for p,v in zip(original_basis,basis_scales))
        if c['positive_basis_scales']!=list(map(str,basis_scales)) or c['chart_basis']!=[wire(p) for p in basis]: return False
        lower,upper,target=map(exact,[c['lower_pi_times_K'],c['upper_pi_times_K'],c['requested_gap']])
        if not 0<lower<=upper or F(c['gap'])!=upper-lower or not F(1,10**18)<=target<=1 or type(c['target_met']) is not bool or c['target_met']!=(upper-lower<=target) or len(c['charts'])!=2: return False
        for j,ch in enumerate(c['charts']):
            m=chart_moments(basis[j],basis[1-j])
            if ch['order']!=[j,1-j] or ch['moments']!={k:list(map(str,v)) for k,v in m.items()}: return False
            if not _verify_sign(ch['mass_positive'],poly.scale(m['mass'],-1),True) or not _verify_sign(ch['energy_positive'],poly.scale(m['energy'],-1),True) or not _verify_sign(ch['upper_nonpositive'],poly.add(m['numerator'],poly.scale(m['denominator'],-upper))): return False
        w=c['witness'];j=w['chart'];parameter=exact(w['parameter'])
        if type(j) is not int or j not in (0,1) or not -1<=parameter<=1: return False
        u=_add(basis[j],_scale(basis[1-j],parameter))
        coefficients=[F(0),F(0)];coefficients[j]=basis_scales[j];coefficients[1-j]=parameter*basis_scales[1-j]
        return w['terms']==wire(u) and w['original_basis_coefficients']==list(map(str,coefficients)) and F(ratio(u)['pi_times_area_K'])==lower
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError): return False


J_NUM=[0,9,15,0,12]
J_DEN=[2,18,54,62,24]


def _shift(p,h=1):
    return [sum((F(p[j])*comb(j,i)*h**(j-i) for j in range(i,len(p))),F(0)) for i in range(len(p))]


def family_law():
    """V223: exact closed rational family, strict integer monotonicity and limit.

    J(n)=3n(2n+1)(2n²-n+3)/[2(4n+1)(3n+1)(n+1)²].
    Positivity of the n=x+1 shifted difference polynomial proves J(n+1)>J(n)
    for every integer n>=1. Leading coefficients give the exact limit 1/2.
    """
    def prod(*ps):
        out=[F(1)]
        for p in ps: out=poly.mul(out,p)
        return out
    A,B,C,D=[1,1],[1,2],[1,3],[1,4]
    N_den=prod(D,C,B,A,A,A,A)
    N_num=poly.add(poly.add(prod(C,B,A,A,A,A),poly.scale(prod(D,B,A,A,A),-4)),
                   poly.add(poly.scale(prod(D,C,A,A),6),poly.scale(prod(D,C,B),-3)))
    # Identity J=N/(4ME), M=n²/[(2n+1)(n+1)²], E=n/[2(2n+1)].
    lhs=prod(N_num,B,B,A,A,J_DEN)
    rhs=prod(J_NUM,[0,0,0,2],N_den)
    if poly.trim(lhs)!=poly.trim(rhs): raise ArithmeticError('Closed family identity failed')
    difference=poly.add(poly.mul(_shift(J_NUM),J_DEN),poly.scale(poly.mul(J_NUM,_shift(J_DEN)),-1))
    shifted=_shift(difference)
    gap=poly.add(J_DEN,poly.scale(J_NUM,-2))
    if not all(c>=0 for c in shifted) or shifted[0]<=0 or not all(c>=0 for c in gap) or gap[0]<=0:
        raise ArithmeticError('Coefficient positivity proof failed')
    return {'format':'compact_family_law_v223','n_domain':'all_integers_n>=1',
        'function':'h_n=y^n-1/(n+1)','ratio_convention':CONVENTION,
        'pi_K_numerator':list(map(str,J_NUM)),'pi_K_denominator':list(map(str,J_DEN)),
        'quartic_numerator':list(map(str,N_num)),'quartic_denominator':list(map(str,N_den)),
        'algebraic_identity_replayed':True,'difference_numerator':list(map(str,difference)),
        'difference_after_n_equals_x_plus_1':list(map(str,shifted)),
        'one_half_gap_twice_numerator':list(map(str,gap)),
        'denominator_positive_for_n_ge_1':True,'strictly_increasing_in_integer_n':True,
        'limit_pi_times_area_K':'1/2','limit_probability_ratio':'2',
        'sharp_constant_claimed':False}


def family_value(n):
    integer(n,1,10**12)
    law=family_law()
    M=F(n*n,(2*n+1)*(n+1)**2);E=F(n,2*(2*n+1))
    N=sum((F(comb(4,j))*(-F(1,n+1))**(4-j)/F(n*j+1) for j in range(5)),F(0))
    J=poly.evaluate(list(map(F,J_NUM)),n)/poly.evaluate(list(map(F,J_DEN)),n)
    if J!=N/(4*M*E): raise ArithmeticError('Closed moment ratio failed')
    return {'format':'compact_family_value_v223','n':n,'law':law,
        'normalized_mean':'0','normalized_mass':str(M),'normalized_energy':str(E),
        'normalized_quartic':str(N),'pi_times_area_K':str(J),'probability_ratio':str(4*J),
        'trial_polynomial_expanded':False,'full_space_sharpness_claimed':False}


def original_cap(n):
    """The original C16 coefficients in t, retained exactly by the bridge."""
    integer(n,1,1024)
    result=[F(comb(n,j)) for j in range(n+1)]
    result[0]-=F(2**n,n+1)
    return result


def evaluate_cap(n,original=True):
    """V224: compact cap evaluator with original 2^n scale fully bound."""
    integer(n,1,1024 if original else MAX_POWER)
    if type(original) is not bool: raise ValueError('original must be boolean')
    scale=F(2**n) if original else F(1)
    function=centered_function([[n,1]],scale)
    rc=ratio(function['terms']);law=family_value(n)
    return {'format':'compact_cap_v224','n':n,'original_scale':original,
        'input_function':'(1+t)^n-2^n/(n+1)' if original else '((1+t)/2)^n-1/(n+1)',
        'function':function,'ratio':rc,'family_value':law,
        'moment_scale_factors':{'mass':str(scale**2),'energy':str(scale**2),'quartic':str(scale**4)},
        'dense_original_expansion_needed_for_moment_evaluation':False,
        'scope':SCOPE}


def bridge_original(raw_t_coefficients,n):
    """V224: check every original monomial coefficient, then evaluate compactly."""
    integer(n,1,1024)
    if not isinstance(raw_t_coefficients,(list,tuple)) or len(raw_t_coefficients)!=n+1:
        raise ValueError('Original polynomial must have exactly n+1 coefficients')
    p=[exact(x) for x in raw_t_coefficients]
    if p!=original_cap(n): raise ValueError('Input is not the exact original C16 polynomial for this n')
    return {'format':'compact_original_bridge_v224','n':n,'original_t_coefficients':list(map(str,p)),
        'coefficient_identity':'binomial expansion of (1+t)^n, then subtract 2^n/(n+1)',
        'verified_coefficient_count':n+1,'evaluation':evaluate_cap(n,True),
        'normalization_substitution_without_scale_proof':False}


def verify(c,expected_terms=None,expected_n=None,expected_original_coefficients=None):
    """Strict replay, not a search rerun; also binds original C16 input."""
    try:
        fmt=c['format']
        if expected_terms is not None:
            bound=c.get('terms')
            if bound is None and 'basis' in c: bound=c['basis'][0]
            if bound is None and 'function' in c: bound=c['function']['terms']
            if bound is None and 'evaluation' in c: bound=c['evaluation']['function']['terms']
            if bound is None or sparse(bound)!=sparse(expected_terms): return False
        if expected_n is not None and (type(expected_n) is not int or c.get('n')!=expected_n): return False
        if expected_original_coefficients is not None:
            if fmt!='compact_original_bridge_v224' or c['original_t_coefficients']!=list(map(str,[exact(x) for x in expected_original_coefficients])): return False
        if fmt=='compact_plane_v222': return verify_plane(c)
        makers={'compact_mean_v215':lambda:mean_certificate(c['terms']),
            'compact_centered_v216':lambda:centered_function(c['weights'],c['scale']),
            'compact_mass_form_v217':lambda:mass_form(c['weights']),
            'compact_energy_form_v218':lambda:energy_form(c['weights']),
            'compact_quartic_v219':lambda:quartic_moment(c['terms']),
            'compact_moments_v220':lambda:moments(c['terms']),
            'compact_ratio_v220':lambda:ratio(c['terms']),
            'compact_residual_v221':lambda:residual(c['terms']),
            'compact_direction_v221':lambda:direction(c['terms'],c['direction']),
            'compact_family_law_v223':family_law,
            'compact_family_value_v223':lambda:family_value(c['n']),
            'compact_cap_v224':lambda:evaluate_cap(c['n'],c['original_scale']),
            'compact_original_bridge_v224':lambda:bridge_original(c['original_t_coefficients'],c['n'])}
        return fmt in makers and same(c,makers[fmt]())
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError,OverflowError): return False


STAGES=[
 {'version':215,'capability':'Sparse powers of y and exact probability mean kernel','callable':'compact_gn.mean_certificate','evidence':['C16/kernel_evidence.json'],'outcome':'implemented_and_evidence_replayed'},
 {'version':216,'capability':'Centered sparse basis and explicit nonzero scale representative','callable':'compact_gn.centered_function','evidence':['C16/kernel_evidence.json'],'outcome':'implemented_and_evidence_replayed'},
 {'version':217,'capability':'Closed pair mass Gram kernel for centered powers','callable':'compact_gn.mass_form','evidence':['C16/kernel_evidence.json'],'outcome':'implemented_and_evidence_replayed'},
 {'version':218,'capability':'Closed exact sphere Dirichlet pair kernel including constant handling','callable':'compact_gn.energy_form','evidence':['C16/kernel_evidence.json'],'outcome':'implemented_and_evidence_replayed'},
 {'version':219,'capability':'Sparse exponent-sum quartic convolution and exact moments','callable':'compact_gn.quartic_moment','evidence':['C16/kernel_evidence.json','C16/benchmark.json'],'outcome':'implemented_and_evidence_replayed'},
 {'version':220,'capability':'Probability GN ratio and original 2^n moment scaling','callable':'compact_gn.ratio','evidence':['C16/original_cases.json'],'outcome':'implemented_and_evidence_replayed'},
 {'version':221,'capability':'Complete sparse Euler-Lagrange residual and exact improvement derivative','callable':'compact_gn.residual','evidence':['C16/residuals.json'],'outcome':'implemented_and_evidence_replayed'},
 {'version':222,'capability':'Global two-chart projective-plane maximum with compressed moment kernels','callable':'compact_gn.plane','evidence':['C16/planes.json'],'outcome':'implemented_and_evidence_replayed'},
 {'version':223,'capability':'Symbolic rational family, strict integer monotonicity and exact asymptotic limit','callable':'compact_gn.family_law','evidence':['C16/family.json'],'outcome':'implemented_and_evidence_replayed'},
 {'version':224,'capability':'Unified evaluator and verified original monomial bridge for n4/16/64','callable':'compact_gn.bridge_original','evidence':['C16/original_cases.json','C16/baseline.json'],'outcome':'implemented_and_evidence_replayed'},
]
