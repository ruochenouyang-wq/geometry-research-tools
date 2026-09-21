"""Exact source-function remainder evidence, independent of spectral transfer.

Analytic sums have uniform Taylor error bounds on unit S2. Axis profiles have
exact probability-L2 orthogonal projection errors. No sampling proves a bound.
"""
from math import factorial, isqrt
import re
from dependencies import F, same, canonical, arithmetic

ANALYTIC_FORMAT='sphere_analytic_function_remainder_v1'
L2_FORMAT='sphere_axis_profile_projection_v1'
SCOPE='unit_S2_probability_measure'
MAX_DEGREE=24
MAX_MONOMIALS=4096
MAX_PAIRS=262144
MAX_EXP_ENVELOPE=4096
ZERO=(0,0,0)


def rational(value):
    if type(value) in (int,F): return F(value)
    if isinstance(value,str) and re.fullmatch(r'-?\d+(?:/[1-9]\d*)?',value): return F(value)
    raise ValueError('Use an exact integer or rational string, never float or boolean')


def _integer(value,lower,upper,name):
    if type(value) is not int or not lower<=value<=upper:
        raise ValueError(name+' must be an integer in '+str(lower)+'..'+str(upper))
    return value


def _fields(value,allowed,name):
    if not isinstance(value,dict) or not set(value)<=set(allowed):
        raise ValueError('Unknown fields or malformed '+name)


def polynomial(raw):
    """Degree24 xyz rational polynomial without the old coefficient-size cap."""
    if not isinstance(raw,dict) or len(raw)>MAX_MONOMIALS:
        raise ValueError('Expected a bounded xyz exponent-to-rational dictionary')
    result={}
    for key,value in raw.items():
        if not isinstance(key,str) or re.fullmatch(r'\d+,\d+,\d+',key) is None:
            raise ValueError('Polynomial keys must be strings i,j,k')
        exponent=tuple(map(int,key.split(',')))
        if sum(exponent)>MAX_DEGREE: raise ValueError('Source polynomial degree exceeds24')
        if exponent in result: raise ValueError('Duplicate normalized xyz exponent')
        coefficient=rational(value)
        result[exponent]=coefficient
    return {e:c for e,c in result.items() if c}


def encode(p): return {','.join(map(str,e)):str(c) for e,c in sorted(p.items()) if c}


def _add(*ps):
    out={}
    for p in ps:
        for e,c in p.items(): out[e]=out.get(e,F(0))+c
    return {e:c for e,c in out.items() if c}


def _scale(p,c): return {e:c*a for e,a in p.items() if c*a}


def _multiply(p,q):
    if len(p)*len(q)>MAX_PAIRS: raise ValueError('Polynomial multiplication pair budget exceeded')
    out={}
    for e,a in p.items():
        for f,b in q.items():
            exponent=tuple(x+y for x,y in zip(e,f))
            if sum(exponent)>MAX_DEGREE: raise ValueError('Expanded approximation degree exceeds24')
            out[exponent]=out.get(exponent,F(0))+a*b
    return {e:c for e,c in out.items() if c}


def normalize_function(function):
    """Canonical exact input with documented defaults and unknown-field rejection.

    analytic_sum defaults: polynomial={}, terms=[]; coefficient=1, argument={}.
    axis_profile defaults: axis=2, amplitude=1, offset=0, abs_power exponent=1.
    The kind and each nonempty term's function / profile remain mandatory.
    """
    if not isinstance(function,dict): raise ValueError('Function model must be an object')
    kind=function.get('kind')
    if kind=='analytic_sum':
        _fields(function,('kind','polynomial','terms'),'analytic sum')
        p=encode(polynomial(function.get('polynomial',{})))
        raw=function.get('terms',[])
        if not isinstance(raw,list) or len(raw)>12: raise ValueError('At most12 analytic terms are supported')
        grouped={}
        for term in raw:
            _fields(term,('function','argument','coefficient'),'analytic term')
            name=term.get('function')
            if name not in ('exp','sin','log1p'): raise ValueError('Unsupported analytic function')
            argument_polynomial=polynomial(term.get('argument',{}))
            # Validate the original expression before coefficients can cancel:
            # even a zero coefficient does not define an undefined logarithm.
            if name=='log1p' and sum(map(abs,argument_polynomial.values()),F(0))>=1:
                raise ValueError('log1p source requires the argument L1 coefficient radius<1 before term combination')
            argument=encode(argument_polynomial)
            coefficient=rational(term.get('coefficient',1))
            key=(name,canonical(argument))
            old=grouped.get(key,(argument,F(0)))
            grouped[key]=(argument,old[1]+coefficient)
        terms=[{'function':name,'argument':argument,'coefficient':str(coefficient)}
            for (name,_),(argument,coefficient) in sorted(grouped.items()) if coefficient]
        return {'kind':kind,'polynomial':p,'terms':terms}
    if kind=='axis_profile':
        profile=function.get('profile')
        if profile not in ('step','abs_power'): raise ValueError('Unsupported axis profile')
        allowed=['kind','axis','profile','amplitude','offset']+(['exponent'] if profile=='abs_power' else [])
        _fields(function,allowed,'axis profile')
        axis=_integer(function.get('axis',2),0,2,'axis')
        result={'kind':kind,'axis':axis,'profile':profile,
            'amplitude':str(rational(function.get('amplitude',1))),
            'offset':str(rational(function.get('offset',0)))}
        if profile=='abs_power':
            exponent=rational(function.get('exponent',1))
            if exponent<=F(-1,2): raise ValueError('abs_power requires exponent>-1/2 for probability L2')
            result['exponent']=str(exponent)
        return result
    raise ValueError('Unsupported function kind')


def sqrt_upper(value_squared,bits=40):
    """Outward dyadic upper bound on sqrt of an exact nonnegative rational."""
    value=rational(value_squared);_integer(bits,0,256,'sqrt_bits')
    if value<0: raise ValueError('Squared quantity must be nonnegative')
    denominator=1<<bits
    floor_square=value.numerator*denominator**2//value.denominator
    numerator=isqrt(floor_square)
    if F(numerator**2,denominator**2)<value: numerator+=1
    return F(numerator,denominator)


def _series_coefficient(name,k):
    if name=='exp': return F(1,factorial(k))
    if name=='sin': return F((-1)**((k-1)//2),factorial(k)) if k%2 else F(0)
    return F((-1)**(k+1),k) if k else F(0)


def _analytic_term(term,order):
    name=term['function'];argument=polynomial(term['argument'])
    radius=sum(map(abs,argument.values()),F(0))
    if name=='log1p' and radius>=1:
        raise ValueError('log1p Taylor evidence requires the argument L1 coefficient radius<1')
    coefficients=[_series_coefficient(name,k) for k in range(order+1)]
    last_nonzero=max((k for k,c in enumerate(coefficients) if c),default=0)
    power={ZERO:F(1)};approximation={};series=list(map(str,coefficients))
    for k in range(last_nonzero+1):
        coefficient=coefficients[k]
        approximation=_add(approximation,_scale(power,coefficient))
        if k<last_nonzero: power=_multiply(power,argument)
    omitted=order+1
    if name=='exp':
        ratio=radius/F(order+2)
        if ratio<1:
            error=radius**omitted/F(factorial(omitted))/(1-ratio)
            tail={'method':'absolute_exponential_series_geometric_tail',
                  'first_omitted_degree':omitted,'first_omitted_absolute_term':str(radius**omitted/F(factorial(omitted))),
                  'successive_absolute_term_ratio_upper':str(ratio),'geometric_denominator':str(1-ratio)}
        else:
            ceiling=(radius.numerator+radius.denominator-1)//radius.denominator
            if ceiling>MAX_EXP_ENVELOPE: raise ValueError('Exponential envelope integer budget exceeds4096')
            envelope=F(3**ceiling)
            error=envelope*radius**omitted/F(factorial(omitted))
            tail={'method':'exponential_lagrange_with_e_less_than3',
                  'first_omitted_degree':omitted,'radius_ceiling':ceiling,
                  'exponential_envelope':str(envelope),
                  'analytic_fact':'exp(radius)<=exp(ceil(radius))<=3^ceil(radius)'}
    elif name=='sin':
        error=radius**omitted/F(factorial(omitted))
        tail={'method':'sine_lagrange_derivative_bound1','first_omitted_degree':omitted,
              'derivative_absolute_upper':'1'}
    else:
        error=radius**omitted/F(omitted)/(1-radius)
        tail={'method':'absolute_log1p_series_geometric_tail','first_omitted_degree':omitted,
              'first_omitted_absolute_term':str(radius**omitted/F(omitted)),
              'successive_geometric_ratio_upper':str(radius),'geometric_denominator':str(1-radius)}
    coefficient=rational(term['coefficient'])
    return {'source_term':term,'argument_radius':str(radius),
        'range_rule':'abs(x_i)<=1 on unit_S2 implies abs(argument)<=sum_abs_coefficients',
        'series_coefficients':series,'unscaled_polynomial':encode(approximation),
        'weighted_polynomial':encode(_scale(approximation,coefficient)),
        'tail':tail,'unweighted_error_upper':str(error),
        'weighted_error_upper':str(abs(coefficient)*error)}


def analytic_model(function,order=4,sqrt_bits=40):
    """Polynomial plus certified complete uniform remainder for analytic sums."""
    function=normalize_function(function)
    if function['kind']!='analytic_sum': raise ValueError('analytic_model requires analytic_sum')
    _integer(order,0,24,'order');_integer(sqrt_bits,0,256,'sqrt_bits')
    term_proofs=[_analytic_term(term,order) for term in function['terms']]
    result=polynomial(function['polynomial'])
    for proof in term_proofs: result=_add(result,polynomial(proof['weighted_polynomial']))
    error=sum((F(proof['weighted_error_upper']) for proof in term_proofs),F(0))
    return {'format':ANALYTIC_FORMAT,'function':function,'polynomial':encode(result),
        'norm':'Linf','error_upper':str(error),'scope':SCOPE,'order':order,
        'sqrt_bits':sqrt_bits,'sqrt_bits_used':False,'term_proofs':term_proofs,
        'combination_rule':'triangle_inequality_with_absolute_signed_coefficients',
        'bound_statement':'sup_unit_S2_abs(source_function-polynomial)<=error_upper',
        'spectral_transfer_claimed':False,'formal_assistant_checked':False}


def _legendre(degree):
    basis=[[F(1)]]
    if degree: basis.append([F(0),F(1)])
    for l in range(1,degree):
        basis.append(arithmetic.scale(arithmetic.add(
            arithmetic.scale(arithmetic.mul([F(0),F(1)],basis[-1]),2*l+1),
            arithmetic.scale(basis[-2],-l)),F(1,l+1)))
    return basis


def _profile_moment(function,k):
    if function['profile']=='step': return F(1,2*(k+1))
    return F(1)/(F(function['exponent'])+k+1) if k%2==0 else F(0)


def _source_moment(function,k):
    amplitude,offset=F(function['amplitude']),F(function['offset'])
    constant_moment=F(1,k+1) if k%2==0 else F(0)
    return offset*constant_moment+amplitude*_profile_moment(function,k)


def l2_model(function,degree=8,sqrt_bits=40):
    """Exact L2 projection error of a step or integrable absolute-power profile."""
    function=normalize_function(function)
    if function['kind']!='axis_profile': raise ValueError('l2_model requires axis_profile')
    _integer(degree,0,24,'degree');_integer(sqrt_bits,0,256,'sqrt_bits')
    amplitude,offset=F(function['amplitude']),F(function['offset'])
    profile_mean=_profile_moment(function,0)
    profile_squared=F(1,2) if function['profile']=='step' else F(1)/(2*F(function['exponent'])+1)
    target_squared=offset**2+2*offset*amplitude*profile_mean+amplitude**2*profile_squared
    monomial_moments=[_source_moment(function,k) for k in range(degree+1)]
    basis=_legendre(degree);coefficients=[];entries=[];result=[F(0)]
    projection_squared=F(0)
    for l,p in enumerate(basis):
        moment=sum((c*monomial_moments[k] for k,c in enumerate(p)),F(0))
        coefficient=(2*l+1)*moment
        coefficients.append(coefficient)
        result=arithmetic.add(result,arithmetic.scale(p,coefficient))
        contribution=(2*l+1)*moment**2;projection_squared+=contribution
        entries.append({'degree':l,'legendre_coefficients':list(map(str,p)),
            'basis_mass':str(F(1,2*l+1)),'source_inner_product':str(moment),
            'projection_coefficient':str(coefficient),'norm_squared_contribution':str(contribution)})
    remainder=target_squared-projection_squared
    if remainder<0: raise ArithmeticError('Negative exact projection error')
    upper=sqrt_upper(remainder,sqrt_bits)
    axis=function['axis'];xyz={}
    for k,c in enumerate(result):
        if c:
            exponent=[0,0,0];exponent[axis]=k;xyz[tuple(exponent)]=c
    profile_at_zero=F(1) if function['profile']=='abs_power' and F(function['exponent'])==0 else F(0)
    return {'format':L2_FORMAT,'function':function,'polynomial':encode(xyz),
        'norm':'L2_probability','error_upper':str(upper),'error_squared':str(remainder),
        'scope':SCOPE,'degree':degree,'sqrt_bits':sqrt_bits,
        'target_norm_squared':str(target_squared),'projection_norm_squared':str(projection_squared),
        'monomial_source_moments':list(map(str,monomial_moments)),
        'legendre_projection':entries,'axis_power_coefficients':list(map(str,result)),
        'point_value_at_axis_zero':str(offset+amplitude*profile_at_zero),
        'point_convention':'step(0)=0;negative_abs_power_at0_is_set_to0_as_an_L2_equivalent_representative',
        'projection_identity':'exact_error_squared=target_norm_squared-sum_l(2l+1)*source_inner_product_l^2',
        'sqrt_upper_squared':str(upper**2),
        'bound_statement':'probability_L2_norm(source_function-polynomial)<=error_upper',
        'spectral_transfer_claimed':False,'formal_assistant_checked':False}


def verify_model(certificate,expected_function=None,expected_polynomial=None):
    """Replay the complete source/remainder algebra, without sampling or spectra."""
    try:
        if not isinstance(certificate,dict): return False
        function=normalize_function(certificate['function'])
        if expected_function is not None and function!=normalize_function(expected_function): return False
        if expected_polynomial is not None and polynomial(certificate['polynomial'])!=polynomial(expected_polynomial): return False
        if certificate.get('format')==ANALYTIC_FORMAT:
            rebuilt=analytic_model(function,certificate['order'],certificate['sqrt_bits'])
        elif certificate.get('format')==L2_FORMAT:
            rebuilt=l2_model(function,certificate['degree'],certificate['sqrt_bits'])
        else: return False
        return same(certificate,rebuilt)
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError,OverflowError,AttributeError): return False
