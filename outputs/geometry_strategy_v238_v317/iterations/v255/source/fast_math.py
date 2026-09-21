"""Exact optimizations with immutable caches and frozen public verification."""
from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from functools import lru_cache
from math import comb
import json
from time import perf_counter

from backend import F, canonical, digest, pieces, direct, enriched

_COUNTERS = Counter()
_CACHES = []


def _count(name, amount=1):
    _COUNTERS[name] += amount


def _cached(maxsize):
    def wrap(function):
        result = lru_cache(maxsize=maxsize)(function)
        _CACHES.append(result)
        return result
    return wrap


def reset_counters(clear_caches=False):
    if clear_caches:
        for cache in _CACHES:
            cache.cache_clear()
    _COUNTERS.clear()


def counters():
    return dict(_COUNTERS)


def _cell(left_root, right_root, degree, amplitude, offset, core=False):
    """Integrate the returned polynomial in its bounded local coordinate."""
    _count('cell_integrations')
    a, b = left_root**4, right_root**4
    width = b-a
    moments = [offset*(b**(k+1)-a**(k+1))/F(k+1)
               + amplitude*pieces._power_moment(k, left_root, right_root)
               for k in range(degree+1)]
    # s=(t-a)/width, expanded exactly once per degree instead of per basis.
    local_moments = [sum((F(comb(k,j))*(-a)**(k-j)*moments[j]
                          for j in range(k+1)), F(0))/width**k
                     for k in range(degree+1)]
    source_squared = (offset**2*width
                      + 2*offset*amplitude*pieces._power_moment(0,left_root,right_root)
                      + 2*amplitude**2*(right_root**2-left_root**2))
    entries, polynomial, projection_norm = [], [F(0)], F(0)
    for l, basis in enumerate(pieces._shifted_legendre(degree)):
        inner = sum((c*local_moments[k] for k,c in enumerate(basis)), F(0))
        mass = width/F(2*l+1)
        coefficient = inner/mass
        contribution = coefficient**2*mass
        projection_norm += contribution
        polynomial = pieces._add(polynomial, pieces._scale(basis,coefficient))
        entries.append({'degree':l,'shifted_legendre_in_s':list(map(str,basis)),
                        'source_inner_product_dt':str(inner),'basis_mass_dt':str(mass),
                        'coefficient':str(coefficient),'norm_squared_dt':str(contribution)})
    squared = pieces._mul(polynomial,polynomial)
    direct_squared = width*sum((c/F(k+1) for k,c in enumerate(squared)), F(0))
    cross = sum((c*local_moments[k] for k,c in enumerate(polynomial)), F(0))
    residual = source_squared-2*cross+direct_squared
    if cross != projection_norm or direct_squared != projection_norm or residual < 0:
        raise ArithmeticError('Exact local-coordinate projection identity failed')
    return {'left':str(a),'right':str(b),'left_fourth_root':str(left_root),
            'right_fourth_root':str(right_root),'core':core,'degree':degree,
            'coordinate':'s=(abs(t)-left)/(right-left)',
            'polynomial_in_s':list(map(str,polynomial)),
            'source_monomial_moments_dt':list(map(str,moments)),
            'projection':entries,'source_norm_squared_dt':str(source_squared),
            'approximation_norm_squared_dt':str(direct_squared),'cross_integral_dt':str(cross),
            'error_squared_dt':str(residual)}


@_cached(32)
def _shape_text(ratio, degree):
    _count('shape_integrations')
    return canonical(_cell(ratio, F(1), degree, F(1), F(0)))


def piecewise_model(function, levels=4, degree=3, root_ratio='1/2', sqrt_bits=40):
    """The legacy wire certificate, with a cached immutable unit shape."""
    function = pieces.normalize_function(function)
    pieces._integer(levels, 1, pieces.MAX_LEVELS, 'levels')
    pieces._integer(degree, 0, pieces.MAX_DEGREE, 'degree')
    pieces._integer(sqrt_bits, 0, 256, 'sqrt_bits')
    ratio = pieces.rational(root_ratio)
    if not 0 < ratio < 1 or ratio.denominator > 4096:
        raise ValueError('root_ratio must lie in (0,1) with denominator<=4096')
    if function['profile'] == 'step':
        return pieces.piecewise_model(function, levels, degree, root_ratio, sqrt_bits)
    amplitude, offset = F(function['amplitude']), F(function['offset'])
    shape = json.loads(_shape_text(ratio, degree))
    cells = [_cell(F(0), ratio**levels, 0, amplitude, offset, core=True)]
    for j in range(levels-1, -1, -1):
        _count('scaled_cells')
        cells.append(pieces._scaled_cell(shape, ratio**j, amplitude, offset))
    error = sum((F(cell['error_squared_dt']) for cell in cells), F(0))
    geometric = amplitude**2*(F(2, 9)*ratio**(2*levels)
        + F(shape['error_squared_dt'])*(1-ratio**(2*levels))/(1-ratio**2))
    source_norm = offset**2+F(8, 3)*offset*amplitude+2*amplitude**2
    if error != geometric or sum((F(c['source_norm_squared_dt']) for c in cells), F(0)) != source_norm:
        raise ArithmeticError('Complete scaled model identities failed')
    upper = pieces.sqrt_upper(error, sqrt_bits)
    return {'format': pieces.FORMAT, 'function': function, 'scope': pieces.SCOPE,
            'norm': 'L2_probability', 'error_squared': str(error), 'error_upper': str(upper),
            'sqrt_upper_squared': str(upper**2), 'source_norm_squared': str(source_norm),
            'levels': levels, 'degree': degree, 'root_ratio': str(ratio), 'sqrt_bits': sqrt_bits,
            'representation': 'even_reflection_of_positive_half_cells', 'cells': cells,
            'shape_proof': shape,
            'shape_rule': 't=R^4*x scales centered profile coefficients by 1/R and squared residual integral by R^2',
            'source_value_at_axis_zero': str(offset),
            'axis_zero_convention': 'step(0)=0; singular_power(0)=0 as an a.e. representative',
            'endpoint_convention': 'at nonzero seams choose the cell starting at abs(t); final endpoint included',
            'expanded_polynomial_coefficient_slots': 1+levels*(degree+1),
            'shape_coefficient_count': degree+1,
            'bound_statement': 'probability_L2_norm(original_function-piecewise_model)<=error_upper',
            'spectral_transfer_claimed': False, 'global_polynomial': False,
            'formal_assistant_checked': False}


verify_piecewise = pieces.verify_piecewise


def evaluate_many(certificate, coordinates):
    """Replay an owned certificate once, then evaluate all exact coordinates."""
    owned = deepcopy(certificate)
    _count('piecewise_batch_verifications')
    if not pieces.verify_piecewise(owned):
        raise ValueError('Invalid piecewise model')
    points = [pieces.rational(t) for t in coordinates]
    if any(not -1 <= t <= 1 for t in points):
        raise ValueError('Axis coordinate must be in [-1,1]')
    q = owned['function']
    if q['profile'] == 'step':
        a,b = F(q['amplitude']),F(q['offset'])
        return [b+(a if t > 0 else 0) for t in points]
    cells = [(F(c['left']),F(c['right']),tuple(map(F,c['polynomial_in_s'])))
             for c in owned['cells']]
    values = []
    for t in points:
        t = abs(t)
        for a,b,coefficients in cells:
            if a <= t < b or t == b == 1:
                s, value = (t-a)/(b-a),F(0)
                for coefficient in reversed(coefficients):
                    value = value*s+coefficient
                values.append(value)
                break
        else:
            raise ArithmeticError('A verified cover must contain this coordinate')
    return values


def evaluate(certificate, t):
    return evaluate_many(certificate, [t])[0]


@_cached(8192)
def _raw_moment(key, degree, power):
    _count('raw_moment_computations')
    profile,a,b,alpha = key
    basic = F(0) if degree % 2 else F(2,degree+1)
    if power == 0: return basic
    if profile == 'step':
        first = second = F(1,degree+1)
    else:
        first = F(0) if degree % 2 else 2/(alpha+degree+1)
        second = F(0) if degree % 2 else 2/(2*alpha+degree+1)
    return a*first+b*basic if power == 1 else a*a*second+2*a*b*first+b*b*basic


class MomentTable:
    """Normalize one owned source and deduplicate exact (degree,power) moments."""
    def __init__(self,function):
        _count('moment_source_normalizations')
        self.function = direct.normalize(function)
        q = self.function
        self.key = (q['profile'],F(q['amplitude']),F(q['offset']),F(q.get('exponent',0)))

    def moment(self,degree,power=1):
        direct.integer(degree,0,256,'moment_degree')
        direct.integer(power,0,2,'power')
        _count('raw_moment_requests')
        return _raw_moment(self.key,degree,power)

    def integral(self,m,l,k,power=1):
        _count('radial_integral_requests')
        profile,a,b,_ = self.key
        # The radial product has parity l+k. A step is not an even function.
        # Its even part is constant, so only its off-diagonal SAME-parity
        # entries vanish. q^2 is even only when a*(a+2*b)==0.
        if ((l+k)%2 and (profile == 'abs_power' or a == 0
                          or power == 0 or (power == 2 and a+2*b == 0))
                or (l != k and (power == 0 or a == 0
                                or (profile == 'step' and (l+k)%2 == 0)))):
            _count('proven_zero_integrals')
            return F(0)
        return sum((v*self.moment(i,power)
                    for i,v in enumerate(direct.radial_product(m,l,k)) if v),F(0))


def moment(function,degree,power=1):
    return MomentTable(function).moment(degree,power)


def matrix_assembly(function,m=0,mean_zero=True,modes=8,near_tail=0):
    table = MomentTable(function)
    direct.integer(m,0,12,'azimuth_m'); direct.integer(modes,1,32,'modes')
    direct.integer(near_tail,0,32,'near_tail')
    if type(mean_zero) is not bool: raise ValueError('mean_zero must be boolean')
    start = 1 if m == 0 and mean_zero else m
    degrees = list(range(start,start+modes))
    mass = [direct.wide.basis_mass(m,l) for l in degrees]
    v,t = [[[F(0)]*modes for _ in degrees] for _ in range(2)]
    for i,l in enumerate(degrees):
        for j in range(i,modes):
            v[i][j] = v[j][i] = table.integral(m,l,degrees[j])
            t[i][j] = t[j][i] = table.integral(m,l,degrees[j],2)
    a = [row[:] for row in v]
    for i,l in enumerate(degrees): a[i][i] += l*(l+1)*mass[i]
    c = [[F(0)]*modes for _ in degrees]
    for i in range(modes):
        for j in range(i,modes):
            c[i][j] = c[j][i] = t[i][j]-sum(
                (v[i][k]*v[j][k]/mass[k] for k in range(modes) if v[i][k] and v[j][k]),F(0))
    removed_constant = None
    if m == 0 and mean_zero:
        column = [table.integral(m,l,0) for l in degrees]
        removed_constant = {'mass':F(2),'column':column}
        for i in range(modes):
            for j in range(i,modes):
                c[i][j] -= column[i]*column[j]/2
                c[j][i] = c[i][j]
    near = []
    for k in range(start+modes,start+modes+near_tail):
        column = [table.integral(m,l,k) for l in degrees]
        mk = direct.wide.basis_mass(m,k)
        near.append({'degree':k,'mass':mk,'column':column})
        for i in range(modes):
            for j in range(i,modes):
                c[i][j] -= column[i]*column[j]/mk
                c[j][i] = c[i][j]
    if direct.base.inertia(c)[0]: raise ArithmeticError('Exact residual Gram must be positive semidefinite')
    return {'degrees':degrees,'mass':mass,'A':a,'V':v,'T':t,'C':c,
            'removed_constant':removed_constant,'near_tail':near}


class Kernel(direct.Kernel):
    """Drop-in exact kernel with fixed Schur numerators prepared once."""
    def __init__(self,function,m=0,mean_zero=True,modes=8,sqrt_bits=40,near_tail=0):
        self.function = direct.normalize(function)
        self.m,self.mean_zero,self.n,self.near_tail = m,mean_zero,modes,near_tail
        self.sqrt_bits = sqrt_bits
        self.data = matrix_assembly(self.function,m,mean_zero,modes,near_tail)
        self.form = direct.form_bound(self.function,sqrt_bits)
        self.start = self.data['degrees'][0]
        self.tail_start = self.start+modes
        self.beta = direct.tail_bound(self.form,self.tail_start)
        self.remainder_beta = direct.tail_bound(self.form,self.tail_start+near_tail)
        self._a = tuple(map(tuple,self.data['A']))
        self._mass = tuple(self.data['mass'])
        self._evidence = canonical(direct.encode(self.data))
        terms = [(self.remainder_beta,tuple((i,j,self.data['C'][i][j])
                 for i in range(modes) for j in range(i,modes) if self.data['C'][i][j]))]
        for near in self.data['near_tail']:
            entries = []
            for i in range(modes):
                for j in range(i,modes):
                    if near['column'][i] and near['column'][j]:
                        _count('schur_fixed_outer_products')
                        entries.append((i,j,near['column'][i]*near['column'][j]/near['mass']))
            terms.append((direct.tail_bound(self.form,near['degree']),tuple(entries)))
        self._terms = tuple(terms)

    def evidence(self):
        return json.loads(self._evidence)

    def matrix(self,x,kind='lower'):
        x = direct.f.rational(x)
        if kind not in ('lower','upper'): raise ValueError('Unknown comparison kind')
        if kind == 'lower' and x >= self.beta:
            raise ValueError('Strictly positive radial tail required')
        result = [list(row) for row in self._a]
        for i in range(self.n): result[i][i] -= x*self._mass[i]
        if kind == 'lower':
            for beta,entries in self._terms:
                denominator = beta-x
                for i,j,value in entries:
                    _count('schur_variable_divisions')
                    result[i][j] -= value/denominator
            for i in range(self.n):
                for j in range(i+1,self.n): result[j][i] = result[i][j]
        _count('kernel_matrix_evaluations')
        return result


def certify_sector(function,m=0,mean_zero=True,modes=8,bits=32,sqrt_bits=40,near_tail=0):
    direct.integer(bits,8,96,'bits')
    k = Kernel(function,m,mean_zero,modes,sqrt_bits,near_tail)
    analytic_low = direct.tail_bound(k.form,k.start)
    high = min(k.data['A'][i][i]/k.data['mass'][i] for i in range(k.n))
    low = analytic_low-1
    if not k.upper_holds(high): raise ArithmeticError('Ritz trial upper bracket failed')
    for _ in range(bits):
        mid = (low+high)/2
        if k.upper_holds(mid): high = mid
        else: low = mid
    upper = high
    distance,low = F(1),analytic_low-1
    for _ in range(64):
        if k.lower_holds(low): break
        distance *= 2; low = analytic_low-distance
    else: raise ValueError('Lower bracket budget exhausted')
    high = min(upper,k.beta)
    for _ in range(bits):
        mid = (low+high)/2
        if k.lower_holds(mid): low = mid
        else: high = mid
    return direct.sector_certificate(k,low,upper)


def full_ground(function,mean_zero=True,modes=8,bits=32,max_m=4,
                tolerance='1/100000000',sqrt_bits=40,near_tail=0):
    q = direct.normalize(function); direct.integer(max_m,0,12,'max_m')
    form,sectors = direct.form_bound(q,sqrt_bits),[]
    for m in range(max_m+1):
        sectors.append(certify_sector(q,m,mean_zero,modes,bits,sqrt_bits,near_tail))
        if direct.tail_bound(form,m+1) >= min(F(c['upper']) for c in sectors): break
    # The frozen aggregator independently reconstructs every sector.
    return direct.full_certificate(q,mean_zero,sectors,tolerance)


@dataclass(frozen=True, init=False)
class SourceContext:
    """Owned immutable source bytes after one complete frozen verification."""
    _text: str

    def __init__(self,source):
        owned = deepcopy(source)
        _count('source_context_verifications')
        if not direct.verify_full(owned,enriched.FUNCTION,True):
            raise ValueError('Source must replay for the original global mean-zero problem')
        if len(owned['sectors']) < 2:
            raise ValueError('Source must explicitly cover m=0 and m=1')
        object.__setattr__(self,'_text',canonical(owned))

    def source(self):
        return json.loads(self._text)


def validate_source(source):
    return source if type(source) is SourceContext else SourceContext(source)


@_cached(4096)
def _weighted_moment(power):
    _count('weighted_moment_computations')
    if power <= -1: raise ValueError('Moment diverges at the equator')
    return F(4)/((power+1)*(power+3))


@_cached(128)
def _action_terms(power):
    return tuple(enriched.action_terms(power).items())


def _power_inner(a,b):
    result = F(0)
    for r,x in a:
        for s,y in b:
            _count('weighted_moment_requests')
            result += x*y*_weighted_moment(r+s)
    return result


@_cached(64)
def _trial_matrices(ss):
    """Immutable exact prefixes; append only the newest Gram border."""
    n = len(ss)
    previous = _trial_matrices(ss[:-1]) if n > 1 else ((),(),())
    M,H,R = [[list(row)+[F(0)] for row in mat]+[[F(0)]*n] for mat in previous]
    last_basis,last_action = ((ss[-1],F(1)),),_action_terms(ss[-1])
    for j,s in enumerate(ss):
        basis,action = ((s,F(1)),),_action_terms(s)
        _count('gram_border_entries',3)
        M[-1][j] = M[j][-1] = _power_inner(last_basis,basis)
        H[-1][j] = H[j][-1] = _power_inner(last_basis,action)
        R[-1][j] = R[j][-1] = _power_inner(last_action,action)
        if H[-1][j] != _power_inner(basis,last_action):
            raise ArithmeticError('Operator integration-by-parts symmetry failed')
    _count('gram_prefix_extensions')
    return tuple(tuple(map(tuple,mat)) for mat in (M,H,R))


def trial_matrices(powers):
    return _trial_matrices(enriched.exponents(powers))


def trial_statistics(powers,coefficients):
    ss = enriched.exponents(powers)
    if not isinstance(coefficients,(list,tuple)) or len(coefficients) != len(ss):
        raise ValueError('A rational coefficient is required for every power')
    v = tuple(direct.f.rational(x) for x in coefficients)
    if not any(v): raise ValueError('Trial function must be nonzero')
    M,H,R = _trial_matrices(ss)
    mass = enriched._quadratic(M,v)
    if mass <= 0: raise ArithmeticError('Nonpositive trial mass')
    energy = enriched._quadratic(H,v)
    operator_norm = enriched._quadratic(R,v)
    mu = energy/mass
    residual = operator_norm/mass-mu*mu
    if residual < 0: raise ArithmeticError('Negative complete residual norm')
    return {'mass':str(mass),'operator_form':str(energy),
            'operator_norm_squared':str(operator_norm),'rayleigh':str(mu),
            'residual_squared':str(residual)}


def propose_trial(powers,precision_bits=112,iterations=14):
    return enriched.propose_trial(powers,precision_bits,iterations)


def trial_certificate(powers,coefficients,source,tolerance='1/100000000'):
    ss = enriched.exponents(powers)
    tol = direct.f.rational(tolerance)
    if not F(1,10**30) <= tol <= 1: raise ValueError('Tolerance must be in 1e-30..1')
    source = validate_source(source).source()
    stats = trial_statistics(ss,coefficients)
    v = [str(direct.f.rational(x)) for x in coefficients]
    mu,residual = F(stats['rayleigh']),F(stats['residual_squared'])
    gap = enriched.second_eigenvalue_bound()
    beta = F(gap['lower'])
    if mu >= beta:
        raise ValueError('Temple requires Rayleigh quotient below the second-sector bound')
    temple = mu-residual/(beta-mu)
    m1_lower = max(temple,F(source['sectors'][1]['lower']))
    other = [{'azimuth_m':c['azimuth_m'],'lower':c['lower']}
             for c in source['sectors'] if c['azimuth_m'] != 1]
    angular = F(source['angular_tail_lower'])
    lower = min([m1_lower,angular]+[F(c['lower']) for c in other])
    upper = min(mu,F(source['upper']))
    if lower > upper: raise ArithmeticError('Global spectral contradiction')
    M,H,R = trial_matrices(ss)
    return {'format':enriched.FORMAT,'function':deepcopy(enriched.FUNCTION),
            'geometry':'unit_S2','measure':source['measure'],
            'scope':'all_real_mean_zero_H1_on_unit_S2','mean_zero':True,'eigenvalue_index':1,
            'trial':{'azimuth_m':1,'real_component':'cos_phi','radial_parity':'even',
                     'formula':'sqrt(1-z^2)*sum(a_s*abs(z)^s)*cos(phi)',
                     'powers':list(map(str,ss)),'coefficients':v,
                     'strong_operator_domain':'s=0_or_s>3/2','azimuth_factor_cancels':True},
            'matrices':{'M':enriched._encoded_matrix(M),'H':enriched._encoded_matrix(H),
                        'R':enriched._encoded_matrix(R)},
            'statistics':stats,'second_eigenvalue_proof':gap,
            'temple_lower':str(temple),'m1_lower':str(m1_lower),
            'other_sector_lowers':other,'angular_tail_lower':str(angular),
            'source':source,'source_digest':digest(source),
            'lower':str(lower),'upper':str(upper),'exact_width':str(upper-lower),
            'tolerance':str(tol),'status':'target_met' if upper-lower <= tol else 'certified_open',
            'full_infinite_space_covered':True,'original_function_moments_exact':True,
            'function_approximation_error':'0','full_strong_residual_not_projected_residual':True,
            'formal_proof_assistant_checked':False}


def enrich(source,max_terms=16,tolerance='1/100000000',precision_bits=112,iterations=14):
    direct.integer(max_terms,1,16,'max_terms')
    context = validate_source(source)
    started = perf_counter()
    counts = sorted(set([n for n in (1,3,6,10,16) if n <= max_terms]+[max_terms]))
    attempts,best = [],None
    for n in counts:
        begin = perf_counter()
        proposal = propose_trial(enriched.DEFAULT_EXPONENTS[:n],precision_bits,iterations)
        cert = trial_certificate(proposal['powers'],proposal['coefficients'],context,tolerance)
        attempts.append({'terms':n,'powers':proposal['powers'],
                         'coefficient_bits':precision_bits,'iterations':iterations,
                         'elapsed_seconds':perf_counter()-begin,
                         'proposal_seconds':proposal['elapsed_seconds'],
                         'rayleigh':cert['statistics']['rayleigh'],
                         'residual_squared':cert['statistics']['residual_squared'],
                         'lower':cert['lower'],'upper':cert['upper'],
                         'exact_width':cert['exact_width'],'status':cert['status']})
        if best is None or F(cert['exact_width']) < F(best['exact_width']): best = cert
        if best['status'] == 'target_met': break
    _count('final_enriched_frozen_replays')
    if not enriched.verify(best,expected_source=context.source(),expected_tolerance=tolerance):
        raise ArithmeticError('Final independent certificate replay failed')
    return {'format':enriched.RUN_FORMAT,'certificate':best,'attempts':attempts,
            'elapsed_seconds_including_final_replay':perf_counter()-started,
            'status':best['status'],'max_terms':max_terms,
            'precision_goal':str(direct.f.rational(tolerance)),
            'search_is_not_a_global_optimality_claim_over_bases':True}
