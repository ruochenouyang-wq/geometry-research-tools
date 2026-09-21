"""V165--174: exact, finite-dimensional improvements of sphere GN trials.

All moments use dmu=dt/2.  J=N/(4*M*E) is pi times the area-measure
GN quotient, and provides a trial LOWER bound only on its full-space supremum.
Search proposals are untrusted; exact moment identities and Sturm sign proofs
are replayed without running the optimization. No frozen input adapter is used.
"""
from fractions import Fraction as F
from math import comb, sqrt
from pathlib import Path
import json
import re
import sys

sys.dont_write_bytecode = True
from common import arithmetic as a, exact_extrema as ex, same

MAX_DEGREE = 36
SCOPE = 'axisymmetric_mean_zero_polynomial_trial_on_unit_S2'
RATIO = 'pi_times_area_GN_K4=N/(4*M*E),probability_moments'
add, mul, scale, derivative, evaluate = a.add, a.mul, a.scale, a.derivative, a.evaluate


def rational(x):
    if type(x) in (int, F):
        return F(x)
    if isinstance(x, str) and re.fullmatch(r'-?\d+(?:/[1-9]\d*)?', x):
        return F(x)
    raise ValueError('Exact integer or rational coefficient required')


def polynomial(raw, max_degree=MAX_DEGREE):
    if not isinstance(raw, (list, tuple)) or not 1 <= len(raw) <= max_degree+1:
        raise ValueError('Polynomial degree budget exceeded')
    p = [rational(x) for x in raw]
    while len(p)>1 and p[-1] == 0:
        p.pop()
    return p


def strings(p):
    return [str(x) for x in p]


def mean(p):
    return sum((F(x, i+1) for i,x in enumerate(p) if i%2 == 0), F(0))


def power(p, n):
    result = [F(1)]
    for _ in range(n):
        result = mul(result, p)
    return result


def inner(p, q):
    return mean(mul(p, q))


def centered(p):
    return add(p, [-mean(p)])


def laplacian(p):
    """Positive -Delta applied to an axisymmetric polynomial."""
    return scale(derivative(mul([F(1),F(0),F(-1)], derivative(p))), -1)


def moments(raw):
    """V165: exact moments, including degree 144 quartics of degree 36 trials."""
    p = polynomial(raw)
    m, e, n = inner(p,p), inner(p,laplacian(p)), mean(power(p,4))
    return {'format':'gn_moments_v165','scope':'axisymmetric_polynomial_on_unit_S2_no_mean_constraint','measure':'dt/2',
            'u':strings(p),'degree':len(p)-1,'mean':str(mean(p)),
            'mass':str(m),'energy':str(e),'quartic':str(n)}


def normalize(raw):
    """V166: center and fix a rational monic representative; reject constants."""
    p = polynomial(raw)
    c = centered(p)
    if c == [0]:
        raise ValueError('A constant trial has zero centered mass and energy')
    multiplier = 1/c[-1]
    u = scale(c, multiplier)
    return {'format':'gn_normalize_v166','original_u':strings(p),
            'removed_mean':str(mean(p)), 'multiplier':str(multiplier),
            'normalized_u':strings(u), 'moments':moments(u)}


def ratio(raw):
    """V167: exact J=pi*K_area certificate for one admissible trial."""
    p = polynomial(raw)
    m = moments(p)
    M,E,N = [F(m[k]) for k in ('mass','energy','quartic')]
    if mean(p) != 0 or M <= 0 or E <= 0:
        raise ValueError('Require a nonzero, mean-zero, positive-energy trial')
    return {'format':'gn_ratio_v167','scope':SCOPE,'ratio_convention':RATIO,
            'u':strings(p),'moments':m,'pi_K_lower':str(N/(4*M*E)),
            'claim':'trial_lower_bound_only_not_universal_upper'}


def residual(raw):
    """V168: the COMPLETE polynomial Euler--Lagrange residual, degree <=108."""
    cert = ratio(raw); p = polynomial(cert['u'])
    M,E,N = [F(cert['moments'][k]) for k in ('mass','energy','quartic')]
    before = add(add(laplacian(p), scale(p,E/M)), scale(power(p,3),-2*E/N))
    r = centered(before)
    return {'format':'gn_residual_v168','u':strings(p),'ratio':cert,
            'residual':strings(r),'removed_mean':str(mean(before)),
            'full_degree':len(r)-1,'norm_squared':str(inner(r,r)),
            'inner_u':str(inner(r,p)),'full_stationary':r == [0]}


def direction(raw, v=None):
    """V169: exact derivative. Default v=-r has J'=2*J*||r||^2/E."""
    rcert = residual(raw); p = polynomial(raw)
    r = polynomial(rcert['residual'],108)
    v = scale(r,-1) if v is None else polynomial(v,108)
    if mean(v) != 0:
        raise ValueError('Direction must preserve the zero-mean constraint')
    M,E,N = [F(rcert['ratio']['moments'][k]) for k in ('mass','energy','quartic')]
    dm,de,dn = 2*inner(p,v), 2*inner(laplacian(p),v), 4*inner(power(p,3),v)
    J = N/(4*M*E)
    direct = J*(dn/N-dm/M-de/E)
    identity = -2*J*inner(r,v)/E
    if direct != identity:
        raise ArithmeticError('Euler derivative identity failed')
    return {'format':'gn_direction_v169','u':strings(p),'v':strings(v),
            'residual':rcert,'d_mass':str(dm),'d_energy':str(de),'d_quartic':str(dn),
            'inner_residual_direction':str(inner(r,v)),
            'derivative':str(direct),'identity_derivative':str(identity),
            'strict_local_improvement':direct > 0,
            'claim':'local_positive_derivative_only_not_a_finite_step_guarantee'}


def _basis(raw):
    if not isinstance(raw,(list,tuple)) or len(raw) != 2:
        raise ValueError('Two independent trial functions required')
    f,g = [polynomial(x) for x in raw]
    ratio(f); ratio(g)
    if inner(f,f)*inner(g,g)-inner(f,g)**2 <= 0:
        raise ValueError('Plane basis is linearly dependent')
    return f,g


def chart_moments(f,g):
    mass = [inner(f,f),2*inner(f,g),inner(g,g)]
    energy = [inner(f,laplacian(f)),2*inner(f,laplacian(g)),inner(g,laplacian(g))]
    numerator = [comb(4,j)*mean(mul(power(f,4-j),power(g,j))) for j in range(5)]
    return {'mass':a.trim(mass),'energy':a.trim(energy),'numerator':a.trim(numerator),
            'denominator':scale(mul(mass,energy),4)}


def _propose(m):
    """Floating suggestions with exact objective selection, including both ends."""
    n,d = m['numerator'],m['denominator']
    def score(x):
        return evaluate(n,x)/evaluate(d,x)
    grid = [F(j,64) for j in range(-64,65)]
    values = [score(x) for x in grid]
    proposals = [F(-1),F(0),F(1),grid[max(range(len(grid)),key=values.__getitem__)]]
    try:
        nf,df = [float(x) for x in n],[float(x) for x in d]
        def fs(x):
            return evaluate(nf,x)/evaluate(df,x)
        phi = (sqrt(5)-1)/2
        for j in range(1,len(grid)-1):
            if values[j]>=values[j-1] and values[j]>=values[j+1]:
                left,right = float(grid[j-1]),float(grid[j+1])
                for _ in range(48):
                    x,y = right-phi*(right-left),left+phi*(right-left)
                    if fs(x)<fs(y): left=x
                    else: right=y
                proposals.append(F((left+right)/2).limit_denominator(10**6))
    except (OverflowError,ZeroDivisionError):
        pass
    s = max(proposals,key=score)
    return s,score(s)


def _sign(p, strict=False):
    normalizer = 1/max(map(abs,p)) if any(p) else F(1)
    proof = ex.maximize(scale(p,normalizer),tolerance=F(1,10**20),max_nodes=128)
    bound = F(proof['maximum_upper'])
    return {'positive_multiplier':str(normalizer),'extrema':proof}, (bound<0 if strict else bound<=0)


def _verify_sign(proof,p,strict):
    c = rational(proof['positive_multiplier']); e = proof['extrema']
    if c<=0 or e['polynomial'] != strings(scale(p,c)) or e['interval'] != ['-1','1'] or not ex.verify(e):
        return False
    b = F(e['maximum_upper'])
    return b<0 if strict else b<=0


def plane(raw_basis,tolerance=F(1,10**8)):
    """V170: certify a GLOBAL maximum on this complete real projective plane.

    The two charts f+s*g and g+s*f, |s|<=1, cover all nonzero real pairs.
    This is not an upper bound on the quotient outside the recorded plane.
    """
    basis = _basis(raw_basis); tolerance = rational(tolerance)
    if not F(1,10**18)<=tolerance<=1:
        raise ValueError('Plane tolerance must lie in [1e-18,1]')
    chart_data = [chart_moments(*basis),chart_moments(*reversed(basis))]
    proposals = [_propose(m) for m in chart_data]
    chart = max(range(2),key=lambda j:proposals[j][1]); s,lower = proposals[chart]
    upper = lower+tolerance
    attempts=[]
    for _ in range(8):
        signs = [_sign(add(m['numerator'],scale(m['denominator'],-upper))) for m in chart_data]
        attempts.append({'upper':str(upper),'accepted_charts':[x[1] for x in signs]})
        if all(x[1] for x in signs): break
        upper = lower+10*(upper-lower)
    else: raise ArithmeticError('Global plane ceiling remained unproved within budget')
    charts=[]
    for j,m in enumerate(chart_data):
        mp,mok = _sign(scale(m['mass'],-1),True)
        ep,eok = _sign(scale(m['energy'],-1),True)
        if not mok or not eok: raise ArithmeticError('Plane denominator positivity not proved')
        charts.append({'order':[j,1-j],'interval':['-1','1'],
            'moments':{k:strings(v) for k,v in m.items()},'mass_positive':mp,
            'energy_positive':ep,'upper_nonpositive':signs[j][0]})
    witness = add(basis[chart],scale(basis[1-chart],s))
    return {'format':'gn_plane_v170','scope':'all_nonzero_real_combinations_of_recorded_basis_only',
        'ratio_convention':RATIO,'basis':[strings(p) for p in basis],
        'lower':str(lower),'upper':str(upper),'gap':str(upper-lower),
        'requested_gap':str(tolerance),'target_met':upper-lower<=tolerance,
        'witness':{'chart':chart,'parameter':str(s),'u':strings(witness)},
        'charts':charts,'ceiling_attempts':attempts,'universal_GN_upper_bound':False}


def verify_plane(c,expected_basis=None):
    try:
        if c['format']!='gn_plane_v170' or c['scope']!='all_nonzero_real_combinations_of_recorded_basis_only' or c['ratio_convention']!=RATIO or c['universal_GN_upper_bound'] is not False:
            return False
        basis = _basis(c['basis'])
        if expected_basis is not None and basis!=_basis(expected_basis): return False
        lower,upper = rational(c['lower']),rational(c['upper'])
        gap = upper-lower; target=rational(c['requested_gap'])
        if not 0<lower<=upper or rational(c['gap'])!=gap or not F(1,10**18)<=target<=1 or c['target_met']!=(gap<=target) or len(c['charts'])!=2: return False
        for j,ch in enumerate(c['charts']):
            m = chart_moments(basis[j],basis[1-j])
            if ch['order']!=[j,1-j] or ch['interval']!=['-1','1'] or ch['moments']!={k:strings(v) for k,v in m.items()}: return False
            for field,p,strict in [('mass_positive',scale(m['mass'],-1),True),('energy_positive',scale(m['energy'],-1),True),('upper_nonpositive',add(m['numerator'],scale(m['denominator'],-upper)),False)]:
                if not _verify_sign(ch[field],p,strict): return False
        w=c['witness']; j=w['chart']; s=rational(w['parameter'])
        if type(j) is not int or j not in (0,1) or not -1<=s<=1: return False
        u=add(basis[j],scale(basis[1-j],s))
        return w['u']==strings(u) and F(ratio(u)['pi_K_lower'])==lower
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError): return False


def legendre(n):
    out=[[F(1)]]
    if n: out.append([F(0),F(1)])
    for l in range(1,n):
        out.append(scale(add(scale(mul([F(0),F(1)],out[-1]),2*l+1),scale(out[-2],-l)),F(1,l+1)))
    return out


def project_residual(raw,degree_budget=12):
    """V172: L2 spectral projection plus the full, untruncated complement."""
    if type(degree_budget) is not int or not 0<=degree_budget<=MAX_DEGREE:
        raise ValueError('Spectral degree budget must lie in 0..36')
    rc=residual(raw); r=polynomial(rc['residual'],108)
    basis=legendre(len(r)-1); inside=[F(0)]; coefficients=[]
    for l,p in enumerate(basis):
        c=(2*l+1)*inner(r,p)
        coefficients.append(c)
        if 1<=l<=degree_budget: inside=add(inside,scale(p,c))
    outside=add(r,scale(inside,-1))
    return {'format':'gn_project_residual_v172','u':rc['u'],'degree_budget':degree_budget,
        'full':rc,'legendre_coefficients':strings(coefficients),
        'inside':strings(inside),'outside':strings(outside),
        'inside_support':[l for l,c in enumerate(coefficients) if c and 1<=l<=degree_budget],
        'outside_support':[l for l,c in enumerate(coefficients) if c and l>degree_budget],
        'inside_norm_squared':str(inner(inside,inside)),
        'outside_norm_squared':str(inner(outside,outside)),
        'cross_inner':str(inner(inside,outside)),
        'projected_stationary':inside==[0],'full_stationary':r==[0]}


def refine(raw,iterations=2,degree_budget=12,tolerance=F(1,10**8)):
    """V171: monotone accepted plane updates, with honest finite-budget status."""
    if type(iterations) is not int or not 0<=iterations<=8:
        raise ValueError('Iteration budget must lie in 0..8')
    initial=normalize(raw); u=polynomial(initial['normalized_u'])
    if type(degree_budget) is not int or not len(u)-1<=degree_budget<=MAX_DEGREE:
        raise ValueError('Degree budget must include the initial trial and be <=36')
    steps=[]; stop='iteration_budget_exhausted'
    for _ in range(iterations):
        p=project_residual(u,degree_budget)
        v=scale(polynomial(p['inside']),-1)
        if v==[0]:
            stop='full_stationary' if p['full_stationary'] else 'projected_stationary_with_nonzero_full_residual'
            break
        v=scale(v,1/max(map(abs,v)))
        d=direction(u,v)
        if not d['strict_local_improvement']: raise ArithmeticError('Projected descent direction did not improve J')
        try: pc=plane([u,v],tolerance)
        except (ValueError,ArithmeticError) as error:
            steps.append({'accepted':False,'before':ratio(u),'projection':p,'direction':d,'failure':str(error)})
            stop='plane_certificate_budget_failed'; break
        normalized=normalize(pc['witness']['u']); candidate=polynomial(normalized['normalized_u'])
        before,after=ratio(u),ratio(candidate)
        accepted=F(after['pi_K_lower'])>F(before['pi_K_lower'])
        steps.append({'accepted':accepted,'before':before,'projection':p,'direction':d,
                      'plane':pc,'candidate_normalization':normalized,'after':after})
        if not accepted: stop='no_strict_improvement_found'; break
        u=candidate
    return {'format':'gn_refine_v171','original_u':strings(polynomial(raw)),
        'normalization':initial,'iteration_budget':iterations,'degree_budget':degree_budget,
        'steps':steps,'final':ratio(u),'final_residual':residual(u),'stop':stop,
        'full_function_space_global_optimum_claimed':False}


def verify_refine(c,expected_u=None):
    try:
        raw=polynomial(c['original_u'])
        if expected_u is not None and raw!=polynomial(expected_u): return False
        if c['format']!='gn_refine_v171' or c['full_function_space_global_optimum_claimed'] is not False: return False
        norm=normalize(raw)
        if not same(norm,c['normalization']): return False
        u=polynomial(norm['normalized_u']); budget=c['degree_budget']; count=c['iteration_budget']
        if type(count) is not int or not 0<=count<=8 or len(c['steps'])>count or type(budget) is not int or not len(u)-1<=budget<=36: return False
        for j,step in enumerate(c['steps']):
            p=project_residual(u,budget)
            if not same(step['before'],ratio(u)) or not same(step['projection'],p): return False
            v=scale(polynomial(p['inside']),-1)
            if v==[0]: return False
            v=scale(v,1/max(map(abs,v)))
            if not same(step['direction'],direction(u,v)): return False
            if 'failure' in step:
                if step['accepted'] is not False or j!=len(c['steps'])-1 or c['stop']!='plane_certificate_budget_failed': return False
                break
            if not verify_plane(step['plane'],[u,v]): return False
            nc=normalize(step['plane']['witness']['u']); next_u=polynomial(nc['normalized_u'])
            after=ratio(next_u)
            if not same(step['candidate_normalization'],nc) or not same(step['after'],after): return False
            accepted=F(after['pi_K_lower'])>F(ratio(u)['pi_K_lower'])
            if type(step['accepted']) is not bool or step['accepted']!=accepted: return False
            if accepted: u=next_u
            elif j!=len(c['steps'])-1 or c['stop']!='no_strict_improvement_found': return False
        stop=c['stop']; projected=project_residual(u,budget)
        if stop=='iteration_budget_exhausted' and len(c['steps'])!=count: return False
        if stop=='full_stationary' and not projected['full_stationary']: return False
        if stop=='projected_stationary_with_nonzero_full_residual' and (not projected['projected_stationary'] or projected['full_stationary']): return False
        if stop=='no_strict_improvement_found' and (not c['steps'] or c['steps'][-1]['accepted'] or 'plane' not in c['steps'][-1]): return False
        if stop=='plane_certificate_budget_failed' and (not c['steps'] or 'failure' not in c['steps'][-1]): return False
        if stop not in ('iteration_budget_exhausted','full_stationary','projected_stationary_with_nonzero_full_residual','plane_certificate_budget_failed','no_strict_improvement_found'): return False
        return same(c['final'],ratio(u)) and same(c['final_residual'],residual(u))
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError): return False


STAGES = [
    {'version':165,'capability':'Exact mass, energy and quartic moments through trial degree 36','callable':'gn_variation.moments','evidence':['round08/results.json','test_gn_variation.py'],'outcome':'degree36_moments_replayed'},
    {'version':166,'capability':'Exact zero-mean centering and nonzero rational scale normalization','callable':'gn_variation.normalize','evidence':['round08/results.json','test_gn_variation.py'],'outcome':'centering_scale_and_constant_rejection_tested'},
    {'version':167,'capability':'Bound trial to exact area-normalized pi*K GN quotient','callable':'gn_variation.ratio','evidence':['round08/results.json','round08/PROOF.md'],'outcome':'trial_lower_bound_replayed'},
    {'version':168,'capability':'Complete projected Euler-Lagrange residual including cubic high degrees','callable':'gn_variation.residual','evidence':['round08/results.json','test_gn_variation.py'],'outcome':'degree108_full_residual_tested'},
    {'version':169,'capability':'Exact directional derivative and certified positive local direction','callable':'gn_variation.direction','evidence':['round08/results.json','test_gn_variation.py'],'outcome':'positive_direction_and_reversed_sign_replayed'},
    {'version':170,'capability':'Global 2D projective-plane GN maximum via two complete charts','callable':'gn_variation.plane','evidence':['round08/high_degree_plane.json','test_gn_variation.py'],'outcome':'two_chart_high_degree_plane_replayed'},
    {'version':171,'capability':'Finite-budget monotone accepted plane refinements','callable':'gn_variation.refine','evidence':['round08/refinement.json','test_gn_variation.py'],'outcome':'fixed_initial_trial_improved_with_finite_budget'},
    {'version':172,'capability':'Legendre degree-budget projection with full inside/outside residual evidence','callable':'gn_variation.project_residual','evidence':['round08/results.json','test_gn_variation.py'],'outcome':'projected_zero_full_nonzero_distinguished'},
    {'version':173,'capability':'Centered-cap exact trial generation and adaptive target/budget screening','callable':'gn_variation.screen_caps','evidence':['round08/caps.json','test_gn_variation.py'],'outcome':'n2_4_8_12_recorded_and_target_stopping_tested'},
    {'version':174,'capability':'Compile GN trial square to degree24 potential and replay full-spectrum diagnostics','callable':'gn_variation.spectral_diagnostic','evidence':['round08/spectral_diagnostics.json','round08/PROOF.md'],'outcome':'fixed_potential_diagnostics_only_no_universal_GN_claim'},
]


def centered_cap(n):
    if type(n) is not int or not 1<=n<=36: raise ValueError('Cap degree must be in 1..36')
    p=[F(comb(n,j),2**n) for j in range(n+1)];p[0]-=F(1,n+1)
    return p


def screen_caps(degrees=(2,4,8,12),budget=4,target=None):
    """V173: adaptive target stopping and explicit retention of untried caps."""
    if not isinstance(degrees,(list,tuple)) or not degrees or len(set(degrees))!=len(degrees): raise ValueError('Distinct cap degrees required')
    for n in degrees: centered_cap(n)
    if type(budget) is not int or not 0<=budget<=len(degrees): raise ValueError('Invalid screening budget')
    target=None if target is None else rational(target)
    rows=[]; best=None; stop='candidate_list_exhausted'
    for n in sorted(degrees):
        if len(rows)==budget: stop='screening_budget_exhausted';break
        rc=ratio(centered_cap(n)); value=F(rc['pi_K_lower'])
        improved=best is None or value>F(best['ratio']['pi_K_lower'])
        row={'degree':n,'ratio':rc,'improves_best':improved};rows.append(row)
        if improved: best=row
        if target is not None and value>=target: stop='requested_trial_lower_bound_met';break
    tested=[r['degree'] for r in rows]
    return {'format':'gn_caps_v173','degrees':list(degrees),'budget':budget,
        'target':None if target is None else str(target),'trials':rows,'best':best,
        'untried':[n for n in sorted(degrees) if n not in tested],'stop':stop,
        'sharp_GN_constant_claimed':False}


def compile_potential(raw):
    """V174: exact trial-to-potential identity, not a universal GN theorem."""
    u=polynomial(raw,12);rc=ratio(u)
    M,E,N=[F(rc['moments'][k]) for k in ('mass','energy','quartic')]
    q=scale(power(u,2),-2*E/N)
    potential_energy=mean(mul(q,power(u,2)))
    return {'format':'gn_compile_v174','u':strings(u),'ratio':rc,
        'q_coefficients':strings(q),'q_degree':len(q)-1,
        'potential_energy':str(potential_energy),'trial_rayleigh':str((E+potential_energy)/M),
        'required_identity':str(-E/M),
        'claim':'fixed_potential_trial_identity_only_not_a_universal_GN_constant'}


def spectral_diagnostic(raw,**options):
    """V174: independently verified wide-polynomial full-sphere diagnostics.

    A separated lower ground state is a necessary-condition diagnostic for a
    putative full-space GN maximizer. Lack of separation proves no optimality.
    """
    import wide_potential as wide
    compiled=compile_potential(raw)
    try:
        c=wide.full_ground(compiled['q_coefficients'],mean_zero=True,**options)
        if not wide.verify_full(c,expected_q=compiled['q_coefficients'],expected_mean_zero=True):
            raise ArithmeticError('Wide spectral certificate failed replay')
        status='lower_spectral_ground_separated' if F(c['upper'])<F(compiled['trial_rayleigh']) else 'no_separated_lower_ground_proved'
        return {'format':'gn_spectral_diagnostic_v174','compiled':compiled,
            'spectral':c,'status':status,'universal_GN_claim':False}
    except (ValueError,ArithmeticError) as error:
        return {'format':'gn_spectral_diagnostic_v174','compiled':compiled,
            'status':'rejected_or_budget_failed','exception':type(error).__name__,
            'message':str(error),'universal_GN_claim':False}


def verify(c,expected_u=None):
    """Replay certificates and bind the original polynomial when requested."""
    try:
        kind=c['format']
        if kind=='gn_plane_v170':
            return (expected_u is None or polynomial(expected_u)==polynomial(c['basis'][0])) and verify_plane(c)
        if kind=='gn_refine_v171': return verify_refine(c,expected_u)
        if kind=='gn_caps_v173':
            return expected_u is None and same(c,screen_caps(c['degrees'],c['budget'],c['target']))
        if kind=='gn_spectral_diagnostic_v174':
            if c['universal_GN_claim'] is not False or not verify(c['compiled'],expected_u): return False
            if c['status']=='rejected_or_budget_failed':
                return isinstance(c.get('message'),str) and c.get('exception') in ('ValueError','ArithmeticError') and 'spectral' not in c
            import wide_potential as wide
            s=c['spectral'];q=c['compiled']['q_coefficients']
            if not wide.verify_full(s,expected_q=q,expected_mean_zero=True): return False
            status='lower_spectral_ground_separated' if F(s['upper'])<F(c['compiled']['trial_rayleigh']) else 'no_separated_lower_ground_proved'
            return c['status']==status
        mapping={'gn_moments_v165':moments,'gn_normalize_v166':normalize,'gn_ratio_v167':ratio,
            'gn_residual_v168':residual,'gn_direction_v169':lambda u:direction(u,c['v']),
            'gn_project_residual_v172':lambda u:project_residual(u,c['degree_budget']),
            'gn_compile_v174':compile_potential}
        if kind not in mapping: return False
        raw=polynomial(c['original_u'] if kind=='gn_normalize_v166' else c['u'])
        if expected_u is not None and raw!=polynomial(expected_u): return False
        return same(c,mapping[kind](raw))
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError): return False
