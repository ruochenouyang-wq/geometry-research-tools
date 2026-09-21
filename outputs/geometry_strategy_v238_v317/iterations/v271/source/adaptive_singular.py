"""Structure-generated rational trials with complete original-sphere proofs.

Candidate policies are human-designed and independently checked by certificate
replay. A small basis is never the mathematical domain of the claimed bound.
"""
from copy import deepcopy
from functools import lru_cache
from heapq import heappop, heappush
from time import perf_counter
import json

from backend import F, FROZEN, direct, enriched, canonical, digest

FORMAT = 'adaptive_singular_temple_v1'
ORIGINAL = deepcopy(enriched.FUNCTION)
rat = direct.f.rational


def normalize(function):
    q = direct.normalize(function)
    if q != ORIGINAL:
        raise ValueError('This stage supports the original singular potential only')
    return q


def generated_powers(function, count=16):
    """The operator shifts suggest the additive semigroup <2, 2+alpha>."""
    q = normalize(function)
    direct.integer(count, 1, 16, 'count')
    increments = (F(2), 2+F(q['exponent']))
    heap, seen, result = [F(0)], {F(0)}, []
    while len(result) < count:
        s = heappop(heap)
        result.append(s)
        for increment in increments:
            t = s+increment
            if t not in seen:
                seen.add(t)
                heappush(heap, t)
    return enriched.exponents(result)


def trial_matrices(function, powers):
    normalize(function)
    return enriched.trial_matrices(powers)


def statistics(function, powers, coefficients):
    normalize(function)
    return enriched.trial_statistics(powers, coefficients)


@lru_cache(maxsize=32)
def _source_json(qkey, source_modes, source_bits):
    q = json.loads(qkey)
    if q == ORIGINAL and source_modes == 4 and source_bits == 28:
        source = json.loads((FROZEN/'certificates'/'singular_n4.json').read_text())
    else:
        source = direct.full_ground(q, mean_zero=True, modes=source_modes, bits=source_bits)
    if not direct.verify_full(source, q, True):
        raise ArithmeticError('Original-function source did not replay')
    return canonical(source)


def source_certificate(function, modes=4, bits=28):
    q = normalize(function)
    direct.integer(modes, 1, 16, 'source_modes')
    direct.integer(bits, 8, 64, 'source_bits')
    return json.loads(_source_json(canonical(q), modes, bits))


def gap_bound(function, policy='full_radial'):
    q = normalize(function)
    if policy != 'full_radial':
        raise ValueError('Unknown spectral gap policy')
    form = direct.form_bound(q, 40)
    return {'policy': policy, 'form': form,
            'space': 'one_real_cos_m1_component_both_radial_parities',
            'free_second_eigenvalue': '6',
            'beta': str(direct.tail_bound(form, 2)),
            'sin_component_is_unitarily_equivalent': True}


def certificate(function, powers, coefficients, source, tolerance='1/100000000',
                gap_policy='full_radial'):
    q = normalize(function)
    ss = enriched.exponents(powers)
    v = [rat(x) for x in coefficients]
    tol = rat(tolerance)
    if not F(1, 10**30) <= tol <= 1:
        raise ValueError('Tolerance must lie in 1e-30..1')
    if not direct.verify_full(source, q, True) or len(source['sectors']) < 2:
        raise ValueError('Require a replayable full mean-zero source with m0 and m1')
    stats = statistics(q, ss, v)
    gap = gap_bound(q, gap_policy)
    mu, variance, beta = F(stats['rayleigh']), F(stats['residual_squared']), F(gap['beta'])
    if mu >= beta:
        raise ValueError('Temple is not applicable: mu >= beta')
    temple = mu-variance/(beta-mu)
    m1 = max(temple, F(source['sectors'][1]['lower']))
    others = [{'azimuth_m': c['azimuth_m'], 'lower': c['lower']}
              for c in source['sectors'] if c['azimuth_m'] != 1]
    angular = F(source['angular_tail_lower'])
    lower = min([m1, angular]+[F(c['lower']) for c in others])
    upper = min(mu, F(source['upper']))
    if lower > upper:
        raise ArithmeticError('Inconsistent global spectral bounds')
    M, A, R = trial_matrices(q, ss)
    return {'format': FORMAT, 'function': q, 'geometry': 'unit_S2',
            'measure': source['measure'], 'mean_zero': True,
            'scope': 'all_real_mean_zero_H1_on_unit_S2', 'eigenvalue_index': 1,
            'powers': list(map(str, ss)), 'coefficients': list(map(str, v)),
            'trial_space': 'one_real_cos_m1_even_radial_trial',
            'strong_domain': 's=0_or_s>3/2', 'statistics': stats,
            'matrix_digest': digest([[[str(x) for x in row] for row in mat] for mat in (M,A,R)]),
            'gap': gap, 'temple_lower': str(temple), 'm1_lower': str(m1),
            'other_sector_lowers': others, 'angular_tail_lower': str(angular),
            'source': deepcopy(source), 'source_digest': digest(source),
            'lower': str(lower), 'upper': str(upper), 'exact_width': str(upper-lower),
            'tolerance': str(tol), 'status': 'target_met' if upper-lower <= tol else 'certified_open',
            'full_infinite_space_covered': True, 'function_approximation_error': '0',
            'full_residual_not_projected_residual': True,
            'formal_proof_assistant_checked': False}


def verify(cert, expected_function=None, expected_mean_zero=None, expected_tolerance=None):
    try:
        if not isinstance(cert, dict):
            return False
        if cert.get('format') == direct.FULL:
            normalize(cert['function'])
            return direct.verify_full(cert, expected_function, expected_mean_zero, expected_tolerance)
        if cert.get('format') != FORMAT:
            return False
        if expected_function is not None and normalize(expected_function) != cert['function']:
            return False
        if expected_mean_zero is not None and expected_mean_zero is not True:
            return False
        if expected_tolerance is not None and rat(expected_tolerance) != F(cert['tolerance']):
            return False
        rebuilt = certificate(cert['function'], cert['powers'], cert['coefficients'],
                              cert['source'], cert['tolerance'], cert['gap']['policy'])
        return canonical(rebuilt) == canonical(cert)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


def _matvec(A, v):
    return [sum((a*b for a,b in zip(row,v)), F(0)) for row in A]


def _congruence(A, C):
    return [[sum((C[k][i]*A[k][l]*C[l][j]
                  for k in range(len(A)) for l in range(len(A))), F(0))
             for j in range(len(C[0]))] for i in range(len(C[0]))]


def residual_diagnosis(function, powers, coefficients, transform=None):
    """Pythagorean split of the complete residual, using the mass inverse."""
    ss = enriched.exponents(powers)
    stats = statistics(function,ss,coefficients)
    v = list(map(rat,coefficients))
    M,A,_ = trial_matrices(function,ss)
    C = transform or [[F(i == j) for j in range(len(ss))] for i in range(len(ss))]
    mu, mass = F(stats['rayleigh']), F(stats['mass'])
    Av,Mv = _matvec(A,v),_matvec(M,v)
    d = [x-mu*y for x,y in zip(Av,Mv)]
    projected = [sum((row[j]*x for row,x in zip(C,d)),F(0)) for j in range(len(C[0]))]
    gram = _congruence(M,C)
    solved = enriched._ldl_solver(gram)(projected)
    inside = sum((x*y for x,y in zip(projected,solved)),F(0))/mass
    outside = F(stats['residual_squared'])-inside
    if inside < 0 or outside < 0:
        raise ArithmeticError('Residual projection decomposition failed')
    return {'inside_squared': str(inside), 'outside_squared': str(outside),
            'complete_squared': stats['residual_squared'],
            'projector_dimension': len(C[0]),
            'action': 'candidate_coefficients' if inside > outside/16 else 'basis_or_gap',
            'action_is_a_heuristic_not_a_proof': True}


def cancellation_transform(function, powers):
    q = normalize(function)
    ss = enriched.exponents(powers)
    delta = 2+F(q['exponent'])
    if F(0) not in ss or delta not in ss:
        raise ValueError('Leading cancellation needs both 0 and 2+alpha in the basis')
    ratio = F(q['amplitude'])/(delta*(delta-1))
    eliminated, constant = ss.index(delta), ss.index(F(0))
    free = [i for i in range(len(ss)) if i != eliminated]
    C = [[F(i == j) for j in free] for i in range(len(ss))]
    C[eliminated][free.index(constant)] = ratio
    return C


def _inverse_proposal(function, powers, precision_bits, iterations, transform=None,
                      adaptive_iterations=False):
    ss = enriched.exponents(powers)
    direct.integer(precision_bits, 32, 256, 'precision_bits')
    direct.integer(iterations, 0, 64, 'iterations')
    M, A, _ = trial_matrices(function, ss)
    C = transform or [[F(i == j) for j in range(len(ss))] for i in range(len(ss))]
    mass, form = _congruence(M,C), _congruence(A,C)
    linear_solve = enriched._ldl_solver(form)
    v = [F(i == 0) for i in range(len(C[0]))]
    denominator = 2**precision_bits
    used = 0
    for _ in range(iterations):
        x = linear_solve(_matvec(mass,v))
        scale = max(abs(t) for t in x)
        v = [F(round(t/scale*denominator), denominator) for t in x]
        used += 1
        if adaptive_iterations:
            diagnostic = residual_diagnosis(function,ss,_matvec(C,v),C)
            if F(diagnostic['inside_squared']) <= F(diagnostic['outside_squared'])/16:
                break
    result = _matvec(C,v)
    return {'powers': list(map(str,ss)), 'coefficients': list(map(str,result)),
            'precision_bits': precision_bits, 'iterations': used,
            'iteration_budget': iterations, 'free_coefficients': len(v),
            'diagnosis': residual_diagnosis(function,ss,result,C)}


def proposal(function, powers, precision_bits=112, iterations=14, **options):
    q = normalize(function)
    cancellation = options.pop('leading_cancellation', False)
    adaptive = options.pop('adaptive_iterations', False)
    if type(cancellation) is not bool or type(adaptive) is not bool or options:
        raise ValueError('Unknown or malformed proposal options')
    ss = enriched.exponents(powers)
    available = 0 in ss and 2+F(q['exponent']) in ss
    if (cancellation and available) or adaptive:
        result = _inverse_proposal(q, ss, precision_bits, iterations,
                                   cancellation_transform(q,ss) if cancellation and available else None,
                                   adaptive)
    else:
        result = enriched.propose_trial(ss, precision_bits, iterations)
    result['leading_cancellation_applied'] = cancellation and available
    if 'diagnosis' not in result:
        result['diagnosis'] = residual_diagnosis(q,ss,result['coefficients'])
    return result


def residual_terms(function, powers, coefficients):
    """Merge the complete residual function, including outside-space powers."""
    stats = statistics(function, powers, coefficients)
    mu, terms = F(stats['rayleigh']), {}
    for s, value in zip(enriched.exponents(powers), coefficients):
        v = rat(value)
        for r, c in enriched.action_terms(s).items():
            terms[r] = terms.get(r, F(0))+v*c
        terms[s] = terms.get(s, F(0))-mu*v
    return {r: c for r, c in terms.items() if c}


def residual_frontier(function, powers, coefficients, count=3):
    direct.integer(count, 1, 8, 'candidate_budget')
    ss = enriched.exponents(powers)
    candidates = []
    for r, c in residual_terms(function, ss, coefficients).items():
        s = r+2
        if s not in ss and F(3,2) < s <= 64:
            # A proposal ranking, not a bound on the residual norm.
            candidates.append((abs(c)/(s*(s-1)), s))
    return tuple(s for _, s in sorted(candidates, key=lambda x: (-x[0],x[1]))[:count])


def greedy_addition(function, powers, coefficients, source, precision_bits=112,
                    iterations=14, candidate_budget=3, tolerance='1/100000000',
                    proposal_options=None):
    """Choose among a bounded residual-driven frontier by certified width."""
    ss = enriched.exponents(powers)
    if len(ss) >= 16:
        raise ValueError('No basis budget remains')
    trials = []
    best = None
    for s in residual_frontier(function, ss, coefficients, candidate_budget):
        new = tuple(sorted(ss+(s,)))
        p = proposal(function, new, precision_bits, iterations, **(proposal_options or {}))
        try:
            c = certificate(function, new, p['coefficients'], source, tolerance)
        except ValueError as error:
            if 'Temple is not applicable' not in str(error):
                raise
            trials.append({'power': str(s), 'status': 'temple_gap_unavailable'})
            continue
        score = (F(c['exact_width']), F(c['statistics']['residual_squared']))
        trials.append({'power': str(s), 'status': c['status'],
                       'exact_width': c['exact_width']})
        if best is None or score < best[0]:
            best = (score, p)
    if best is None:
        raise ValueError('Residual frontier has no certifiable candidate')
    return dict(best[1], frontier_trials=trials)


def solve(function, tolerance='1/100000000', max_terms=16,
          precision_bits=112, iterations=14, **budget):
    q = normalize(function)
    direct.integer(max_terms, 1, 16, 'max_terms')
    allowed = {'source_modes', 'source_bits', 'basis_policy', 'candidate_budget',
               'leading_cancellation', 'adaptive_iterations'}
    if set(budget)-allowed:
        raise ValueError('Unknown solve budget fields')
    start = perf_counter()
    source = source_certificate(q, budget.get('source_modes', 4), budget.get('source_bits', 28))
    best = direct.full_certificate(q, True, source['sectors'], tolerance)
    counts = sorted(set([n for n in (1, 3, 6, 10, 16) if n <= max_terms]+[max_terms]))
    policy = budget.get('basis_policy', 'nested')
    if policy not in ('nested', 'residual'):
        raise ValueError('Unknown basis policy')
    if policy == 'residual':
        counts = list(range(1, max_terms+1))
    attempts = []
    previous = None
    proposal_options = {'leading_cancellation': budget.get('leading_cancellation', False),
                        'adaptive_iterations': budget.get('adaptive_iterations', False)}
    for n in counts:
        begin = perf_counter()
        if policy == 'residual' and previous is not None:
            p = greedy_addition(q, previous['powers'], previous['coefficients'], source,
                                precision_bits, iterations, budget.get('candidate_budget', 3), tolerance,
                                proposal_options)
            powers = p['powers']
        else:
            powers = generated_powers(q, n)
            p = proposal(q, powers, precision_bits, iterations, **proposal_options)
        previous = p
        try:
            c = certificate(q, powers, p['coefficients'], source, tolerance)
        except ValueError as error:
            if 'Temple is not applicable' not in str(error):
                raise
            attempts.append({'terms': n, 'powers': list(map(str,powers)),
                             'status': 'temple_gap_unavailable', 'reason': str(error)})
            continue
        attempts.append({'terms': n, 'powers': list(map(str,powers)),
                         'statistics': c['statistics'], 'exact_width': c['exact_width'],
                         'frontier_trials': p.get('frontier_trials', []),
                         'diagnosis': p['diagnosis'], 'iterations_used': p['iterations'],
                         'elapsed_seconds': perf_counter()-begin,
                         'status': c['status']})
        if F(c['exact_width']) < F(best['exact_width']):
            best = c
        if best['status'] == 'target_met':
            break
    if not verify(best, q, True, tolerance):
        raise ArithmeticError('Final exact certificate replay failed')
    return {'certificate': best, 'attempts': attempts, 'status': best['status'],
            'elapsed_seconds': perf_counter()-start,
            'policy': policy,
            'human_designed_policy': True,
            'no_claim_of_universal_or_basis_global_optimality': True}
