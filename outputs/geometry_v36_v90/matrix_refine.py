"""Budgeted refinement of V33 matrix envelopes; standard library only.

Search is untrusted. Every accepted bound is replayed by dual_v33/family_v29.
The objective is the fixed Poisson matrix envelope, not spectral sharpness.
"""
from fractions import Fraction as F
from time import perf_counter
import math
import algebra as a
import error_bounds as e
import dual_v33 as dual
import family_v29 as family
import float_sdp


def _eye(d):
    return [[F(i == j) for j in range(d)] for i in range(d)]


def _add(x, y):
    return [[u + v for u, v in zip(r, s)] for r, s in zip(x, y)]


def _scale(x, value):
    return [[value * z for z in row] for row in x]


def congruence(t, x):
    return a.matmul(a.matmul(t, x), a.transpose(t))


def transform_problem(directions, cost, transform):
    """v_new=T v, H=T G T^T, C_new=T^-T C T^-1, all exactly.

    An arbitrary invertible rational T is supported; it need not be symmetric.
    Return transformed directions/cost and the inverse used for pullback.
    """
    dirs, _, c = dual.setup(directions, cost)
    d = len(dirs)
    t = [list(map(F, row)) for row in transform]
    if len(t) != d or any(len(row) != d for row in t):
        raise ValueError('Coordinate transform size')
    inverse = a.inverse(t)
    new_dirs = []
    for row in t:
        p = [F(0)]
        for coefficient, direction in zip(row, dirs):
            p = e.add(p, e.scale(direction, coefficient))
        new_dirs.append(p)
    new_cost = congruence(a.transpose(inverse), c)
    # Internal coordinates can have larger rational denominators than the
    # user-facing polynomial parser permits. They are exact search data, not a
    # new input problem; only the ORIGINAL directions enter final certificates.
    dual.matrix(new_cost, d, True)
    return new_dirs, new_cost, inverse


def psd_factors(matrix):
    """Exact LDL rank-one factors (pivot, vector), including singular PSDs."""
    n = len(matrix)
    rest = dual.matrix(matrix, n)
    answer = []
    for k in range(n):
        pivot = rest[k][k]
        if not pivot:
            continue  # PSD and zero diagonal imply a zero row/column.
        vector = [F(0)] * n
        vector[k] = F(1)
        for i in range(k + 1, n):
            vector[i] = rest[i][k] / pivot
        answer.append((pivot, vector))
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                rest[i][j] -= pivot * vector[i] * vector[j]
    return answer


def preconditioner(derivatives):
    """Exact coefficient decorrelation, then rational powers-of-two scaling."""
    d = len(derivatives)
    width = max(map(len, derivatives))
    rows = [p + [F(0)] * (width - len(p)) for p in derivatives]
    gram = a.matmul(rows, a.transpose(rows))
    factors = psd_factors(gram)
    if len(factors) != d:
        raise ValueError('Degenerate derivative coordinates')
    lower = a.transpose([vector for _, vector in factors])
    diagonal = []
    for pivot, _ in factors:
        # Integer bit lengths avoid converting huge rational coefficients to float.
        exponent = (pivot.numerator.bit_length() - pivot.denominator.bit_length()) // 2
        diagonal.append(F(1, 2 ** exponent) if exponent >= 0 else F(2 ** (-exponent)))
    inverse = a.inverse(lower)
    return [[diagonal[i] * z for z in row] for i, row in enumerate(inverse)]


def effective_samples(derivatives, samples, dominance=True):
    """Discard duplicates and EXACT Loewner-dominated finite constraints.

    This only simplifies the finite proposal. Full-interval proof is unchanged.
    """
    samples = sorted(set(map(F, samples)))
    if not samples or any(abs(t) > 1 for t in samples):
        raise ValueError('Finite sample domain')
    unique = {}
    for t in samples:
        matrix = dual.at(derivatives, t)
        key = tuple(tuple(row) for row in matrix)
        unique.setdefault(key, (t, matrix))
    values = list(unique.values())
    if dominance:
        values = [(t, matrix) for index, (t, matrix) in enumerate(values)
                  if not any(index != j and a.inertia(_add(other, _scale(matrix, -1)))[0] == 0
                             for j, (_, other) in enumerate(values))]
    # All equal constraints became one before the strict-index domination pass.
    return values


def merge_atoms(directions, atoms):
    """Merge atoms at identical A(t), preserving both matrix sum and objective."""
    dirs, derivatives, _ = dual.setup(directions)
    d = len(dirs)
    grouped = {}
    for atom in atoms:
        t = F(atom['t'])
        matrix = dual.matrix(atom['matrix'], d)
        key = tuple(tuple(row) for row in dual.at(derivatives, t))
        if key not in grouped:
            grouped[key] = [t, [[F(0)] * d for _ in range(d)]]
        grouped[key][1] = _add(grouped[key][1], matrix)
    return [{'t': str(t), 'matrix': dual.encode(matrix)} for t, matrix in grouped.values()]


def _total(atoms, d):
    answer = [[F(0)] * d for _ in range(d)]
    for atom in atoms:
        answer = _add(answer, [list(map(F, row)) for row in atom['matrix']])
    return answer


def fill_slack(directions, cost, certificate, samples):
    """Assign each exact PSD rank-one slack component to its best finite node.

    The sum of dual atoms becomes C. Each added objective term is nonnegative.
    """
    dirs, derivatives, c = dual.setup(directions, cost)
    if not dual.verify_lower(certificate) or certificate['directions'] != [list(map(str, p)) for p in dirs] or certificate['cost'] != dual.encode(c):
        raise ValueError('Slack fill certificate binding')
    atoms = merge_atoms(dirs, certificate['atoms'])
    nodes = sorted(set(F(t) for t in samples) | {F(atom['t']) for atom in atoms})
    if not nodes or any(abs(t) > 1 for t in nodes):
        raise ValueError('Slack fill sample domain')
    if len(atoms) + len(c) > 128:
        nodes = sorted({F(atom['t']) for atom in atoms})
    candidates = [(t, dual.at(derivatives, t)) for t in nodes]
    for pivot, vector in psd_factors(certificate['cost_slack']):
        matrix = [[pivot * x * y for y in vector] for x in vector]
        t, _ = max(candidates, key=lambda pair: dual.inner(matrix, pair[1]))
        atoms.append({'t': str(t), 'matrix': dual.encode(matrix)})
    return dual.lower_certificate(dirs, c, merge_atoms(dirs, atoms))


def repair_dual(directions, cost, atoms, error=F(1, 10**10), samples=None, fill=True):
    """Maximize a common atom multiplier in BOTH directions with exact PSD tests.

    The bisection stops on an absolute dual-objective error, not on an
    unrelated multiplier error. No floating eigenvalue is accepted as proof.
    """
    dirs, derivatives, c = dual.setup(directions, cost)
    error = F(error)
    if error <= 0:
        raise ValueError('Positive dual repair error required')
    atoms = merge_atoms(dirs, atoms)
    if not 1 <= len(atoms) <= 128:
        raise ValueError('Dual atom budget')
    d = len(dirs)
    total = _total(atoms, d)
    low_value = sum((dual.inner([list(map(F, row)) for row in atom['matrix']],
                                dual.at(derivatives, atom['t'])) for atom in atoms), F(0))
    if not any(z for row in total for z in row):
        certificate = dual.lower_certificate(dirs, c, atoms)
    else:
        def feasible(multiplier):
            return a.inertia(_add(c, _scale(total, -multiplier)))[0] == 0
        # Taking traces gives beta <= trace(C)/trace(total), an exact bracket
        # invariant under arbitrarily large/small common proposal scaling.
        lo = F(0)
        hi = sum(c[i][i] for i in range(d)) / sum(total[i][i] for i in range(d))
        if feasible(hi):
            lo = hi
        for _ in range(256):
            if (hi - lo) * low_value <= error:
                break
            mid = (lo + hi) / 2
            if feasible(mid):
                lo = mid
            else:
                hi = mid
        if (hi - lo) * low_value > error:
            raise ValueError('Dual repair objective-error budget not reached')
        scaled = [{'t': atom['t'], 'matrix': dual.encode(_scale([list(map(F, row)) for row in atom['matrix']], lo))} for atom in atoms]
        certificate = dual.lower_certificate(dirs, c, scaled)
    return fill_slack(dirs, c, certificate, samples or [atom['t'] for atom in atoms]) if fill else certificate


def sparsify_dual(directions, cost, certificate, error, samples=None):
    """Drop PSD atoms only within an exact objective-loss budget, then refill."""
    dirs, derivatives, c = dual.setup(directions, cost)
    if not dual.verify_lower(certificate) or certificate['directions'] != [list(map(str, p)) for p in dirs] or certificate['cost'] != dual.encode(c):
        raise ValueError('Sparse certificate binding')
    error = F(error)
    if error < 0:
        raise ValueError('Nonnegative sparse error required')
    atoms = merge_atoms(dirs, certificate['atoms'])
    scored = sorted(((dual.inner([list(map(F, row)) for row in atom['matrix']], dual.at(derivatives, atom['t'])), i) for i, atom in enumerate(atoms)))
    removed, loss = set(), F(0)
    for contribution, i in scored:
        if len(removed) + 1 < len(atoms) and loss + contribution <= error:
            removed.add(i)
            loss += contribution
    kept = [atom for i, atom in enumerate(atoms) if i not in removed]
    candidate = dual.lower_certificate(dirs, c, kept)
    return fill_slack(dirs, c, candidate, samples or [atom['t'] for atom in kept]), loss


def reconstruct_psd(matrix, denominator):
    """Round LDL square-root factors, so quantization remains PSD by construction."""
    d = len(matrix)
    result = [[F(0)] * d for _ in range(d)]
    for pivot, vector in psd_factors(matrix):
        root = math.sqrt(float(pivot))
        rounded = [F(round(root * float(x) * denominator), denominator) for x in vector]
        result = _add(result, [[x * y for y in rounded] for x in rounded])
    return result


def error_allocation(tolerance, gap, objective_scale, round_index=0):
    """Absolute search errors and a relative range error, in original cost units."""
    tolerance, gap, objective_scale = map(F, (tolerance, gap, objective_scale))
    if tolerance <= 0 or gap < 0 or objective_scale <= 0:
        raise ValueError('Invalid error allocation')
    work = max(tolerance, min(16 * tolerance, gap / 8))
    work = max(tolerance, work / (2 ** min(round_index, 8)))
    return {'search': work / 8, 'repair': tolerance / 128,
            'quantize': tolerance / 128, 'sparse': tolerance / 256,
            'range': max(F(1, 10**12), min(F(1, 10**4), tolerance / (32 * objective_scale)))}


def warm_solve(samples, cost, final_mu, initial=None, max_steps=280, initial_mu=None):
    """Untrusted normalized log-barrier search with feasible warm starts.

    Only the small linear/Cholesky kernels are reused from float_sdp; the old
    solver is unchanged. A return flag never controls certificate acceptance.
    """
    d = len(cost)
    indices = [(i, j) for i in range(d) for j in range(i, d)]
    n = len(indices)
    def unpack(x):
        g = [[0.0] * d for _ in range(d)]
        for value, (i, j) in zip(x, indices):
            g[i][j] = g[j][i] = value
        return g
    objective = [cost[i][j] * (1 if i == j else 2) for i, j in indices]
    alpha = max(1.0, 2 * max(sum(aa[i][i] for i in range(d)) for aa in samples))
    g = [[float(initial[i][j]) for j in range(d)] for i in range(d)] if initial is not None else [[alpha * (i == j) for j in range(d)] for i in range(d)]
    used_warm = initial is not None
    # An inherited optimum lies on the boundary. A ridge on the final barrier
    # scale can make the first Newton Hessian nearly rank deficient, so start
    # farther inside before following the new finite problem's central path.
    ridge = max(math.sqrt(final_mu), initial_mu or 0.0, 1e-8) if used_warm else max(final_mu, 1e-12)
    for _ in range(64):
        candidate = [[g[i][j] + (ridge if i == j else 0.0) for j in range(d)] for i in range(d)]
        try:
            for aa in samples:
                float_sdp.inv_logdet([[candidate[i][j] - aa[i][j] for j in range(d)] for i in range(d)])
            g = candidate
            break
        except ValueError:
            ridge *= 4
    else:
        raise ValueError('Warm candidate has no finite positive repair')
    x = [g[i][j] for i, j in indices]
    def evaluate(x, mu, derivatives=False):
        g = unpack(x)
        value = sum(u * v for u, v in zip(objective, x))
        gradient = objective[:]
        hessian = [[0.0] * n for _ in range(n)]
        for aa in samples:
            inv, logdet = float_sdp.inv_logdet([[g[i][j] - aa[i][j] for j in range(d)] for i in range(d)])
            value -= mu * logdet
            if derivatives:
                for k, (i, j) in enumerate(indices):
                    gradient[k] -= mu * inv[i][j] * (1 if i == j else 2)
                    for l, (u, v) in enumerate(indices):
                        if i == j and u == v:
                            z = inv[i][u] ** 2
                        elif i == j:
                            z = 2 * inv[i][u] * inv[i][v]
                        elif u == v:
                            z = 2 * inv[i][u] * inv[j][u]
                        else:
                            z = 2 * (inv[i][u] * inv[j][v] + inv[i][v] * inv[j][u])
                        hessian[k][l] += mu * z
        return value, gradient, hessian
    mu = max(final_mu, initial_mu if initial_mu is not None else 0.05)
    steps, stages = 0, 0
    converged = False
    while steps < max_steps:
        stages += 1
        for _ in range(48):
            if steps >= max_steps:
                break
            value, gradient, hessian = evaluate(x, mu, True)
            direction = float_sdp.linear(hessian, [-z for z in gradient])
            slope = sum(u * v for u, v in zip(gradient, direction))
            steps += 1
            if not math.isfinite(slope) or slope > 1e-15:
                raise ValueError('Nonfinite or uphill Newton candidate')
            if -slope / 2 < max(1e-22, mu * 1e-9):
                break
            rate = 1.0
            for _ in range(64):
                candidate = [z + rate * v for z, v in zip(x, direction)]
                try:
                    new_value = evaluate(candidate, mu)[0]
                except ValueError:
                    new_value = float('inf')
                if new_value <= value + 0.01 * rate * slope:
                    x = candidate
                    break
                rate *= 0.5
            else:
                break
        if mu <= final_mu:
            converged = True
            break
        mu = max(final_mu, mu / 6)
    return unpack(x), mu, {'newton_steps': steps, 'barrier_stages': stages,
                          'reached_final_barrier': converged, 'warm_start': used_warm}


def _template(directions, weight, epsilon, leaves, evaluator):
    dirs, w, _, _, polynomial = family.assemble(directions, weight)
    proof = evaluator(polynomial, epsilon, leaves)
    # This binding also rejects an oracle returning a certificate for another q.
    return family.build(dirs, w, proof, 'poisson_matrix_v30')


def _upper(template, cost):
    return dual.inner(cost, [list(map(F, row)) for row in template['quadratic_bound']])


def _denominator(cost, error, maximum):
    # Entry error <= 1/(2D); a strict finite repair adds at most O(d/D).
    target = max(F(2**20), 8 * len(cost) * sum(abs(z) for row in cost for z in row) / error)
    d = 1
    while d < target and d < maximum:
        d *= 2
    return min(d, maximum)


def _quantize_candidate(raw, constraints, denominator):
    d = len(raw)
    if any(not math.isfinite(z) for row in raw for z in row):
        raise ValueError('Nonfinite proposed matrix')
    g = [[F(round((raw[i][j] + raw[j][i]) * denominator / 2), denominator) for j in range(d)] for i in range(d)]
    shift = F(0)
    for _ in range(64):
        candidate = _add(g, _scale(_eye(d), shift))
        if all(a.inertia(_add(candidate, _scale(aa, -1))) == [0, 0, d] for _, aa in constraints):
            return candidate, shift
        shift = F(1, denominator) if not shift else 2 * shift
    raise ValueError('Finite rational candidate repair budget')


def _separation_nodes(template, limit=6):
    """Use witness plus local polynomial stationary proposals to avoid stalling.

    Numerical root locations are only finite constraints; the old all-interval
    certificate remains the sole primal authority.
    """
    proof = template['range_proof']
    q, start, end, coordinate = family.reduced_polynomial(proof['polynomial'])
    derivative = e.derivative(q)
    locations = [F(proof['witness']), start, end]
    # Bracket sign changes; include both signs after the t^2 lift.
    grid = [float(start + (end - start) * F(i, 96)) for i in range(97)]
    poly = list(map(float, derivative))
    def value(x):
        result = 0.0
        for coefficient in reversed(poly):
            result = result * x + coefficient
        return result
    for lo, hi in zip(grid, grid[1:]):
        left, right = value(lo), value(hi)
        if left * right < 0:
            for _ in range(40):
                mid = (lo + hi) / 2
                if value(lo) * value(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            locations.append(F(round((lo + hi) / 2 * 2**32), 2**32))
    locations = sorted(set(locations), key=lambda t: e.evaluate(q, t), reverse=True)[:limit]
    nodes = set()
    for location in locations:
        if coordinate == 't_squared':
            location = F(round(math.sqrt(max(0.0, float(location))) * 2**32), 2**32)
        location = max(F(-1), min(F(1), location))
        nodes.update((location, -location))
    return nodes


def update_samples(samples, proposed, derivatives, lower, maximum):
    """Keep new separators even when full; replace least productive old nodes.

    Deleting a finite proposal constraint does not delete any retained global
    certificate. Dual objective mass ranks old support; the exact lower proof
    continues to include all its atoms, whether or not searched next round.
    """
    proposed = set(proposed)
    score = {F(atom['t']): dual.inner([list(map(F, row)) for row in atom['matrix']],
                                     dual.at(derivatives, atom['t'])) for atom in lower['atoms']}
    chosen = set(sorted(proposed)[:maximum])
    old = sorted(set(samples) - chosen, key=lambda t: (score.get(t, F(0)), -abs(t)), reverse=True)
    chosen.update(old[:maximum - len(chosen)])
    return chosen


_DEFAULTS = {'rounds': 8, 'newton_steps': 280, 'range_leaves': 128,
             'max_samples': 48, 'max_denominator': 2**44,
             'precondition': True, 'dominance': True, 'sparsify': True}


def _budget(raw):
    if raw is None:
        return dict(_DEFAULTS)
    if type(raw) is int:
        raw = {'rounds': raw}
    if not isinstance(raw, dict) or set(raw) - set(_DEFAULTS):
        raise ValueError('Unknown matrix refinement budget')
    result = dict(_DEFAULTS, **raw)
    for key, lo, hi in [('rounds', 0, 24), ('newton_steps', 1, 1200),
                        ('range_leaves', 1, 256), ('max_samples', 9, 96),
                        ('max_denominator', 2**10, 2**50)]:
        if type(result[key]) is not int or not lo <= result[key] <= hi:
            raise ValueError('Invalid ' + key + ' budget')
    if any(type(result[key]) is not bool for key in ('precondition', 'dominance', 'sparsify')):
        raise ValueError('Refinement switches must be bool')
    return result


def refine(directions, cost=None, tolerance=F(1, 10**6), seed_result=None,
           budget=None, *, range_evaluator=None, candidate_solver=None):
    """Return exact global bounds, retaining every verified seed bound.

    budget is a round count or a dict of _DEFAULTS overrides. With rounds=0,
    only safe initialization/seed replay occurs. range_evaluator(q,epsilon,
    max_leaves) must return a legacy family-compatible range certificate.
    candidate_solver follows warm_solve's signature; its output is untrusted.
    """
    started = perf_counter()
    tolerance = F(tolerance)
    if not F(1, 10**10) <= tolerance <= F(1, 10):
        raise ValueError('Refinement tolerance: 1e-10 to 1e-1')
    limits = _budget(budget)
    evaluator = family.maximum if range_evaluator is None else range_evaluator
    solver = warm_solve if candidate_solver is None else candidate_solver
    dirs, derivatives, c = dual.setup(directions, cost)
    d = len(dirs)
    transform = preconditioner(derivatives) if limits['precondition'] else _eye(d)
    _, transformed_cost, inverse = transform_problem(dirs, c, transform)
    transformed_derivatives = []
    for row in transform:
        polynomial = [F(0)]
        for coefficient, derivative in zip(row, derivatives):
            polynomial = e.add(polynomial, e.scale(derivative, coefficient))
        transformed_derivatives.append(polynomial)
    # A scalar normalization keeps float gradients near one without changing C.
    objective_scale = max(sum(transformed_cost[i][i] for i in range(d)), F(1, 10**100))
    normalized_cost = _scale(transformed_cost, 1 / objective_scale)
    samples = {F(i, 4) for i in range(-4, 5)}
    trace = []
    if seed_result is not None:
        if not isinstance(seed_result, dict):
            raise ValueError('Seed must contain a verified certificate')
        seed = seed_result.get('certificate', seed_result)
        if not dual.verify(seed) or seed['dual']['directions'] != [list(map(str, p)) for p in dirs] or seed['dual']['cost'] != dual.encode(c):
            raise ValueError('Seed certificate does not match the requested problem')
        best_template, best_dual = seed['template'], seed['dual']
        samples.update(F(atom['t']) for atom in best_dual['atoms'])
        # Untrusted cached samples never acquire authority from the seed proof.
        warm_state = seed_result.get('warm_state', {})
        cached = warm_state.get('samples', []) if isinstance(warm_state, dict) else []
        if not isinstance(cached, list):
            cached = []
        for value in cached[:limits['max_samples']]:
            try:
                t = F(value)
                if -1 <= t <= 1:
                    samples.add(t)
            except (ValueError, TypeError, ZeroDivisionError):
                pass
    else:
        fallback_weight = congruence(inverse, _eye(d))
        try:
            best_template = _template(dirs, fallback_weight, F(1, 10**8), limits['range_leaves'], evaluator)
        except (ValueError, TypeError, KeyError, IndexError, ArithmeticError):
            best_template = _template(dirs, fallback_weight, F(1, 10**8), limits['range_leaves'], family.maximum)
            trace.append({'outcome': 'range_evaluator_fallback'})
        zero = {'t': '0', 'matrix': dual.encode([[F(0)] * d for _ in range(d)])}
        best_dual = dual.lower_certificate(dirs, c, [zero])
        best_dual = fill_slack(dirs, c, best_dual, samples)
    best_upper = _upper(best_template, c)
    total_steps = 0
    last_mu = None
    for index in range(limits['rounds']):
        gap = best_upper - F(best_dual['lower_bound'])
        if gap <= tolerance:
            break
        if len(samples) > limits['max_samples']:
            # Keep seed support first; finite pruning never discards its proof.
            support = {F(atom['t']) for atom in best_dual['atoms']}
            samples = set(sorted(support)[:limits['max_samples']])
            samples.update(t for t in sorted({F(i, 4) for i in range(-4, 5)} - samples)[:max(0, limits['max_samples'] - len(samples))])
        constraints = effective_samples(transformed_derivatives, samples, limits['dominance'])
        allocation = error_allocation(tolerance, gap, max(best_upper, F(1)), index)
        denominator = _denominator(transformed_cost, allocation['quantize'], limits['max_denominator'])
        final_mu = float(allocation['search'] / (d * len(constraints) * objective_scale))
        previous_g = [list(map(F, row)) for row in best_template['quadratic_bound']]
        warm = congruence(transform, previous_g)
        stats = {}
        try:
            raw, actual_mu, stats = solver(
                [[[float(z) for z in row] for row in aa] for _, aa in constraints],
                [[float(z) for z in row] for row in normalized_cost], final_mu,
                initial=warm, max_steps=limits['newton_steps'],
                initial_mu=max(final_mu, min(0.05, (last_mu or float(gap / objective_scale)) * 4)))
            if not math.isfinite(actual_mu) or actual_mu <= 0:
                raise ValueError('Invalid candidate barrier scale')
            total_steps += stats.get('newton_steps', 0)
            g, shift = _quantize_candidate(raw, constraints, denominator)
            # Rational PSD atoms are reconstructed in normalized coordinates,
            # then mapped by T^T Z T to the ORIGINAL objective/directions.
            atoms = []
            rational_mu = F(str(actual_mu)) * objective_scale
            for t, aa in constraints:
                inv = a.inverse(_add(g, _scale(aa, -1)))
                y = reconstruct_psd(_scale(inv, rational_mu), denominator)
                atoms.append({'t': str(t), 'matrix': dual.encode(congruence(a.transpose(transform), y))})
            lower = repair_dual(dirs, c, atoms, allocation['repair'], samples)
            sparse_loss = F(0)
            if limits['sparsify']:
                sparse, sparse_loss = sparsify_dual(dirs, c, lower, allocation['sparse'], samples)
                if F(sparse['lower_bound']) >= F(lower['lower_bound']):
                    lower = sparse
            if F(lower['lower_bound']) > F(best_dual['lower_bound']):
                best_dual = lower
            last_mu = actual_mu
            # This is the only place a numerical candidate can become an upper
            # bound: the full interval proof must pass the old family verifier.
            weight = congruence(inverse, g)
            objective = max(F(1), dual.inner(c, weight))
            range_error = error_allocation(tolerance, gap, objective, index)['range']
            template = _template(dirs, weight, range_error, limits['range_leaves'], evaluator)
            upper = _upper(template, c)
            if upper < best_upper:
                best_upper, best_template = upper, template
            proposed_nodes = _separation_nodes(template, min(6, limits['max_samples'] // 2))
            before = len(samples)
            new_count = len(proposed_nodes - samples)
            samples = update_samples(samples, proposed_nodes, derivatives, lower, limits['max_samples'])
            current = dual.certificate(best_template, best_dual, tolerance)
            trace.append({'round': index + 1, 'outcome': 'certified',
                          'constraint_count': len(constraints), 'sample_count': len(samples),
                          'new_samples': new_count, 'dropped_samples': max(0, before + new_count - len(samples)), 'lower': current['lower_bound'],
                          'upper': current['upper_bound'], 'gap': current['gap'],
                          'range_status': template['range_proof']['status'],
                          'range_error': str(range_error), 'finite_diagonal_shift': str(shift),
                          'denominator': denominator, 'dual_atoms': len(best_dual['atoms']),
                          'sparse_loss_budget_used': str(sparse_loss),
                          'error_allocation': {k: str(v) for k, v in allocation.items()}, **stats})
        except (ValueError, OverflowError, ZeroDivisionError, TypeError, KeyError, IndexError) as exc:
            trace.append({'round': index + 1, 'outcome': 'proposal_failed', 'reason': str(exc), **stats})
            break
    certificate = dual.certificate(best_template, best_dual, tolerance)
    return {'algorithm': 'matrix_refine', 'certificate': certificate,
            'statistics': {'exchange_rounds': sum('round' in row for row in trace),
                           'newton_steps': total_steps, 'sample_count': len(samples),
                           'dual_atoms': len(best_dual['atoms']), 'model_calls': 0,
                           'seed_reused': seed_result is not None},
            'conditioning': {'transform': dual.encode(transform),
                             'transformed_cost': dual.encode(transformed_cost),
                             'objective_scale': str(objective_scale)},
            'warm_state': {'samples': list(map(str, sorted(samples)))},
            'trace': trace, 'seconds': perf_counter() - started}


def verify(result):
    """Verify just the mathematical certificate; search metadata is untrusted."""
    try:
        return dual.verify(result.get('certificate', result))
    except (AttributeError, TypeError):
        return False


if __name__ == '__main__':
    import argparse
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directions', default='[[0,1],[0,0,1]]', help='JSON polynomial coefficient arrays')
    parser.add_argument('--cost', default=None, help='JSON rational positive definite cost')
    parser.add_argument('--tolerance', default='1/1000000')
    parser.add_argument('--rounds', type=int, default=8)
    parser.add_argument('--seed', help='Previous result JSON for verified continuation')
    args = parser.parse_args()
    seed = None
    if args.seed:
        with open(args.seed, encoding='utf-8') as handle:
            seed = json.load(handle)
    result = refine(json.loads(args.directions), json.loads(args.cost) if args.cost else None,
                    F(args.tolerance), seed, args.rounds)
    print(json.dumps(result, indent=2))
