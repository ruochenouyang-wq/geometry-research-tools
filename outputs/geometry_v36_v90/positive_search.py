"""Certified fixed-potential optimization over positive exponential functions.

All acceptance uses Fractions. Floating Newton/soft-min steps only propose.
An atom certificate bounds the BEST BARTA LOWER in the specified ansatz;
it is deliberately not a spectral upper bound. See PROOF_POSITIVE_SEARCH.md.
"""
from fractions import Fraction as F
from math import exp, isfinite, log, sqrt
import algebra as a
import error_bounds as e
import family_v29 as family


MODEL = 'unit_sphere_positive_exponential_fixed_potential_v1'
SCOPE = 'ground energy of -Delta+q(t) on the full unit two-sphere, t=x_3'
ANSATZ_SCOPE = 'sup_theta inf_{-1<=t<=1} H exp(sum theta_i b_i)/exp(sum theta_i b_i), for the stated fixed q and basis'
LIMITATION = 'ansatz_upper bounds the best Barta lower in this basis; it is not a spectral upper bound'
MAX_DIMENSION = 6
MAX_DEGREE = 6


def rational(x):
    if isinstance(x, F):
        z = x
    else:
        z = e.base.rational(x)
    if abs(z) > 10**6 or z.denominator > 2**60:
        raise ValueError('Rational input budget exceeded')
    return z


def polynomial(raw, max_degree=MAX_DEGREE):
    if not isinstance(raw, list) or not 1 <= len(raw) <= max_degree+1:
        raise ValueError('Polynomial degree budget exceeded')
    return e.trim([rational(x) for x in raw])


def potential(raw):
    q = polynomial(raw)
    return e.base.potential(list(map(str, q)))


def encode(x):
    return list(map(str, x))


def setup(q, basis):
    q = potential(q)
    if not isinstance(basis, list) or not 1 <= len(basis) <= MAX_DIMENSION:
        raise ValueError('One to six basis polynomials required')
    out = []
    for raw in basis:
        b = polynomial(raw)
        b[0] = F(0)  # Multiplicative normalization of psi has no effect.
        out.append(e.trim(b))
    rows = [b+[F(0)]*(MAX_DEGREE+1-len(b)) for b in out]
    if len(a.independent_rows(rows)) != len(out):
        raise ValueError('Basis must be independent modulo constants')
    return q, out


def monomial_basis(q, degree=6, symmetry=True):
    q = potential(q)
    if type(degree) is not int or not 1 <= degree <= MAX_DEGREE:
        raise ValueError('Logarithmic degree must lie in 1..6')
    even = symmetry and not any(q[1::2])
    powers = list(range(2 if even else 1, degree+1, 2 if even else 1))
    if not powers:
        raise ValueError('Even symmetry needs degree at least two')
    return [[F(0)]*j+[F(1)] for j in powers]


def logarithm(basis, theta):
    if not isinstance(theta, list) or len(theta) != len(basis):
        raise ValueError('Parameter dimension mismatch')
    theta = [rational(x) for x in theta]
    s = [F(0)]
    for b, value in zip(basis, theta):
        s = e.add(s, e.scale(b, value))
    return theta, s


def local_polynomial(q, basis, theta):
    q, basis = setup(q, basis)
    theta, s = logarithm(basis, theta)
    return e.local_energy(q, s)


class BernsteinBackend:
    """Hook contract: name, minimum(p,tolerance,max_leaves), verify_minimum(p,proof).

    verify_minimum must independently return exact (inf_lower, inf_upper),
    or raise ValueError. The caller must supply a trusted backend verifier
    again when verifying a certificate made with a non-default backend.
    """
    name = 'rational_adaptive_bernstein_minimum_v1'

    def minimum(self, p, tolerance, max_leaves):
        return family.maximum(e.scale(p, -1), tolerance, max_leaves)

    def verify_minimum(self, p, proof):
        if not family.verify_range(proof) or list(map(F, proof['polynomial'])) != e.scale(p, -1):
            raise ValueError('Range proof is not bound to the local energy')
        return -F(proof['maximum_upper']), -F(proof['maximum_lower'])


DEFAULT_RANGE_BACKEND = BernsteinBackend()


def _backend(backend):
    backend = DEFAULT_RANGE_BACKEND if backend is None else backend
    if not isinstance(backend.name, str) or not backend.name:
        raise ValueError('Range backend requires a stable name')
    return backend


def lower_certificate(q, basis, theta, tolerance=F(1, 10**6), max_leaves=128,
                      range_proof=None, range_backend=None):
    q, basis = setup(q, basis)
    theta, s = logarithm(basis, theta)
    p = e.local_energy(q, s)
    backend = _backend(range_backend)
    if range_proof is None:
        range_proof = backend.minimum(p, tolerance, max_leaves)
    lo, hi = backend.verify_minimum(p, range_proof)
    if not isinstance(lo, F) or not isinstance(hi, F) or lo > hi:
        raise ValueError('Backend must return ordered exact Fraction bounds')
    return {'method': 'positive_exponential_lower_v1', 'model': MODEL, 'scope': SCOPE,
            'q_coefficients': encode(q), 'basis': [encode(b) for b in basis],
            'theta': encode(theta), 's_coefficients': encode(s),
            'local_energy_coefficients': encode(p), 'range_backend': backend.name,
            'range_proof': range_proof, 'spectral_lower': str(lo),
            'candidate_infimum_upper': str(hi), 'range_gap': str(hi-lo),
            'formal_assistant_checked': False}


def verify_lower(cert, range_backend=None):
    try:
        backend = _backend(range_backend)
        if cert['range_backend'] != backend.name:
            return False
        return cert == lower_certificate(cert['q_coefficients'], cert['basis'], cert['theta'],
                                         range_proof=cert['range_proof'], range_backend=backend)
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError, AttributeError):
        return False


def atom_statistics(q, basis, atoms):
    q, basis = setup(q, basis)
    d = len(basis)
    if not isinstance(atoms, list) or not 1 <= len(atoms) <= 128:
        raise ValueError('One to 128 probability atoms required')
    deriv = [e.derivative(b) for b in basis]
    lap = [family.laplacian(b) for b in basis]
    mass, c = F(0), F(0)
    linear = [F(0)]*d
    gram = [[F(0)]*d for _ in range(d)]
    canonical = []
    for atom in atoms:
        if not isinstance(atom, dict) or set(atom) != {'t', 'weight'}:
            raise ValueError('Probability atom requires t and weight only')
        t, weight = rational(atom['t']), rational(atom['weight'])
        if not -1 <= t <= 1 or weight < 0:
            raise ValueError('Invalid probability atom')
        mass += weight
        c += weight*e.evaluate(q, t)
        v = [e.evaluate(b, t) for b in deriv]
        for i in range(d):
            linear[i] += weight*e.evaluate(lap[i], t)
            for j in range(d):
                gram[i][j] += weight*(1-t*t)*v[i]*v[j]
        canonical.append({'t': str(t), 'weight': str(weight)})
    if mass != 1:
        raise ValueError('Probability mass must equal one exactly')
    return q, basis, canonical, c, linear, gram


def ansatz_upper_certificate(q, basis, atoms):
    """max_theta sum mu_k E_theta(t_k), an ansatz upper only.

    Singular positive semidefinite Gram matrices are accepted iff 2 A z=l
    is consistent. Otherwise the quadratic is unbounded above and no finite
    upper is emitted. No numerical inverse or epsilon pivot is accepted.
    """
    q, basis, atoms, c, linear, gram = atom_statistics(q, basis, atoms)
    if a.inertia(gram)[0]:
        raise ArithmeticError('Probability Gram matrix cannot be indefinite')
    try:
        z, pivots = a.rectangular([[2*x for x in row] for row in gram], linear, len(basis))
    except ValueError as exc:
        raise ValueError('Atom quadratic unbounded above: linear term in Gram nullspace') from exc
    upper = c + sum((x*y for x, y in zip(linear, z)), F(0))/2
    return {'method': 'positive_exponential_ansatz_upper_v1', 'model': MODEL,
            'scope': ANSATZ_SCOPE, 'limitation': LIMITATION,
            'q_coefficients': encode(q), 'basis': [encode(b) for b in basis],
            'atoms': atoms, 'constant': str(c), 'linear': encode(linear),
            'gram': [encode(row) for row in gram], 'gram_rank': len(pivots),
            'quadratic_maximizer': encode(z), 'ansatz_upper': str(upper),
            'formal_assistant_checked': False}


def verify_ansatz_upper(cert):
    try:
        return cert == ansatz_upper_certificate(cert['q_coefficients'], cert['basis'], cert['atoms'])
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError):
        return False


def mean(p):
    return sum((x/F(j+1) for j, x in enumerate(p) if j % 2 == 0), F(0))


def rayleigh_certificate(q, trial):
    """A genuine full-sphere spectral upper, using a rational polynomial trial."""
    q, trial = potential(q), polynomial(trial)
    norm = mean(e.mul(trial, trial))
    if norm <= 0:
        raise ValueError('Nonzero polynomial trial required')
    derivative = e.derivative(trial)
    kinetic = mean(e.mul([1, 0, -1], e.mul(derivative, derivative)))
    potential_energy = mean(e.mul(q, e.mul(trial, trial)))
    return {'method': 'positive_search_rayleigh_upper_v1', 'model': MODEL, 'scope': SCOPE,
            'q_coefficients': encode(q), 'trial_power_coefficients': encode(trial),
            'norm_squared': str(norm), 'kinetic_numerator': str(kinetic),
            'potential_numerator': str(potential_energy),
            'spectral_upper': str((kinetic+potential_energy)/norm),
            'formal_assistant_checked': False}


def verify_rayleigh(cert):
    try:
        return cert == rayleigh_certificate(cert['q_coefficients'], cert['trial_power_coefficients'])
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError):
        return False


def parity_certificate(q, s):
    """Exact symmetrization witness: E_even - (E_s(t)+E_s(-t))/2 >= 0."""
    q, s = potential(q), polynomial(s)
    if any(q[1::2]):
        raise ValueError('Parity reduction requires an even potential')
    s[0] = F(0)
    s = e.trim(s)
    even = e.trim([x if j % 2 == 0 else F(0) for j, x in enumerate(s)])
    odd = e.trim([x if j % 2 else F(0) for j, x in enumerate(s)])
    original = e.local_energy(q, s)
    averaged = e.trim([x if j % 2 == 0 else F(0) for j, x in enumerate(original)])
    derivative = e.derivative(odd)
    gain = e.mul([1, 0, -1], e.mul(derivative, derivative))
    if e.add(e.local_energy(q, even), e.scale(averaged, -1)) != gain:
        raise ArithmeticError('Parity identity failed')
    return {'method': 'positive_exponential_even_projection_v1',
            'q_coefficients': encode(q), 's_coefficients': encode(s),
            'even_s_coefficients': encode(even), 'odd_s_derivative': encode(derivative),
            'gain_coefficients': encode(gain),
            'identity': 'E_even(t)-(E_s(t)+E_s(-t))/2=(1-t^2)*(s_odd_prime(t))^2',
            'conclusion': 'For even q, projecting log psi to its even part cannot decrease the exact Barta lower.',
            'formal_assistant_checked': False}


def verify_parity(cert):
    try:
        return cert == parity_certificate(cert['q_coefficients'], cert['s_coefficients'])
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError):
        return False


def propose_rayleigh(q, degree=6):
    q = potential(q)
    if type(degree) is not int or not 1 <= degree <= 6:
        raise ValueError('Rayleigh degree must lie in 1..6')
    coefficients = e.polynomial_trial(q, degree+1, even=not any(q[1::2]))
    legendre = [[F(1)], [F(0), F(1)]]
    for n in range(1, degree):
        legendre.append(e.scale(e.add(e.scale(e.mul([0, 1], legendre[-1]), 2*n+1),
                                     e.scale(legendre[-2], -n)), F(1, n+1)))
    p = [F(0)]
    for c, b in zip(coefficients, legendre):
        p = e.add(p, e.scale(b, c))
    # Keep the public rational input budget after the basis conversion.
    p = [F(round(float(x)*2**40), 2**40) for x in p]
    return rayleigh_certificate(q, p)


def optimization_certificate(lower, upper, tolerance=F(1, 10**4), rayleigh=None,
                             range_backend=None):
    tolerance = rational(tolerance)
    if not F(1, 10**12) <= tolerance <= 1:
        raise ValueError('Optimization tolerance must lie in [1e-12,1]')
    if not verify_lower(lower, range_backend) or not verify_ansatz_upper(upper):
        raise ValueError('Invalid lower or ansatz upper certificate')
    if any(lower[key] != upper[key] for key in ('q_coefficients', 'basis')):
        raise ValueError('Lower and ansatz upper use different fixed problems')
    lo, hi = F(lower['spectral_lower']), F(upper['ansatz_upper'])
    if lo > hi:
        raise ValueError('Inconsistent ansatz bounds')
    theta = list(map(F, lower['theta']))
    linear, gram = list(map(F, upper['linear'])), [list(map(F, row)) for row in upper['gram']]
    stationarity = [linear[i]-2*sum((gram[i][j]*theta[j] for j in range(len(theta))), F(0))
                    for i in range(len(theta))]
    average = F(upper['constant'])+sum((x*y for x, y in zip(linear, theta)), F(0))
    average -= sum((theta[i]*gram[i][j]*theta[j] for i in range(len(theta))
                    for j in range(len(theta))), F(0))
    if average < lo or average > hi:
        raise ValueError('Invalid finite-average weak duality decomposition')
    result = {'method': 'positive_exponential_global_optimization_v1', 'model': MODEL,
              'scope': ANSATZ_SCOPE, 'limitation': LIMITATION,
              'lower_certificate': lower, 'ansatz_upper_certificate': upper,
              'ansatz_lower': str(lo), 'ansatz_upper': str(hi), 'ansatz_gap': str(hi-lo),
              'finite_average_at_candidate': str(average),
              'averaged_gradient': encode(stationarity),
              'quadratic_stationarity_gap': str(hi-average),
              'average_contact_and_range_gap': str(average-lo),
              'tolerance': str(tolerance),
              'status': 'ansatz_gap_closed' if hi-lo <= tolerance else 'ansatz_gap_open',
              'formal_assistant_checked': False}
    if rayleigh is not None:
        if not verify_rayleigh(rayleigh) or rayleigh['q_coefficients'] != lower['q_coefficients']:
            raise ValueError('Rayleigh certificate uses a different potential')
        spectral_upper = F(rayleigh['spectral_upper'])
        if spectral_upper < lo:
            raise ValueError('Inconsistent spectral bounds')
        result.update(rayleigh_certificate=rayleigh, spectral_lower=str(lo),
                      spectral_upper=str(spectral_upper), spectral_gap=str(spectral_upper-lo))
    return result


def verify(cert, range_backend=None):
    try:
        return cert == optimization_certificate(cert['lower_certificate'], cert['ansatz_upper_certificate'],
                                               cert['tolerance'], cert.get('rayleigh_certificate'), range_backend)
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError):
        return False


def poisson_initialization(q, basis):
    """Exact coefficient least-squares projection; full Poisson when representable."""
    q, basis = setup(q, basis)
    r, _ = family.poisson(q)
    rows = [b+[F(0)]*(7-len(b)) for b in basis]
    target = r+[F(0)]*(7-len(r))
    gram = [[sum((x*y for x, y in zip(u, v)), F(0)) for v in rows] for u in rows]
    rhs = [sum((x*y for x, y in zip(u, target)), F(0)) for u in rows]
    return a.solve(gram, rhs)


def _float_system(q, basis, points):
    deriv = [e.derivative(b) for b in basis]
    lap = [family.laplacian(b) for b in basis]
    q0 = q[0]
    centered = e.add(q, [-q0])
    return [(float(e.evaluate(centered, t)),
             [float(e.evaluate(b, t)) for b in lap],
             [float(e.evaluate(b, t)) for b in deriv], float(1-t*t)) for t in points]


def _softmin(theta, rows, temperature, derivatives=True):
    values, gradients = [], []
    d = len(theta)
    for c, linear, deriv, w in rows:
        ds = sum(x*y for x, y in zip(theta, deriv))
        values.append(c+sum(x*y for x, y in zip(theta, linear))-w*ds*ds)
        gradients.append([linear[i]-2*w*ds*deriv[i] for i in range(d)])
    minimum = min(values)
    weights = [exp(-(x-minimum)/temperature) for x in values]
    total = sum(weights)
    weights = [x/total for x in weights]
    value = minimum-temperature*log(total)
    if not derivatives:
        return value, weights
    gradient = [sum(mu*g[i] for mu, g in zip(weights, gradients)) for i in range(d)]
    # Positive negative-Hessian: 2 E[w vv^T] + Cov(grad)/temperature.
    hessian = [[sum(mu*(2*row[3]*row[2][i]*row[2][j]+
                            (g[i]-gradient[i])*(g[j]-gradient[j])/temperature)
                         for mu, row, g in zip(weights, rows, gradients))
                for j in range(d)] for i in range(d)]
    return value, weights, gradient, hessian


def _propose(q, basis, points, initial, max_newton, final_temperature):
    rows = _float_system(q, basis, points)
    theta = list(map(float, initial))
    scale = max(1.0, sum(abs(float(x)) for x in q[1:]))
    temperature = scale/4
    stages = []
    while temperature > final_temperature:
        stages.append(temperature)
        temperature /= 4
    stages.append(final_temperature)
    for temperature in stages:
        for _ in range(max_newton):
            value, weights, gradient, hessian = _softmin(theta, rows, temperature)
            if max(map(abs, gradient)) < 1e-10:
                break
            regularizer = max(1e-12, max(hessian[i][i] for i in range(len(theta)))*1e-13)
            for i in range(len(theta)):
                hessian[i][i] += regularizer
            try:
                step = e.solve_linear(hessian, gradient)
            except ValueError:
                break
            directional = sum(x*y for x, y in zip(gradient, step))
            accepted = False
            length = 1.0
            for _ in range(30):
                candidate = [x+length*y for x, y in zip(theta, step)]
                if all(isfinite(x) and abs(x) <= 10**5 for x in candidate):
                    new_value, _ = _softmin(candidate, rows, temperature, False)
                    if new_value >= value+1e-4*length*directional:
                        theta, accepted = candidate, True
                        break
                length /= 2
            if not accepted or length*max(map(abs, step)) < 1e-13:
                break
    _, weights = _softmin(theta, rows, stages[-1], False)
    return [F(round(x*2**32), 2**32) for x in theta], weights


def _probability_atoms(points, weights):
    # Positive dyadic numerators preserve a nondegenerate grid Gram matrix;
    # exact normalization is essential for weak duality.
    numerators = [max(1, round(x*2**38)) for x in weights]
    total = sum(numerators)
    return [{'t': str(t), 'weight': str(F(n, total))} for t, n in zip(points, numerators)]


def _exchange_points(lower):
    """Bernstein minimum witnesses only propose new rational collocation nodes."""
    if lower['range_backend'] != DEFAULT_RANGE_BACKEND.name:
        return []
    proof = lower['range_proof']
    x = F(proof['witness'])
    if proof['coordinate'] == 't_squared':
        t = F(round(sqrt(float(x))*2**30), 2**30)
        return [-t, t]
    return [x]


def search(q, basis=None, degree=6, symmetry=True, tolerance=F(1, 10**4),
           range_tolerance=F(1, 10**6), max_leaves=128, grid_size=33,
           max_exchanges=4, max_newton=30, range_backend=None):
    """Bounded deterministic proposal + exact acceptance, no external packages.

    Retains zero and projected-Poisson starts, every accepted range proof, and
    rejected proposals. A budget exhaustion is an open gap, never optimality.
    """
    q = potential(q)
    basis = monomial_basis(q, degree, symmetry) if basis is None else basis
    q, basis = setup(q, basis)
    tolerance, range_tolerance = rational(tolerance), rational(range_tolerance)
    if not F(1, 10**12) <= tolerance <= 1 or not F(1, 10**12) <= range_tolerance <= 1:
        raise ValueError('Tolerance budget')
    if type(grid_size) is not int or not 9 <= grid_size <= 97:
        raise ValueError('Grid size must lie in 9..97')
    if type(max_exchanges) is not int or not 1 <= max_exchanges <= 8:
        raise ValueError('Exchange budget must lie in 1..8')
    if type(max_newton) is not int or not 1 <= max_newton <= 60:
        raise ValueError('Newton budget must lie in 1..60')
    backend = _backend(range_backend)
    points = [F(-1)+F(2*i, grid_size-1) for i in range(grid_size)]
    theta = poisson_initialization(q, basis)
    attempts, failures = [], []
    # Keep a moderate exact guard even if a numerical atom proposal becomes
    # almost singular. Every interior grid weight is strictly positive.
    uppers = [ansatz_upper_certificate(q, basis,
              [{'t': str(t), 'weight': str(F(1, len(points)))} for t in points])]
    for name, candidate in [('zero', [F(0)]*len(basis)), ('projected_poisson', theta)]:
        cert = lower_certificate(q, basis, candidate, range_tolerance, max_leaves, range_backend=backend)
        attempts.append({'proposal': name, 'certificate': cert})
    baseline = attempts[-1]['certificate']
    best = max((row['certificate'] for row in attempts), key=lambda c: F(c['spectral_lower']))
    for iteration in range(max_exchanges):
        try:
            theta, weights = _propose(q, basis, points, theta, max_newton, max(1e-10, float(tolerance)/200))
            candidate = lower_certificate(q, basis, theta, range_tolerance, max_leaves, range_backend=backend)
            attempts.append({'proposal': 'softmin_exchange_'+str(iteration), 'certificate': candidate})
            if F(candidate['spectral_lower']) > F(best['spectral_lower']):
                best = candidate
            upper = ansatz_upper_certificate(q, basis, _probability_atoms(points, weights))
            uppers.append(upper)
            if F(upper['ansatz_upper'])-F(best['spectral_lower']) <= tolerance:
                break
            added = False
            for point in _exchange_points(candidate):
                if point not in points and len(points) < 128:
                    points.append(point)
                    added = True
            if not added:
                break
        except (ValueError, ArithmeticError, OverflowError) as exc:
            failures.append({'iteration': iteration, 'reason': str(exc)})
            break
    upper = min(uppers, key=lambda c: F(c['ansatz_upper']))
    rayleigh = propose_rayleigh(q)
    cert = optimization_certificate(best, upper, tolerance, rayleigh, backend)
    if not verify(cert, backend):
        raise ArithmeticError('Exact acceptance failed')
    return {'method': 'positive_exponential_search_result_v1', 'certificate': cert,
            'poisson_baseline': baseline, 'attempts': attempts, 'failed_proposals': failures,
            'ansatz_upper_attempts': uppers,
            'improvement_over_projected_poisson': str(F(best['spectral_lower'])-F(baseline['spectral_lower'])),
            'proposal_warning': 'Floating grids and Newton termination are untrusted; only attached exact certificates establish bounds.'}


def shift_certificate(cert, shift, range_backend=None):
    """Energy translation q -> q+c; the same theta, basis, atoms and trial apply."""
    if not verify(cert, range_backend):
        raise ValueError('Invalid source optimization certificate')
    shift = rational(shift)
    old = cert['lower_certificate']
    q = list(map(F, old['q_coefficients']))
    q[0] += shift
    backend = _backend(range_backend)
    range_proof = None
    # Only the built-in backend object has the schema used by this fast path.
    # A subclass may override minimum/verify_minimum with a different proof.
    if backend is DEFAULT_RANGE_BACKEND:
        proof = old['range_proof']
        p = local_polynomial(q, old['basis'], old['theta'])
        range_proof = family.range_certificate(e.scale(p, -1), proof['intervals'],
                                                proof['witness'], proof['tolerance'])
    lower = lower_certificate(q, old['basis'], old['theta'], range_proof=range_proof,
                              range_backend=backend)
    upper = ansatz_upper_certificate(q, old['basis'], cert['ansatz_upper_certificate']['atoms'])
    rayleigh = None
    if 'rayleigh_certificate' in cert:
        rayleigh = rayleigh_certificate(q, cert['rayleigh_certificate']['trial_power_coefficients'])
    return optimization_certificate(lower, upper, cert['tolerance'], rayleigh, range_backend)


def demo():
    """Small reproducible research cases; prints exact certificates as JSON."""
    return {name: search(q, degree=degree) for name, q, degree in
            [('quadratic_6_degree2', [0, 0, 6], 2),
             ('quadratic_6_degree6', [0, 0, 6], 6),
             ('quadratic_20_degree6', [0, 0, 20], 6),
             ('mixed_quartic_degree4', [0, 1, 3, 0, 2], 4)]}


if __name__ == '__main__':
    import json
    print(json.dumps(demo(), indent=2))
