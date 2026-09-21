"""Exact error certificates for the axisymmetric unit-sphere ground state.

Floating-point routines propose functions. Only rational arithmetic accepts
certificates. No optional dependency. See ERROR_ANALYSIS.md for the proofs.
"""
from fractions import Fraction as F
from math import comb, cos, pi, sqrt, copysign, isfinite
import spectral_certifier as base


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(a, b):
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
                 for i in range(max(len(a), len(b)))])


def scale(a, c):
    return trim([c*x for x in a])


def mul(a, b):
    out = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out)


def derivative(a):
    return [i*a[i] for i in range(1, len(a))] or [F(0)]


def evaluate(a, t):
    out = 0
    for c in reversed(a):
        out = out*t + c
    return out


def bernstein_coefficients(p):
    """p(2*x-1) in Bernstein basis on [0,1], all exact."""
    p = trim([F(x) for x in p])
    n = len(p)-1
    power = [sum((p[j]*comb(j,k)*2**k*(-1)**(j-k)
                  for j in range(k, n+1)), F(0)) for k in range(n+1)]
    return [sum((power[j]*F(comb(i,j), comb(n,j))
                 for j in range(i+1)), F(0)) for i in range(n+1)]


def split_bernstein(b):
    rows = [b]
    while len(rows[-1]) > 1:
        rows.append([(a+c)/2 for a,c in zip(rows[-1], rows[-1][1:])])
    return [r[0] for r in rows], [r[-1] for r in reversed(rows)]


def polynomial_range(p, depth=3):
    if type(depth) is not int or not 0 <= depth <= 7:
        raise ValueError('Bernstein subdivision depth must be in 0..7')
    leaves = [bernstein_coefficients(p)]
    for _ in range(depth):
        leaves = [child for b in leaves for child in split_bernstein(b)]
    return min(min(b) for b in leaves), max(max(b) for b in leaves)


def local_energy(q, s):
    """H exp(s) / exp(s); s is a rational polynomial in t."""
    ds = derivative(s)
    return add(add(q, mul([0, 2], ds)),
               scale(mul([1, 0, -1], add(derivative(ds), mul(ds, ds))), -1))


def vector_input(raw, maximum=32):
    if not isinstance(raw, list) or not 1 <= len(raw) <= maximum:
        raise ValueError('Invalid coefficient count')
    values = [base.rational(x) for x in raw]
    if any(abs(x) > 10**6 or x.denominator > 2**60 for x in values):
        raise ValueError('Candidate coefficient is too large')
    return values


def enclosure_fields(lo, hi):
    return {'lower': str(lo), 'upper': str(hi), 'exact_width': str(hi-lo),
            'decimal_enclosure': [base.outward_decimal(lo), base.outward_decimal(hi, upper=True)]}


def barta_certificate(q, s, depth=3):
    q = base.potential([str(x) for x in q])
    s = vector_input([str(x) for x in s], 17)
    # The constant in s only normalizes exp(s); canonicalize it away.
    s[0] = F(0)
    s = trim(s)
    energy = local_energy(q, s)
    lo, hi = polynomial_range(energy, depth)
    return {'method': 'positive_exponential_local_energy_v1', 'model': base.MODEL,
            'scope': base.SCOPE, 'eigenvalue_index': 1,
            'q_coefficients': [str(x) for x in q], 's_power_coefficients': [str(x) for x in s],
            'subdivision_depth': depth, 'local_energy_coefficients': [str(x) for x in energy],
            'formal_assistant_checked': False, **enclosure_fields(lo, hi)}


def verify_barta(cert):
    try:
        if not isinstance(cert, dict):
            return False
        expected = barta_certificate(cert['q_coefficients'], cert['s_power_coefficients'], cert['subdivision_depth'])
        return cert == expected
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return False


def banded_inertia(matrix):
    """Exact LDL in O(N*w^2) rational operations for fixed bandwidth w.

    A zero pivot falls back to the general exact 1x1/2x2 implementation.
    No pivot is replaced by a numerical epsilon.
    """
    a = [[F(x) for x in row] for row in matrix]
    n = len(a)
    if any(len(row) != n for row in a) or any(a[i][j] != a[j][i] for i in range(n) for j in range(n)):
        raise ValueError('Symmetric square matrix required')
    width = max((j-i for i in range(n) for j in range(i+1,n) if a[i][j]), default=0)
    neg, pos = 0, 0
    for k in range(n):
        d = a[k][k]
        if not d:
            return base.inertia(matrix)
        neg += d < 0
        pos += d > 0
        end = min(n, k+width+1)
        for i in range(k+1, end):
            for j in range(i, end):
                a[i][j] -= a[i][k]*a[k][j]/d
                a[j][i] = a[i][j]
    return [neg, 0, pos]


def certify_banded(q, k=1, modes=8, bits=36):
    """Same brackets, arithmetic and certificate as the dense baseline."""
    base.check_sizes(k, modes, bits)
    q = base.potential([str(x) for x in q])
    data = base.assemble(q, modes)
    def upper(x):
        count = banded_inertia(base.shifted(data, x))
        return count[0]+count[1] >= k
    def lower(x):
        return x < data['beta'] and banded_inertia(base.shifted(data,x,tail=True))[0] < k
    lo, hi = data['qmin']-1, data['qmax']+(k-1)*k+1
    for _ in range(bits):
        mid = (lo+hi)/2
        if upper(mid): hi = mid
        else: lo = mid
    up = hi
    lo, hi, distance = data['qmin']-1, up, F(1)
    while not lower(lo):
        distance *= 2
        lo = data['qmin']-distance
    for _ in range(bits):
        mid = (lo+hi)/2
        if lower(mid): lo = mid
        else: hi = mid
    cert = {'model': base.MODEL, 'q_coefficients': [str(x) for x in q],
            'eigenvalue_index': k, 'modes': modes, 'bisection_steps': bits,
            'evidence': {'analytic_potential_lower': str(data['qmin']),
                         'tail_lower': str(data['beta']),
                         'lower_schur_inertia': banded_inertia(base.shifted(data,lo,tail=True)),
                         'upper_finite_inertia': banded_inertia(base.shifted(data,up))},
            'scope': base.SCOPE, 'verification': 'exact_rational_with_analytic_tail',
            'formal_assistant_checked': False, **enclosure_fields(lo,up)}
    if not base.verify(cert):
        raise ArithmeticError('Dense independent inertia check rejected certificate')
    return cert


def jacobi_lowest(matrix):
    """Untrusted floating proposal; cyclic symmetric Jacobi rotations."""
    a = [list(map(float,row)) for row in matrix]
    n = len(a)
    v = [[float(i==j) for j in range(n)] for i in range(n)]
    for _ in range(50):
        off = max((abs(a[i][j]) for i in range(n) for j in range(i+1,n)), default=0)
        if off < 1e-14:
            break
        for p in range(n):
            for q in range(p+1,n):
                apq = a[p][q]
                if abs(apq) < 1e-16:
                    continue
                tau = (a[q][q]-a[p][p])/(2*apq)
                t = copysign(1.0,tau)/(abs(tau)+sqrt(1+tau*tau))
                c = 1/sqrt(1+t*t)
                s = t*c
                app, aqq = a[p][p], a[q][q]
                for i in range(n):
                    if i not in (p,q):
                        aip, aiq = a[i][p], a[i][q]
                        a[i][p] = a[p][i] = c*aip-s*aiq
                        a[i][q] = a[q][i] = s*aip+c*aiq
                    vip, viq = v[i][p], v[i][q]
                    v[i][p], v[i][q] = c*vip-s*viq, s*vip+c*viq
                a[p][p], a[q][q], a[p][q], a[q][p] = app-t*apq, aqq+t*apq, 0.0, 0.0
    index = min(range(n), key=lambda i:a[i][i])
    return [v[i][index] for i in range(n)]


def polynomial_trial(q, modes, even=False):
    base.check_sizes(1,modes)
    data = base.assemble(q,modes)
    mass = [float(x) for x in data['M']]
    indices = list(range(0,modes,2)) if even else list(range(modes))
    if even and any(q[i] for i in range(1,len(q),2)):
        raise ValueError('Even symmetry requires an even potential')
    h = [[float(data['A'][i][j])/sqrt(mass[i]*mass[j]) for j in indices] for i in indices]
    v = jacobi_lowest(h)
    c = [0.0]*modes
    for j,i in enumerate(indices): c[i] = v[j]/sqrt(mass[i])
    norm = max(abs(x) for x in c)
    return [F(round(x/norm*2**44),2**44) for x in c]


def residual_statistics(q, c):
    """Full residual, including every mode created by q*c."""
    n, degree = len(c), len(q)-1
    h = [F(0)]*(n+degree)
    for j,value in enumerate(c):
        h[j] += j*(j+1)*value
        for i,w in base.q_times_p(q,j).items():
            h[i] += value*w
    mass = [F(2,2*i+1) for i in range(len(h))]
    norm = sum((mass[i]*c[i]**2 for i in range(n)),F(0))
    if not norm:
        raise ValueError('Zero trial function')
    mu = sum((mass[i]*c[i]*h[i] for i in range(n)),F(0))/norm
    finite = sum((mass[i]*(h[i]-mu*c[i])**2 for i in range(n)),F(0))/norm
    tail = sum((mass[i]*h[i]**2 for i in range(n,len(h))),F(0))/norm
    return mu, finite, tail


def even_data(q,modes):
    base.check_sizes(1,modes)
    if modes < 3 or any(q[i] for i in range(1,len(q),2)):
        raise ValueError('Even-sector separation requires an even potential and at least 3 ambient modes')
    full = base.assemble(q,modes)
    indices = list(range(0,modes,2))
    data = {key:[[full[key][i][j] for j in indices] for i in indices] for key in ['A','G']}
    data['M'] = [full['M'][i] for i in indices]
    data['qmin'],data['qmax'] = full['qmin'],full['qmax']
    first_omitted = modes if modes%2 == 0 else modes+1
    data['beta'] = first_omitted*(first_omitted+1)+full['qmin']
    return data


def even_separation(q,modes,bits=32):
    data = even_data(q,modes)
    def holds(x):
        return x < data['beta'] and banded_inertia(base.shifted(data,x,tail=True))[0] < 2
    distance = F(1)
    lo,hi = data['qmin']-distance, data['qmax']+7
    while not holds(lo):
        distance *= 2
        lo = data['qmin']-distance
    for _ in range(bits):
        mid = (lo+hi)/2
        if holds(mid): lo = mid
        else: hi = mid
    return {'type':'even_schur_second_eigenvalue','ambient_modes':modes,'lower':str(lo),
            'tail_lower':str(data['beta']),
            'schur_inertia':banded_inertia(base.shifted(data,lo,tail=True))}


def gap_value(q, proof):
    if not isinstance(proof,dict):
        raise ValueError('Missing spectral separation proof')
    if proof.get('type') in ('potential_minmax','even_potential_minmax'):
        even = proof['type'] == 'even_potential_minmax'
        if even and any(q[i] for i in range(1,len(q),2)):
            raise ValueError('Even potential required')
        depth = proof['subdivision_depth']
        qlo, _ = polynomial_range(q,depth)
        expected = {'type':proof['type'],'subdivision_depth':depth,'potential_lower':str(qlo)}
        if proof != expected:
            raise ValueError('Invalid potential range proof')
        return (6 if even else 2)+qlo
    if proof.get('type') == 'even_schur_second_eigenvalue':
        data = even_data(q,proof['ambient_modes'])
        gamma = base.rational(proof['lower'])
        if gamma >= data['beta']: raise ValueError('Even tail is not positive')
        count = base.inertia(base.shifted(data,gamma,tail=True))
        if count[0] >= 2: raise ValueError('Even second eigenvalue is not separated')
        expected = {'type':'even_schur_second_eigenvalue','ambient_modes':proof['ambient_modes'],
                    'lower':str(gamma),'tail_lower':str(data['beta']),'schur_inertia':count}
        if proof != expected: raise ValueError('Invalid even-sector proof')
        return gamma
    if proof.get('type') == 'schur_second_eigenvalue':
        cert = proof['certificate']
        if not base.verify(cert) or cert['eigenvalue_index'] != 2 or base.potential(cert['q_coefficients']) != q:
            raise ValueError('Invalid second-eigenvalue certificate')
        return F(cert['lower'])
    raise ValueError('Unknown separation proof')


def temple_certificate(q, c, gap_proof):
    q = base.potential([str(x) for x in q])
    c = vector_input([str(x) for x in c])
    if isinstance(gap_proof,dict) and str(gap_proof.get('type','')).startswith('even_'):
        if any(c[i] for i in range(1,len(c),2)) or any(q[i] for i in range(1,len(q),2)):
            raise ValueError('An even-sector certificate requires both potential and trial function to be even')
    mu, finite, tail = residual_statistics(q,c)
    gamma = gap_value(q,gap_proof)
    if gamma <= mu:
        raise ValueError('No certified separation: gamma must exceed the Rayleigh quotient')
    correction = (finite+tail)/(gamma-mu)
    return {'method':'full_residual_temple_v1','model':base.MODEL,'scope':base.SCOPE,
            'eigenvalue_index':1,'q_coefficients':[str(x) for x in q],
            'trial_legendre_coefficients':[str(x) for x in c], 'gap_proof':gap_proof,
            'gap_lower':str(gamma),'rayleigh_quotient':str(mu),
            'finite_residual_squared':str(finite),'tail_residual_squared':str(tail),
            'residual_squared':str(finite+tail),
            'excited_weight_upper':str(min(F(1),(finite+tail)/(gamma-mu)**2)),
            'formal_assistant_checked':False,
            **enclosure_fields(mu-correction,mu)}


def verify_temple(cert):
    try:
        if not isinstance(cert,dict): return False
        expected = temple_certificate(cert['q_coefficients'],cert['trial_legendre_coefficients'],cert['gap_proof'])
        return expected == cert
    except (KeyError,TypeError,ValueError,ZeroDivisionError):
        return False


def certify_temple(q, modes=8, symmetry=False):
    q = base.potential([str(x) for x in q])
    even = symmetry and not any(q[i] for i in range(1,len(q),2))
    c = polynomial_trial(q,modes,even=even)
    mu,_,_ = residual_statistics(q,c)
    depth = 3
    qlo,_ = polynomial_range(q,depth)
    gap = {'type':'even_potential_minmax' if even else 'potential_minmax',
           'subdivision_depth':depth,'potential_lower':str(qlo)}
    if even and 6+qlo <= mu:
        gap = even_separation(q,max(3,modes))
    elif not even and 2+qlo <= mu:
        # A modestly accurate bound suffices when the true gap is appreciable.
        second = certify_banded(q,2,max(2,modes),16)
        if F(second['lower']) <= mu:
            # Failure may be a search-precision problem, not a true cluster.
            second = certify_banded(q,2,max(2,modes),40)
        gap = {'type':'schur_second_eigenvalue','certificate':second}
    cert = temple_certificate(q,c,gap)
    if not verify_temple(cert): raise ArithmeticError('Temple certificate rejected')
    return cert


def solve_linear(a,b):
    a = [list(row)+[value] for row,value in zip(a,b)]
    n = len(b)
    for k in range(n):
        pivot = max(range(k,n),key=lambda i:abs(a[i][k]))
        a[k],a[pivot] = a[pivot],a[k]
        if abs(a[k][k]) < 1e-18: raise ValueError('Singular candidate system')
        for i in range(k+1,n):
            factor = a[i][k]/a[k][k]
            for j in range(k,n+1): a[i][j] -= factor*a[k][j]
    x = [0.0]*n
    for i in reversed(range(n)):
        x[i] = (a[i][n]-sum(a[i][j]*x[j] for j in range(i+1,n)))/a[i][i]
    return x


def chebyshev_table(x, degree):
    t, d, dd = [1.0,x],[0.0,1.0],[0.0,0.0]
    for j in range(2,degree+1):
        t.append(2*x*t[-1]-t[-2])
        d.append(2*t[-2]+2*x*d[-1]-d[-2])
        dd.append(4*d[-2]+2*x*dd[-1]-dd[-2])
    return t,d,dd


def exponential_trial(q, degree, continuation=8):
    """Untrusted damped Newton collocation for s=sum b_j*T_j, psi=exp(s)."""
    if type(degree) is not int or not 1 <= degree <= 16:
        raise ValueError('Exponential log degree must be in 1..16')
    if type(continuation) is not int or not 1 <= continuation <= 32:
        raise ValueError('Continuation steps must be in 1..32')
    qfloat = list(map(float,q))
    points = [cos(pi*j/degree) for j in range(degree+1)]
    tables = [chebyshev_table(x,degree) for x in points]
    state = [0.0]*(degree+1)  # b_1,...,b_degree, energy
    def system(y, strength):
        rows, values = [],[]
        for x,(_,d,dd) in zip(points,tables):
            ds = sum(y[j-1]*d[j] for j in range(1,degree+1))
            dds = sum(y[j-1]*dd[j] for j in range(1,degree+1))
            values.append(-(1-x*x)*(dds+ds*ds)+2*x*ds+strength*evaluate(qfloat,x)-y[-1])
            rows.append([-(1-x*x)*(dd[j]+2*ds*d[j])+2*x*d[j] for j in range(1,degree+1)]+[-1.0])
        return rows,values
    for stage in range(1,continuation+1):
        strength = stage/continuation
        for _ in range(35):
            jac,res = system(state,strength)
            score = max(map(abs,res))
            if score < 1e-12: break
            delta = solve_linear(jac,[-x for x in res])
            step = 1.0
            accepted = False
            for _ in range(20):
                candidate = [x+step*dx for x,dx in zip(state,delta)]
                if all(isfinite(x) for x in candidate) and max(map(abs,system(candidate,strength)[1])) < score:
                    state, accepted = candidate, True
                    break
                step /= 2
            if not accepted: break
    cheb = [[F(1)],[F(0),F(1)]]
    for j in range(2,degree+1):
        cheb.append(add(mul([0,2],cheb[-1]),scale(cheb[-2],-1)))
    power = [F(0)]
    for j in range(1,degree+1):
        # Rationalization belongs to the proposal. The verifier rechecks it.
        coefficient = F(round(state[j-1]*2**44),2**44)
        power = add(power,scale(cheb[j],coefficient))
    power[0] = F(0)
    return trim(power)


def certify_exponential(q, degree=6, depth=3):
    q = base.potential([str(x) for x in q])
    s = exponential_trial(q,degree)
    cert = barta_certificate(q,s,depth)
    if not verify_barta(cert): raise ArithmeticError('Local-energy certificate rejected')
    return cert


def verify(cert):
    if not isinstance(cert,dict): return False
    if cert.get('method') == 'positive_exponential_local_energy_v1': return verify_barta(cert)
    if cert.get('method') == 'full_residual_temple_v1': return verify_temple(cert)
    return base.verify(cert)
