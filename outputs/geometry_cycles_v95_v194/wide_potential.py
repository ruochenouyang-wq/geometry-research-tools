"""V115–124: exact broad-band polynomial and bounded exponential sphere spectra.

The search/Schur structure is adapted from frozen V91–94. Input validation,
associated-Legendre multiplication, range evidence, certificate construction and
verification here are independent of its degree-six certificate validator.
Only pure rational inertia, enclosure formatting and Sturm evidence are reused.
"""
from copy import deepcopy
from math import factorial
from common import F, exact, same, base, arithmetic, exact_extrema

CHEAP_FORMAT = 'wide_cheap_sphere_ground_v115'
SECTOR_FORMAT = 'wide_sphere_sector_v121'
FULL_FORMAT = 'wide_full_sphere_ground_v122'
RANGE_FORMAT = 'wide_potential_range_v117_v119'
KERNEL_FORMAT = 'wide_schur_kernel_v120'
TAYLOR_FORMAT = 'exp_linear_taylor_v123'
EXP_FORMAT = 'exp_linear_sphere_ground_v124'
canonical_equal = same


def potential(q):
    """Degree <=24, |coefficient|<=10^6, denominator<=2^2048; exact only."""
    if not isinstance(q, (list, tuple)) or not 1 <= len(q) <= 25:
        raise ValueError('Require 1..25 rational power coefficients')
    q = [exact(x) for x in q]
    if any(abs(x) > 10**6 or x.denominator > 2**2048 for x in q):
        raise ValueError('Coefficient magnitude or denominator exceeds budget')
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q


def parameters(m, mean_zero, k, modes):
    if type(m) is not int or not 0 <= m <= 64 or type(mean_zero) is not bool:
        raise ValueError('Require integer 0<=m<=64 and boolean mean_zero')
    if type(k) is not int or type(modes) is not int or not 1 <= k <= modes <= 64:
        raise ValueError('Require integer 1<=k<=modes<=64')


def scope(mean_zero):
    return ('all_mean_zero_H1_functions_on_unit_S2' if mean_zero
            else 'all_H1_functions_on_unit_S2')


def basis_mass(m, l):
    if type(m) is not int or not 0 <= m <= 64 or type(l) is not int or not m <= l <= 256:
        raise ValueError('Invalid associated-Legendre basis indices')
    return F(2*factorial(l+m), (2*l+1)*factorial(l-m))


def _times_t(series, m):
    result = {}
    for l, c in series.items():
        result[l+1] = result.get(l+1, F(0)) + c*F(l-m+1, 2*l+1)
        if l > m:
            result[l-1] = result.get(l-1, F(0)) + c*F(l+m, 2*l+1)
    return {l:c for l,c in result.items() if c}


def q_times_basis(q, m, l):
    q = potential(q)
    basis_mass(m, l)
    if l+len(q)-1 > 256:
        raise ValueError('Multiplication exceeds supported basis budget')
    result, term = {}, {l:F(1)}
    for power, c in enumerate(q):
        for degree, value in term.items():
            result[degree] = result.get(degree, F(0))+c*value
        if power+1 < len(q):
            term = _times_t(term, m)
    return {l:c for l,c in result.items() if c}


def assemble(q, m=0, start=1, modes=8):
    parameters(m, False, 1, modes)
    basis_mass(m, start)
    if start+modes+24 > 256:
        raise ValueError('Assembly exceeds supported basis budget')
    q = potential(q)
    degrees = list(range(start, start+modes))
    mass = [basis_mass(m, l) for l in degrees]
    matrix = [[F(0) for _ in degrees] for _ in degrees]
    for j, l in enumerate(degrees):
        for degree, c in q_times_basis(q, m, l).items():
            if start <= degree < start+modes:
                matrix[degree-start][j] += mass[degree-start]*c
        matrix[j][j] += l*(l+1)*mass[j]
    if any(matrix[i][j] != matrix[j][i] for i in range(modes) for j in range(i)):
        raise ArithmeticError('Associated-Legendre matrix is not symmetric')
    return {'A':matrix, 'M':mass, 'degrees':degrees}


def _range(q, method, proof):
    if method == 'sturm':
        lo, hi = -F(proof['negative']['maximum_upper']), F(proof['positive']['maximum_upper'])
    elif method == 'bernstein':
        lo = min(F(cell['lower']) for cell in proof['cells'])
        hi = max(F(cell['upper']) for cell in proof['cells'])
    else:
        raise ValueError('Unknown range method')
    return {'format': RANGE_FORMAT, 'q_coefficients':[str(x) for x in q],
            'interval':['-1','1'], 'method':method, 'proof':proof,
            'lower':str(lo), 'upper':str(hi)}


def sturm_range(q, tolerance=F(1,10**10), max_nodes=64):
    q = potential(q)
    tolerance = exact(tolerance)
    positive = exact_extrema.maximize(q, (-1,1), tolerance, max_nodes)
    negative = exact_extrema.maximize([-x for x in q], (-1,1), tolerance, max_nodes)
    return _range(q, 'sturm', {'positive':positive, 'negative':negative})


def bernstein_range(q, depth=1):
    q = potential(q)
    if type(depth) is not int or not 0 <= depth <= 5:
        raise ValueError('Require integer 0<=depth<=5')
    cells = []
    for j in range(2**depth):
        left, right = -1+F(2*j, 2**depth), -1+F(2*(j+1), 2**depth)
        lo, hi = exact_extrema.bernstein_bounds(q, left, right)
        cells.append({'interval':[str(left), str(right)], 'lower':str(lo), 'upper':str(hi)})
    return _range(q, 'bernstein', {'depth':depth, 'cells':cells})


def potential_range(q, strategy='auto', tolerance=F(1,10**10), max_nodes=64):
    """Choose bounded-cost Bernstein for degree>8; Sturm otherwise.

    This deterministic cost heuristic makes no claim that the cheap bound is as
    tight as Sturm. Explicit strategies remain available for both methods.
    """
    q = potential(q)
    if strategy == 'auto':
        strategy = 'bernstein' if len(q)>9 else 'sturm'
    if strategy == 'bernstein':
        return bernstein_range(q)
    if strategy == 'sturm':
        return sturm_range(q, tolerance, max_nodes)
    raise ValueError('Unknown range strategy')


def verify_range(cert, expected_q=None):
    try:
        q = potential(cert['q_coefficients'])
        if expected_q is not None and q != potential(expected_q):
            return False
        proof, method = cert['proof'], cert['method']
        if method == 'sturm':
            for name, p in [('positive',q), ('negative',[-x for x in q])]:
                evidence = proof[name]
                if not exact_extrema.verify(evidence):
                    return False
                if potential(evidence['polynomial']) != p or evidence['interval'] != ['-1','1']:
                    return False
            expected = _range(q, method, {'positive':proof['positive'], 'negative':proof['negative']})
        elif method == 'bernstein':
            expected = bernstein_range(q, proof['depth'])
        else:
            return False
        return same(cert, expected)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


def reuse_range(q, evidence):
    """Replay supplied evidence, without extrema root isolation or range search."""
    if not verify_range(evidence, q):
        raise ValueError('Range evidence is invalid or bound to another potential')
    return deepcopy(evidence)


def range_bounds(q, evidence):
    evidence = reuse_range(q, evidence)
    return F(evidence['lower']), F(evidence['upper'])


def _cheap_ground_certificate(q, mean_zero, range_proof):
    q = potential(q)
    if type(mean_zero) is not bool:
        raise ValueError('mean_zero must be boolean')
    proof = reuse_range(q, range_proof)
    # These moments integrate the ORIGINAL q, rather than replacing it with a
    # sampled potential or the looser range bounds in a candidate numerator.
    integral_q = sum((c*F(2,j+1) for j,c in enumerate(q) if j%2==0),F(0))
    integral_q_t2 = sum((c*F(2,j+3) for j,c in enumerate(q) if j%2==0),F(0))
    data = [('z',1,0,[0,0,1],F(2,3),integral_q_t2),
            ('x',1,1,[1,0,-1],F(4,3),integral_q-integral_q_t2)]
    if not mean_zero:
        data.insert(0,('1',0,0,[1],F(2),integral_q))
    trials=[]
    for function,l,m,radial_squared,mass,numerator in data:
        trials.append({'sphere_function':function, 'degree_l':l, 'azimuth_m':m,
                       'radial_squared_coefficients':[str(v) for v in radial_squared],
                       'radial_mass':str(mass), 'kinetic_eigenvalue':l*(l+1),
                       'potential_integral':str(numerator),
                       'rayleigh_quotient':str(l*(l+1)+numerator/mass)})
    lower = (2 if mean_zero else 0)+F(proof['lower'])
    upper = min(F(trial['rayleigh_quotient']) for trial in trials)
    if lower > upper:
        raise ArithmeticError('Verified potential range contradicts a Rayleigh bound')
    return {'format':CHEAP_FORMAT, 'geometry':'unit_S2', 'scope':scope(mean_zero),
            'q_coefficients':[str(c) for c in q], 'mean_zero':mean_zero,
            'eigenvalue_index_in_constrained_space':1, 'range_proof':proof,
            'laplacian_ground_lower':2 if mean_zero else 0,
            'trial_coordinate_convention':'t=z; x=sqrt(1-t^2)*cos(phi)',
            'azimuthal_normalization':'cancels separately within each Rayleigh quotient',
            'trials':trials, 'status':'certified_exact' if lower==upper else 'certified_bound_open_gap',
            'analytic_dependencies':['unit sphere Laplacian ground 0 and mean-zero ground 2',
              'Rayleigh variational principle on the indicated full H1 space'],
            'formal_assistant_checked':False, **arithmetic.enclosure_fields(lower,upper)}


def cheap_ground_enclosure(q, mean_zero=True, range_proof=None):
    """Full-sphere bound from qlo and real trials 1 (if allowed), z and x.

    Uses only exact polynomial moments and a whole-domain range certificate.
    There is no finite spectral matrix, eigenvalue search or omitted-tail claim.
    """
    q=potential(q)
    if type(mean_zero) is not bool:
        raise ValueError('mean_zero must be boolean')
    proof=bernstein_range(q) if range_proof is None else range_proof
    return _cheap_ground_certificate(q,mean_zero,proof)


def verify_cheap_ground(cert, expected_q=None, expected_mean_zero=None):
    """Replay range, original moments, admissible real trials and both bounds."""
    try:
        if not isinstance(cert,dict):
            return False
        q=potential(cert['q_coefficients'])
        if expected_q is not None and q!=potential(expected_q):
            return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool
                                               or cert['mean_zero'] is not expected_mean_zero):
            return False
        expected=_cheap_ground_certificate(q,cert['mean_zero'],cert['range_proof'])
        return same(cert,expected)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


class Kernel:
    def __init__(self, q, m=0, mean_zero=True, modes=12, range_proof=None):
        parameters(m, mean_zero, 1, modes)
        self.q = potential(q)
        self.m, self.mean_zero, self.n = m, mean_zero, modes
        self.start = 1 if m == 0 and mean_zero else m
        self.tail_start = self.start+modes
        self.range_proof = potential_range(self.q) if range_proof is None else range_proof
        self.qlo, self.qhi = range_bounds(self.q, self.range_proof)
        self.beta = self.tail_start*(self.tail_start+1)+self.qlo
        data = assemble(self.q, m, self.start, modes)
        self.a, self.mass = data['A'], data['M']
        self.couplings = []
        # A degree-d multiplier couples at most d omitted radial degrees.
        # Projection removes l=0 only AFTER multiplication by the whole q.
        for l in range(self.tail_start, self.tail_start+len(self.q)-1):
            column = [(i-self.start, self.mass[i-self.start]*v)
                      for i, v in q_times_basis(self.q, m, l).items()
                      if self.start <= i < self.tail_start and v]
            self.couplings.append((l, basis_mass(m, l), column))

    def evidence(self):
        return {'A': [[str(x) for x in row] for row in self.a],
                'mass': [str(x) for x in self.mass],
                'tail_couplings': [{'degree': l, 'mass': str(mass),
                    'column': [[i, str(v)] for i, v in column]}
                    for l, mass, column in self.couplings]}

    def matrix(self, x, kind):
        x = F(x)
        if kind not in ('lower', 'upper_tail', 'upper_ritz'):
            raise ValueError('Unknown Schur comparison kind')
        if kind != 'upper_ritz' and x >= self.beta:
            raise ValueError('Infinite radial tail must be strictly positive')
        a = [row[:] for row in self.a]
        for i in range(self.n):
            a[i][i] -= x*self.mass[i]
        if kind == 'upper_ritz':
            return a
        for l, mass, column in self.couplings:
            diagonal = l*(l+1)+(self.qlo if kind == 'lower' else self.qhi)
            denominator = mass*(diagonal-x)
            for i, c in column:
                for j, d in column:
                    a[i][j] -= c*d/denominator
        return a

    def count(self, x, kind, independent=False):
        matrix = self.matrix(x, kind)
        return base.inertia(matrix) if independent else arithmetic.banded_inertia(matrix)

    def lower_holds(self, k, x):
        return x < self.beta and self.count(x, 'lower')[0] < k

    def upper_holds(self, k, x, kind):
        if kind == 'upper_tail' and x >= self.beta:
            return False
        negative, zero, _ = self.count(x, kind)
        return negative+zero >= k


def sector_certificate(kernel, k, lower, upper, upper_kind, independent=False):
    parameters(kernel.m, kernel.mean_zero, k, kernel.n)
    lo, hi = F(lower), F(upper)
    if lo > hi or upper_kind not in ('upper_tail', 'upper_ritz'):
        raise ValueError('Invalid spectral endpoints or upper comparison')
    low_count = kernel.count(lo, 'lower', independent)
    up_count = kernel.count(hi, upper_kind, independent)
    if low_count[0] >= k or up_count[0]+up_count[1] < k:
        raise ValueError('Exact spectral counts do not establish the endpoints')
    return {'format': SECTOR_FORMAT, 'geometry': 'unit_S2',
            'q_coefficients': [str(x) for x in kernel.q], 'azimuth_m': kernel.m,
            'mean_zero': kernel.mean_zero,
            'projection': 'remove_l0' if kernel.m == 0 and kernel.mean_zero else 'none',
            'scope': 'single_azimuth_sector_only',
            'degree_start': kernel.start, 'modes': kernel.n, 'eigenvalue_index': k,
            'range_proof': kernel.range_proof, 'kernel_evidence': kernel.evidence(),
            'radial_tail_start': kernel.tail_start, 'radial_tail_lower': str(kernel.beta),
            'upper_kind': upper_kind, 'lower_inertia': low_count, 'upper_inertia': up_count,
            'formal_assistant_checked': False, **arithmetic.enclosure_fields(lo, hi)}


def verify_sector(cert, expected_q=None, expected_m=None, expected_mean_zero=None,
                  expected_k=None, independent=True):
    try:
        if not isinstance(cert, dict):
            return False
        if expected_q is not None and potential(expected_q) != potential(cert['q_coefficients']):
            return False
        if expected_m is not None and (type(expected_m) is not int or cert['azimuth_m'] != expected_m):
            return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool
                                               or cert['mean_zero'] is not expected_mean_zero):
            return False
        if expected_k is not None and (type(expected_k) is not int or cert['eigenvalue_index'] != expected_k):
            return False
        kernel = Kernel(cert['q_coefficients'], cert['azimuth_m'], cert['mean_zero'],
                        cert['modes'], cert['range_proof'])
        expected = sector_certificate(kernel, cert['eigenvalue_index'], base.rational(cert['lower']),
                    base.rational(cert['upper']), cert['upper_kind'], independent=independent)
        return canonical_equal(cert, expected)
    except (ValueError, KeyError, TypeError, ArithmeticError, IndexError):
        return False


def certify_sector(q, m=0, mean_zero=True, k=1, modes=12, bits=44, range_proof=None):
    parameters(m, mean_zero, k, modes)
    if type(bits) is not int or not 8 <= bits <= 160:
        raise ValueError('Require integer bits in 8..160')
    kernel = Kernel(q, m, mean_zero, modes, range_proof)
    if len(kernel.q) == 1:
        l = kernel.start+k-1
        exact = l*(l+1)+kernel.q[0]
        cert = sector_certificate(kernel, k, exact, exact, 'upper_ritz')
    else:
        unperturbed_low = kernel.start*(kernel.start+1)
        l = kernel.start+k-1
        low, high = unperturbed_low+kernel.qlo-1, l*(l+1)+kernel.qhi+1
        kind = 'upper_tail' if high < kernel.beta else 'upper_ritz'
        if not kernel.upper_holds(k, high, kind):
            raise ArithmeticError('Analytic Ritz upper bracket failed')
        for _ in range(bits):
            mid = (low+high)/2
            if kernel.upper_holds(k, mid, kind):
                high = mid
            else:
                low = mid
        if kind == 'upper_ritz' and high < kernel.beta:
            kind = 'upper_tail'
            low = unperturbed_low+kernel.qlo-1
            for _ in range(bits):
                mid = (low+high)/2
                if kernel.upper_holds(k, mid, kind):
                    high = mid
                else:
                    low = mid
        upper = high
        distance = F(1)
        low = unperturbed_low+kernel.qlo-distance
        # Termination follows from the mass term dominating as low -> -infty.
        while not kernel.lower_holds(k, low):
            distance *= 2
            low = unperturbed_low+kernel.qlo-distance
        high = upper
        for _ in range(bits):
            mid = (low+high)/2
            if kernel.lower_holds(k, mid):
                low = mid
            else:
                high = mid
        cert = sector_certificate(kernel, k, low, upper, kind)
    if not verify_sector(cert, independent=True):
        raise ArithmeticError('Independent dense sector verifier rejected certificate')
    return cert


def full_certificate(q, mean_zero, sectors, range_proof, tolerance, independent=True):
    q = potential(q)
    if type(mean_zero) is not bool or not isinstance(sectors, list) or not 1 <= len(sectors) <= 65:
        raise ValueError('At least one consecutive azimuth sector is required')
    tolerance = exact(tolerance)
    if not F(1, 10**30) <= tolerance <= 1:
        raise ValueError('Require 1e-30 <= tolerance <= 1')
    qlo, _ = range_bounds(q, range_proof)
    for m, cert in enumerate(sectors):
        if not verify_sector(cert, expected_q=q, expected_m=m,
                             expected_mean_zero=mean_zero, expected_k=1, independent=independent):
            raise ValueError('Sector proof has wrong input, index, projection, or evidence')
    next_m = len(sectors)
    angular_tail = next_m*(next_m+1)+qlo
    finite_low = min(F(c['lower']) for c in sectors)
    upper = min(F(c['upper']) for c in sectors)
    lower = min(finite_low, angular_tail)
    covered = angular_tail >= upper
    return {'format': FULL_FORMAT, 'geometry': 'unit_S2', 'scope': scope(mean_zero),
            'q_coefficients': [str(x) for x in q], 'mean_zero': mean_zero,
            'potential_symmetry': 'axisymmetric_polynomial_degree_at_most_24',
            'eigenvalue_index_in_constrained_space': 1,
            'sectors': sectors, 'range_proof': range_proof,
            'omitted_azimuth_m_start': next_m, 'angular_tail_lower': str(angular_tail),
            'angular_tail_cannot_improve_best_upper': covered,
            'tolerance': str(tolerance),
            'status': 'certified_target_met' if upper-lower <= tolerance else 'certified_bound_open_gap',
            'formal_assistant_checked': False, **arithmetic.enclosure_fields(lower, upper)}


def verify_full(cert, expected_q=None, expected_mean_zero=None, independent=True):
    try:
        if not isinstance(cert, dict):
            return False
        if expected_q is not None and potential(expected_q) != potential(cert['q_coefficients']):
            return False
        if expected_mean_zero is not None and (type(expected_mean_zero) is not bool
                                               or cert['mean_zero'] is not expected_mean_zero):
            return False
        expected = full_certificate(cert['q_coefficients'], cert['mean_zero'], cert['sectors'],
                                     cert['range_proof'], base.rational(cert['tolerance']), independent)
        return canonical_equal(cert, expected)
    except (ValueError, KeyError, TypeError, ArithmeticError, IndexError):
        return False


def full_ground(q, mean_zero=True, modes=8, max_modes=32, bits=44, max_m=16,
                tolerance=F(1, 10**10), range_proof=None):
    parameters(0, mean_zero, 1, modes)
    if type(max_modes) is not int or not modes <= max_modes <= 64:
        raise ValueError('Require modes <= max_modes <= 64')
    if type(max_m) is not int or not 0 <= max_m <= 64:
        raise ValueError('Require 0 <= max_m <= 64')
    tolerance = exact(tolerance)
    if not F(1, 10**30) <= tolerance <= 1:
        raise ValueError('Require 1e-30 <= tolerance <= 1')
    q = potential(q)
    proof = potential_range(q) if range_proof is None else reuse_range(q, range_proof)
    qlo, _ = range_bounds(q, proof)
    sectors = []
    for m in range(max_m+1):
        size = modes
        while True:
            cert = certify_sector(q, m=m, mean_zero=mean_zero, modes=size, bits=bits, range_proof=proof)
            if F(cert['exact_width']) <= tolerance or size == max_modes:
                break
            size = min(max_modes, size*2)
        sectors.append(cert)
        next_m = m+1
        if next_m*(next_m+1)+qlo >= min(F(c['upper']) for c in sectors):
            break
    return full_certificate(q, mean_zero, sectors, proof, tolerance)


def kernel_certificate(q, m=0, mean_zero=True, modes=8, range_proof=None):
    kernel = Kernel(q, m, mean_zero, modes, range_proof)
    return {'format':KERNEL_FORMAT, 'q_coefficients':[str(x) for x in kernel.q],
            'azimuth_m':m, 'mean_zero':mean_zero, 'modes':modes,
            'degree_start':kernel.start, 'radial_tail_start':kernel.tail_start,
            'radial_tail_lower':str(kernel.beta), 'range_proof':kernel.range_proof,
            **kernel.evidence()}


def verify_kernel(cert, expected_q=None):
    try:
        if expected_q is not None and potential(expected_q) != potential(cert['q_coefficients']):
            return False
        expected = kernel_certificate(cert['q_coefficients'], cert['azimuth_m'],
                     cert['mean_zero'], cert['modes'], cert['range_proof'])
        return same(cert, expected)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


def exponential_approximation(a, degree=12):
    """Taylor polynomial with a rational uniform Lagrange remainder on [-1,1]."""
    a = exact(a)
    if abs(a)>3 or type(degree) is not int or not 0<=degree<=24:
        raise ValueError('Require |a|<=3 and integer 0<=degree<=24')
    polynomial = potential([a**k/F(factorial(k)) for k in range(degree+1)])
    exponent = (abs(a).numerator+abs(a).denominator-1)//abs(a).denominator
    envelope = F(3**exponent)
    delta = envelope*abs(a)**(degree+1)/factorial(degree+1)
    return {'format':TAYLOR_FORMAT, 'original_potential':{'kind':'exp_linear', 'a':str(a)},
            'interval':['-1','1'], 'degree':degree,
            'polynomial_coefficients':[str(x) for x in polynomial],
            'absolute_error_bound':str(delta), 'exponential_envelope':str(envelope),
            'envelope_exponent':exponent,
            'bound_statement':'sup_on_interval_abs(exp(a*t)-polynomial)<=absolute_error_bound',
            'analytic_dependencies':['Taylor theorem with Lagrange remainder for real exp',
              'exp monotonicity and exp(x+y)=exp(x)*exp(y)',
              'e<3: sum_(k>=2) 1/k! < sum_(k>=2) 1/2^(k-1)=1'],
            'formal_assistant_checked':False}


def verify_approximation(cert, expected_a=None):
    try:
        a = exact(cert['original_potential']['a'])
        if expected_a is not None and a != exact(expected_a):
            return False
        return same(cert, exponential_approximation(a, cert['degree']))
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


def perturb_certificate(approximation, spectral_certificate, tolerance=F(1,10**8)):
    if not verify_approximation(approximation):
        raise ValueError('Invalid uniform approximation certificate')
    q = approximation['polynomial_coefficients']
    if not verify_full(spectral_certificate, expected_q=q, expected_mean_zero=True):
        raise ValueError('Require a full mean-zero spectrum proof for the approximating polynomial')
    tolerance = exact(tolerance)
    if not F(1,10**30)<=tolerance<=1:
        raise ValueError('Require 1e-30<=tolerance<=1')
    delta = F(approximation['absolute_error_bound'])
    lo, hi = F(spectral_certificate['lower'])-delta, F(spectral_certificate['upper'])+delta
    return {'format':EXP_FORMAT, 'geometry':'unit_S2', 'scope':scope(True),
            'original_potential':deepcopy(approximation['original_potential']),
            'mean_zero':True, 'eigenvalue_index_in_constrained_space':1,
            'approximation':deepcopy(approximation), 'polynomial_spectrum':deepcopy(spectral_certificate),
            'perturbation_bound':str(delta), 'tolerance':str(tolerance),
            'analytic_dependencies':['minmax principle on the same mean-zero H1 space',
              'uniform potential perturbation changes each eigenvalue by at most delta'],
            'status':'certified_target_met' if hi-lo<=tolerance else 'certified_bound_open_gap',
            'formal_assistant_checked':False, **arithmetic.enclosure_fields(lo,hi)}


def verify_exponential(cert, expected_a=None):
    try:
        if not verify_approximation(cert['approximation'], expected_a):
            return False
        expected = perturb_certificate(cert['approximation'], cert['polynomial_spectrum'],
                                       exact(cert['tolerance']))
        return same(cert, expected)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


def exponential_ground(a, degree=12, modes=8, max_modes=12, bits=28, max_m=8,
                       tolerance=F(1,10**8), range_strategy='auto'):
    approximation = exponential_approximation(a, degree)
    q = approximation['polynomial_coefficients']
    proof = potential_range(q, range_strategy)
    spectrum = full_ground(q, True, modes, max_modes, bits, max_m, tolerance,
                           range_proof=proof)
    return perturb_certificate(approximation, spectrum, tolerance)
