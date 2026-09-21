"""External q=4 sphere GN trials using the frozen V90 polynomial verifier.

Only a two-dimensional function space is optimized, globally over real
coefficients. These certificates do NOT bound the ratio above on all H^1(S^2).
Run with python3 -B to keep the frozen release free of new bytecode files.
"""

from fractions import Fraction as F
from math import comb, pi, sqrt
from pathlib import Path
import json
import sys
from time import perf_counter

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'geometry_v36_v90'))
import exact_extrema as extrema
from error_bounds import add, mul, scale, derivative, evaluate, trim

FORMAT = 'sphere_gn_q4_two_plane_v1'
SCOPE = 'all_nonzero_real_combinations_of_the_two_recorded_basis_functions'
RATIO = 'pi_times_K4_for_standard_area_measure_on_unit_S2'
INTERVAL = ['-1', '1']


def mean(p):
    """Integral against dt/2, the axisymmetric unit sphere probability law."""
    return sum((F(c, i+1) for i, c in enumerate(p) if i % 2 == 0), F(0))


def power(p, n):
    result = [F(1)]
    for _ in range(n):
        result = mul(result, p)
    return result


def strings(p):
    return [str(v) for v in p]


def basis_input(raw):
    if not isinstance(raw, list) or len(raw) != 2:
        raise ValueError('Exactly two basis functions required')
    out = []
    for p in raw:
        if not isinstance(p, list) or not 1 <= len(p) <= 7:
            raise ValueError('Trial basis degree must be at most six')
        if any(isinstance(v, bool) or not isinstance(v, (int, str, F)) for v in p):
            raise ValueError('Exact rational basis coefficients required')
        out.append(trim([F(v) for v in p]))
    f, g = out
    if mean(f) != 0 or mean(g) != 0:
        raise ValueError('Both basis functions must have zero sphere mean')
    if mean(mul(f, f))*mean(mul(g, g)) - mean(mul(f, g))**2 <= 0:
        raise ValueError('Basis must be linearly independent')
    return out


def chart_moments(f, g):
    """Polynomials in s for u(t)=f(t)+s*g(t), with exact integration in t."""
    mass = [mean(mul(f, f)), 2*mean(mul(f, g)), mean(mul(g, g))]
    df, dg = derivative(f), derivative(g)
    w = [F(1), F(0), F(-1)]
    energy = [mean(mul(w, mul(df, df))),
              2*mean(mul(w, mul(df, dg))), mean(mul(w, mul(dg, dg)))]
    numerator = [comb(4, j)*mean(mul(power(f, 4-j), power(g, j))) for j in range(5)]
    # Standard sphere area = 4*pi. Thus pi*K4 = N/(4*M*E).
    denominator = scale(mul(mass, energy), F(4))
    return {'mass': trim(mass), 'energy': trim(energy),
            'numerator': trim(numerator), 'denominator': denominator}


def propose(moments):
    """Untrusted floating search; neither upper bound nor coverage uses it."""
    n, d = ([float(x) for x in moments[k]] for k in ('numerator', 'denominator'))
    def score(x):
        return evaluate(n, x)/evaluate(d, x)
    grid = [-1 + j/512 for j in range(1025)]
    values = [score(x) for x in grid]
    proposals = [F(-1), F(1)]
    r = (sqrt(5)-1)/2
    for j in range(1, len(grid)-1):
        if values[j] >= values[j-1] and values[j] >= values[j+1]:
            left, right = grid[j-1], grid[j+1]
            for _ in range(75):
                a, b = right-r*(right-left), left+r*(right-left)
                if score(a) < score(b):
                    left = a
                else:
                    right = b
            proposals.append(F((left+right)/2).limit_denominator(10**10))
    def exact_score(s):
        return evaluate(moments['numerator'], s)/evaluate(moments['denominator'], s)
    s = max(proposals, key=exact_score)
    return s, exact_score(s)


def sign_proof(p, strict=False):
    cert = extrema.maximize(p, tolerance=F(1, 10**22), max_nodes=256)
    if not extrema.verify(cert):
        raise ArithmeticError('Frozen polynomial verifier rejected its output')
    bound = F(cert['maximum_upper'])
    accepted = bound < 0 if strict else bound <= 0
    return cert, accepted


def certify(raw_basis, gap_target=F(1, 10**10)):
    basis = basis_input(raw_basis)
    chart_data = [chart_moments(*basis), chart_moments(*reversed(basis))]
    proposals = [propose(m) for m in chart_data]
    witness_chart = max(range(2), key=lambda j: proposals[j][1])
    s, lower = proposals[witness_chart]
    # A proposed ceiling is accepted only by polynomial sign certificates.
    # Failure enlarges the ceiling, never silently accepts floating evidence.
    upper = lower + gap_target
    proof_attempts = 0
    for _ in range(8):
        proofs = []
        for m in chart_data:
            proof_attempts += 1
            proof, accepted = sign_proof(add(m['numerator'], scale(m['denominator'], -upper)))
            proofs.append(proof)
            if not accepted:
                break
        if len(proofs) == 2 and accepted:
            break
        upper = lower + 10*(upper-lower)
    else:
        raise ArithmeticError('Unable to certify the proposed global ceiling')
    charts = []
    for j, m in enumerate(chart_data):
        mass_proof, mass_ok = sign_proof(scale(m['mass'], -1), strict=True)
        energy_proof, energy_ok = sign_proof(scale(m['energy'], -1), strict=True)
        if not mass_ok or not energy_ok:
            raise ArithmeticError('Denominator positivity not certified')
        charts.append({'order': [j, 1-j],
                       'moments': {key: strings(p) for key, p in m.items()},
                       'mass_positive': mass_proof, 'energy_positive': energy_proof,
                       'upper_residual_nonpositive': proofs[j]})
    cert = {'format': FORMAT, 'scope': SCOPE, 'ratio': RATIO,
            'basis': [strings(p) for p in basis],
            'lower': str(lower), 'upper': str(upper), 'gap': str(upper-lower),
            'witness': {'chart': witness_chart, 'parameter': str(s)},
            'charts': charts}
    if not verify(cert, expected_basis=raw_basis):
        raise ArithmeticError('New reduction verifier rejected the certificate')
    return cert, {'upper_sign_proof_attempts': proof_attempts,
                  'requested_gap': str(gap_target), 'target_met': upper-lower <= gap_target}


def verify(cert, expected_basis=None):
    """Replay exact evidence and the two-chart reduction without search."""
    try:
        if cert['format'] != FORMAT or cert['scope'] != SCOPE or cert['ratio'] != RATIO:
            return False
        basis = basis_input(cert['basis'])
        if expected_basis is not None and basis != basis_input(expected_basis):
            return False
        lower, upper = F(cert['lower']), F(cert['upper'])
        if not 0 < lower <= upper or F(cert['gap']) != upper-lower:
            return False
        if len(cert['charts']) != 2:
            return False
        all_moments = []
        for j, chart in enumerate(cert['charts']):
            if chart['order'] != [j, 1-j]:
                return False
            m = chart_moments(basis[j], basis[1-j])
            all_moments.append(m)
            if chart['moments'] != {key: strings(p) for key, p in m.items()}:
                return False
            claims = [('mass_positive', scale(m['mass'], -1), True),
                      ('energy_positive', scale(m['energy'], -1), True),
                      ('upper_residual_nonpositive',
                       add(m['numerator'], scale(m['denominator'], -upper)), False)]
            for field, polynomial, strict in claims:
                proof = chart[field]
                if (proof['polynomial'] != strings(polynomial)
                        or proof['interval'] != INTERVAL or not extrema.verify(proof)):
                    return False
                bound = F(proof['maximum_upper'])
                if bound >= 0 if strict else bound > 0:
                    return False
        w = cert['witness']
        if type(w['chart']) is not int or w['chart'] not in (0, 1):
            return False
        s = F(w['parameter'])
        if not -1 <= s <= 1:
            return False
        m = all_moments[w['chart']]
        den = evaluate(m['denominator'], s)
        return den > 0 and lower == evaluate(m['numerator'], s)/den
    except (ArithmeticError, ValueError, KeyError, TypeError, IndexError):
        return False


def centered_cap(n):
    # ((1+t)/2)^n - 1/(n+1), a smooth mean-zero polynomial on S^2.
    p = [F(comb(n, j), 2**n) for j in range(n+1)]
    p[0] -= F(1, n+1)
    return p


def run():
    cases = [
        ('legendre_P1_P3', [[F(0), F(1)], [F(0), F(-3, 2), F(0), F(5, 2)]]),
        ('centered_caps_5_6', [centered_cap(5), centered_cap(6)]),
    ]
    result = {'problem': 'zero-mean unit S2 q=4 Gagliardo-Nirenberg ratio',
              'source': 'https://arxiv.org/abs/2204.12414',
              'analytic_dependencies': ['area measure d_sigma = 4*pi*d_mu',
                                        'axisymmetric energy density (1-t^2)*(u_prime)^2',
                                        'homogeneous ratio: two charts cover every real projective direction'],
              'global_scope': SCOPE,
              'no_claim_of_universal_upper_bound_or_open_problem_solution': True,
              'no_change_to_frozen_v90': True, 'cases': []}
    target = HERE / 'certificates' / 'interpolation'
    target.mkdir(parents=True, exist_ok=True)
    for name, basis in cases:
        started = perf_counter()
        cert, details = certify(basis)
        cert = json.loads(json.dumps(cert))
        accepted = verify(cert, expected_basis=basis)
        elapsed = perf_counter()-started
        path = target / (name+'.json')
        path.write_text(json.dumps(cert, indent=2, sort_keys=True)+'\n')
        lower, upper = F(cert['lower']), F(cert['upper'])
        result['cases'].append({'case': name, 'certificate': str(path.relative_to(HERE)),
            'basis': cert['basis'], 'verified': accepted,
            'pi_K4_lower': str(lower), 'pi_K4_upper': str(upper), 'pi_K4_gap': str(upper-lower),
            'K4_area_approximation_only': float(lower)/pi,
            'decimal_conversion_is_not_an_outward_rounded_certificate': True,
            'seconds_including_json_replay': elapsed, 'search': details,
            'witness': cert['witness']})
    path = HERE / 'interpolation_results.json'
    path.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


if __name__ == '__main__':
    run()
