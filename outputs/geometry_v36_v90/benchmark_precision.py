"""Small sequential mathematical benchmarks with saved, replayable certificates.

Default: three sequential repetitions; elapsed time includes generation and
independent verification, excludes imports and output-file I/O. These compare
algorithms/proof routes on stated mathematical problems, never model ability.
Run --verify to replay every saved certificate without repeating any search.
"""
import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from math import comb, isfinite
from pathlib import Path
import platform
from statistics import median
import sys
from time import perf_counter

import dual_v33 as dual
import exact_extrema as extrema
import family_v29 as family
import matrix_refine
import optimize_v34 as matrix_old
import positive_search as positive
import uniform_refine as uniform


ROOT = Path(__file__).resolve().parent
FORMAT = 'bounded_mathematical_benchmark_v1'
VERIFIERS = {
    'family_range': family.verify_range,
    'exact_extrema': extrema.verify,
    'matrix_optimization': dual.verify,
    'positive_lower': positive.verify_lower,
    'positive_optimization': positive.verify,
    'rayleigh': positive.verify_rayleigh,
    'poisson_family': family.verify,
    'uniform': uniform.verify,
}


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(k): encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    if isinstance(value, float) and not isfinite(value):
        raise ValueError('Nonfinite floating metadata is not a JSON benchmark value')
    return value


def json_bytes(value):
    return (json.dumps(encode(value), ensure_ascii=False, sort_keys=True, indent=2,
                       allow_nan=False)+'\n').encode('utf-8')


def save(path, value):
    data = json_bytes(value)
    path.write_bytes(data)
    return {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(data).hexdigest(),
            'bytes': len(data)}


def fraction_width(lower, upper):
    return str(F(upper)-F(lower))


@dataclass
class Case:
    identifier: str
    group: str
    method: str
    problem: str
    inputs: dict
    budget: dict
    scope: str
    produce: object
    certificates: object
    readout: object


def _range_readout(cert, exact):
    return {'maximum_lower': cert['maximum_lower'], 'maximum_upper': cert['maximum_upper'],
            'range_width': fraction_width(cert['maximum_lower'], cert['maximum_upper']),
            'status': cert['status'],
            'target_met': cert['status'] in ('exact', 'enclosed') if exact else cert['status']=='scale_enclosed',
            'refinement_operations': cert['refinement_nodes'] if exact else
                    (len(cert['intervals'])-1 if cert['intervals'] is not None else 0),
            'isolation_complete': cert.get('isolation_complete'),
            'claim_proved': True}


def range_cases():
    # Same four polynomials as test_exact_extrema.comparison_rows().
    roots = [F(-4, 5), F(-1, 4), F(2, 3)]
    derivative = [F(-1)]
    for r in roots:
        derivative = positive.e.mul(derivative, [-r, F(1)])
    polynomials = [
        ('cubic_irrational', [0, 1, 0, -1]),
        ('flat_peak_one_third', [-F(comb(8, k))*(-F(1, 3))**(8-k) for k in range(9)]),
        ('two_unequal_peaks', [F(0)]+[v/F(i+1) for i, v in enumerate(derivative)]),
        ('even_octic_flat_peak', [F(-1, 81), 0, F(4, 27), 0, F(-2, 3), 0, F(4, 3), 0, -1]),
    ]
    result = []
    for name, q in polynomials:
        for method, exact in [('legacy_bernstein', False), ('exact_sturm', True)]:
            result.append(Case('range_'+name+'_'+method, 'extrema', method, name,
                {'polynomial': q, 'interval': [-1, 1]},
                {'tolerance': F(1, 10**12), 'allowed_refinements': 32,
                 'max_leaves': None if exact else 33, 'max_nodes': 32 if exact else None,
                 'comparability': 'Both allow 32 refinement operations; work per operation differs.'},
                'Maximum of this fixed polynomial on the entire interval [-1,1]; no spectral claim.',
                (lambda q=q: extrema.maximize(q, tolerance=F(1, 10**12), max_nodes=32)) if exact else
                (lambda q=q: family.maximum(q, F(1, 10**12), max_leaves=33)),
                lambda cert, exact=exact: [('maximum', 'exact_extrema' if exact else 'family_range', cert)],
                lambda cert, exact=exact: _range_readout(cert, exact)))
    return result


def _matrix_readout(result):
    cert = result['certificate']
    return {'optimization_lower': cert['lower_bound'], 'optimization_upper': cert['upper_bound'],
            'optimization_gap': cert['gap'], 'status': cert['status'],
            'target_met': cert['status']=='global_gap_closed', 'claim_proved': True,
            'statistics': result['statistics'],
            'proposal_failures': [row for row in result['trace'] if row.get('outcome')=='proposal_failed'],
            'limitation': cert['limitation']}


def matrix_cases():
    problems = [
        ('matrix_standard', [[0, 1], [0, 0, 1]], [[1, 0], [0, 1]]),
        ('matrix_nondiagonal', [[0, 1], [0, 0, 1]], [[2, F(1, 2)], [F(1, 2), 1]]),
        ('matrix_scaled_same_objective', [[0, 1000], [0, 0, F(1, 1000)]],
         [[F(1, 10**6), 0], [0, 10**6]]),
    ]
    result = []
    for name, directions, cost in problems:
        for method in ['v34', 'matrix_refine']:
            budget = {'tolerance': F(1, 10**5), 'maximum_exchange_rounds': 5,
                      'range_leaves': 128,
                      'comparability': 'Same original directions, cost, tolerance and round cap; each round may do different work. No equal-Newton-step or equal-time claim.'}
            if method == 'matrix_refine':
                budget['newton_steps_per_round_cap'] = 280
            else:
                budget['newton_steps_per_round_cap'] = 'legacy solver internal schedule'
            result.append(Case(name+'_'+method, 'matrix', method, name,
                {'directions': directions, 'cost': cost}, budget,
                'Minimum trace(C G) for the fixed Poisson matrix envelope on all t in [-1,1]; not the sharp spectral inequality.',
                (lambda directions=directions, cost=cost: matrix_old.solve(directions, cost,
                    F(1, 10**5), max_rounds=5, max_range_leaves=128)) if method=='v34' else
                (lambda directions=directions, cost=cost: matrix_refine.refine(directions, cost,
                    F(1, 10**5), budget={'rounds': 5, 'range_leaves': 128, 'newton_steps': 280})),
                lambda result: [('matrix_optimization', 'matrix_optimization', result['certificate'])],
                _matrix_readout))
    return result


def _fixed_poisson(q):
    basis = positive.monomial_basis(q, degree=6)
    theta = positive.poisson_initialization(q, basis)
    return {'lower': positive.lower_certificate(q, basis, theta),
            'rayleigh': positive.propose_rayleigh(q)}


def _positive_readout(result, fixed=False):
    if fixed:
        lo, hi = result['lower']['spectral_lower'], result['rayleigh']['spectral_upper']
        return {'spectral_lower': lo, 'spectral_upper': hi, 'spectral_gap': fraction_width(lo, hi),
                'status': 'fixed_candidate_verified', 'target_met': None, 'claim_proved': True,
                'ansatz_gap': None, 'range_gap': result['lower']['range_gap'],
                's_coefficients': result['lower']['s_coefficients']}
    cert = result['certificate']
    keys = ['spectral_lower', 'spectral_upper', 'spectral_gap', 'ansatz_lower',
            'ansatz_upper', 'ansatz_gap', 'status', 'quadratic_stationarity_gap',
            'average_contact_and_range_gap']
    return {**{k: cert[k] for k in keys}, 'target_met': cert['status']=='ansatz_gap_closed',
            'claim_proved': True, 'range_gap': cert['lower_certificate']['range_gap'],
            's_coefficients': cert['lower_certificate']['s_coefficients'],
            'proposal_failures': result['failed_proposals'],
            'limitation': cert['limitation'], 'accepted_candidate_count': len(result['attempts'])}


def positive_cases():
    result = []
    for name, q in [('quadratic_6', [0, 0, 6]), ('quadratic_20', [0, 0, 20])]:
        for method, degree in [('fixed_poisson', None), ('exponential_degree2', 2), ('exponential_degree6', 6)]:
            fixed = degree is None
            result.append(Case('positive_'+name+'_'+method, 'positive', method, name,
                {'q_coefficients': q},
                {'log_degree': 2 if fixed else degree, 'rayleigh_degree': 6,
                 'range_tolerance': F(1, 10**6), 'range_leaves': 128,
                 'optimization_tolerance': None if fixed else F(1, 10**4),
                 'exchange_rounds': 0 if fixed else 4,
                 'newton_steps_per_temperature': 0 if fixed else 30,
                 'comparability': 'Same fixed q and independently computed degree-6 Rayleigh trial; log ansatz dimension differs intentionally.'},
                'Full unit two-sphere ground energy for this fixed q. Ansatz upper concerns only the best Barta lower in the specified log basis; spectral upper comes solely from Rayleigh.',
                (lambda q=q: _fixed_poisson(q)) if fixed else
                (lambda q=q, degree=degree: positive.search(q, degree=degree)),
                (lambda result: [('positive_lower', 'positive_lower', result['lower']),
                                 ('rayleigh_upper', 'rayleigh', result['rayleigh'])]) if fixed else
                (lambda result: [('positive_optimization', 'positive_optimization', result['certificate']),
                                 ('rayleigh_upper', 'rayleigh', result['certificate']['rayleigh_certificate'])]),
                lambda result, fixed=fixed: _positive_readout(result, fixed)))
    # Honest under-budget/near-singular probability proposal retained as a
    # separate diagnostic, not substituted into the main fixed-budget pair.
    result.append(Case('positive_underbudget_degree6', 'positive', 'underbudget_degree6',
        'quadratic_6_underbudget', {'q_coefficients': [0, 0, 6]},
        {'log_degree': 6, 'optimization_tolerance': F(1, 10**10), 'exchange_rounds': 1,
         'newton_steps_per_temperature': 1},
        'Same fixed q=6t^2; deliberately inadequate optimization budget, requiring an honest open gap.',
        lambda: positive.search([0, 0, 6], degree=6, max_newton=1, max_exchanges=1,
                                 tolerance=F(1, 10**10)),
        lambda result: [('positive_optimization', 'positive_optimization', result['certificate'])],
        _positive_readout))
    return result


def _uniform_counterexample():
    cert = uniform.synthesize_log_sobolev(coefficient=F(1, 7))
    witness = cert['counterexample']
    amplitude = F(witness['amplitude'])
    rayleigh = positive.rayleigh_certificate([0, amplitude], witness['trial_power_coefficients'])
    objective_upper = F(rayleigh['spectral_upper'])+F(cert['coefficient'])*amplitude**2
    if objective_upper != F(witness['spectral_objective_upper']) or objective_upper >= 0:
        raise ArithmeticError('Independent Rayleigh counterexample check failed')
    return {'certificate': cert, 'independent_rayleigh': rayleigh,
            'independent_spectral_objective_upper': str(objective_upper)}


def _uniform_readout(result, legacy=False, counterexample=False):
    cert = result['certificate'] if counterexample else result
    if legacy:
        return {'coefficient': cert['scale'], 'status': 'proved', 'target_met': True,
                'claim_proved': True, 'all_real': True,
                'proof_route': 'exact fixed Poisson all-amplitude envelope'}
    record = {'coefficient': cert['coefficient'], 'status': cert['status'],
              'target_met': cert['status']=='proved', 'claim_proved': cert['status']=='proved',
              'all_real': cert.get('all_real', True), 'scope': cert['scope'],
              'proof_route': cert['method'], 'leaf_count': cert.get('leaf_count'),
              'failure_witness': cert.get('failure_witness'),
              'trusted_analytic_theorem': cert.get('trusted_analytic_theorem')}
    if counterexample:
        record.update(counterexample_valid=True, counterexample=cert['counterexample'],
                      independent_rayleigh_upper=result['independent_rayleigh']['spectral_upper'],
                      spectral_objective_upper=result['independent_spectral_objective_upper'])
    return record


def uniform_cases():
    scope = 'lambda_1(-Delta+a*t) >= -c*a^2 for every real a on the full unit two-sphere.'
    definitions = [
        ('uniform_poisson_1_4', 'fixed_poisson_1_4', F(1, 4), {'order': 1},
         lambda: family.synthesize([0, 1]), True, False),
        ('uniform_barta_1_5_all_real', 'barta_cubic_1_5', F(1, 5), {'order': 3, 'max_leaves': 256},
         lambda: uniform.synthesize(coefficient=F(1, 5), order=3, all_real=True), False, False),
        ('uniform_sharp_1_6', 'log_sobolev_sharp_1_6', F(1, 6), {'analytic_rule': 'sphere log-Sobolev + exact moment-series rule'},
         lambda: uniform.synthesize_log_sobolev(), False, False),
        ('uniform_false_1_7', 'rayleigh_counterexample_1_7', F(1, 7), {'rayleigh_degree': 1},
         _uniform_counterexample, False, True),
        ('uniform_linear_1_5_obstruction', 'barta_linear_1_5_obstruction', F(1, 5), {'order': 1, 'max_leaves': 256},
         lambda: uniform.synthesize(coefficient=F(1, 5), order=1, all_real=True), False, False),
        ('uniform_quartic_1_6_obstruction', 'barta_quartic_1_6_obstruction', F(1, 6), {'order': 4, 'max_leaves': 256},
         lambda: uniform.synthesize(coefficient=F(1, 6), order=4, all_real=True), False, False),
        ('uniform_cubic_1_5_unresolved', 'barta_cubic_1_5_underbudget', F(1, 5), {'order': 3, 'max_leaves': 1},
         lambda: uniform.synthesize(coefficient=F(1, 5), order=3, all_real=True, max_leaves=1), False, False),
    ]
    result = []
    for identifier, method, coefficient, budget, producer, legacy, counterexample in definitions:
        result.append(Case(identifier, 'uniform', method, 'all_real_linear_spectral_bound',
            {'direction': [0, 1], 'requested_coefficient': coefficient, 'domain': 'all real amplitudes'},
            {**budget, 'comparability': 'Proof routes and coefficient strengths differ; this is not a matched-work timing race.'},
            scope, producer,
            (lambda raw: [('uniform_disproof', 'uniform', raw['certificate']),
                          ('independent_rayleigh', 'rayleigh', raw['independent_rayleigh'])]) if counterexample else
            (lambda raw, legacy=legacy: [('uniform_claim', 'poisson_family' if legacy else 'uniform', raw)]),
            lambda raw, legacy=legacy, counterexample=counterexample:
                _uniform_readout(raw, legacy, counterexample)))
    return result


def all_cases():
    return range_cases()+matrix_cases()+positive_cases()+uniform_cases()


def source_hashes():
    names = ['benchmark_precision.py', 'exact_extrema.py', 'matrix_refine.py',
             'optimize_v34.py', 'positive_search.py', 'uniform_refine.py', 'family_v29.py',
             'dual_v33.py', 'error_bounds.py', 'algebra.py', 'float_sdp.py']
    return {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names}


def bound_to_input(group, inputs, budget, method, verifier, cert):
    """A valid certificate for a different problem must not pass a benchmark."""
    def vector(values):
        return positive.e.trim(list(map(F, values)))
    try:
        if group == 'extrema':
            return (vector(cert['polynomial']) == vector(inputs['polynomial']) and
                    (list(map(F, cert['interval'])) == list(map(F, inputs['interval']))
                     if verifier=='exact_extrema' else list(map(F, inputs['interval']))==[F(-1), F(1)]))
        if group == 'matrix':
            return ([vector(p) for p in cert['dual']['directions']] == [vector(p) for p in inputs['directions']]
                    and [list(map(F, row)) for row in cert['dual']['cost']] ==
                        [list(map(F, row)) for row in inputs['cost']])
        if group == 'positive':
            lower = cert['lower_certificate'] if verifier=='positive_optimization' else cert
            if vector(lower['q_coefficients']) != vector(inputs['q_coefficients']):
                return False
            if verifier != 'rayleigh':
                if len(lower['s_coefficients'])-1 > budget['log_degree']:
                    return False
                if method=='fixed_poisson' and vector(lower['s_coefficients']) != family.poisson(inputs['q_coefficients'])[0]:
                    return False
            return True
        if group == 'uniform':
            if verifier=='rayleigh':
                return vector(cert['q_coefficients']) == [F(0), F(1)] and vector(cert['trial_power_coefficients']) == [F(1), F(-1, 2)]
            if verifier=='poisson_family':
                return ([vector(v) for v in cert['directions']] == [vector(inputs['direction'])]
                        and F(cert['scale']) == F(inputs['requested_coefficient']))
            return (vector(cert['direction']) == vector(inputs['direction'])
                    and F(cert['coefficient']) == F(inputs['requested_coefficient'])
                    and 'all real amplitudes' in cert['scope'])
        return False
    except (TypeError, KeyError, ValueError, IndexError, ZeroDivisionError):
        return False


def run_once(case, repetition, certificate_dir):
    record = {'case_id': case.identifier, 'repetition': repetition,
              'certificates': [], 'exception': None}
    raw, extracted = None, []
    start = perf_counter()
    try:
        raw = case.produce()
        generated = perf_counter()
        extracted = case.certificates(raw)
        checks = []
        for label, verifier, certificate in extracted:
            mathematical_valid = bool(VERIFIERS[verifier](certificate))
            input_valid = bound_to_input(case.group, case.inputs, case.budget, case.method, verifier, certificate)
            checks.append({'label': label, 'verifier': verifier,
                           'mathematically_valid': mathematical_valid,
                           'bound_to_benchmark_input': input_valid,
                           'valid': mathematical_valid and input_valid})
        readout = case.readout(raw)
        stopped = perf_counter()
        valid = all(row['valid'] for row in checks) and bool(checks)
        record.update(generation_seconds=generated-start,
                      external_verification_and_readout_seconds=stopped-generated,
                      seconds_including_verification=stopped-start,
                      certificate_valid=valid, mathematical=encode(readout), checks=checks)
    except Exception as exc:
        stopped = perf_counter()
        record.update(seconds_including_verification=stopped-start,
                      certificate_valid=False, mathematical=None,
                      exception={'type': type(exc).__name__, 'message': str(exc)})
    # Save the original certificates unchanged, one file per certificate per
    # repetition. Raw search output separately retains all candidate traces.
    stem = case.identifier+'_repeat'+str(repetition)
    for label, verifier, certificate in extracted:
        reference = save(certificate_dir/(stem+'_'+label+'.certificate.json'), certificate)
        record['certificates'].append({**reference, 'label': label, 'verifier': verifier})
    record['raw_result'] = save(certificate_dir/(stem+'.run.json'),
                               raw if raw is not None else {'exception': record['exception']})
    return record


def aggregate(case, runs):
    valid_runs = [r for r in runs if r['certificate_valid']]
    fingerprints = [tuple(item['sha256'] for item in r['certificates']) for r in valid_runs]
    return {'case_id': case.identifier, 'group': case.group, 'method': case.method,
            'problem': case.problem, 'inputs': encode(case.inputs), 'budget': encode(case.budget),
            'scope': case.scope, 'repetitions': len(runs), 'valid_runs': len(valid_runs),
            'median_seconds_including_verification': median(r['seconds_including_verification'] for r in runs),
            'minimum_seconds_including_verification': min(r['seconds_including_verification'] for r in runs),
            'maximum_seconds_including_verification': max(r['seconds_including_verification'] for r in runs),
            'exact_certificate_fingerprints_stable': len(set(fingerprints))==1 if fingerprints else False,
            'mathematical': valid_runs[0]['mathematical'] if valid_runs else None,
            'runs': runs}


def comparisons(rows):
    result = []
    baselines = {'extrema': ('legacy_bernstein', 'range_width'),
                 'matrix': ('v34', 'optimization_gap'),
                 'positive': ('fixed_poisson', 'spectral_gap'),
                 'uniform': ('fixed_poisson_1_4', 'coefficient')}
    for row in rows:
        if row['group'] not in baselines:
            continue
        baseline_method, metric = baselines[row['group']]
        if row['method']==baseline_method or row['mathematical'] is None:
            continue
        old = next((r for r in rows if r['group']==row['group'] and r['problem']==row['problem']
                    and r['method']==baseline_method and r['mathematical'] is not None), None)
        if old is None:
            continue
        before, after = old['mathematical'], row['mathematical']
        proved = before.get('claim_proved') and after.get('claim_proved')
        change = F(before[metric])-F(after[metric])
        oldtime, newtime = old['median_seconds_including_verification'], row['median_seconds_including_verification']
        result.append({'group': row['group'], 'problem': row['problem'],
            'baseline': old['case_id'], 'candidate': row['case_id'], 'metric': metric,
            'baseline_exact': before[metric], 'candidate_exact': after[metric],
            'requested_metric_reduction': str(change),
            'verified_metric_improvement': str(change) if proved else None,
            'candidate_claim_proved': bool(after.get('claim_proved')),
            'quality_negative_gain': change < 0 if proved else None,
            'baseline_median_seconds': oldtime, 'candidate_median_seconds': newtime,
            'median_seconds_saved': oldtime-newtime,
            'runtime_negative_gain': newtime > oldtime,
            'observed_time_ratio_baseline_over_candidate': oldtime/newtime if newtime else None,
            'timing_scope': 'Three sequential local wall-clock repetitions including verification; work differs, no universal speed guarantee or model comparison.'})
    return result


def run_benchmark(output=None, groups=None, repetitions=3, progress=True):
    if repetitions != 3:
        raise ValueError('This benchmark protocol requires exactly three sequential repetitions')
    output = ROOT/'results'/'math_benchmark.json' if output is None else Path(output).resolve()
    output.relative_to(ROOT)
    certificate_dir = output.parent/'math_certificates'
    certificate_dir.mkdir(parents=True, exist_ok=True)
    cases = [c for c in all_cases() if groups is None or c.group in groups]
    if not cases:
        raise ValueError('No selected benchmark cases')
    case_runs = {case.identifier: [] for case in cases}
    benchmark_start = perf_counter()
    for repetition in range(1, repetitions+1):
        previous_group = None
        for case in cases:
            if progress and case.group != previous_group:
                print('Repetition '+str(repetition)+'/3: '+case.group, flush=True)
                previous_group = case.group
            record = run_once(case, repetition, certificate_dir)
            case_runs[case.identifier].append(record)
    rows = [aggregate(case, case_runs[case.identifier]) for case in cases]
    result = {'format': FORMAT, 'description': 'Bounded mathematical algorithm/proof-route comparison; not a model capability benchmark.',
        'protocol': {'repetitions': repetitions, 'execution': 'strictly sequential; complete case order repeated three times',
                     'timed_work': 'generation including any internal checks, plus an external exact certificate verification and result readout',
                     'excluded_from_case_timings': ['module imports', 'certificate serialization and file I/O'],
                     'correctness': 'Every original certificate is saved and independently replayed; verified unresolved, ansatz obstruction and disproved outcomes remain distinct.',
                     'limitations': ['Round/node caps do not imply equal work per round/node.',
                                     'Local wall-clock observations include uncontrolled system load and interpreter effects.',
                                     'Negative precision or runtime gains and failed claims are retained.',
                                     'A fixed-potential ansatz upper is not a spectral upper.',
                                     'The sharp log-Sobolev route depends on an explicitly named analytic theorem.']},
        'environment': {'python': sys.version, 'platform': platform.platform(), 'machine': platform.machine()},
        'source_sha256': source_hashes(), 'total_seconds_including_file_io': perf_counter()-benchmark_start,
        'case_count': len(rows), 'run_count': sum(len(r['runs']) for r in rows),
        'all_certificates_valid': all(r['valid_runs']==repetitions for r in rows),
        'cases': rows, 'comparisons': comparisons(rows)}
    save(output, result)
    return result


def replay(path):
    """Integrity + mathematical replay, never search; report valid failed claims."""
    path = Path(path).resolve()
    result = json.loads(path.read_text())
    if result.get('format') != FORMAT:
        raise ValueError('Unknown benchmark format')
    checks, failures = [], []
    for case in result['cases']:
        for run in case['runs']:
            for item in run['certificates']:
                target = (ROOT/item['path']).resolve()
                target.relative_to(ROOT)
                try:
                    data = target.read_bytes()
                    integrity = hashlib.sha256(data).hexdigest()==item['sha256']
                    certificate = json.loads(data)
                    valid = integrity and bool(VERIFIERS[item['verifier']](certificate))
                    valid = valid and bound_to_input(case['group'], case['inputs'], case['budget'],
                                                    case['method'], item['verifier'], certificate)
                except (ValueError, KeyError, TypeError, OSError) as exc:
                    valid = False
                record = {'case_id': case['case_id'], 'repetition': run['repetition'],
                          'label': item['label'], 'valid': valid}
                checks.append(record)
                if not valid:
                    failures.append(record)
    current = source_hashes()
    return {'certificate_count': len(checks), 'all_saved_certificates_valid': bool(checks) and not failures,
            'failed_certificate_checks': failures,
            'generation_failures_recorded': sum(run['exception'] is not None for case in result['cases'] for run in case['runs']),
            'source_files_changed_since_run': [name for name, digest in result['source_sha256'].items()
                                               if current.get(name)!=digest],
            'note': 'Replayed certificate validity includes honest open, obstruction, and disproved records; it does not turn those into proved inequalities.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', help='Output JSON within this module directory; default results/math_benchmark.json')
    parser.add_argument('--groups', nargs='+', choices=['extrema', 'matrix', 'positive', 'uniform'])
    parser.add_argument('--verify', metavar='JSON', help='Replay every saved certificate from a benchmark JSON')
    args = parser.parse_args()
    if args.verify:
        result = replay(args.verify)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0 if result['all_saved_certificates_valid'] else 1
    result = run_benchmark(args.output, args.groups)
    print(json.dumps({'output': str(Path(args.output).resolve() if args.output else ROOT/'results'/'math_benchmark.json'),
                      'cases': result['case_count'], 'runs': result['run_count'],
                      'all_certificates_valid': result['all_certificates_valid'],
                      'total_seconds': result['total_seconds_including_file_io']}, indent=2))
    return 0 if result['all_certificates_valid'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
