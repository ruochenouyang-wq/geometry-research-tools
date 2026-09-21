"""Run the frozen scope/reference set without deleting failures or editing code.

Standard library only. Published decimals are compatibility controls, never
mathematical certificates. Each run retains its full certificates and JSON.
"""
import sys
sys.dont_write_bytecode = True
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime, timezone
from time import perf_counter
import copy
import hashlib
import json

import projected_spectrum as projected

HERE = Path(__file__).resolve().parent
LEGACY = HERE.parent / 'geometry_v36_v90'
base = projected.base


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def overlaps(first, second):
    return max(F(first['lower']), F(second['lower'])) <= min(F(first['upper']), F(second['upper']))


def reference_check(case, cert):
    ref = case['reference']
    low, high = F(cert['lower']), F(cert['upper'])
    kind = ref['kind']
    if kind == 'analytic_exact':
        value = F(ref['value'])
        return {'kind': kind, 'compatible': low <= value <= high, 'exact_reference': str(value)}
    if 'value' in ref:
        value, radius = F(ref['value']), F(ref['absolute_compatibility_radius'])
        return {'kind': kind, 'compatible': low <= value + radius and high >= value - radius,
                'reference_band': [str(value-radius), str(value+radius)],
                'comparison_is_not_proof': True}
    if kind == 'analytic_coarse_enclosure':
        return {'kind': kind, 'compatible': low <= F(ref['upper']) and high >= F(ref['lower']),
                'analytic_band': [ref['lower'], ref['upper']]}
    return {'kind': kind, 'compatible': None, 'requires_cross_case_relation': 'R5'}


def invoke(case):
    q = case['potential_coefficients']
    space = case['requested_space']
    if space == 'full_sphere_mean_zero':
        options = dict(mean_zero=True, modes=8, max_modes=32, bits=44, max_m=16,
                       tolerance=F(1, 10**10))
        cert = projected.full_ground(q, **options)
        verify = lambda c: projected.verify_full(c, expected_q=q, expected_mean_zero=True)
    else:
        options = dict(m=case['m'], mean_zero=space == 'axisymmetric_mean_zero',
                       k=case['eigenvalue_index'], modes=12, bits=44)
        cert = projected.certify_sector(q, **options)
        verify = lambda c: projected.verify_sector(c, expected_q=q, expected_m=options['m'],
                        expected_mean_zero=options['mean_zero'], expected_k=options['k'])
    return cert, verify, {key: str(value) if isinstance(value, F) else value for key, value in options.items()}


def relations(records):
    by_id = {row['id']: row for row in records}
    outcomes = []
    def evaluate(label, ids, predicate, detail):
        rows = [by_id[x] for x in ids]
        if any('lower' not in row for row in rows):
            result = None
        else:
            result = bool(predicate(*rows))
        outcomes.append({'id': label, 'cases': ids, 'passed': result, 'criterion': detail})
    for label, ids in [('R1', ['E07', 'E08', 'E09']), ('R2', ['E10', 'E11', 'E12'])]:
        evaluate(label, ids, lambda ordinary, compressed, full:
                 F(compressed['upper']) < F(ordinary['lower']) and overlaps(compressed, full),
                 'Projected and ordinary intervals strictly separate; two distinct mean-zero scopes agree.')
    evaluate('R3', ['E13', 'E20'], lambda sector, full:
             overlaps(sector, full) and full.get('winning_absolute_m') == 1 and full.get('angular_tail_covered') is True,
             'Full mean-zero ground includes all sectors and agrees with the m=1 branch.')
    evaluate('R4', ['E01', 'E02'], lambda zero, shifted:
             F(shifted['lower'])-F(zero['lower']) == F(3, 2)
             and F(shifted['upper'])-F(zero['upper']) == F(3, 2),
             'Exact constant potential shift is 3/2.')
    evaluate('R5', ['E21', 'E22', 'E23'], lambda ordinary, compressed, full:
             F(full['upper']) < F(compressed['lower']) and F(compressed['upper']) < F(ordinary['lower']),
             'For q=t+t^2, full mean-zero < m0 compression < ordinary m0 second eigenvalue.')
    return outcomes


def behavioral_controls(destination):
    controls = []
    for identifier, threshold, expected in [('C01', F(119, 50), 'proved'),
                                            ('C02', F(12, 5), 'refuted_by_spectral_existence')]:
        started = perf_counter()
        try:
            cert = projected.weighted_poincare([0, 0, 2], threshold, modes=8, max_modes=32,
                                               bits=44, max_m=16, tolerance=F(1, 10**10))
            valid = projected.verify_inequality(cert, expected_q=[0, 0, 2], expected_threshold=threshold)
            path = destination / (identifier + '.json')
            write(path, cert)
            controls.append({'id': identifier, 'operation': 'weighted_poincare', 'threshold': str(threshold),
                 'expected_status': expected, 'status': cert['status'], 'verified': valid,
                 'passed': valid and cert['status'] == expected,
                 'refutation_has_explicit_function': False if identifier == 'C02' else None,
                 'certificate': str(path.relative_to(HERE)), 'seconds': perf_counter()-started})
        except Exception as error:
            controls.append({'id': identifier, 'passed': False, 'error_type': type(error).__name__, 'error': str(error)})
    try:
        cert = projected.full_ground([0, 0, 2], mean_zero=True, modes=8, max_modes=8,
                      bits=44, max_m=0, tolerance=F(1, 10**10))
        path = destination / 'C03.json'
        write(path, cert)
        valid = projected.verify_full(cert, expected_q=[0, 0, 2], expected_mean_zero=True)
        controls.append({'id': 'C03', 'operation': 'full_ground', 'resource_limit': 'max_m=0',
              'status': cert['status'], 'lower': cert['lower'], 'upper': cert['upper'],
              'verified': valid, 'angular_tail_covered': cert['angular_tail_cannot_improve_best_upper'],
              'passed': valid and cert['status'] == 'certified_bound_open_gap'
                        and cert['angular_tail_cannot_improve_best_upper'] is False,
              'certificate': str(path.relative_to(HERE))})
    except Exception as error:
        controls.append({'id': 'C03', 'passed': False, 'error_type': type(error).__name__, 'error': str(error)})
    return controls


def mutation_audit(certificates):
    checks = []
    sector = certificates['E08']
    for field, value in [('projection', 'none'), ('degree_start', 0), ('azimuth_m', 1),
                         ('eigenvalue_index', 2), ('mean_zero', False), ('lower', '100')]:
        changed = copy.deepcopy(sector)
        changed[field] = value
        checks.append({'name': 'sector_' + field, 'rejected': not projected.verify_sector(changed)})
    full = certificates['E20']
    for name, change in [
        ('full_scope', lambda c: c.update(scope='all_H1_functions_on_unit_S2')),
        ('remove_middle_sector', lambda c: c['sectors'].pop(0)),
        ('false_angular_tail', lambda c: c.update(angular_tail_lower='1000')),
        ('wrong_q', lambda c: c.update(q_coefficients=['0', '0', '3'])),
    ]:
        changed = copy.deepcopy(full)
        change(changed)
        checks.append({'name': name, 'rejected': not projected.verify_full(changed)})
    return checks


def main():
    manifest = json.loads((HERE / 'external_cases.json').read_text())
    digest = hashlib.sha256(json.dumps({'cases': manifest['cases'],
                        'cross_case_relations': manifest['cross_case_relations']}, ensure_ascii=False,
                        sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if digest != manifest['frozen_cases_sha256']:
        raise ValueError('Frozen case hash mismatch; explicit amendment is required')
    run_id = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    destination = HERE / 'results' / 'certificates' / run_id
    source_paths = [HERE / 'projected_spectrum.py', HERE / 'sphere_sectors.py', LEGACY / 'spectral_certifier.py']
    hashes = {str(p.relative_to(HERE.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    records, certificates, legacy_records = [], {}, []
    for case in manifest['cases']:
        started = perf_counter()
        record = {'id': case['id'], 'name': case['name'], 'requested_space': case['requested_space'],
                  'potential_coefficients': case['potential_coefficients'], 'eigenvalue_index': case['eigenvalue_index'],
                  'm': case['m']}
        try:
            cert, verifier, options = invoke(case)
            generated = perf_counter()
            valid = verifier(cert)
            verified = perf_counter()
            path = destination / (case['id'] + '.json')
            write(path, cert)
            check = reference_check(case, cert)
            width_pass = F(cert['exact_width']) <= F(case['acceptance']['maximum_interval_width'])
            record.update(status=cert.get('status', 'certified_sector_enclosure'), verified=valid,
                          lower=cert['lower'], upper=cert['upper'], exact_width=cert['exact_width'],
                          requested_width=case['acceptance']['maximum_interval_width'], width_met=width_pass,
                          reference_check=check, options=options,
                          acceptance='passed' if valid and width_pass and check['compatible'] is not False else 'failed',
                          certificate=str(path.relative_to(HERE)),
                          generation_seconds=generated-started, independent_replay_seconds=verified-generated)
            if 'sectors' in cert:
                record.update(winning_absolute_m=min(cert['sectors'], key=lambda c: F(c['upper']))['azimuth_m'],
                              included_sectors=[c['azimuth_m'] for c in cert['sectors']],
                              angular_tail_covered=cert['angular_tail_cannot_improve_best_upper'],
                              sector_modes=[c['modes'] for c in cert['sectors']])
            certificates[case['id']] = cert
        except Exception as error:
            record.update(acceptance='failed', status='exception_retained',
                          error_type=type(error).__name__, error=str(error))
        records.append(record)
        print(json.dumps({'case': case['id'], 'acceptance': record['acceptance'],
                          'lower_approx': float(F(record['lower'])) if 'lower' in record else None,
                          'width_approx': float(F(record['exact_width'])) if 'exact_width' in record else None,
                          'reference_compatible': record.get('reference_check', {}).get('compatible')}, ensure_ascii=False), flush=True)
        if case['requested_space'] == 'axisymmetric_full':
            comparison = {'case': case['id'], 'old_scope': 'ordinary_m0_without_mean_zero_constraint'}
            try:
                started = perf_counter()
                old = base.certify(case['potential_coefficients'], k=case['eigenvalue_index'], modes=12, bits=44)
                comparison.update(verified=base.verify(old), lower=old['lower'], upper=old['upper'],
                                  exact_width=old['exact_width'], seconds=perf_counter()-started,
                                  intervals_overlap=overlaps(old, record) if 'lower' in record else None)
                path = destination / (case['id'] + '_frozen_v90.json')
                write(path, old)
                comparison['certificate'] = str(path.relative_to(HERE))
            except Exception as error:
                comparison.update(error_type=type(error).__name__, error=str(error), intervals_overlap=None)
            legacy_records.append(comparison)
    cross = relations(records)
    controls = behavioral_controls(destination)
    try:
        mutations = mutation_audit(certificates)
    except Exception as error:
        mutations = [{'name': 'mutation_audit_incomplete', 'rejected': False, 'error': str(error)}]
    ending_hashes = {str(p.relative_to(HERE.parent)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths}
    report = {'run_id': run_id, 'frozen_cases_sha256': digest, 'case_count': len(records),
              'source_sha256': hashes, 'source_changed_during_run': ending_hashes != hashes,
              'cases': records, 'legacy_scope_comparisons': legacy_records, 'cross_case_relations': cross,
              'behavioral_controls': controls, 'mutation_audit': mutations,
              'summary': {'cases_passed': sum(r['acceptance'] == 'passed' for r in records),
                          'cases_failed': sum(r['acceptance'] != 'passed' for r in records),
                          'all_cross_relations_passed': all(r['passed'] is True for r in cross),
                          'all_legacy_comparisons_overlap': all(r.get('verified') and r.get('intervals_overlap') for r in legacy_records),
                          'all_controls_passed': all(r['passed'] for r in controls),
                          'all_mutations_rejected': all(r['rejected'] for r in mutations)},
              'limitations': ['Published decimals are finite-precision compatibility references, not proofs.',
                              'Local 80-mode float linear references are not externally published or certified.',
                              'One serial run is not a performance guarantee or monotonic speedup study.',
                              'A valid existence-based spectral refutation need not include an explicit function.',
                              'No failed, open, unsupported, or incompatible case was removed.']}
    path = HERE / 'results' / ('benchmark_projected_' + run_id + '.json')
    write(path, report)
    write(HERE / 'results' / 'benchmark_projected.json', report)
    print(json.dumps({'report': str(path), 'summary': report['summary']}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
