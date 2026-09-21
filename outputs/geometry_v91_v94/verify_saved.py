"""Replay all saved examples and check frozen dependencies, without re-search."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys
from time import perf_counter

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
import projected_spectrum as projected
import research_tool
import benchmark_projected as benchmark


def require(value, message):
    if not value:
        raise AssertionError(message)


def read(relative):
    return json.loads((HERE/relative).read_text())


def main():
    started = perf_counter()
    cases = read('external_cases.json')
    digest = sha256(json.dumps({'cases': cases['cases'],
        'cross_case_relations': cases['cross_case_relations']}, ensure_ascii=False,
        sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    require(digest == cases['frozen_cases_sha256'], 'Frozen external cases changed')
    report = read('results/benchmark_projected.json')
    require(report['frozen_cases_sha256'] == digest, 'Benchmark used a different case set')
    require(len(report['cases']) == len(cases['cases']) == 23, 'External case inventory')
    require([r['id'] for r in report['cases']] == [c['id'] for c in cases['cases']],
            'Case identity/order mismatch')
    count = 0
    for case, record in zip(cases['cases'], report['cases']):
        cert = read(record['certificate'])
        if case['requested_space'] == 'full_sphere_mean_zero':
            valid = projected.verify_full(cert, expected_q=case['potential_coefficients'],
                                           expected_mean_zero=True)
        else:
            valid = projected.verify_sector(cert, expected_q=case['potential_coefficients'],
                expected_m=case['m'], expected_mean_zero=case['requested_space'] == 'axisymmetric_mean_zero',
                expected_k=case['eigenvalue_index'])
        require(valid, 'Certificate/input verification failed: '+case['id'])
        for field in ('lower', 'upper', 'exact_width'):
            require(cert[field] == record[field], 'Recorded endpoints changed: '+case['id'])
        require(F(cert['exact_width']) <= F(case['acceptance']['maximum_interval_width']),
                'External requested width was not reached')
        check = benchmark.reference_check(case, cert)
        require(check['compatible'] is not False, 'External reference incompatible: '+case['id'])
        require(check == record['reference_check'], 'Stored reference check differs')
        count += 1
    old_count = 0
    by_id = {case['id']: case for case in cases['cases']}
    for item in report['legacy_scope_comparisons']:
        cert = read(item['certificate'])
        case = by_id[item['case']]
        require(case['requested_space'] == 'axisymmetric_full', 'Invalid legacy comparison scope')
        require(projected.base.verify(cert), 'Legacy reference certificate failed')
        require(projected.potential(cert['q_coefficients']) == projected.potential(case['potential_coefficients'])
                and cert['eigenvalue_index'] == case['eigenvalue_index'], 'Legacy proof input mismatch')
        corresponding = next(r for r in report['cases'] if r['id'] == item['case'])
        require(benchmark.overlaps(cert, corresponding), 'Old and new ordinary spectra disagree')
        old_count += 1
    require(old_count == 5, 'Legacy comparison inventory')
    controls = report['behavioral_controls']
    require([r['id'] for r in controls] == ['C01', 'C02', 'C03'], 'Control inventory')
    for row, threshold, expected in zip(controls[:2], ['119/50', '12/5'],
                                       ['proved', 'refuted_by_spectral_existence']):
        cert = read(row['certificate'])
        require(projected.verify_inequality(cert, expected_q=[0, 0, 2], expected_threshold=threshold),
                'Control inequality verification')
        require(cert['status'] == expected, 'Wrong inequality control verdict')
    limited = read(controls[2]['certificate'])
    require(projected.verify_full(limited, expected_q=[0, 0, 2], expected_mean_zero=True),
            'Limited-budget proof failed')
    require(len(limited['sectors']) == 1 and limited['lower'] == '2'
            and limited['status'] == 'certified_bound_open_gap'
            and not limited['angular_tail_cannot_improve_best_upper'], 'Limited-budget scope lost')
    relations = benchmark.relations(report['cases'])
    require(all(row['passed'] for row in relations), 'Cross-case mathematical relation failed')
    require(relations == report['cross_case_relations'], 'Stored relations changed')
    service_count = 0
    store = research_tool.CertificateStore()
    for path in sorted(store.directory.glob('*.json')):
        store.get(path.stem)
        service_count += 1
    frozen = HERE.parent/'geometry_v36_v90'
    frozen_manifest = json.loads((frozen/'MANIFEST.json').read_text())
    for item in frozen_manifest['entries']:
        data = (frozen/item['path']).read_bytes()
        require(len(data) == item['bytes'] and sha256(data).hexdigest() == item['sha256'],
                'Frozen V90 dependency changed: '+item['path'])
    result = {'all_passed': True, 'external_case_certificates': count,
              'legacy_comparison_certificates': old_count, 'behavior_control_certificates': 3,
              'service_certificates': service_count,
              'all_top_level_certificates_replayed': count+old_count+3+service_count,
              'frozen_V90_manifest_files_checked': len(frozen_manifest['entries']),
              'external_case_set_sha256': digest, 'cross_case_relations_passed': len(relations),
              'seconds_without_search': perf_counter()-started,
              'formal_proof_assistant_checked': False}
    (HERE/'results/verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    main()
