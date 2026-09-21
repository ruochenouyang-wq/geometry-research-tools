"""Replay saved external-task proofs, binding every proof to its intended input."""
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json
import sys
from time import perf_counter

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
FROZEN = HERE.parent / 'geometry_v36_v90'
sys.path.insert(0, str(FROZEN))
import v03_tail
import positive_search
import compiler_v23
import interpolation_trials as gn


def read(relative):
    return json.loads((HERE/relative).read_text())


def check(condition, label):
    if not condition:
        raise AssertionError(label)


def run():
    started = perf_counter()
    sph = read('spheroidal_results.json')
    expected_cases = {(c, n) for c in (1, 4, 10) for n in (0, 1, 2)}
    check(len(sph['core_cases']) == 9, 'Nine distinct external spectral cases required')
    check({(r['c_squared'], r['n']) for r in sph['core_cases']} == expected_cases,
          'External spectral case inventory')
    count = 0
    for row in sph['core_cases']:
        cert = read(row['selected']['certificate'])
        check(cert['q_coefficients'] == ['0', '0', str(row['c_squared'])], 'Wrong potential')
        check(cert['eigenvalue_index'] == row['n']+1, 'Wrong spectral index')
        check(v03_tail.verify(cert, independent=True), 'Independent dense inertia replay')
        for key in ('lower', 'upper'):
            check(row['selected'][key] == cert[key], 'Summary must match certificate')
        check(F(cert['upper'])-F(cert['lower']) <= F(1, 10**8), 'Spectral target width')
        count += 1
    check(len(sph['positive_function_ground_trials']) == 2, 'Two positive-function trials')
    check({r['c_squared'] for r in sph['positive_function_ground_trials']} == {4, 10},
          'Positive-function input inventory')
    for row in sph['positive_function_ground_trials']:
        cert, baseline = read(row['certificate']), read(row['baseline_certificate'])
        expected_q = ['0', '0', str(row['c_squared'])]
        check(cert['lower_certificate']['q_coefficients'] == expected_q, 'Positive proof input')
        check(baseline['q_coefficients'] == expected_q, 'Baseline input')
        check(positive_search.verify(cert), 'Positive-function replay')
        check(positive_search.verify_lower(baseline), 'Poisson baseline replay')
        check(row['optimized_barta_lower'] == cert['spectral_lower'], 'Positive summary lower')
        check(row['spectral_upper_from_rayleigh'] == cert['spectral_upper'], 'Positive summary upper')
        check(row['projected_poisson_lower'] == baseline['spectral_lower'], 'Baseline summary')
        count += 2
    expected_bases = {
        'legendre_P1_P3': [[F(0), F(1)], [F(0), F(-3, 2), F(0), F(5, 2)]],
        'centered_caps_5_6': [gn.centered_cap(5), gn.centered_cap(6)]}
    interpolation = read('interpolation_results.json')
    check(len(interpolation['cases']) == 2, 'Two interpolation trials')
    check({r['case'] for r in interpolation['cases']} == set(expected_bases), 'GN inventory')
    for row in interpolation['cases']:
        cert = read(row['certificate'])
        check(gn.verify(cert, expected_basis=expected_bases[row['case']]), 'GN exact replay')
        check(F(cert['gap']) <= F(1, 10**10), 'GN family maximum gap')
        check(row['pi_K4_lower'] == cert['lower'] and row['pi_K4_upper'] == cert['upper'],
              'GN summary bounds')
        count += 1
    # These are schema/scope probes, not a fake reduction of an external theorem.
    control = {'geometry': 'unit_sphere', 'function_space': 'full_sphere',
               'eigenvalue_index': 1, 'parameters': ['a'], 'potential': 'a*t',
               'penalty': 'a**2/6', 'domain': {'kind': 'all_real'}, 'threshold': '0'}
    check(compiler_v23.verify(compiler_v23.compile_goal(control)), 'Supported compiler control')
    probes = [
        ('zero_mean_constraint_not_expressible', dict(control, constraints=['integral(u)=0'])),
        ('exponential_barycenter_constraint_not_expressible',
         dict(control, constraints=['integral(exp(u)*x_i)=0 for i=1,2,3'])),
        ('arbitrary_metric_not_expressible', dict(control, geometry='arbitrary_closed_surface'))]
    scope_results = []
    for name, raw in probes:
        try:
            compiler_v23.compile_goal(raw)
        except ValueError as exc:
            scope_results.append({'probe': name, 'rejected': True, 'reason': str(exc)})
        else:
            raise AssertionError('Unsupported probe unexpectedly accepted: '+name)
    manifest = json.loads((FROZEN/'MANIFEST.json').read_text())
    for entry in manifest['entries']:
        data = (FROZEN/entry['path']).read_bytes()
        check(len(data) == entry['bytes'] and sha256(data).hexdigest() == entry['sha256'],
              'Frozen release file differs: '+entry['path'])
    result = {'all_passed': True, 'saved_certificates_replayed': count,
              'binding_to_expected_external_inputs_checked': True,
              'frozen_manifest_files_checked': len(manifest['entries']),
              'frozen_manifest_unchanged': True, 'scope_probes': scope_results,
              'seconds': perf_counter()-started,
              'formal_proof_assistant_checked': False,
              'reference_table_comparison_is_separate_from_proof_verification': True}
    (HERE/'verification_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))
    return result


if __name__ == '__main__':
    run()
