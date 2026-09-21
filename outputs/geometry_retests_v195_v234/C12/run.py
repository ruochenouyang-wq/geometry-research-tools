"""Reproduce original C12 retest, ten capabilities, and full-scope certificates."""
from fractions import Fraction as F
from pathlib import Path
from datetime import datetime, timezone
import io
import json
import sys
import tempfile
import time
import unittest

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from support import save, digest, original_case, old_portal, old_constraints, source_hash, PREVIOUS
import global_constraints as g
import test_global_constraints as tests

C, X, Y, Z = tests.C, tests.X, tests.Y, tests.Z


def main():
    original = original_case('C12')
    if original['input'] != {'q': ['0'], 'orthogonality': ['1', 'z']}:
        raise AssertionError('Original C12 changed; do not silently reuse a stale mapping')
    old_request = {'op': 'constraints', 'q': ['0'], 'threshold': '2', 'cutoff': -1,
                   'rows': [['1'], ['0', '1']], 'm': 0,
                   'scope': 'all_real_H1_functions_on_unit_S2'}
    with tempfile.TemporaryDirectory() as directory:
        try:
            old_response = old_portal.Service(directory).call(old_request)
        except (ValueError, TypeError) as error:
            old_response = {'status': 'rejected', 'exception': type(error).__name__, 'message': str(error)}
    axisymmetric = old_constraints.certify_moment_sector([0], [[1], [0, 1]], m=0, modes=6, bits=40)
    if not old_constraints.verify(axisymmetric):
        raise AssertionError('Old axisymmetric certificate failed its own verifier')
    old_path = save('C12/baseline/axisymmetric_only.json', axisymmetric)
    save('C12/baseline.json', {'original_case': original, 'old_full_request': old_request,
        'old_full_response': old_response, 'old_axisymmetric_certificate': old_path,
        'old_axisymmetric_decimal_enclosure': axisymmetric['decimal_enclosure'],
        'old_axisymmetric_scope': axisymmetric['scope'],
        'failure_cause': 'General moments previously acted on one real m component; adding a global quantifier was unsupported. The old value near 6 is correct only in the separately declared axisymmetric subspace.',
        'retested_utc': datetime.now(timezone.utc).isoformat()})
    stages, cases, certificate_paths = [], [], []

    def certificate(name, thunk):
        start = time.perf_counter()
        result = thunk()
        elapsed = time.perf_counter()-start
        replayed = g.verify(json.loads(json.dumps(result)))
        if not replayed:
            raise AssertionError('New certificate replay failed: '+name)
        path = save('C12/certificates/'+name+'.json', result)
        certificate_paths.append(path)
        row = {'case': name, 'path': path, 'verified': replayed, 'generation_seconds': elapsed,
               'format': result['format'], 'scope': result['scope']}
        for key in ('lower', 'upper', 'exact_width', 'decimal_enclosure', 'status'):
            if key in result:
                row[key] = result[key]
        cases.append(row)
        return result

    def stage(version, capability, call, thunk, outcome):
        start = time.perf_counter()
        result = thunk()
        path = save('C12/stages/V%d.json' % version,
                    {'version': version, 'callable': call, 'wall_seconds': time.perf_counter()-start,
                     'output': result, 'outcome': outcome})
        stages.append({'version': version, 'capability': capability, 'callable': call,
                       'evidence': [path], 'outcome': outcome})
        return result

    stage(195, 'Validate all global moments of a real trial and a justified Poincare lower bound',
          'global_constraints.trial_enclosure',
          lambda: certificate('C12_admissible_x', lambda: g.trial_enclosure({}, [C, Z], X, 1)),
          'The original full-sphere problem is sharply enclosed by [2,2]; x is orthogonal to both 1 and z.')
    stage(196, 'Assemble every real harmonic including constant and sine/cosine modes',
          'global_constraints.full_harmonics', lambda: g.full_harmonics(2),
          'All 9 real harmonics of degrees 0..2 have exact Gram and independently integrated gradient matrices.')
    stage(197, 'Mass-weighted nullspace of complete cross-m and cross-component moments',
          'global_constraints.constraint_nullspace',
          lambda: g.constraint_nullspace([C, {'1,0,0': 1, '0,0,1': 1}], 2),
          'The x+z moment couples different m blocks; no constraint support is dropped.')
    stage(198, 'Congruence of the full coupled energy and non-diagonal mass matrix',
          'global_constraints.reduced_forms',
          lambda: g.reduced_forms({'1,1,0': '1/2'}, [{'0,0,0': 1, '1,0,0': 1, '0,0,1': 1}], 2),
          'Dense T^T M T retains every off-diagonal entry; mixed xyz potential is integrated exactly.')
    stage(199, 'Transform every finite-to-tail harmonic coupling through the global nullspace',
          'global_constraints.transformed_tail',
          lambda: g.transformed_tail({'1,0,0': 1, '0,0,1': 1}, [C, Z], 2),
          'All 7 real degree-3 tail modes are represented, including sine components.')
    nonconstant = stage(200, 'Full-sphere constrained Schur spectral enclosure with entire infinite tail',
          'global_constraints.certify', lambda: certificate('nonconstant_x_plus_z',
              lambda: g.certify({'1,0,0': 1, '0,0,1': 1}, [C, Z], L=2, bits=36, snap=False)),
          'A nonconstant, nonaxisymmetric problem is certified on the complete constrained H1 space.')
    c12 = stage(201, 'Snap to an exact Laplace level only when both exact inertia inequalities hold',
          'global_constraints.exact_level', lambda: certificate('C12_original_full_exact',
              lambda: g.exact_level({}, [C, Z], L=2)),
          'Original C12 full-space sharp value is exactly 2; no axisymmetric certificate is reused as global evidence.')
    stage(202, 'Transfer verified spectra between identical global constraint row spaces',
          'global_constraints.equivalent_constraints', lambda: certificate('C12_scaled_redundant_transfer',
              lambda: g.equivalent_constraints(c12, [{'0,0,0': -3}, {'0,0,1': 2}, {'0,0,1': 4}, {}])),
          'Scaling, duplicates and zero rows preserve the complete constrained function space and reuse the spectral proof.')
    stage(203, 'Increase the complete harmonic head while preserving all constraint support',
          'global_constraints.adaptive', lambda: certificate('adaptive_all_linear_removed',
              lambda: g.adaptive({}, [C, X, Y, Z], start_L=1, max_L=2)),
          'L=1 has no admissible head direction; L=2 resolves the unchanged problem and certifies exact 6.')
    stage(204, 'Prove a global weighted inequality or provide an explicit rational violating function',
          'global_constraints.inequality', lambda: certificate('C12_refute_threshold3',
              lambda: g.inequality({}, [C, Z], 3, L=2)),
          'The overstrong threshold 3 is refuted with an explicitly verified admissible real polynomial trial.')

    certificate('C12_admissible_y', lambda: g.trial_enclosure({}, [C, Z], Y, 1))
    certificate('C12_prove_threshold2', lambda: g.inequality({}, [C, Z], 2, L=2))
    certificate('cross_m_exact2', lambda: g.certify({}, [C, {'1,0,0': 1, '0,0,1': 1}], 2))
    certificate('all_linear_exact6', lambda: g.certify({}, [C, X, Y, Z], 2))
    certificate('no_constraints_exact0', lambda: g.certify({}, [], 0))
    certificate('mixed_constant_not_mean_zero', lambda: g.certify({}, [{'0,0,0': 1, '1,0,0': 1}], 2, bits=40))
    certificate('mixed_all_components_six_over17', lambda: g.certify({},
        [{'0,0,0': 1, '1,0,0': 1, '0,1,0': 2, '0,0,1': 3}], 2, bits=40))
    certificate('nonconstant_xy', lambda: g.certify({'1,1,0': '1/2'}, [C, Z], 3, bits=32))
    certificate('schur_tail_explicit_refutation', lambda: g.inequality(
        {'1,0,0': 1, '0,0,1': 1}, [C, Z], F(39, 20), L=1, bits=32))
    certificate('budget_failure_retained', lambda: g.adaptive({}, [C, X, Y, Z], start_L=1, max_L=1))
    certificate('open_gap_retained', lambda: g.adaptive({'1,0,0': 1, '0,0,1': 1}, [C, Z],
        tolerance=F(1, 10**12), start_L=1, max_L=1, bits=16))
    certificate('undetermined_threshold', lambda: g.inequality_record(nonconstant,
        (F(nonconstant['lower'])+F(nonconstant['upper']))/2))
    save('C12/original_binding.json', {'original_case': original,
        'source_challenge_book': str(PREVIOUS/'challenges.json'),
        'source_challenge_book_sha256': source_hash(PREVIOUS/'challenges.json'),
        'original_potential_to_xyz': {'original_q': ['0'], 'xyz_potential': {}},
        'original_moment_mapping': [
            {'original_name': '1', 'xyz_polynomial': {'0,0,0': '1'},
             'harmonic': {'degree': 0, 'm': 0, 'part': 'real', 'coefficient': '1'}},
            {'original_name': 'z', 'xyz_polynomial': {'0,0,1': '1'},
             'harmonic': {'degree': 1, 'm': 0, 'part': 'real', 'coefficient': '1'}}],
        'full_scope': g.SCOPE, 'certificate_path': 'C12/certificates/C12_original_full_exact.json',
        'certificate_content_sha256': digest(c12),
        'original_acceptance_met': c12['lower'] == c12['upper'] == '2',
        'axisymmetric_acceptance_separately_verified': F(axisymmetric['lower']) <= 6 <= F(axisymmetric['upper'])})
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(
        unittest.defaultTestLoader.loadTestsFromModule(tests))
    (ROOT/'C12/tests.log').write_text(stream.getvalue())
    validation = {'all_passed': result.wasSuccessful() and all(row['verified'] for row in cases),
        'tests_run': result.testsRun, 'failures': len(result.failures), 'errors': len(result.errors),
        'certificates_replayed': len(certificate_paths), 'original_C12_full_scope_exact2': c12['lower'] == c12['upper'] == '2',
        'prior_axisymmetric_scope_kept_separate': True,
        'formal_assistant_checked': False}
    save('C12/STAGES.json', stages)
    save('C12/RESULTS.json', {'original_case_id': 'C12', 'cases': cases,
        'certificates': certificate_paths, 'validation': validation,
        'preserved_failures': ['budget_failure_retained', 'open_gap_retained', 'undetermined_threshold']})
    save('C12/VALIDATION.json', validation)
    print(json.dumps(validation, indent=2))
    if not validation['all_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
