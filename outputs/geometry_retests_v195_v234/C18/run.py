"""Reproduce the unchanged C18 stress case, ten increments and their evidence."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fractions import Fraction as F
import copy
import io
import json
import time
import unittest
import statistics
import parity_spectrum as p
from support import save, original_case, old_anisotropic as old


Q = {'2,0,0': 1, '0,2,0': 2, '0,0,2': 3}
TARGET = F(1, 10**8)
WALL_BUDGET = 90


def record(version, capability, call, evidence, outcome):
    path = save(f'C18/v{version}_evidence.json', p._json(evidence))
    return {'version': version, 'capability': capability, 'callable': call,
            'evidence': [path], 'outcome': outcome}


def save_certificate(name, certificate):
    return save(f'C18/certificates/{name}.json', certificate)


def main():
    began = time.perf_counter()
    original = original_case('C18')
    from_original = {','.join(map(str, entry['powers'])): entry['coefficient'] for entry in original['input']['q_xyz']}
    assert old.polynomial(from_original) == old.polynomial(Q)
    assert original['input']['target_width'] == str(TARGET)
    save('C18/ORIGINAL_TASK.json', original)
    baseline_rows = []
    for L, bits in ((3, 40), (4, 40), (4, 48)):
        start = time.perf_counter()
        cert = old.certify(Q, L=L, bits=bits)
        elapsed = time.perf_counter()-start
        verified = old.verify(json.loads(json.dumps(cert)), expected_q=Q)
        path = save_certificate(f'baseline_L{L}_bits{bits}', cert)
        baseline_rows.append({'solver': 'frozen_V182', 'L': L, 'dimension': len(cert['mass']), 'bits': bits,
                              'lower': cert['lower'], 'upper': cert['upper'], 'width': cert['exact_width'],
                              'width_decimal': float(F(cert['exact_width'])), 'seconds': elapsed,
                              'target': str(TARGET), 'target_met': F(cert['exact_width']) <= TARGET,
                              'wall_budget_seconds': WALL_BUDGET, 'within_wall_budget': elapsed <= WALL_BUDGET,
                              'verified': verified, 'certificate': path})
    try:
        old.certify(Q, L=4, bits=80)
        bit80 = {'rejected': False}
    except ValueError as error:
        bit80 = {'rejected': True, 'exception': 'ValueError', 'reason': str(error),
                 'interpretation': 'unsupported precision request, not a measured accuracy failure'}
    baseline = {'case': 'C18', 'original_input_preserved': True, 'target_width': str(TARGET),
                'allowance_seconds_per_attempt': WALL_BUDGET, 'attempts': baseline_rows,
                'bits80': bit80,
                'diagnosis': 'L3 and L4 retain the same competing odd-parity harmonics; 40 to48 bits does not remove the finite-tail comparison gap',
                'old_degree_cap': 4}
    save('C18/BASELINE.json', baseline)

    stages = []
    bound = p.potential_range(Q)
    stages.append(record(225, 'Exact simplex Bernstein potential range on the sphere',
                         'parity_spectrum.potential_range / verify_range', bound,
                         'Original three-axis potential has globally certified range [1,3]'))
    example = {'0,0,0': 1, '1,0,0': 2, '0,1,1': 3, '1,1,1': 4, '2,2,0': 5}
    stages.append(record(226, 'Exact sign-flip symmetry and eight character projections',
                         'parity_spectrum.reflection_symmetry / parity_project',
                         {'C18': p.reflection_symmetry(Q), 'invalid_xy': p.reflection_symmetry({'1,1,0': 1}),
                          'projectors': [{'parity': list(c), 'polynomial': p.parity_project(example, c)} for c in p.PARITIES]},
                         'C18 accepts all eight blocks; xy exposes only partial symmetry and cannot be split this way'))
    dimensions = []
    for degree in range(7):
        counts = [len(p.degree_basis(c, degree)) for c in p.PARITIES]
        dimensions.append({'degree': degree, 'character_dimensions': counts, 'sum': sum(counts)})
    stages.append(record(227, 'Complete character harmonic bases and exact positive masses',
                         'parity_spectrum.degree_basis / parity_basis',
                         {'degree_counts': dimensions, 'winner_L5': p.parity_basis((1, 0, 0), 5),
                          '000_L4': p.parity_basis((0, 0, 0), 4)},
                         'All degree counts sum to 2l+1; the 000 block removes the constant only'))
    stages.append(record(228, 'Coupled forms within a character and exact cross-block checks',
                         'parity_spectrum.block_form / cross_block_form',
                         {'winner_L3': p.block_form(Q, (1, 0, 0), 3),
                          'C18_cross_100_010': p.cross_block_form(Q, (1, 0, 0), (0, 1, 0)),
                          'xy_cross_100_010': p.cross_block_form({'1,1,0': 1}, (1, 0, 0), (0, 1, 0), 1),
                          'trial_x_upper': old.rayleigh(Q, {'1,0,0': 1})},
                         'The original x trial gives18/5; C18 cross blocks vanish and xy has exact nonzero1/15 coupling'))
    stages.append(record(229, 'Complete same-character tail columns and correct first omitted degree',
                         'parity_spectrum.tail_couplings / next_omitted_degree',
                         {'winner_L3': p.tail_couplings(Q, (1, 0, 0), 3),
                          'winner_L4': p.tail_couplings(Q, (1, 0, 0), 4),
                          'zero_character_L3': p.tail_couplings(Q, (0, 0, 0), 3)},
                         'Odd winner omits degree5 after either L3 orL4; even000 omits degree4 afterL3'))
    block = p.block_ground(Q, (1, 0, 0), L=3, bits=40)
    block_path = save_certificate('winner_L3', block)
    floors = [p.analytic_floor({}, c) for c in p.PARITIES]
    stages.append(record(230, 'Infinite-tail Schur enclosure for a reflection character',
                         'parity_spectrum.block_ground / verify_block / analytic_floor',
                         {'block_certificate': block_path, 'block_lower': block['lower'], 'block_upper': block['upper'],
                          'unperturbed_all_character_floors': floors},
                         'All block endpoints use exact inertia with same-character graded tails and non-unit masses'))
    initial = p.full_ground(Q, L=3, bits=40, target=TARGET)
    initial_path = save_certificate('initial_full_L3', initial)
    stages.append(record(231, 'All-eight-character global ground merge with valid pruning',
                         'parity_spectrum.full_ground / verify_full',
                         {'certificate': initial_path, 'lower': initial['lower'], 'upper': initial['upper'],
                          'width': initial['exact_width'], 'status': initial['status'],
                          'characters': [{'parity': c['parity'], 'format': c['format'], 'lower': c['lower']}
                                         for c in initial['characters']]},
                         'Three competitive odd blocks computed; five entire-character analytic floors preserved'))
    diagnosis = p.error_diagnosis(Q, (1, 0, 0), L=3, bits=32, extra_bits=24)
    for label in ('precision_low', 'precision_high', 'larger_block'):
        save_certificate('diagnosis_'+label, diagnosis[label])
    stages.append(record(232, 'Measured precision-versus-tail diagnosis and targeted same-parity refinement',
                         'parity_spectrum.error_diagnosis', diagnosis,
                         '32 to56 bits at fixed degree scarcely changes width; degree3 to5 closes the target'))
    cache = p.FormCache()
    reuse = []
    for label, potential, character, bits in [('first', Q, (1, 0, 0), 32),
             ('precision_reuse', Q, (1, 0, 0), 48), ('other_character', Q, (0, 1, 0), 32),
             ('other_potential', {'2,0,0': 1}, (1, 0, 0), 32)]:
        start = time.perf_counter()
        certificate = p.block_ground(potential, character, L=3, bits=bits, cache=cache)
        reuse.append({'case': label, 'seconds': time.perf_counter()-start,
                      'verified': p.verify_block(certificate, potential, character),
                      'work': copy.deepcopy(cache.stats), 'lower': certificate['lower'], 'upper': certificate['upper']})
    stages.append(record(233, 'Input-bound exact finite-form and potential-range reuse',
                         'parity_spectrum.FormCache / block_ground(cache=...)', reuse,
                         'A second precision run reuses the matrix; changing parity or potential builds another matrix'))
    drivers = []
    for repeat in range(3):
        driver = p.adaptive_ground(Q, target=TARGET, start_L=3, max_L=9, bits=40,
                                   wall_budget_seconds=WALL_BUDGET)
        verified = p.verify_result(json.loads(json.dumps(driver)), Q, TARGET)
        path = save_certificate(f'C18_adaptive_repeat{repeat+1}', driver)
        drivers.append({'repeat': repeat+1, 'certificate': path, 'verified': verified,
                        'status': driver['status'], 'width': driver['certificate']['exact_width'],
                        'seconds': driver['elapsed_seconds'], 'within_wall_budget': driver['within_wall_budget'],
                        'work': driver['work'], 'trace': driver['trace']})
    # Stable name for the exact original acceptance proof.
    final = json.loads((Path(__file__).resolve().parents[1]/drivers[0]['certificate']).read_text())
    final_path = save_certificate('C18_original_target', final)
    quartic_q = {'4,0,0': '1/2', '2,2,0': '-1/3', '0,0,2': 1}
    quartic = p.adaptive_ground(quartic_q, target=F(1, 10**7), start_L=3, max_L=7)
    quartic_path = save_certificate('general_even_quartic', quartic)
    open_budget = p.adaptive_ground({'4,0,0': -100}, start_L=1, max_L=1, bits=16)
    open_path = save_certificate('insufficient_budget_open', open_budget)
    stages.append(record(234, 'Requested-width full-space driver with competitor-only refinement and honest budgets',
                         'parity_spectrum.adaptive_ground / verify_result',
                         {'original_case': 'C18', 'original_certificate': final_path, 'repeats': drivers,
                          'quartic_certificate': quartic_path, 'budget_open_certificate': open_path,
                          'status': final['status'], 'requested_width': final['requested_width'],
                          'lower': final['certificate']['lower'], 'upper': final['certificate']['upper']},
                         'Unchanged C18 target1e-8 met in three runs; insufficient max_L remains openly unclosed'))

    suite = unittest.defaultTestLoader.loadTestsFromName('test_parity_spectrum')
    stream = io.StringIO()
    tested = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    (Path(__file__).resolve().parent/'tests.log').write_text(stream.getvalue())
    saved_replays = []
    for path in sorted((Path(__file__).resolve().parent/'certificates').glob('*.json')):
        c = json.loads(path.read_text())
        form = c['format']
        if form == old.FORMAT:
            ok = old.verify(c)
        elif form == p.BLOCK_FORMAT:
            ok = p.verify_block(c)
        elif form == p.FULL_FORMAT:
            ok = p.verify_full(c)
        elif form == p.DRIVER_FORMAT:
            ok = p.verify_result(c)
        else:
            ok = False
        saved_replays.append({'path': str(path.relative_to(Path(__file__).resolve().parents[1])), 'verified': ok})
    validation = {'all_passed': tested.wasSuccessful() and all(row['verified'] for row in saved_replays)
                  and all(row['status'] == 'target_met' and row['within_wall_budget'] for row in drivers)
                  and bit80['rejected'], 'tests_run': tested.testsRun, 'failures': len(tested.failures),
                  'errors': len(tested.errors), 'saved_certificates_replayed': len(saved_replays),
                  'original_input_preserved': True, 'original_target_width_preserved': str(TARGET),
                  'formal_assistant_checked': False}
    elapsed_total = time.perf_counter()-began
    results = {'case': 'C18', 'potential': p.a.encode(p.a.polynomial(Q)), 'scope': p.SCOPE,
               'target_width': str(TARGET), 'baseline': baseline_rows, 'unsupported_bits80': bit80,
               'new_runs': drivers, 'median_new_seconds': statistics.median(row['seconds'] for row in drivers),
               'new_lower': final['certificate']['lower'], 'new_upper': final['certificate']['upper'],
               'new_width': final['certificate']['exact_width'], 'new_width_decimal': float(F(final['certificate']['exact_width'])),
               'width_improvement_vs_old_L4_bits48': float(F(baseline_rows[-1]['width'])/F(final['certificate']['exact_width'])),
               'original_certificate': final_path, 'validation': validation, 'saved_replays': saved_replays,
               'whole_reproduction_elapsed_seconds': elapsed_total,
               'timing_note': 'One Python process; exact moments may be warm. Each attempt has the same90-second allowance; old methods stop at their degree cap. Driver times include final replay; no model speed claim.'}
    save('C18/STAGES.json', stages)
    save('C18/RESULTS.json', results)
    save('C18/VALIDATION.json', validation)
    print(json.dumps({'validation': validation, 'new_width': results['new_width_decimal'],
                      'median_new_seconds': results['median_new_seconds'],
                      'whole_reproduction_seconds': elapsed_total}, indent=2))
    if not validation['all_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
