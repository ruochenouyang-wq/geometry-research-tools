"""Reproduce the four unchanged challenges and replay the frozen twenty-case bank.

Only this release directory is written. Each new original case runs in its own
process with the original 90-second outer limit. Mathematical comparisons are
exact rationals; elapsed time is an observation, not a verified theorem.
"""
import importlib
import json
from math import comb
from pathlib import Path
import subprocess
import sys
import time

from support import ROOT, PREVIOUS, F, save, source_hash, original_case, old_portal
import research_tool as tool

CASES = ('C12', 'C14', 'C16', 'C18')


def previous_integrity():
    snapshot = json.loads((ROOT/'PREVIOUS_AT_START.json').read_text())
    current = {str(p.relative_to(PREVIOUS)): source_hash(p)
               for p in PREVIOUS.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    differences = [p for p in sorted(set(current) | set(snapshot['files']))
                   if current.get(p) != snapshot['files'].get(p)]
    assert not snapshot['preexisting_manifest_differences']
    assert not differences, differences
    return {'unchanged_files': len(current), 'differences': differences,
            'manifest_sha256': source_hash(PREVIOUS/'MANIFEST.json')}


def replay_previous():
    folder = PREVIOUS/'cross_challenge'
    prior = json.loads((folder/'result.json').read_text())
    assert prior['challenge_sha256'] == source_hash(PREVIOUS/'challenges.json')
    import wide_potential as wide
    import global_parameter as glob
    import anisotropic as xyz
    special = {wide.CHEAP_FORMAT: wide.verify_cheap_ground,
               'cheap_global_affine_sphere_v135': glob.verify_cheap_objective,
               'anisotropic_cheap_ground_v175': xyz.verify_cheap}
    rows = []
    for case in prior['results']:
        original = original_case(case['id'])
        assert case['original_input'] == original['input']
        assert case['original_acceptance'] == original['acceptance']
        for item in case['certificates']:
            p = folder/item['path']
            assert source_hash(p) == item['sha256']
            proof = json.loads(p.read_text())
            assert special.get(proof['format'], old_portal.verify)(proof) is True
            rows.append({'case': case['id'], 'path': str(p.relative_to(PREVIOUS)),
                         'sha256': item['sha256'], 'verified': True})
    assert len(rows) == 56
    return {'cases_bound_to_original': 20, 'certificates_verified': len(rows),
            'previous_complete_cases': sum(r['status'] in ('solved', 'narrow_target_met')
                                           for r in prior['results']),
            'proofs': rows}


def run_case(identifier):
    start = time.perf_counter()
    original = original_case(identifier)
    inp = original['input']
    service = tool.Service(ROOT/'service_certificates')
    rows = []

    def call(request):
        brief = service.call(request)
        proof = service.store.get(brief['certificate_id'])
        assert tool.verify(proof)
        rows.append({'request': request, 'result': brief})
        return proof['certificate']

    if identifier == 'C12':
        assert inp == {'q': ['0'], 'orthogonality': ['1', 'z']}
        c = call({'op': 'constraints', 'q': {},
                  'constraints': [{'0,0,0': '1'}, {'0,0,1': '1'}], 'L': 2})
        assert F(c['lower']) == F(c['upper']) == 2
        from support import old_constraints
        axis = json.loads((PREVIOUS/'cross_challenge/certs/C12_axisymmetric_exact_6.json').read_text())
        assert old_constraints.verify_sector(axis)
        assert F(axis['lower']) == F(axis['upper']) == 6
        for coordinate in ('1,0,0', '0,1,0'):
            trial = tool.constraints.trial_enclosure({}, c['constraints'], {coordinate: '1'})
            assert tool.constraints.verify(trial)
            assert F(trial['lower']) == F(trial['upper']) == 2
        detail = {'full_ground': '2', 'separately_axisymmetric_ground': '6',
                  'admissible_trials': ['x', 'y']}
    elif identifier == 'C14':
        c = call({'op': 'spectrum', 'q': inp['q'], 'spaces': inp['spaces'],
                  'k': 9, 'quantities': ['ordered', 'count', 'trace']})
        assert c['all_requested_targets_met']
        detail = {}
        for row in c['results']:
            space = row['request']['space']
            count, trace = (8, 20) if space == 'full_mean_zero' else (9, 27)
            a = row['answers']
            assert a['count']['count_lower'] == a['count']['count_upper'] == count
            assert F(a['trace']['lower']) == F(a['trace']['upper']) == trace
            detail[space] = {'negative_count': count, 'negative_trace': trace}
    elif identifier == 'C16':
        assert inp['n_values'] == [4, 16, 64]
        detail = []
        for n in inp['n_values']:
            raw = [F(comb(n, j)) for j in range(n+1)]
            raw[0] -= F(2**n, n+1)
            c = call({'op': 'cap_original', 'n': n,
                      'raw_t_coefficients': list(map(str, raw))})
            moments = c['evaluation']['ratio']['moments']
            mass = 2**(2*n)*(F(1, 2*n+1)-F(1, (n+1)**2))
            energy = F(n*2**(2*n-1), 2*n+1)
            quartic = sum(F(comb(4,j))*(-F(2**n,n+1))**(4-j)*F(2**(j*n),j*n+1)
                          for j in range(5))
            assert F(moments['mean']) == 0
            assert F(moments['mass']) == mass
            assert F(moments['energy']) == energy
            assert F(moments['quartic']) == quartic
            ratio = quartic/(mass*energy)
            assert F(c['evaluation']['ratio']['probability_ratio']) == ratio
            assert F(c['evaluation']['ratio']['pi_times_area_K']) == ratio/4
            detail.append({'n': n, 'original_coefficients_checked': len(raw),
                           'probability_ratio': str(ratio), 'pi_times_area_K': str(ratio/4),
                           'all_four_moments_independently_checked': True})
    else:
        assert inp['mean_zero'] is True
        q = {','.join(map(str, t['powers'])): t['coefficient'] for t in inp['q_xyz']}
        c = call({'op': 'parity_ground', 'q': q, 'target': inp['target_width'],
                  'start_L': 3, 'max_L': 9, 'bits': 40, 'wall_budget_seconds': 90})
        full = c['certificate']
        assert len(full['characters']) == 8
        assert F(full['exact_width']) <= F(inp['target_width'])
        assert F(3) <= F(full['lower']) <= F(full['upper']) <= F(18,5)
        detail = {k: full[k] for k in ('lower', 'upper', 'exact_width')}
        detail['complete_reflection_blocks'] = 8
        detail['max_computed_dimension'] = max(c['work']['built_dimensions'])
    elapsed = time.perf_counter()-start
    assert elapsed < 90
    result = {'case': identifier, 'original': original, 'original_target_met': True,
              'unchanged_input_scope_and_target': True, 'detail': detail,
              'elapsed_including_store_and_replay_seconds': elapsed,
              'outer_process_limit_seconds': 90, 'service_calls': rows}
    save(f'acceptance/{identifier}.json', result)
    return result


def stages():
    all_rows = []
    for identifier in CASES:
        rows = json.loads((ROOT/identifier/'STAGES.json').read_text())
        assert len(rows) == 10
        for row in rows:
            assert row['outcome'] and row['outcome'] != 'pending_run'
            names = row['callable'].split(' / ')
            module_name = names[0].split('.')[0]
            module = importlib.import_module(module_name)
            for name in names:
                symbol = name.split('(', 1)[0].split('.')[-1]
                assert callable(getattr(module, symbol))
            for evidence in row['evidence']:
                assert (ROOT/evidence).is_file(), evidence
            all_rows.append(dict(row, case=identifier))
    assert [r['version'] for r in all_rows] == list(range(195, 235))
    return all_rows


def main():
    previous_integrity()
    for identifier in CASES:
        subprocess.run([sys.executable, '-B', str(Path(__file__).resolve()), '--case', identifier],
                       check=True, timeout=90)
    prior = replay_previous()
    stage_rows = stages()
    result = {'original_problem_bank_sha256': source_hash(PREVIOUS/'challenges.json'),
              'four_retested_targets_met': True,
              'cases': [json.loads((ROOT/f'acceptance/{c}.json').read_text()) for c in CASES],
              'previous_replay': prior, 'previous_integrity': previous_integrity(),
              'executable_increments': len(stage_rows), 'versions': [195, 234],
              'twenty_case_status': {'previously_complete_replayed': 16,
                                    'newly_completed_original_cases': 4,
                                    'completed_original_targets': 20},
              'scope_note': 'A fixed calibration/stress bank; not 20 open problems solved.'}
    save('ROOT_ACCEPTANCE.json', result)
    save('ITERATIONS.json', stage_rows)
    print(json.dumps({'original_retests': '4/4', 'prior_proofs_replayed': 56,
                      'executable_increments': 40, 'old_files_unchanged': 440}))


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == '--case': run_case(sys.argv[2])
    else: main()
