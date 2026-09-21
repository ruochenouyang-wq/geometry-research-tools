"""Independent cycle-5 review: global quantifiers, exact extrema and replay.

Run: python3 -B -m unittest -v test_global_review
The analytical references use C(c)=2+c and elementary square completions.
"""
import copy
import json
from pathlib import Path
import unittest
from unittest.mock import patch
from fractions import Fraction as F

import global_parameter as G
import witness as W

HERE = Path(__file__).resolve().parent
QUARTIC = {'q0': [0], 'direction': [0], 'penalty': [0, 0, -1, 0, 1], 'domain': [-1, 1]}
SPHERE = {'q0': [0, 0, 2], 'direction': [0, 1], 'penalty': [0, 0, F(1, 10)], 'domain': [-2, 2]}


def read(name):
    return json.loads((HERE/'round05'/name).read_text())


def objective(q0=(0,), direction=(0,), penalty=(0,), domain=(-1, 1)):
    return dict(q0=list(q0), direction=list(direction), penalty=list(penalty), domain=list(domain))


class GlobalIndependentReview(unittest.TestCase):
    def test_three_point_sampling_is_not_a_global_lower_bound(self):
        # 2+s^4-s² = 7/4+(s²-1/2)², so the exact min is 7/4.
        self.assertTrue(all(2+s**4-s**2 == 2 for s in [-1, 0, 1]))
        self.assertEqual(2+F(1, 2)**4-F(1, 2)**2, F(29, 16))
        cert = G.branch_and_bound(QUARTIC, max_leaves=1, tolerance=F(1, 10**10), modes=2)
        self.assertEqual(F(cert['lower']), F(7, 4))
        self.assertEqual(F(cert['upper']), 2)
        self.assertEqual(cert['status'], 'certified_bound_open_gap')
        bad = copy.deepcopy(cert)
        bad.update(lower='2', gap='0', status='certified_target_met')
        self.assertFalse(G.verify(bad))

    def test_chord_then_penalty_minimization_has_correct_sign(self):
        # C(s)=2+s, penalty=(s-1/4)²; sum=2+(s+1/4)².
        obj = objective(direction=[1], penalty=[F(1, 16), -F(1, 2), 1])
        left, right = [G.feasible_point(obj, s, modes=2) for s in [-1, 1]]
        proof = G.concavity_cell(obj, [-1, 1], left, right)
        self.assertEqual(list(map(F, proof['chord_plus_penalty'])), [F(33, 16), F(1, 2), 1])
        self.assertEqual(F(proof['lower']), 2)
        global_cert = G.branch_and_bound(obj, max_leaves=6, tolerance=F(1, 10**10), modes=2)
        self.assertEqual((F(global_cert['lower']), F(global_cert['upper'])), (2, 2))
        self.assertEqual(F(global_cert['incumbent_parameter']), -F(1, 4))
        swapped = copy.deepcopy(left)
        swapped['parameter'] = '1'
        with self.assertRaises(ValueError):
            G.concavity_cell(obj, [-1, 1], swapped, right)

    def test_lipschitz_lower_uses_farthest_distance_and_global_penalty_floor(self):
        obj = objective(direction=[1], penalty=[0, 0, 1])
        point = G.feasible_point(obj, 0, modes=2)
        proof = G.lipschitz_cell(obj, [-1, 1], [point])
        self.assertEqual(F(proof['lipschitz_constant']), 1)
        self.assertEqual(F(proof['lower']), 1)
        # True global min of 2+s+s² is 7/4. The coarse bound is valid and loose.
        self.assertLess(F(proof['lower']), F(7, 4))
        wrong = copy.deepcopy(point)
        wrong['parameter'] = '1/2'
        with self.assertRaises(ValueError):
            G.lipschitz_cell(obj, [-1, 1], [wrong])

    def test_complete_cover_rejects_omitted_and_duplicated_intervals(self):
        cert = read('global_quartic.json')
        self.assertGreater(len(cert['cells']), 1)
        intervals = [list(map(F, c['interval'])) for c in cert['cells']]
        self.assertEqual(intervals[0][0], -1)
        self.assertEqual(intervals[-1][1], 1)
        self.assertTrue(all(b == c for (_, b), (c, _) in zip(intervals, intervals[1:])))
        for mutation in ('omit', 'duplicate', 'tiny_gap'):
            bad = copy.deepcopy(cert)
            if mutation == 'omit':
                bad['cells'].pop(len(bad['cells'])//2)
            elif mutation == 'duplicate':
                bad['cells'].append(copy.deepcopy(bad['cells'][0]))
            else:
                bad['cells'][0]['interval'][1] = str(F(bad['cells'][0]['interval'][1])-F(1, 10**15))
            self.assertFalse(G.verify(bad), mutation)

    def test_two_irrational_minimizers_are_both_enclosed(self):
        cert = read('global_quartic.json')
        result = G.minimizer_enclosure(cert)
        intervals = [list(map(F, v)) for v in result['all_minimizers_in']]
        # Exact tests for containment of +/-sqrt(1/2); no floating roots used.
        plus = any(b >= 0 and b*b >= F(1, 2) and (a <= 0 or a*a <= F(1, 2)) for a, b in intervals)
        minus = any(a <= 0 and a*a >= F(1, 2) and (b >= 0 or b*b <= F(1, 2)) for a, b in intervals)
        self.assertTrue(plus and minus)
        self.assertFalse(result['uniqueness_claimed'])
        bad = copy.deepcopy(cert)
        bad['uniqueness_claimed'] = True
        self.assertFalse(G.verify(bad))
        bad = copy.deepcopy(cert)
        bad['minimizer_intervals'] = [v for v in bad['minimizer_intervals'] if F(v[1]) < 0]
        self.assertFalse(G.verify(bad))

    def test_constant_objective_preserves_every_minimizer_at_equality(self):
        obj = objective(q0=[3], direction=[0], penalty=[-1], domain=[-2, 5])
        cert = G.branch_and_bound(obj, max_leaves=4, modes=2)
        self.assertEqual((F(cert['lower']), F(cert['upper'])), (4, 4))
        self.assertTrue(all(c['status'] == 'live' for c in cert['cells']))
        self.assertEqual(G.minimizer_enclosure(cert)['all_minimizers_in'], [['-2', '5']])
        bad = copy.deepcopy(cert)
        bad['cells'][0]['status'] = 'excluded'
        bad['minimizer_intervals'] = []
        self.assertFalse(G.verify(bad))

    def test_closed_boundary_minimum_and_open_domain_rejection(self):
        obj = objective(direction=[1])
        cert = G.branch_and_bound(obj, max_leaves=1, modes=2)
        self.assertEqual((F(cert['lower']), F(cert['upper'])), (1, 1))
        self.assertEqual(cert['incumbent_parameter'], '-1')
        self.assertTrue(any(F(a) <= -1 <= F(b) for a, b in cert['minimizer_intervals']))
        wrong = copy.deepcopy(obj)
        wrong['left_open'] = True
        with self.assertRaises(ValueError):
            G.canonical_objective(wrong)
        bad = copy.deepcopy(cert)
        bad['objective']['domain'] = ['-1/2', '1']
        self.assertFalse(G.verify(bad))

    def test_polynomial_proofs_bind_both_domain_and_penalty(self):
        p = [0, 0, -1, 0, 1]
        evidence = G.penalty_enclosure(p, [-1, 0])
        self.assertEqual(F(evidence['lower']), -F(1, 4))
        # Symmetry gives equal numerical extrema on [0,1], but not permission
        # to reuse a certificate proving a different interval.
        with self.assertRaises(ValueError):
            G.penalty_enclosure(p, [0, 1], evidence)
        with self.assertRaises(ValueError):
            G.penalty_enclosure([1, 0, -1, 0, 1], [-1, 0], evidence)

    def test_cache_accepts_equivalent_point_but_rejects_stale_penalty_and_scope(self):
        plain = objective()
        square = objective(penalty=[0, 0, 1])
        cache = {}
        point = G.cached_point(plain, 0, cache, modes=2)
        with patch.object(G.projected, 'full_ground', side_effect=AssertionError('cache miss')):
            self.assertEqual(G.cached_point(square, 0, cache), point)
        G.cached_point(plain, 1, cache, modes=2)
        with self.assertRaises(ValueError):
            G.cached_point(square, 1, cache)
        poisoned = {'0': copy.deepcopy(point)}
        poisoned['0']['spectrum'] = G.projected.full_ground([0], mean_zero=False, modes=2, max_modes=2)
        with self.assertRaises(ValueError):
            G.cached_point(plain, 0, poisoned)

    def test_resume_binds_whole_objective_and_preserves_certified_bounds(self):
        seed = G.branch_and_bound(QUARTIC, max_leaves=2, tolerance=F(1, 10**12), modes=2)
        resumed = G.resume(QUARTIC, seed, max_leaves=4, tolerance=F(1, 10**12), modes=2)
        self.assertGreaterEqual(F(resumed['lower']), F(seed['lower']))
        self.assertLessEqual(F(resumed['upper']), F(seed['upper']))
        for key, value in seed['points'].items():
            self.assertEqual(resumed['points'][key], value)
        changed = copy.deepcopy(QUARTIC)
        changed['penalty'] = [0]
        with self.assertRaises(ValueError):
            G.resume(changed, seed, max_leaves=4)
        self.assertFalse(G.verify(resumed, expected_objective=changed))
        bad = copy.deepcopy(resumed)
        bad['objective_sha256'] = '0'*64
        self.assertFalse(G.verify(bad))

    def test_threshold_equality_gap_and_spectral_existence_scope(self):
        cert = read('open_budget.json')
        lo, hi = F(cert['lower']), F(cert['upper'])
        self.assertLess(lo, hi)
        self.assertEqual(G.threshold_certificate(cert, lo)['status'], 'proved')
        self.assertEqual(G.threshold_certificate(cert, hi)['status'], 'undetermined')
        refuted = G.threshold_certificate(cert, hi+F(1, 1000))
        self.assertEqual(refuted['status'], 'refuted_by_feasible_parameter_spectral_existence')
        self.assertFalse(refuted['explicit_function_witness_included'])
        self.assertFalse(G.verify_threshold(refuted, expected_threshold=hi))
        bad = copy.deepcopy(refuted)
        bad['explicit_function_witness_included'] = True
        self.assertFalse(G.verify_threshold(bad))

    def test_cheap_objective_has_actual_trial_upper_and_search_free_lower(self):
        with (patch.object(G.projected, 'full_ground', side_effect=AssertionError('spectrum search')),
              patch.object(G.E, 'maximize', side_effect=AssertionError('extrema search'))):
            cert = G.cheap_objective_enclosure(SPHERE)
            self.assertTrue(G.verify_cheap_objective(cert, expected_objective=SPHERE))
        self.assertEqual((F(cert['lower']), F(cert['upper'])), (0, F(12, 5)))
        self.assertEqual(cert['trial_function'], 'x1')
        # Independent associated-Legendre integral: q=2t², u=P1^1 cos(phi).
        trial = W.rayleigh_certificate([0, 0, 2], [1], m=1, degrees=[1])
        self.assertEqual(F(trial['rayleigh_quotient']), F(cert['upper']))
        for field, value in [('upper', '1'), ('function_space', 'all_H1_functions'),
                             ('feasible_parameter', '1'), ('trial_function', 'x3')]:
            bad = copy.deepcopy(cert)
            bad[field] = value
            self.assertFalse(G.verify_cheap_objective(bad), field)

    def test_cheap_interval_ranges_on_negative_domain_and_negative_even_terms(self):
        obj = objective(q0=[1, 2, -3], direction=[-2, 0, 1], penalty=[2, -1, -2, 0, 1], domain=[-3, -1])
        cert = G.cheap_objective_enclosure(obj)
        self.assertTrue(G.verify_cheap_objective(cert))
        # At the midpoint s=-2, q=5+2t-5t². Coordinate x3 gives R=2+5-3=4.
        # penalty(-2)=2+2-8+16=12, hence the actual upper is 16.
        self.assertEqual(F(cert['feasible_parameter']), -2)
        self.assertEqual(F(cert['upper']), 16)
        # Term bounds: q >= 3-2-6=-5; p >=2+1-18+1=-14.
        self.assertEqual(F(cert['potential_lower']), -5)
        self.assertEqual(F(cert['penalty_lower']), -14)
        self.assertEqual(F(cert['lower']), -17)

    def test_replay_every_saved_global_and_threshold_without_search(self):
        records = []
        for path in sorted((HERE/'round05').glob('*.json')):
            raw = json.loads(path.read_text())
            if isinstance(raw, dict) and raw.get('format') == G.FORMAT:
                records.append((path.name, raw))
        self.assertGreaterEqual(len(records), 9)
        with (patch.object(G.projected, 'full_ground', side_effect=AssertionError('spectrum search')),
              patch.object(G.E, 'maximize', side_effect=AssertionError('extrema search')),
              patch.object(G, 'branch_and_bound', side_effect=AssertionError('optimization search'))):
            for name, cert in records:
                with self.subTest(name=name):
                    self.assertTrue(G.verify(cert, expected_objective=cert['objective']))
            thresholds = read('stage_144.json')['certificates']
            self.assertEqual(len(thresholds), 3)
            for cert in thresholds:
                self.assertTrue(G.verify_threshold(cert, expected_threshold=cert['threshold']))

    def test_replay_saved_individual_stage_evidence_bound_to_original_inputs(self):
        stages = read('STAGES.json')
        self.assertEqual([s['version'] for s in stages], list(range(135, 145)))
        for stage in stages:
            self.assertTrue(callable(getattr(G, stage['callable'].split('.')[-1])))
            for name in stage['evidence']:
                self.assertTrue((HERE/name).is_file(), name)
        with (patch.object(G.projected, 'full_ground', side_effect=AssertionError('spectrum search')),
              patch.object(G.E, 'maximize', side_effect=AssertionError('extrema search'))):
            self.assertTrue(G.verify_cheap_objective(read('stage_135.json'), expected_objective=SPHERE))
            penalty = read('stage_136.json')
            self.assertEqual(G.penalty_enclosure(QUARTIC['penalty'], QUARTIC['domain'], penalty), penalty)
            saved = read('global_sphere.json')
            points = [saved['points'][str(s)] for s in [-2, 0, 2]]
            lip = read('stage_137.json')
            self.assertEqual(G.lipschitz_cell(SPHERE, [-2, 2], points, lip['penalty_bounds']), lip)
            self.assertEqual(points[1], read('stage_138.json'))
            chord = read('stage_140.json')
            self.assertEqual(G.concavity_cell(SPHERE, [-2, 2], points[0], points[2], chord['negative_maximum']), chord)
            cache_record = read('stage_141.json')
            self.assertEqual(cache_record['cache_before'], cache_record['cache_after'])
            cached = {'0': cache_record['reused_point']}
            self.assertEqual(G.cached_point(SPHERE, 0, cached), points[1])
            location = read('stage_142.json')
            self.assertTrue(G.verify(location['certificate']))
            self.assertEqual(G.minimizer_enclosure(location['certificate']), location['result'])


if __name__ == '__main__':
    unittest.main()
