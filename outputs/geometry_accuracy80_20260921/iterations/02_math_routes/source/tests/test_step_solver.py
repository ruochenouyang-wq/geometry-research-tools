"""Public development tests; no hidden or holdout inputs are read."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import json
import sys
import unittest

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import step_solver as s

Q = {'kind': 'axis_profile', 'axis': 2, 'profile': 'step', 'amplitude': '-3/5', 'offset': '0'}


def integral(p, side):
    return sum((a*F(1 if side == 1 or k % 2 == 0 else -1, k+1) for k, a in p.items()), F(0))


def mul(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            out[i+j] = out.get(i+j, F(0))+a*b
    return out


def energy_entry(a, b, side, m, q):
    """Independent gradient-form calculation, not H applied to the trial."""
    da = {k-1: k*v for k, v in a.items() if k}
    db = {k-1: k*v for k, v in b.items() if k}
    w = {0: F(1), 2: F(-1)}
    mass_weight = w if m else {0: F(1)}
    gradient_weight = mul(mass_weight, w)
    return integral(mul(mul(da, db), gradient_weight), side)+(m*(m+1)+q)*integral(mul(mul(a, b), mass_weight), side)


class StepTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.solution = s.solve(Q, degree=8, tolerance='1/10000000000')
        cls.cert = cls.solution['certificate']
        cls.source = cls.cert['source']

    def test_h14_meets_1e10_and_replays(self):
        self.assertEqual(self.cert['status'], 'target_met')
        self.assertLessEqual(F(self.cert['exact_width']), F(1, 10**10))
        self.assertTrue(s.verify(self.cert, Q, self.source, '1/10000000000', True))

    def test_h14_cold_budget_report(self):
        # Wall-clock is reported, not a brittle timing assertion in the tests.
        self.assertGreater(self.solution['elapsed_seconds_including_verify'], 0)
        self.assertEqual(len(self.cert['trial']['basis']), 16)

    def test_full_space_m0_meets_target(self):
        r = s.solve(Q, mean_zero=False, tolerance='1/10000000000')
        self.assertEqual(r['certificate']['trial']['azimuth_m'], 0)
        self.assertEqual(r['status'], 'target_met')
        self.assertTrue(s.verify(r['certificate'], expected_mean_zero=False))

    def test_task_convenience_entry(self):
        r = s.solve({'function': Q, 'mean_zero': False, 'tolerance': '1/100000000'})
        self.assertTrue(s.verify(r['certificate'], Q, expected_mean_zero=False))

    def test_public_changed_amplitudes(self):
        for amp in ('-1/3', '-1', '3/5'):
            for mz in (False, True):
                r = s.solve(dict(Q, amplitude=amp), degree=10, mean_zero=mz, tolerance='1/10000000000')
                self.assertEqual(r['status'], 'target_met')
                self.assertTrue(s.verify(r['certificate'], dict(Q, amplitude=amp), expected_mean_zero=mz))

    def test_zero_amplitude_exact_constant(self):
        for mz, eigenvalue in ((True, F(7, 3)), (False, F(1, 3))):
            r = s.solve(dict(Q, amplitude='0', offset='1/3'), mean_zero=mz)
            self.assertEqual(F(r['certificate']['lower']), eigenvalue)
            self.assertEqual(F(r['certificate']['upper']), eigenvalue)

    def test_constant_offset_translation(self):
        q = dict(Q, offset='7/5')
        r = s.solve(q, degree=8, tolerance='1/10000000000')
        self.assertEqual(r['proposal']['coefficients'], self.solution['proposal']['coefficients'])
        self.assertEqual(F(r['certificate']['statistics']['rayleigh'])-F(self.cert['statistics']['rayleigh']), F(7, 5))
        self.assertEqual(F(r['certificate']['exact_width']), F(self.cert['exact_width']))

    def test_axis_rotation_binding(self):
        r = s.solve(dict(Q, axis=0))
        self.assertTrue(s.verify(r['certificate'], dict(Q, axis=0)))
        self.assertFalse(s.verify(r['certificate'], Q))

    def test_constant_basis_exact_moments(self):
        a, b = F(-3, 5), F(2, 7)
        for m in (0, 1):
            M, H, R = s.trial_matrices(dict(Q, offset=str(b)), [(0, 0), (0, 1)], m)
            half = F(1) if m == 0 else F(2, 3)
            h = m*(m+1)+b
            self.assertEqual(M[0][0], 2*half)
            self.assertEqual(H[0][0], half*(2*h+a))
            self.assertEqual(R[0][0], half*(h*h+(h+a)**2))

    def test_matrices_independent_gradient_form(self):
        terms = s.basis(4)
        for m in (0, 1):
            M, H, R = s.trial_matrices(Q, terms, m)
            for i, ti in enumerate(terms):
                for j, tj in enumerate(terms):
                    expected = F(0)
                    for side in (-1, 1):
                        a = {ti[1]: F(1)} if ti[0] in (0, side) else {}
                        b = {tj[1]: F(1)} if tj[0] in (0, side) else {}
                        expected += energy_entry(a, b, side, m, F(-3, 5) if side == 1 else F(0))
                    self.assertEqual(H[i][j], expected)
                    self.assertEqual(M[i][j], M[j][i])
                    self.assertEqual(R[i][j], R[j][i])

    def test_strong_residual_independent_merged_polynomials(self):
        terms = self.cert['trial']['basis']; v = self.cert['trial']['coefficients']
        self.assertEqual(s.trial_statistics(Q, terms, v), s._direct_statistics(Q, s._basis(terms), v, 1))
        self.assertGreater(F(self.cert['statistics']['residual_squared']), 0)

    def test_piecewise_one_sided_support_orthogonal(self):
        M, H, R = s.trial_matrices(Q, [(-1, 2), (1, 2)])
        self.assertEqual((M[0][1], H[0][1], R[0][1]), (0, 0, 0))

    def test_non_c1_basis_rejected(self):
        for terms in ([(0, 0), (1, 0)], [(0, 0), (-1, 1)]):
            with self.assertRaises(ValueError): s.trial_matrices(Q, terms)

    def test_zero_and_dependent_trials_rejected(self):
        with self.assertRaises(ValueError): s.trial_statistics(Q, s.basis(2), [0]*4)
        with self.assertRaises(ValueError): s.trial_matrices(Q, [(0, 2), (-1, 2), (1, 2)])

    def test_complete_gap_plugin_and_wrong_mode_rejected(self):
        import spectral_gap
        proof = spectral_gap.certified_gap(Q, m=1, mean_zero=True, modes=4, near_tail=4, bits=12)
        c = s.trial_certificate(Q, self.cert['trial']['basis'], self.cert['trial']['coefficients'], self.source, proof)
        self.assertTrue(s.verify(c, Q))
        proof['azimuth_m'] = 0
        with self.assertRaises(ValueError): s.trial_certificate(Q, self.cert['trial']['basis'], self.cert['trial']['coefficients'], self.source, proof)

    def test_strict_temple_gap_failure(self):
        # r=z has energy 6-3/10, above the analytic beta=6-3/5.
        with self.assertRaisesRegex(ValueError, 'strict'):
            s.trial_certificate(Q, [(0, 0), (0, 1)], [0, 1], self.source)

    def test_underresolved_returns_honest_open(self):
        r = s.solve(Q, source=self.source, degree=2, tolerance='1/10000000000')
        self.assertEqual(r['status'], 'certified_open')
        self.assertTrue(s.verify(r['certificate']))

    def test_tamper_trial_and_full_residual(self):
        for field in ('rayleigh', 'residual_squared', 'operator_norm_squared'):
            c = deepcopy(self.cert); c['statistics'][field] = '0'
            self.assertFalse(s.verify(c))
        c = deepcopy(self.cert); c['trial']['coefficients'][0] = '2'
        self.assertFalse(s.verify(c))

    def test_tamper_source_angular_or_missing_m0(self):
        for field in ('angular_tail_lower', 'lower'):
            c = deepcopy(self.cert); c['source'][field] = '1000'
            self.assertFalse(s.verify(c))
        c = deepcopy(self.cert); c['source']['sectors'].pop(0)
        self.assertFalse(s.verify(c))

    def test_tamper_scope_target_and_status(self):
        for field, value in (('mean_zero', False), ('scope', 'finite_trial_plane'), ('status', 'certified_open'), ('full_infinite_space_covered', False)):
            c = deepcopy(self.cert); c[field] = value
            self.assertFalse(s.verify(c))
        self.assertFalse(s.verify(self.cert, expected_tolerance='1/100'))
        self.assertFalse(s.verify(self.cert, expected_function=dict(Q, amplitude='-2/3')))

    def test_pointwise_gap_tamper_rejected(self):
        for field, value in (('lower', '100'), ('eigenvalue_index', 1), ('mean_zero', False)):
            c = deepcopy(self.cert); c['gap'][field] = value
            self.assertFalse(s.verify(c))

    def test_json_roundtrip_and_bad_inputs(self):
        self.assertTrue(s.verify(json.loads(json.dumps(self.cert))))
        for x in (None, [], {}, 1, 'certificate'): self.assertFalse(s.verify(x))
        with self.assertRaises(ValueError): s.solve(dict(Q, profile='abs_power', exponent='-1/4'))
        with self.assertRaises(ValueError): s.solve(Q, mean_zero=1)
        with self.assertRaises(ValueError): s.basis(True)
        with self.assertRaises(ValueError): s.solve({'kind': 'approximation', 'function': Q})


if __name__ == '__main__':
    unittest.main()
