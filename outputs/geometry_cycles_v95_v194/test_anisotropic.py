import copy
import json
import unittest
from fractions import Fraction as F
import anisotropic as a


XY = {'1,1,0': '1'}
X = {'1,0,0': '1'}
RADIUS = {'2,0,0': '1', '0,2,0': '1', '0,0,2': '1'}
RX = [[0, 0, -1], [0, 1, 0], [1, 0, 0]]
R345 = [[0, 0, 1], ['4/5', '-3/5', 0], ['3/5', '4/5', 0]]


class AnisotropicTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.xy2 = a.certify(XY, 2, 24)
        cls.xy3 = a.certify(XY, 3, 24)
        cls.rotated_x = a.transfer_rotation(X, [0, 1], RX)

    def test_v175_sparse_algebra(self):
        p = a.polynomial({'1,0,0': '1/2', '0,1,0': '-1/3'})
        self.assertEqual(a.encode(a.multiply(p, p)),
                         {'0,2,0': '1/9', '1,1,0': '-1/3', '2,0,0': '1/4'})
        self.assertEqual(a.derivative(p, 0), {a.ZERO: F(1, 2)})

    def test_v175_cheap_bound_has_actual_mathematical_output(self):
        c = a.cheap_ground_enclosure(XY)
        self.assertTrue(a.verify_cheap(c, expected_q=XY))
        self.assertEqual((c['lower'], c['upper']), ('1', '9/5'))
        self.assertEqual(len(c['trials']), 9)

    def test_v175_cheap_input_scope_and_trial_bound(self):
        c = a.cheap_ground_enclosure(XY)
        self.assertFalse(a.verify_cheap(c, expected_q=X))
        bad = copy.deepcopy(c)
        bad['scope'] = 'all_H1_functions'
        self.assertFalse(a.verify_cheap(bad))
        bad = copy.deepcopy(c)
        bad['trials'][bad['selected_trial']]['upper'] = '1'
        self.assertFalse(a.verify_cheap(bad))

    def test_v175_reject_ambiguous_and_out_of_scope_inputs(self):
        for q in ({'1,1,0': 0.5}, {'1,1,0': True}, {'-1,0,0': 1},
                  {'5,0,0': 1}, {'1,0': 1}, {(True, 0, 0): 1}):
            with self.assertRaises(ValueError):
                a.polynomial(q)

    def test_v176_known_moments(self):
        cases = {(0, 0, 0): F(1), (2, 0, 0): F(1, 3), (4, 0, 0): F(1, 5),
                 (2, 2, 0): F(1, 15), (2, 2, 2): F(1, 105), (1, 2, 0): F(0)}
        for e, value in cases.items():
            self.assertEqual(a.sphere_moment(e), value)
        self.assertEqual(a.integrate(a.power(a.R2, 4)), 1)

    def test_v177_harmonics_homogeneous_harmonic_and_mean_zero(self):
        for l in range(1, 9):
            basis = a.harmonic_basis(l, l)
            self.assertEqual(len(basis), 2*l+1)
            for b in basis:
                p = b['polynomial']
                self.assertTrue(all(sum(e) == l for e in p))
                self.assertEqual(a.integrate(p), 0)
                self.assertEqual(a.add(*(a.derivative(a.derivative(p, i), i) for i in range(3))), {})

    def test_v178_complete_gram_and_independent_energy(self):
        data = a.gram_energy(4)
        self.assertEqual(len(data['basis']), 24)
        self.assertEqual(data['mass'][:3], [F(1, 3)]*3)
        for i, b in enumerate(data['basis']):
            self.assertEqual(data['energy'][i][i]/data['mass'][i], b['degree']*(b['degree']+1))

    def test_v179_cross_azimuth_couplings_are_present(self):
        data = a.coupled_form(XY, 2)
        self.assertTrue(any(data['potential_form'][i][j] and p['m'] != q['m']
                            for i, p in enumerate(data['basis']) for j, q in enumerate(data['basis'])))
        self.assertEqual(a.rayleigh(XY, {'1,0,0': 1, '0,1,0': -1})['upper'], '9/5')

    def test_v180_every_tail_harmonic_and_mass(self):
        columns = a.tail_couplings(XY, 3)
        self.assertEqual(len(columns), 9+11)
        self.assertEqual({c['degree'] for c in columns}, {4, 5})
        self.assertTrue(all(c['mass'] > 0 and len(c['entries']) == 15 for c in columns))
        self.assertTrue(any(c['mass'] != 1 for c in columns))

    def test_v180_beyond_band_vanishes(self):
        kept = a.harmonic_basis(1, 2)
        beyond = a.harmonic_basis(5, 5)
        q = a.polynomial(XY)
        self.assertTrue(all(a.inner(p['polynomial'], a.multiply(q, b['polynomial'])) == 0
                            for p in beyond for b in kept))

    def test_v181_schur_uses_mass_and_has_consistent_sign(self):
        k = a.CoupledKernel(XY, 2)
        lo, up, ritz = (k.matrix(F(9, 5), kind) for kind in ('lower', 'upper_tail', 'upper_ritz'))
        self.assertTrue(all(lo[i][i] <= up[i][i] <= ritz[i][i] for i in range(len(lo))))
        with self.assertRaises(ValueError):
            k.matrix(k.beta, 'lower')

    def test_v182_zero_and_constant_potentials_exact(self):
        for q, value in (({}, '2'), ({'0,0,0': '7/3'}, '13/3')):
            c = a.certify(q, 1)
            self.assertEqual(c['lower'], value)
            self.assertEqual(c['upper'], value)
            self.assertTrue(a.verify(c, expected_q=q))

    def test_v182_true_anisotropy_and_improved_truncation(self):
        self.assertTrue(a.verify(self.xy3, expected_q=XY))
        self.assertLess(F(self.xy3['upper']), F(9, 5))
        self.assertLess(F(self.xy3['exact_width']), F(self.xy2['exact_width'])/100)
        self.assertLessEqual(F(self.xy2['lower']), F(self.xy3['upper']))
        self.assertLessEqual(F(self.xy3['lower']), F(self.xy2['upper']))

    def test_v182_json_roundtrip(self):
        self.assertTrue(a.verify(json.loads(json.dumps(self.xy3))))

    def test_v182_wrong_input_and_scope_rejected(self):
        self.assertFalse(a.verify(self.xy3, expected_q=X))
        for field, value in [('scope', 'all_H1_functions'), ('mean_zero', False),
                             ('mean_zero', 1), ('geometry', 'general_surface'),
                             ('eigenvalue_index_in_constrained_space', 2)]:
            bad = copy.deepcopy(self.xy3)
            bad[field] = value
            self.assertFalse(a.verify(bad))

    def test_v182_matrix_mass_and_tail_mutations_rejected(self):
        bad = copy.deepcopy(self.xy3)
        bad['mass'][0] = '1'
        self.assertFalse(a.verify(bad))
        bad = copy.deepcopy(self.xy3)
        bad['finite_form'][0][1] = '999'
        self.assertFalse(a.verify(bad))
        bad = copy.deepcopy(self.xy3)
        bad['complete_tail_couplings'].pop()
        self.assertFalse(a.verify(bad))
        bad = copy.deepcopy(self.xy3)
        bad['complete_tail_couplings'][0]['mass'] = '1'
        self.assertFalse(a.verify(bad))
        bad = copy.deepcopy(self.xy3)
        bad['complete_tail_couplings'][0]['entries'][0] = '1'
        self.assertFalse(a.verify(bad))

    def test_v182_false_endpoint_rejected(self):
        bad = copy.deepcopy(self.xy3)
        bad['lower'] = '2'
        self.assertFalse(a.verify(bad))

    def test_v183_sphere_identity_gives_exact_constant(self):
        reduction = a.sphere_reduce(RADIUS)
        self.assertEqual(reduction['reduced'], {'0,0,0': '1'})
        self.assertEqual(reduction['quotient'], {'0,0,0': '1'})
        bound = a.potential_range(RADIUS)
        self.assertEqual((bound['lower'], bound['upper']), ('1', '1'))
        c = a.certify(RADIUS, 2)
        self.assertEqual((c['lower'], c['upper']), ('3', '3'))
        self.assertTrue(a.verify(c))

    def test_v183_reduction_never_worsens_selected_coefficient_range(self):
        for q in (XY, {'0,0,2': 1}, {'0,0,4': 1}, RADIUS):
            old, new = a.potential_range(q, False), a.potential_range(q)
            self.assertLessEqual(F(new['upper'])-F(new['lower']), F(old['upper'])-F(old['lower']))

    def test_v183_corrupt_identity_rejected(self):
        c = a.certify(RADIUS, 1)
        c['range_proof']['sphere_identity']['quotient']['0,0,0'] = '2'
        self.assertFalse(a.verify(c))

    def test_v184_rotation_agrees_with_coupled_solver(self):
        coupled = a.certify(X, 3, 28)
        rotated = self.rotated_x
        self.assertTrue(a.verify_rotation(rotated, expected_q=X))
        self.assertLessEqual(F(coupled['lower']), F(rotated['upper']))
        self.assertLessEqual(F(rotated['lower']), F(coupled['upper']))

    def test_v184_oblique_rational_axis_square(self):
        q = {'2,0,0': '9/25', '1,1,0': '24/25', '0,2,0': '16/25'}
        c = a.transfer_rotation(q, [0, 0, 1], R345)
        self.assertTrue(a.verify_rotation(c, expected_q=q))

    def test_v184_wrong_rotation_and_child_rejected(self):
        bad = copy.deepcopy(self.rotated_x)
        bad['rotation'][2][0] = '2'
        self.assertFalse(a.verify_rotation(bad))
        bad = copy.deepcopy(self.rotated_x)
        bad['source_certificate']['mean_zero'] = False
        self.assertFalse(a.verify_rotation(bad))
        with self.assertRaises(ValueError):
            a.transfer_rotation(XY, [0, 1], RX)

    def test_trial_rejects_nonzero_mean_and_sphere_zero(self):
        for u in ({'0,0,0': 1}, {'0,0,0': -1, **RADIUS}):
            with self.assertRaises(ValueError):
                a.rayleigh(XY, u)


if __name__ == '__main__':
    unittest.main()
