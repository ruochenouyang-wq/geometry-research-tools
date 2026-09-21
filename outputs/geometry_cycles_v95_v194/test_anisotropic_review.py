"""Independent Cartesian moments, coupled-tail norms and rotation audit."""
import copy
import json
from pathlib import Path
from fractions import Fraction as F
import unittest
import anisotropic as a

HERE = Path(__file__).resolve().parent
XY = {'1,1,0': '1'}
X = {'1,0,0': '1'}
RX = [[0, 0, -1], [0, 1, 0], [1, 0, 0]]
R345 = [[0, 0, 1], ['4/5', '-3/5', 0], ['3/5', '4/5', 0]]


def certificate_nodes(value):
    if isinstance(value, dict):
        if value.get('format') in (a.FORMAT, a.ROTATION_FORMAT, a.CHEAP_FORMAT):
            yield value
        for child in value.values():
            yield from certificate_nodes(child)
    elif isinstance(value, list):
        for child in value:
            yield from certificate_nodes(child)


class IndependentAnisotropicReview(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cert = a.certify(XY, L=2, bits=20)
        cls.rotated = a.transfer_rotation(X, [0, 1], RX)

    def test_xyz_normalization_and_cross_m_form_closed_formula(self):
        data = a.coupled_form(XY, L=1)
        labels = [(b['m'], b['part']) for b in data['basis']]
        self.assertEqual(labels, [(0, 'real'), (1, 'real'), (1, 'imag')])
        expected = [[F(2, 3), 0, 0], [0, F(2, 3), F(1, 15)],
                    [0, F(1, 15), F(2, 3)]]
        self.assertEqual(data['A'], expected)
        self.assertEqual(data['mass'], [F(1, 3)]*3)
        # x-y yields the exact finite trial 2-(1/15)/(1/3)=9/5.
        self.assertEqual(F(a.rayleigh(XY, {'1,0,0': 1, '0,1,0': -1})['upper']), F(9, 5))

    def test_cross_m_zero_to_two_exact_moment(self):
        data = a.coupled_form(XY, L=2)
        # P2(z)=(3z²-1)/2 and imaginary m=2 harmonic=6xy.
        # Integral P2(z)*xy*(6xy)=3*(3/105-1/15)=-4/35.
        self.assertEqual(data['basis'][3]['m'], 0)
        self.assertEqual(data['basis'][7]['m'], 2)
        self.assertEqual(data['potential_form'][3][7], F(-4, 35))
        self.assertEqual(data['potential_form'][7][3], F(-4, 35))

    def test_cheap_enclosure_is_mathematical_and_bound_to_input(self):
        c = a.cheap_ground_enclosure(XY)
        self.assertEqual((c['lower'], c['upper']), ('1', '9/5'))
        self.assertTrue(a.verify_cheap(c, expected_q=XY))
        self.assertFalse(a.verify_cheap(c, expected_q=X))
        bad = copy.deepcopy(c)
        bad['lower'] = '9/5'
        self.assertFalse(a.verify_cheap(bad))

    def test_exact_complete_tail_norm_xy_times_x(self):
        columns = a.tail_couplings(XY, L=1)
        # xy*x=x^2*y. Its degree-one projection is y/5. Independent
        # normalized moments give ||x^2*y-y/5||²=1/35-1/75=8/525.
        norm = sum((c['entries'][1]**2/c['mass'] for c in columns), F(0))
        self.assertEqual(norm, F(8, 525))
        self.assertEqual(len(columns), 5+7)
        self.assertEqual({c['degree'] for c in columns}, {2, 3})

    def test_projection_removes_constant_before_tail_norm(self):
        columns = a.tail_couplings(X, L=1)
        # q*x=x² -> projected operator removes 1/3. The remaining
        # degree-two norm is 1/5-1/9=4/45, not unprojected 1/5.
        norm = sum((c['entries'][1]**2/c['mass'] for c in columns), F(0))
        self.assertEqual(norm, F(4, 45))
        self.assertEqual(len(columns), 5)
        self.assertTrue(all(c['degree'] == 2 for c in columns))

    def test_schur_correction_independent_mass_formula(self):
        kernel = a.CoupledKernel(XY, L=1)
        x = F(0)
        finite = kernel.matrix(x, 'upper_ritz')
        lower = kernel.matrix(x, 'lower')
        correction = sum((c['entries'][1]**2/(c['mass']*(c['degree']*(c['degree']+1)-1))
                          for c in kernel.columns), F(0))
        self.assertEqual(finite[1][1]-lower[1][1], correction)
        self.assertEqual(correction, F(8, 5775))  # every active tail degree is 3.
        with self.assertRaises(ValueError):
            kernel.matrix(kernel.beta, 'upper_tail')

    def test_far_tail_all_components_vanish_after_degree_band(self):
        q = a.polynomial({'1,0,1': 2, '0,2,0': -1})
        kept = a.harmonic_basis(1, 2)
        far = a.harmonic_basis(5, 5)
        for h in far:
            for b in kept:
                self.assertEqual(a.inner(h['polynomial'], a.multiply(q, b['polynomial'])), 0)

    def test_tail_omission_mass_projection_and_range_mutations_rejected(self):
        for mutation in ('omit_column', 'mass', 'zero_mean', 'range', 'cross_m'):
            bad = copy.deepcopy(self.cert)
            if mutation == 'omit_column':
                bad['complete_tail_couplings'].pop()
            elif mutation == 'mass':
                bad['complete_tail_couplings'][0]['mass'] = '1'
            elif mutation == 'zero_mean':
                bad['mean_zero'] = False
            elif mutation == 'range':
                bad['range_proof']['lower'] = '0'
            else:
                bad['finite_form'][1][2] = '0'
            self.assertFalse(a.verify(bad), mutation)

    def test_sphere_ideal_reduction_identity_and_range(self):
        q = a.add(a.polynomial(X), a.scale(a.add(a.R2, {a.ZERO: -1}), 7))
        reduction = a.sphere_reduce(q)
        self.assertEqual(a.polynomial(reduction['reduced']), a.polynomial(X))
        proof = a.potential_range(q)
        self.assertEqual((proof['lower'], proof['upper']), ('-1', '1'))
        left = a.add(q, a.scale(a.polynomial(reduction['reduced']), -1))
        right = a.multiply(a.add(a.R2, {a.ZERO: -1}), a.polynomial(reduction['quotient']))
        self.assertEqual(left, right)

    def test_rotation_exact_orthogonality_and_oblique_axis(self):
        q = a.rotated_potential([0, 0, 1], R345)
        self.assertEqual(a.encode(q), {'2,0,0': '9/25', '1,1,0': '24/25', '0,2,0': '16/25'})
        rotated = a.transfer_rotation(q, [0, 0, 1], R345)
        self.assertTrue(a.verify_rotation(rotated))
        self.assertTrue(a.verify_rotation(self.rotated, expected_q=X))
        for badrotation in ([[1, 0, 0], [0, 1, 0], [0, 0, -1]],
                            [[0, 0, -1], [0, 1, 0], [2, 0, 0]]):
            with self.assertRaises(ValueError):
                a.transfer_rotation(X, [0, 1], badrotation)

    def test_rotation_mod_sphere_identity_is_safe_current_limitation(self):
        equivalent = a.add(a.polynomial(X), a.add(a.R2, {a.ZERO: -1}))
        self.assertEqual(a.polynomial(a.sphere_reduce(equivalent)['reduced']), a.polynomial(X))
        # Current transfer API requires exact polynomial equality, a stronger
        # sufficient condition than equality modulo r²-1. It safely rejects.
        with self.assertRaises(ValueError):
            a.transfer_rotation(equivalent, [0, 1], RX)
        bad = copy.deepcopy(self.rotated)
        bad['potential'] = a.encode(equivalent)
        self.assertFalse(a.verify_rotation(bad))

    def test_rotation_child_binding_and_radius_zero_trial(self):
        bad = copy.deepcopy(self.rotated)
        bad['source_certificate']['q_coefficients'] = ['0', '2']
        self.assertFalse(a.verify_rotation(bad))
        bad = copy.deepcopy(self.rotated)
        bad['rotation'][0][0] = '1/1000'
        self.assertFalse(a.verify_rotation(bad))
        with self.assertRaises(ValueError):
            a.rayleigh(XY, a.add(a.R2, {a.ZERO: -1}))

    def test_saved_certificate_roundtrips(self):
        found = 0
        for path in sorted((HERE/'round09').rglob('*.json')):
            if path.name == 'INDEPENDENT_TESTS.json':
                continue
            for c in certificate_nodes(json.loads(path.read_text())):
                valid = (a.verify_rotation(c) if c['format'] == a.ROTATION_FORMAT else
                         a.verify_cheap(c) if c['format'] == a.CHEAP_FORMAT else a.verify(c))
                self.assertTrue(valid, str(path)+':'+c['format'])
                found += 1
        self.assertGreaterEqual(found, 4, 'Need actual saved mathematical certificates')


if __name__ == '__main__':
    unittest.main()
