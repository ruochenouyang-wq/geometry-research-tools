import copy
from pathlib import Path
import sys
import unittest

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import spectral_gap as g


def function(amplitude='-2/5', exponent='-4/9'):
    return {'kind': 'axis_profile', 'axis': 2, 'profile': 'abs_power',
            'amplitude': amplitude, 'offset': '0', 'exponent': exponent}


class SpectralGapTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.q = function()
        cls.cert = g.certified_gap(cls.q, modes=8, near_tail=8, bits=16)

    def test_h07_gap_replays_and_separates(self):
        self.assertTrue(g.verify_gap(self.cert, self.q, 1, True))
        self.assertGreater(g.F(self.cert['beta']), g.F(19, 10))
        self.assertLess(g.F(self.cert['beta']), g.F(self.cert['radial_tail_lower']))

    def test_constant_spectrum_exact_checks_all_spaces(self):
        q = dict(self.q, amplitude='0', offset='-1/3')
        for m, mean, exact in ((0, False, g.F(5, 3)), (0, True, g.F(17, 3)), (1, True, g.F(17, 3))):
            with self.subTest(m=m, mean=mean):
                c = g.certified_gap(q, m, mean, modes=3, near_tail=2, bits=20)
                self.assertTrue(g.verify_gap(c, q, m, mean))
                self.assertLessEqual(g.F(c['beta']), exact)
                self.assertLess(exact-g.F(c['beta']), g.F(1, 10000))

    def test_binding_rejects_wrong_function_space_and_false_bools(self):
        self.assertFalse(g.verify_gap(self.cert, function('-1/3'), 1, True))
        self.assertFalse(g.verify_gap(self.cert, self.q, 0, True))
        self.assertFalse(g.verify_gap(self.cert, self.q, True, True))
        self.assertFalse(g.verify_gap(self.cert, self.q, 1, 1))
        self.assertFalse(g.verify_gap(self.cert, self.q, 1, False))

    def test_all_proof_components_reject_tampering(self):
        edits = (
            lambda c: c.update(beta='6'),
            lambda c: c.update(lower='6', beta='6'),
            lambda c: c.update(eigenvalue_index=True),
            lambda c: c.update(full_radial_space_covered=1),
            lambda c: c['comparison_inertia'].__setitem__(0, 0),
            lambda c: c['comparison_matrix'][0].__setitem__(0, '0'),
            lambda c: c['kernel_evidence']['C'][0].__setitem__(0, '0'),
            lambda c: c['kernel_evidence']['near_tail'].pop(),
            lambda c: c['form_bound'].update(energy_factor='1'),
            lambda c: c.update(remainder_tail_lower='1000'),
        )
        for edit in edits:
            c = copy.deepcopy(self.cert)
            edit(c)
            with self.subTest(edit=edit):
                self.assertFalse(g.verify_gap(c, self.q, 1, True))

    def test_custom_form_requires_replay_and_function_binding(self):
        form = {'function': self.q, 'form': g.direct.form_bound(self.q, 40)}
        def check(c, expected_function=None):
            return c == {'function': g.direct.normalize(expected_function),
                         'form': g.direct.form_bound(expected_function, 40)}
        c = g.certified_gap(self.q, modes=4, bits=8, form_certificate=form, form_verifier=check)
        self.assertTrue(g.verify_gap(c, self.q, 1, True, form_verifier=check))
        self.assertFalse(g.verify_gap(c, self.q, 1, True))
        c['form_certificate']['form']['mass_offset'] = '100'
        self.assertFalse(g.verify_gap(c, self.q, 1, True, form_verifier=check))
        with self.assertRaises(ValueError):
            g.certified_gap(self.q, form_certificate=form, form_verifier=lambda *a, **kw: 1)

    def test_strong_residual_domain_and_native_types(self):
        for exponent in ('-1/2', '-2/3'):
            with self.assertRaises(ValueError):
                g.certified_gap(function(exponent=exponent))
        c = copy.deepcopy(self.cert)
        c['comparison_inertia'] = tuple(c['comparison_inertia'])
        self.assertFalse(g.verify_gap(c))
        c = copy.deepcopy(self.cert)
        c['cycle'] = c
        self.assertFalse(g.verify_gap(c))


if __name__ == '__main__':
    unittest.main()
