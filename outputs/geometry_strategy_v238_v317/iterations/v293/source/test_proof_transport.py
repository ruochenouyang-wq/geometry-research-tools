import unittest
from copy import deepcopy
from fractions import Fraction as F
import proof_transport as pt


class TransportTests(unittest.TestCase):
    def test_constant_spectrum_both_spaces(self):
        for mean_zero, expected in [(False, '-7/3'), (True, '-1/3')]:
            cert = pt.constant_spectrum('-7/3', mean_zero)
            self.assertTrue(pt.verify(cert))
            self.assertEqual(F(cert['lower']), F(expected))
            self.assertEqual(cert['lower'], cert['upper'])

    def test_constant_rejects_false_claims(self):
        cert = pt.constant_spectrum(1)
        cert['lower'] = '4'
        self.assertFalse(pt.verify(cert))
        cert = pt.constant_spectrum(1)
        cert['parameters']['function']['amplitude'] = '1'
        self.assertFalse(pt.verify(cert))
        with self.assertRaises(ValueError):
            pt.constant_spectrum(1, 1)

    def test_spectral_axis_reuses_real_source(self):
        from pathlib import Path
        import json
        source = json.loads((Path(pt.backend.FROZEN)/'certificates/singular_n4.json').read_text())
        cert = pt.spectral_axis(source, 0)
        self.assertTrue(pt.verify(cert))
        self.assertEqual(cert['function']['axis'], 0)
        self.assertEqual(cert['lower'], source['lower'])
        cert['sources'][0]['lower'] = '99'
        self.assertFalse(pt.verify(cert))
        with self.assertRaises(ValueError):
            pt.spectral_axis(source, 3)

    def test_spectral_shift_preserves_width_and_space(self):
        for mean_zero in (False, True):
            source = pt.constant_spectrum('-1/3', mean_zero)
            cert = pt.spectral_shift(source, '-7/5')
            self.assertTrue(pt.verify(cert))
            self.assertEqual(F(cert['upper']), F(source['upper'])-F(7, 5))
            self.assertEqual(cert['exact_width'], source['exact_width'])
            self.assertIs(cert['mean_zero'], mean_zero)
            cert['mean_zero'] = not mean_zero
            self.assertFalse(pt.verify(cert))

    def test_approximation_rotation_keeps_L2_scope(self):
        source = pt.backend.pieces.piecewise_model({'kind': 'axis_profile', 'profile': 'abs_power',
              'exponent': '-1/4'}, levels=2, degree=1)
        cert = pt.approximation_axis(source, 1)
        self.assertTrue(pt.verify(cert))
        self.assertEqual(cert['error_upper'], source['error_upper'])
        self.assertFalse(cert['spectral_transfer_claimed'])
        with self.assertRaises(ValueError):
            pt.spectral_axis(source, 1)
        with self.assertRaises(ValueError):
            pt.approximation_axis(pt.constant_spectrum(), 0)

    def test_approximation_shift_transforms_actual_representative(self):
        source = pt.backend.pieces.piecewise_model({'kind': 'axis_profile', 'profile': 'step',
              'amplitude': '3/2', 'offset': '-1/2'})
        cert = pt.approximation_shift(source, '7/3')
        self.assertTrue(pt.verify(cert))
        self.assertEqual(cert['error_upper'], '0')
        for t in ['-1', '0', '1/4', '1']:
            self.assertEqual(pt.evaluate(cert, t), pt.backend.pieces.evaluate(source, t)+F(7, 3))
        with self.assertRaises(ValueError):
            pt.evaluate(cert, '2')

    def test_signed_and_zero_approximation_scaling(self):
        source = pt.backend.pieces.piecewise_model({'kind': 'axis_profile', 'profile': 'abs_power',
              'exponent': '-1/4', 'offset': '2'}, levels=2, degree=1)
        for factor in ('-3/2', '0', '2/3'):
            cert = pt.approximation_scale(source, factor)
            self.assertTrue(pt.verify(cert))
            self.assertEqual(F(cert['error_upper']), abs(F(factor))*F(source['error_upper']))
            self.assertEqual(pt.evaluate(cert, '1/2'), F(factor)*pt.backend.pieces.evaluate(source, '1/2'))
        with self.assertRaises(ValueError):
            pt.approximation_scale(pt.constant_spectrum(), 2)


if __name__ == '__main__':
    unittest.main()
