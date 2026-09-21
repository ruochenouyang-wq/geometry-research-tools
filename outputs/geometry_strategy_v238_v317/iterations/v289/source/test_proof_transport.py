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


if __name__ == '__main__':
    unittest.main()
