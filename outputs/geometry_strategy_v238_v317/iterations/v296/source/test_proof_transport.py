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

    def test_interval_intersection_verifies_every_premise(self):
        import json
        sources = [json.loads((pt.backend.FROZEN/'certificates'/name).read_text())
                   for name in ('step_n4.json', 'step_n8.json')]
        cert = pt.spectral_intersection(sources)
        self.assertTrue(pt.verify(cert))
        self.assertEqual(F(cert['lower']), max(F(s['lower']) for s in sources))
        self.assertEqual(F(cert['upper']), min(F(s['upper']) for s in sources))
        cert['sources'][0]['upper'] = '999'
        self.assertFalse(pt.verify(cert))
        for pair in ([pt.constant_spectrum(0), pt.constant_spectrum(1)],
                     [pt.constant_spectrum(0, True), pt.constant_spectrum(0, False)]):
            with self.assertRaises(ValueError):
                pt.spectral_intersection(pair)

    def test_best_L2_evidence_preserves_its_witness(self):
        q = {'kind': 'axis_profile', 'profile': 'abs_power', 'exponent': '-1/4'}
        sources = [pt.backend.pieces.piecewise_model(q, levels=n, degree=2) for n in (1, 4)]
        cert = pt.approximation_select(sources)
        self.assertTrue(pt.verify(cert))
        self.assertEqual(cert['selected_source_index'], 1)
        self.assertEqual(cert['error_upper'], sources[1]['error_upper'])
        self.assertEqual(pt.evaluate(cert, '1/8'), pt.backend.pieces.evaluate(sources[1], '1/8'))
        cert['selected_source_index'] = 0
        self.assertFalse(pt.verify(cert))
        with self.assertRaises(ValueError):
            pt.approximation_select([sources[0], pt.approximation_axis(sources[1], 0)])

    def test_affine_chain_normalization_preserves_order(self):
        source = pt.backend.pieces.piecewise_model({'kind': 'axis_profile', 'profile': 'step',
                                                  'amplitude': '3/2', 'offset': '1/3'})
        cert = pt.approximation_axis(source, 0)
        cert = pt.approximation_shift(cert, '2/3')
        cert = pt.approximation_scale(cert, '-3/2')
        cert = pt.approximation_shift(cert, '7/4')
        cert = pt.approximation_axis(cert, 1)
        flat = pt.normalize_transport(cert)
        self.assertTrue(pt.verify(flat))
        self.assertEqual(flat['function'], cert['function'])
        self.assertEqual(flat['error_upper'], cert['error_upper'])
        self.assertEqual(flat['sources'][0], source)
        self.assertLess(len(pt.backend.canonical(flat)), len(pt.backend.canonical(cert)))
        for t in ('-1', '0', '1/4', '1'):
            self.assertEqual(pt.evaluate(flat, t), pt.evaluate(cert, t))
        self.assertEqual(pt.normalize_transport(flat), flat)
        spec = pt.spectral_shift(pt.spectral_axis(pt.spectral_shift(pt.constant_spectrum('1/2'), '2/3'), 0), '-4/5')
        simple = pt.normalize_transport(spec)
        self.assertEqual(simple['lower'], spec['lower'])
        simple['parameters']['factor'] = '2'
        self.assertFalse(pt.verify(simple))


if __name__ == '__main__':
    unittest.main()
