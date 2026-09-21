"""Independent exact audit cases; no optimizer runs, network, or model calls."""
import copy
from fractions import Fraction as F
import unittest

import algebra as a
import error_bounds as e
import positive_search as p


Q = [0, 0, F(20, 9)]
B = [[0, 0, 1]]
ATOMS = [{'t': '1/2', 'weight': '1'}]


class WrappedSubclass(p.BernsteinBackend):
    """A valid explicit backend with deliberately different proof serialization."""
    name = 'independent_audit_wrapped_bernstein_v1'

    def minimum(self, polynomial, tolerance, max_leaves):
        return {'wrapped': p.DEFAULT_RANGE_BACKEND.minimum(polynomial, tolerance, max_leaves)}

    def verify_minimum(self, polynomial, proof):
        if not isinstance(proof, dict) or set(proof) != {'wrapped'}:
            raise ValueError('Expected independently verified wrapped proof')
        return p.DEFAULT_RANGE_BACKEND.verify_minimum(polynomial, proof['wrapped'])


def combined(backend=None, rayleigh=True):
    lower = p.lower_certificate(Q, B, [F(-1, 3)], range_backend=backend)
    upper = p.ansatz_upper_certificate(Q, B, ATOMS)
    trial = p.rayleigh_certificate(Q, [1]) if rayleigh else None
    return p.optimization_certificate(lower, upper, rayleigh=trial, range_backend=backend)


class BoundMeaningAudit(unittest.TestCase):
    def test_exact_local_energy_identity_and_infimum(self):
        lower = p.lower_certificate(Q, B, [F(-1, 3)])
        expected = [F(2, 3), F(0), F(-2, 9), F(0), F(4, 9)]
        self.assertEqual(list(map(F, lower['local_energy_coefficients'])), expected)
        for t in map(F, ['-1', '-1/2', '0', '1/3', '1/2', '1']):
            self.assertEqual(e.evaluate(expected, t), F(23, 36) + F(4, 9)*(t*t-F(1, 4))**2)
        self.assertEqual(F(lower['spectral_lower']), F(23, 36))

    def test_certified_spectral_lower_can_exceed_smaller_basis_ansatz_upper(self):
        stronger = p.lower_certificate(Q, [[0, 0, 1], [0, 0, 0, 0, 1]],
                                       [F(-1, 3), F(-1, 100)], max_leaves=16)
        scalar = p.ansatz_upper_certificate(Q, B, ATOMS)
        self.assertTrue(p.verify_lower(stronger))
        self.assertTrue(p.verify_ansatz_upper(scalar))
        self.assertGreater(F(stronger['spectral_lower']), F(scalar['ansatz_upper']))
        self.assertNotIn('spectral_upper', scalar)

    def test_closed_ansatz_gap_does_not_close_spectral_gap(self):
        cert = combined()
        self.assertTrue(p.verify(cert))
        self.assertEqual(cert['status'], 'ansatz_gap_closed')
        self.assertEqual(F(cert['ansatz_gap']), 0)
        self.assertEqual(F(cert['spectral_upper']), F(20, 27))
        self.assertEqual(F(cert['spectral_gap']), F(11, 108))

    def test_without_rayleigh_no_spectral_upper_is_emitted(self):
        cert = combined(rayleigh=False)
        self.assertTrue(p.verify(cert))
        for key in ('spectral_upper', 'spectral_gap', 'rayleigh_certificate'):
            self.assertNotIn(key, cert)
        self.assertEqual(cert['limitation'], p.LIMITATION)

    def test_rayleigh_kinetic_potential_and_norm_match_manual_integrals(self):
        cert = p.rayleigh_certificate([1, 1, 2], [1, 1])
        self.assertEqual(F(cert['norm_squared']), F(4, 3))
        self.assertEqual(F(cert['kinetic_numerator']), F(2, 3))
        self.assertEqual(F(cert['potential_numerator']), F(46, 15))
        self.assertEqual(F(cert['spectral_upper']), F(14, 5))
        self.assertTrue(p.verify_rayleigh(cert))

    def test_rayleigh_allows_signed_trial_and_is_scale_invariant(self):
        left = p.rayleigh_certificate([0], [0, 1])
        right = p.rayleigh_certificate([0], [0, -7])
        self.assertEqual(F(left['spectral_upper']), 2)
        self.assertEqual(left['spectral_upper'], right['spectral_upper'])
        self.assertTrue(p.verify_rayleigh(right))

    def test_nonorthogonal_basis_change_preserves_actual_bounds(self):
        q = [1, 1, 2]
        basis = [[0, 1], [0, 0, 1]]
        transformed = [[0, 1, 2], [0, 3, -1]]
        transform = [[F(1), F(2)], [F(3), F(-1)]]
        theta = [F(-1, 5), F(-1, 4)]
        theta_new = a.matvec(a.transpose(a.inverse(transform)), theta)
        left = p.lower_certificate(q, basis, theta)
        right = p.lower_certificate(q, transformed, theta_new)
        self.assertEqual(left['local_energy_coefficients'], right['local_energy_coefficients'])
        self.assertEqual(left['spectral_lower'], right['spectral_lower'])
        atoms = [{'t': '-1/2', 'weight': '1/3'}, {'t': '1/3', 'weight': '2/3'}]
        self.assertEqual(p.ansatz_upper_certificate(q, basis, atoms)['ansatz_upper'],
                         p.ansatz_upper_certificate(q, transformed, atoms)['ansatz_upper'])


class ProbabilityAndTamperAudit(unittest.TestCase):
    def test_atom_splitting_reordering_and_zero_mass_preserve_moments(self):
        original = p.ansatz_upper_certificate(Q, B, ATOMS)
        split = p.ansatz_upper_certificate(Q, B, [
            {'t': '-1', 'weight': '0'}, {'t': '1/2', 'weight': '2/3'},
            {'t': '1/2', 'weight': '1/3'}])
        for key in ('constant', 'linear', 'gram', 'gram_rank', 'ansatz_upper'):
            self.assertEqual(original[key], split[key])
        self.assertTrue(p.verify_ansatz_upper(split))

    def test_rank_one_gram_with_nullspace_is_certified_exactly(self):
        upper = p.ansatz_upper_certificate([0, 0, 1], [[0, 1], [0, 0, 0, 1]],
                                           [{'t': '0', 'weight': '1'}])
        self.assertEqual(upper['gram'], [['1', '0'], ['0', '0']])
        self.assertEqual(upper['gram_rank'], 1)
        self.assertEqual(upper['linear'], ['0', '0'])
        self.assertEqual(upper['ansatz_upper'], '0')
        self.assertTrue(p.verify_ansatz_upper(upper))

    def test_almost_normalized_probability_is_not_accepted(self):
        for delta in (F(1, 10**12), F(-1, 10**12)):
            with self.assertRaises(ValueError):
                p.ansatz_upper_certificate(Q, B, [{'t': '1/2', 'weight': str(1+delta)}])

    def test_probability_gram_and_optimizer_tampering_are_rejected(self):
        source = p.ansatz_upper_certificate(Q, B, ATOMS)
        mutations = [lambda c: c['gram'][0].__setitem__(0, '-1'),
                     lambda c: c['quadratic_maximizer'].__setitem__(0, '0'),
                     lambda c: c['linear'].__setitem__(0, '0'),
                     lambda c: c['atoms'][0].__setitem__('weight', '1/2'),
                     lambda c: c.update(spectral_upper=c['ansatz_upper'])]
        for mutate in mutations:
            bad = copy.deepcopy(source); mutate(bad)
            self.assertFalse(p.verify_ansatz_upper(bad))

    def test_certificates_for_different_potential_or_basis_cannot_be_spliced(self):
        source = combined()
        bad = copy.deepcopy(source)
        bad['rayleigh_certificate'] = p.rayleigh_certificate([0], [1])
        self.assertFalse(p.verify(bad))
        for upper in (p.ansatz_upper_certificate([1, 0, F(20, 9)], B, ATOMS),
                      p.ansatz_upper_certificate(Q, [[0, 0, 2]], ATOMS)):
            with self.assertRaises(ValueError):
                p.optimization_certificate(source['lower_certificate'], upper)

    def test_scope_relabeling_and_ansatz_to_spectral_copy_are_rejected(self):
        source = combined()
        for changes in ({'spectral_upper': source['ansatz_upper']},
                        {'scope': p.SCOPE}, {'limitation': ''}, {'status': 'spectral_exact'}):
            bad = copy.deepcopy(source); bad.update(changes)
            self.assertFalse(p.verify(bad))


class RangeBackendAudit(unittest.TestCase):
    def test_subclass_wrapped_backend_shift_uses_its_explicit_verifier(self):
        backend = WrappedSubclass(); source = combined(backend)
        self.assertTrue(p.verify(source, backend)); self.assertFalse(p.verify(source))
        shifted = p.shift_certificate(source, F(1, 7), backend)
        self.assertTrue(p.verify(shifted, backend)); self.assertFalse(p.verify(shifted))
        self.assertEqual(set(shifted['lower_certificate']['range_proof']), {'wrapped'})
        for key in ('ansatz_lower', 'ansatz_upper', 'spectral_lower', 'spectral_upper'):
            self.assertEqual(F(shifted[key])-F(source[key]), F(1, 7))
        self.assertEqual(shifted['ansatz_gap'], source['ansatz_gap'])
        self.assertEqual(shifted['spectral_gap'], source['spectral_gap'])

    def test_default_shift_preserves_full_partition_and_proof_gap(self):
        lower = p.lower_certificate(Q, [[0, 0, 1], [0, 0, 0, 0, 1]],
                                    [F(-1, 3), F(-1, 100)], max_leaves=4)
        basis = lower['basis']
        atoms = [{'t': '-1/2', 'weight': '1/2'}, {'t': '1/3', 'weight': '1/2'}]
        source = p.optimization_certificate(lower, p.ansatz_upper_certificate(Q, basis, atoms))
        shifted = p.shift_certificate(source, F(-1, 7))
        self.assertEqual(shifted['lower_certificate']['range_proof']['intervals'],
                         source['lower_certificate']['range_proof']['intervals'])
        self.assertEqual(shifted['lower_certificate']['range_gap'], lower['range_gap'])
        self.assertEqual(shifted['ansatz_gap'], source['ansatz_gap'])
        self.assertTrue(p.verify(shifted))

    def test_backend_name_alone_does_not_confer_proof_trust(self):
        cert = combined(WrappedSubclass())
        cert['lower_certificate']['range_backend'] = p.DEFAULT_RANGE_BACKEND.name
        self.assertFalse(p.verify(cert))

    def test_wrapped_range_proof_for_another_local_energy_is_rejected(self):
        backend = WrappedSubclass(); source = combined(backend)
        other = p.lower_certificate(Q, B, [F(0)], range_backend=backend)
        source['lower_certificate']['range_proof'] = other['range_proof']
        self.assertFalse(p.verify(source, backend))

    def test_backend_bounds_must_be_ordered_exact_fractions(self):
        for bounds in ((0, F(1)), (F(1), F(0)), (0.0, 1.0)):
            class BadReturn(WrappedSubclass):
                def verify_minimum(self, polynomial, proof):
                    return bounds
            with self.subTest(bounds=bounds), self.assertRaises(ValueError):
                p.lower_certificate(Q, B, [F(-1, 3)], range_backend=BadReturn())

    def test_default_range_partition_and_infimum_claims_are_replayed(self):
        lower = p.lower_certificate(Q, [[0, 0, 1], [0, 0, 0, 0, 1]],
                                    [F(-1, 3), F(-1, 100)], max_leaves=4)
        bad = copy.deepcopy(lower); bad['range_proof']['intervals'].pop()
        self.assertFalse(p.verify_lower(bad))
        bad = copy.deepcopy(lower); bad['candidate_infimum_upper'] = '1000'
        self.assertFalse(p.verify_lower(bad))


if __name__ == '__main__':
    unittest.main()
