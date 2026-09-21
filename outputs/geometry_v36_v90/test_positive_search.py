"""Mathematical regression, adversarial binding, and proposal/acceptance tests."""
import copy
from fractions import Fraction as F
import unittest
from unittest.mock import patch
import error_bounds as e
import family_v29 as family
import positive_search as p


B2 = [[0, 0, 1]]


class ExactPositiveTests(unittest.TestCase):
    def test_known_local_energy_and_exact_improvement(self):
        cert = p.lower_certificate([0, 0, 6], B2, [F(-4, 5)])
        self.assertTrue(p.verify_lower(cert))
        self.assertEqual(cert['local_energy_coefficients'], ['8/5', '0', '-34/25', '0', '64/25'])
        self.assertEqual(F(cert['spectral_lower']), F(2271, 1600))
        poisson = p.lower_certificate([0, 0, 6], B2, [-1])
        self.assertEqual(poisson['spectral_lower'], '1')
        self.assertGreater(F(cert['spectral_lower'])-F(poisson['spectral_lower']), F(2, 5))

    def test_exact_global_scalar_ansatz_optimum(self):
        q = [0, 0, F(20, 9)]
        lower = p.lower_certificate(q, B2, [F(-1, 3)])
        upper = p.ansatz_upper_certificate(q, B2, [{'t': '1/2', 'weight': '1'}])
        cert = p.optimization_certificate(lower, upper)
        self.assertTrue(p.verify(cert))
        self.assertEqual(cert['ansatz_lower'], '23/36')
        self.assertEqual(cert['ansatz_upper'], '23/36')
        self.assertEqual(cert['ansatz_gap'], '0')
        self.assertEqual(cert['averaged_gradient'], ['0'])
        self.assertEqual(cert['quadratic_stationarity_gap'], '0')
        self.assertEqual(cert['average_contact_and_range_gap'], '0')
        self.assertNotIn('spectral_upper', cert)
        self.assertNotIn('spectral_upper', upper)

    def test_gram_completion_of_square(self):
        upper = p.ansatz_upper_certificate([0, 1, 2], [[0, 1], [0, 0, 1]],
                    [{'t': '-1/2', 'weight': '1/3'}, {'t': '1/3', 'weight': '2/3'}])
        self.assertTrue(p.verify_ansatz_upper(upper))
        gram = [list(map(F, row)) for row in upper['gram']]
        linear = list(map(F, upper['linear']))
        z = list(map(F, upper['quadratic_maximizer']))
        theta = [F(2, 7), F(-3, 4)]
        difference = [x-y for x, y in zip(theta, z)]
        average = F(upper['constant']) + sum(x*y for x, y in zip(linear, theta))
        average -= sum(theta[i]*gram[i][j]*theta[j] for i in range(2) for j in range(2))
        penalty = sum(difference[i]*gram[i][j]*difference[j] for i in range(2) for j in range(2))
        self.assertEqual(F(upper['ansatz_upper'])-average, penalty)

    def test_singular_consistent_probability_gram(self):
        cert = p.ansatz_upper_certificate([0, 0, 6], [[0, 1]],
                        [{'t': '-1', 'weight': '1/2'}, {'t': '1', 'weight': '1/2'}])
        self.assertTrue(p.verify_ansatz_upper(cert))
        self.assertEqual(cert['gram_rank'], 0)
        self.assertEqual(cert['ansatz_upper'], '6')

    def test_singular_unbounded_probability_gram_rejected(self):
        with self.assertRaisesRegex(ValueError, 'unbounded above'):
            p.ansatz_upper_certificate([0, 0, 6], B2, [{'t': '1', 'weight': '1'}])
        with self.assertRaisesRegex(ValueError, 'unbounded above'):
            p.ansatz_upper_certificate([0, 0, 6], [[0, 1], [0, 0, 1]],
                                      [{'t': '0', 'weight': '1'}])

    def test_probability_constraints_exact(self):
        bad_atoms = [[{'t': '0', 'weight': '999999/1000000'}],
                     [{'t': '0', 'weight': '-1'}, {'t': '1/2', 'weight': '2'}],
                     [{'t': '1000001/1000000', 'weight': '1'}],
                     [{'t': '0', 'weight': 1.0}],
                     [{'t': '0', 'weight': '1', 'extra': 'ignored?'}]]
        for atoms in bad_atoms:
            with self.subTest(atoms=atoms), self.assertRaises(ValueError):
                p.ansatz_upper_certificate([0, 1], [[0, 1]], atoms)

    def test_basis_quotient_by_constants_and_dimension_validation(self):
        left = p.lower_certificate([0, 0, 6], [[13, 0, 1]], [F(-4, 5)])
        right = p.lower_certificate([0, 0, 6], B2, [F(-4, 5)])
        self.assertEqual(left, right)
        for basis in ([[1]], [[0, 1], [7, 2]], [[0]*7+[1]]):
            with self.assertRaises(ValueError):
                p.setup([0], basis)
        with self.assertRaises(ValueError):
            p.setup([0], [[0, 1]]*7)
        with self.assertRaises(ValueError):
            p.lower_certificate([0, 0, 6], B2, [-0.8])

    def test_poisson_initialization_reproduces_full_solution(self):
        q = [3, 1, -2, 0, 4, 0, 1]
        basis = p.monomial_basis(q)
        theta = p.poisson_initialization(q, basis)
        self.assertEqual(p.logarithm(basis, theta)[1], family.poisson(q)[0])
        self.assertEqual(p.poisson_initialization([5, 0, 6], B2), [F(-1)])

    def test_even_symmetrization_exact_nonnegative_identity(self):
        q, s = [1, 0, 3, 0, 2], [8, F(2, 3), -1, F(3, 5), F(-1, 9)]
        cert = p.parity_certificate(q, s)
        self.assertTrue(p.verify_parity(cert))
        derivative = list(map(F, cert['odd_s_derivative']))
        gain = list(map(F, cert['gain_coefficients']))
        for t in [F(-1), F(-3, 5), F(0), F(2, 7), F(1)]:
            self.assertEqual(e.evaluate(gain, t), (1-t*t)*e.evaluate(derivative, t)**2)
            self.assertGreaterEqual(e.evaluate(gain, t), 0)
        with self.assertRaises(ValueError):
            p.parity_certificate([0, 1], s)
        cert['gain_coefficients'][0] = '0'
        self.assertFalse(p.verify_parity(cert))

    def test_rayleigh_is_independent_exact_spectral_upper(self):
        constant = p.rayleigh_certificate([0, 0, 6], [1])
        self.assertTrue(p.verify_rayleigh(constant))
        self.assertEqual(constant['spectral_upper'], '2')
        self.assertEqual(p.rayleigh_certificate([7], [3])['spectral_upper'], '7')
        with self.assertRaises(ValueError):
            p.rayleigh_certificate([0], [0])
        trial = p.propose_rayleigh([0, 0, 6])
        self.assertTrue(p.verify_rayleigh(trial))
        self.assertLess(F(trial['spectral_upper']), F(158, 100))
        trial['spectral_upper'] = '0'
        self.assertFalse(p.verify_rayleigh(trial))

    def test_degree_six_local_energy_uses_full_degree_twelve(self):
        cert = p.lower_certificate([0, 0, 1], [[0]*6+[1]], [F(-1, 30)])
        self.assertEqual(len(cert['local_energy_coefficients']), 13)
        self.assertTrue(p.verify_lower(cert))


class SearchAndAdversarialTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.scalar = p.search([0, 0, 6], degree=2)
        cls.sixth = p.search([0, 0, 6], degree=6)

    def test_verified_improvement_and_global_ansatz_gap(self):
        for result in [self.scalar, self.sixth]:
            cert = result['certificate']
            self.assertTrue(p.verify(cert))
            self.assertEqual(cert['status'], 'ansatz_gap_closed')
            self.assertLessEqual(F(cert['ansatz_gap']), F(1, 10000))
            self.assertGreater(F(result['improvement_over_projected_poisson']), F(2, 5))
            self.assertEqual(result['failed_proposals'], [])
            self.assertTrue(all(p.verify_lower(row['certificate']) for row in result['attempts']))
            self.assertEqual(F(cert['ansatz_gap']), F(cert['quadratic_stationarity_gap'])+
                             F(cert['average_contact_and_range_gap']))

    def test_ansatz_upper_provably_is_not_spectral_upper(self):
        # This is a direct certified counterexample to the dangerous relabeling:
        # the sixth-degree spectral lower EXCEEDS the scalar ansatz upper.
        self.assertGreater(F(self.sixth['certificate']['spectral_lower']),
                           F(self.scalar['certificate']['ansatz_upper']))
        self.assertGreater(F(self.scalar['certificate']['spectral_upper']),
                           F(self.sixth['certificate']['spectral_lower']))

    def test_strong_potential_no_false_failure_or_claim(self):
        result = p.search([0, 0, 20], degree=6)
        self.assertTrue(p.verify(result['certificate']))
        self.assertEqual(F(result['poisson_baseline']['spectral_lower']), F(-40, 9))
        self.assertGreater(F(result['improvement_over_projected_poisson']), F(8))
        self.assertGreater(F(result['certificate']['spectral_lower']), F(363, 100))

    def test_general_noneven_potential_search(self):
        result = p.search([0, 1, 3, 0, 2], degree=4)
        self.assertTrue(p.verify(result['certificate']))
        self.assertEqual(len(result['certificate']['lower_certificate']['basis']), 4)
        self.assertGreater(F(result['certificate']['spectral_lower']), F(103, 100))

    def test_exhausted_budget_retains_open_gap(self):
        result = p.search([0, 0, 6], degree=6, max_newton=1, max_exchanges=1,
                          tolerance=F(1, 10**10))
        self.assertTrue(p.verify(result['certificate']))
        self.assertEqual(result['certificate']['status'], 'ansatz_gap_open')
        self.assertGreater(F(result['certificate']['ansatz_gap']), F(1, 10**10))
        self.assertLess(F(result['certificate']['ansatz_upper']), 3)
        self.assertGreater(max(F(c['ansatz_upper']) for c in result['ansatz_upper_attempts']), 10**10)

    def test_failed_proposal_is_retained_and_exact_fallback_used(self):
        with patch.object(p, '_propose', side_effect=ValueError('deliberate numerical failure')):
            result = p.search([0, 0, 6], degree=2)
        self.assertTrue(p.verify(result['certificate']))
        self.assertEqual(result['failed_proposals'], [{'iteration': 0, 'reason': 'deliberate numerical failure'}])
        self.assertEqual(result['certificate']['spectral_lower'], '1')
        self.assertEqual(result['certificate']['status'], 'ansatz_gap_open')

    def test_energy_shift_preserves_all_certified_gaps(self):
        source = self.sixth['certificate']
        translated = p.shift_certificate(source, F(7, 3))
        self.assertTrue(p.verify(translated))
        for key in ['ansatz_lower', 'ansatz_upper', 'spectral_lower', 'spectral_upper']:
            self.assertEqual(F(translated[key])-F(source[key]), F(7, 3))
        self.assertEqual(translated['ansatz_gap'], source['ansatz_gap'])
        self.assertEqual(translated['spectral_gap'], source['spectral_gap'])
        self.assertEqual(translated['lower_certificate']['theta'], source['lower_certificate']['theta'])

    def test_constant_potential(self):
        result = p.search([7], degree=2)
        cert = result['certificate']
        self.assertTrue(p.verify(cert))
        self.assertEqual(cert['spectral_lower'], '7')
        self.assertEqual(cert['spectral_upper'], '7')

    def test_tampering_all_derived_lower_claims_is_rejected(self):
        lower = self.sixth['certificate']['lower_certificate']
        for key in ['spectral_lower', 'candidate_infimum_upper', 'range_gap', 'scope']:
            bad = copy.deepcopy(lower)
            bad[key] = '0'
            self.assertFalse(p.verify_lower(bad), key)
        bad = copy.deepcopy(lower)
        bad['local_energy_coefficients'][0] = '0'
        self.assertFalse(p.verify_lower(bad))
        bad = copy.deepcopy(lower)
        bad['range_proof']['polynomial'][0] = '0'
        self.assertFalse(p.verify_lower(bad))
        bad = copy.deepcopy(lower)
        if bad['range_proof']['intervals'] is not None:
            bad['range_proof']['intervals'].pop()
            self.assertFalse(p.verify_lower(bad))

    def test_tampering_upper_and_combined_claims_is_rejected(self):
        cert = self.sixth['certificate']
        for key in ['ansatz_upper', 'ansatz_gap', 'status', 'spectral_upper', 'spectral_gap',
                    'quadratic_stationarity_gap', 'average_contact_and_range_gap']:
            bad = copy.deepcopy(cert)
            bad[key] = 'wrong'
            self.assertFalse(p.verify(bad), key)
        for key in ['ansatz_upper', 'gram_rank', 'scope', 'limitation']:
            bad = copy.deepcopy(cert['ansatz_upper_certificate'])
            bad[key] = 'wrong'
            self.assertFalse(p.verify_ansatz_upper(bad), key)

    def test_fixed_potential_and_basis_binding(self):
        cert = self.scalar['certificate']
        with self.assertRaises(ValueError):
            p.optimization_certificate(cert['lower_certificate'], self.sixth['certificate']['ansatz_upper_certificate'])
        wrong_q = p.ansatz_upper_certificate([1, 0, 6], B2,
                         cert['ansatz_upper_certificate']['atoms'])
        with self.assertRaises(ValueError):
            p.optimization_certificate(cert['lower_certificate'], wrong_q)
        with self.assertRaises(ValueError):
            p.optimization_certificate(cert['lower_certificate'], cert['ansatz_upper_certificate'],
                                       rayleigh=p.rayleigh_certificate([1], [1]))

    def test_external_range_hook_requires_explicit_trusted_verifier(self):
        class WrappedBackend:
            name = 'test_independent_backend_v1'
            def minimum(self, polynomial, tolerance, max_leaves):
                return {'wrapped': p.DEFAULT_RANGE_BACKEND.minimum(polynomial, tolerance, max_leaves)}
            def verify_minimum(self, polynomial, proof):
                if set(proof) != {'wrapped'}:
                    raise ValueError('Incorrect wrapper')
                return p.DEFAULT_RANGE_BACKEND.verify_minimum(polynomial, proof['wrapped'])
        backend = WrappedBackend()
        lower = p.lower_certificate([0, 0, 6], B2, [F(-4, 5)], range_backend=backend)
        self.assertFalse(p.verify_lower(lower))
        self.assertTrue(p.verify_lower(lower, backend))
        upper = self.scalar['certificate']['ansatz_upper_certificate']
        cert = p.optimization_certificate(lower, upper, range_backend=backend)
        self.assertFalse(p.verify(cert))
        self.assertTrue(p.verify(cert, backend))
        lower['range_proof']['wrapped']['maximum_upper'] = '0'
        self.assertFalse(p.verify_lower(lower, backend))


if __name__ == '__main__':
    unittest.main()
