"""Integration attacks on goal status and mathematical conclusion binding."""
from copy import deepcopy
import unittest
import backend
import verification


class IntegrationBoundaryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.q={'kind':'axis_profile','profile':'step'}
        cls.spectrum=backend.direct.full_ground(cls.q,modes=2,bits=12,
                                               tolerance='1/100000000')
        cls.approximation=backend.pieces.piecewise_model(cls.q)

    def test_valid_open_is_not_target_met(self):
        result=verification.assess(self.spectrum,{'kind':'spectrum','function':self.q})
        self.assertTrue(result['certificate_valid'])
        self.assertFalse(result['target_met'])
        self.assertEqual(result['status'],'certified_open')

    def test_exact_approximation_does_not_solve_spectrum(self):
        result=verification.assess(self.approximation,{'kind':'approximation','function':self.q})
        self.assertTrue(result['target_met'])
        self.assertFalse(verification.verify_any(self.approximation,expected_kind='spectrum'))

    def test_claim_edit_is_not_a_mathematical_proof(self):
        bad=deepcopy(self.spectrum)
        bad.update(lower=bad['upper'],exact_width='0',status='target_met')
        self.assertFalse(verification.verify_any(bad))

    def test_same_type_different_question_rejected(self):
        self.assertFalse(verification.verify_any(self.spectrum,
            expected_function=dict(self.q,amplitude='-1')))
        self.assertFalse(verification.verify_any(self.spectrum,expected_mean_zero=False))
        self.assertFalse(verification.verify_any(self.spectrum,expected_mean_zero=1))


if __name__=='__main__':unittest.main()
