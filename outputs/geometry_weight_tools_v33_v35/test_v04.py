from fractions import Fraction as F
from unittest.mock import patch
import unittest
import v04_search as v

class GuidedTests(unittest.TestCase):
    def test_same_math_enclosure_as_exact_search(self):
        for q,k,n in [([0,1],1,8),([1,-2,3],2,8),([0,'1/1000000',-100],2,20)]:
            a=v.v3.certify(q,k,n,48)
            b=v.guided(q,k,n)['certificate']
            self.assertLessEqual(max(F(a['lower']),F(b['lower'])),min(F(a['upper']),F(b['upper'])))
            self.assertTrue(v.v3.verify(b,independent=True))

    def test_small_number_of_exact_counts(self):
        r=v.guided([0,1],1,8)
        self.assertTrue(r['target_met'])
        self.assertLess(r['exact_count_calls'],12)

    def test_high_precision_refinement(self):
        r=v.guided([0,1],1,14,F(1,10**30))
        self.assertTrue(r['target_met'])
        self.assertEqual(r['path'],'exact_local_refinement')
        self.assertTrue(v.v3.verify(r['certificate']))

    def test_bad_float_proposal_cannot_be_accepted(self):
        with patch.object(v,'proposals',return_value=(10.0,10.0,'upper_ritz')):
            r=v.guided([0,1],1,8)
        self.assertEqual(r['path'],'rational_fallback')
        self.assertLess(F(r['certificate']['upper']),0)
        self.assertTrue(v.v3.verify(r['certificate']))

    def test_true_tail_floor_does_not_trigger_useless_refinement(self):
        r=v.guided([0,1],1,2)
        self.assertFalse(r['target_met'])
        self.assertEqual(r['path'],'tail_bound_limited')

    def test_explicit_budget_failure(self):
        r=v.search([0,20,-100,0,1],max_modes=4)
        self.assertEqual(r['status'],'target_not_met')
        self.assertTrue(v.v3.verify(r['certificate']))

if __name__=='__main__':unittest.main()
