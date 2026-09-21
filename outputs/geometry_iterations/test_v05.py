from copy import deepcopy
from fractions import Fraction as F
import unittest
import v05_range as v

class RangeTests(unittest.TestCase):
    def test_quadratic_extrema_exact(self):
        for q,expected in [([1,-2,3],(F(2,3),F(6))),([0,0,100],(F(0),F(100))),
                           ([0,0,-100],(F(-100),F(0))),([1,1],(F(0),F(2)))]:
            c=v.make_range(q)
            self.assertEqual(v.verify_bounds(list(map(F,q)),c),expected)

    def test_irrational_extrema_are_enclosed(self):
        q=[F(1,9),F(0),F(-2,3),F(0),F(1)]
        c=v.make_range(q,F(1,10**6),64)
        self.assertTrue(c['range_target_met'])
        lo,hi=v.verify_bounds(q,c)
        self.assertLessEqual(lo,0);self.assertGreaterEqual(hi,F(4,9))
        self.assertGreaterEqual(lo,-F(1,10**6))
        self.assertLess(len(c['intervals']),64)

    def test_missing_interval_is_rejected(self):
        q=list(map(F,[0,1,0,1,1]))
        c=v.make_range(q)
        c['intervals'].pop()
        with self.assertRaises(ValueError):v.verify_bounds(q,c)

    def test_forged_extremum_and_success_rejected(self):
        q=[F(1,9),F(0),F(-2,3),F(0),F(1)]
        original=v.make_range(q,F(1,10**12),1)
        self.assertFalse(original['range_target_met'])
        for key,value in [('lower','100'),('range_target_met',True)]:
            c=deepcopy(original);c[key]=value
            with self.assertRaises(ValueError):v.verify_bounds(q,c)

    def test_subdivision_matches_direct_affine_expansion(self):
        q=list(map(F,[1,2,3,-4,5]))
        left,right=v.e.split_bernstein(v.e.bernstein_coefficients(q))
        self.assertEqual(left,v.interval_coefficients(q,F(-1),F(0)))
        self.assertEqual(right,v.interval_coefficients(q,F(0),F(1)))

    def test_range_proof_is_bound_to_operator(self):
        import v03_tail as v3
        c=v3.certify([0,0,100],1,8,30,range_proof=v.make_range([0,0,100]))
        self.assertTrue(v3.verify(c))
        c['q_coefficients']=['0','0','-100']
        self.assertFalse(v3.verify(c))

if __name__=='__main__':unittest.main()
