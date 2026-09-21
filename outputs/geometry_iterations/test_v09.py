from fractions import Fraction as F
from copy import deepcopy
import unittest
import v09_family as v

class FamilyTests(unittest.TestCase):
    def test_uniform_true_claim(self):
        r=v.family([0],[0,1],threshold=F(-3,10),max_cells=32)
        self.assertEqual(r['certificate']['status'],'proved');self.assertTrue(v.verify(r['certificate']))
        self.assertLess(r['point_evaluations'],32)

    def test_verified_counterexample(self):
        c=v.family([0],[0,1],threshold=F(-1,10))['certificate']
        self.assertEqual(c['status'],'disproved');self.assertTrue(v.verify(c))
        self.assertLess(F(c['point_certificate']['upper']),F(-1,10))

    def test_small_budget_is_unresolved(self):
        c=v.family([0],[0,1],max_cells=1)['certificate']
        self.assertEqual(c['status'],'unresolved');self.assertTrue(v.verify(c))

    def test_cover_hole_overlap_and_false_success_rejected(self):
        c=v.family([0],[0,1])['certificate']
        for mutate in [lambda d:d['leaves'].pop(),
                       lambda d:d['leaves'][0]['interval'].__setitem__(1,'1'),
                       lambda d:d.__setitem__('lipschitz_upper','0'),
                       lambda d:d['leaves'][0].__setitem__('uniform_lower','100')]:
            bad=deepcopy(c);mutate(bad);self.assertFalse(v.verify(bad))
        bad=v.family([0],[0,1],max_cells=1)['certificate'];bad['status']='proved';self.assertFalse(v.verify(bad))

    def test_anchor_operator_and_parameter_binding(self):
        c=v.family([0],[0,1],threshold=F(-1,10))['certificate']
        for key,value in [('parameter','0'),('parameter','2'),('direction_coefficients',['0','2']),('eigenvalue_index',2)]:
            bad=deepcopy(c);bad[key]=value;self.assertFalse(v.verify(bad))

    def test_zero_perturbation_is_exact(self):
        c=v.family([3],[0],threshold=2,max_cells=1)['certificate']
        self.assertEqual(c['status'],'proved');self.assertEqual(c['lipschitz_upper'],'0')

if __name__=='__main__':unittest.main()
