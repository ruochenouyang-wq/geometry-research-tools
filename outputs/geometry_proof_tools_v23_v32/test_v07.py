from fractions import Fraction as F
from copy import deepcopy
import unittest
import v07_precision as v

class PrecisionTests(unittest.TestCase):
    def test_exact_linear_solve(self):
        self.assertEqual(v.solve([[2,1],[1,3]],[1,2]),[F(1,5),F(3,5)])

    def test_breaks_fixed_coefficient_plateau(self):
        r=v.refine([0,1],14,96,1)
        before,after=r['before_certificate'],r['certificate']
        self.assertGreater(F(before['exact_width']),F(1,10**30))
        self.assertLess(F(after['exact_width']),F(1,10**40))
        self.assertTrue(v.verify(after))

    def test_known_constant_function_remains_exact(self):
        r=v.refine([F(-2,3)],4)
        self.assertEqual(r['certificate']['lower'],'-2/3')
        self.assertEqual(r['certificate']['exact_width'],'0')

    def test_refinement_does_not_remove_spatial_error(self):
        r=v.refine([0,1],4)
        self.assertGreater(F(r['certificate']['exact_width']),F(1,10**8))
        self.assertGreater(F(r['certificate']['outside_residual_squared']),0)

    def test_corrupted_coefficients_and_gap_rejected(self):
        original=v.refine([0,1],12)['certificate']
        for key,value in [('trial_legendre_coefficients',['1']),('gap_lower','100'),('eigenvalue_index',2)]:
            c=deepcopy(original);c[key]=value;self.assertFalse(v.verify(c))

if __name__=='__main__':unittest.main()
