from fractions import Fraction as F
from copy import deepcopy
import unittest
import v03_tail as v


class TailTests(unittest.TestCase):
    def test_known_sphere_and_shifts(self):
        for k in [1,2,5]:
            c=v.certify([F(-2,3)],k,6,40)
            exact=(k-1)*k-F(2,3)
            self.assertLessEqual(F(c['lower']),exact)
            self.assertGreaterEqual(F(c['upper']),exact)
            self.assertTrue(v.verify(c,independent=True))

    def test_single_tail_column_by_hand(self):
        ker=v.Kernel([0,1],2)
        a=ker.matrix(F(-1),'lower')
        b=ker.matrix(F(-1),'upper_tail')
        finite=ker.matrix(F(-1),'upper_ritz')
        self.assertEqual(finite[1][1]-a[1][1],F(8,45)/6)
        self.assertEqual(finite[1][1]-b[1][1],F(8,45)/8)

    def test_graded_improves_scalar_in_matrix_order(self):
        a=v.Kernel([1,-2,3],4,'scalar_ritz').matrix(F(-1),'lower')
        b=v.Kernel([1,-2,3],4,'graded_ritz').matrix(F(-1),'lower')
        diff=[[b[i][j]-a[i][j] for j in range(4)] for i in range(4)]
        self.assertEqual(v.base.inertia(diff)[0],0)

    def test_both_sides_improve_coarse_bound(self):
        a=v.certify([0,1],1,2,48,'scalar_ritz')
        b=v.certify([0,1],1,2,48,'sandwich')
        self.assertLess(F(b['exact_width']),F(a['exact_width'])/2)
        self.assertLess(F(b['upper']),F(a['upper']))

    def test_matches_independent_previous_enclosure(self):
        for k in [1,2,3]:
            a=v.old.certify_banded([0,1],k,8,40)
            b=v.certify([0,1],k,8,40)
            self.assertLessEqual(max(F(a['lower']),F(b['lower'])),min(F(a['upper']),F(b['upper'])))

    def test_boundary_positivity_is_mandatory(self):
        ker=v.Kernel([0,1000],2)
        for kind in ['lower','upper_tail']:
            with self.assertRaises(ValueError):ker.matrix(ker.beta,kind)
        ker.matrix(ker.beta,'upper_ritz')

    def test_forged_fields_rejected(self):
        c=v.certify([0,1],2,6,30)
        for key,value in [('eigenvalue_index',1),('tail_lower','10000'),('upper_kind','invented'),
                          ('scope','arbitrary manifolds'),('lower','100')]:
            bad=deepcopy(c);bad[key]=value
            self.assertFalse(v.verify(bad),key)
        bad=deepcopy(c);bad['range_proof']['lower']='10'
        self.assertFalse(v.verify(bad))

    def test_invalid_sizes(self):
        for k,n in [(0,2),(3,2),(True,5),(1,65)]:
            with self.assertRaises(ValueError):v.sizes(k,n)

    def test_lower_comparison_cannot_certify_an_upper_endpoint(self):
        kernel=v.Kernel([0,1],2)
        # This is above the lower comparison root, but below the true eigenvalue.
        with self.assertRaises(ValueError):v.certificate(kernel,1,F(-1),F(-159,1000),'lower')


if __name__=='__main__':unittest.main()
