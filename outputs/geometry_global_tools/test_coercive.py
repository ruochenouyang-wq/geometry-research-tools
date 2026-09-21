from fractions import Fraction as F
from copy import deepcopy
import unittest
import coercive as v


class CoercivityTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=v.search([0],[[0,1]],['1/50'],[['1/5']],modes=10,tolerance=F(1,1000))

    def test_global_over_all_reals(self):
        c=self.result['certificate'];self.assertTrue(v.verify(c));self.assertEqual(c['status'],'epsilon_global')
        self.assertGreater(F(c['exterior_lower']),F(c['candidate_upper']))
        self.assertEqual(c['boxed_certificate']['problem']['box'],[['-11','11']])

    def test_unbounded_growth_condition_is_required(self):
        with self.assertRaises(ValueError):v.search([0],[[0,1]],[0],[[0]],modes=6)
        with self.assertRaises(ValueError):v.search([0],[[0,1]],[0],[[-1]],modes=6)

    def test_radius_cannot_be_invented_or_insufficient(self):
        boxed=v.g.search(v.g.problem([0],[[0,1]],[[-6,6]],linear=['1/50'],hessian=[['1/5']]),modes=8)['certificate']
        with self.assertRaises(ValueError):v.certificate(boxed,v.g.v5.make_range([0]))
        for key,value in [('exterior_lower','100'),('hessian_lower','1'),('optimization_scope','all manifolds')]:
            bad=deepcopy(self.result['certificate']);bad[key]=value;self.assertFalse(v.verify(bad))

    def test_two_dimensional_exterior_bound(self):
        p=v.g.problem([0],[[0,1],[0,0,1]],[[-12,12]]*2,linear=[0,0],hessian=[[1,'1/10'],['1/10',1]])
        proof=v.g.v5.make_range([0]);directions=[v.g.v5.make_range([0,1]),v.g.v5.make_range([0,0,1])]
        h,b,offset=v.exterior_data(p,proof,directions)
        self.assertEqual(h,F(9,10));self.assertGreater(v.exterior_lower(p,h,b,offset),0)


if __name__=='__main__':unittest.main()
