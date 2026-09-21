from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import unittest
import global_search as g
import global_quadratic as q


def tilted(box=None):
    return g.problem([0],[[0,1]],box or [[-6,6]],linear=['1/50'],hessian=[['1/5']])


class QuadraticTests(unittest.TestCase):
    def test_integer_inputs_do_not_trigger_float_division(self):
        self.assertIsInstance(q.value(0,[0],[[1]],[1]),F)
        self.assertEqual(q.value(0,[0],[[1]],[1]),F(1,2))

    def test_interior_minimum(self):
        val,x=q.minimum(F(0),[F(-1),F(-1)],[[F(2),F(1)],[F(1),F(2)]],[(F(0),F(1))]*2)
        self.assertEqual(val,F(-1,3));self.assertEqual(x,(F(1,3),F(1,3)))

    def test_indefinite_quadratic_requires_boundary(self):
        val,x=q.minimum(0,[0,0],[[0,1],[1,0]],[(-1,1)]*2)
        self.assertEqual(val,-1)

    def test_singular_minimizer_set_reaches_face(self):
        val,x=q.minimum(0,[0,0],[[2,-2],[-2,2]],[(0,1)]*2)
        self.assertEqual(val,0)

    def test_minimum_on_edge_not_corner(self):
        val,x=q.minimum(F(0),[F(-1),F(-1)],[[F(2),F(0)],[F(0),F(0)]],[(F(0),F(1))]*2)
        self.assertEqual(x,(F(1,2),F(1)));self.assertEqual(val,F(-5,4))

    def test_exact_affine_change(self):
        c,l,h=q.on_unit_box(F(2),[F(1),F(-2)],[[F(3),F(2)],[F(2),F(-1)]],[(F(-3),F(1)),(F(1),F(2))])
        self.assertEqual(q.value(c,l,h,[F(1,3),F(2,5)]),
                         q.value(F(2),[F(1),F(-2)],[[F(3),F(2)],[F(2),F(-1)]],[F(-5,3),F(7,5)]))


class GlobalTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.full=g.search(tilted(),F(1,1000),modes=10,max_leaves=128,seeds=[[3]])['certificate']
        cls.right=g.search(tilted([[0,6]]),F(1,1000),modes=10,max_leaves=128)['certificate']

    def test_global_gap_and_wrong_basin(self):
        self.assertTrue(g.verify(self.full));self.assertTrue(g.verify(self.right))
        self.assertEqual(self.full['status'],'epsilon_global')
        self.assertLess(F(self.full['candidate_parameters'][0]),0)
        self.assertLess(F(self.full['candidate_upper']),F(self.right['global_lower'])-F(1,10))
        self.assertLessEqual(F(self.full['global_gap']),F(1,1000))

    def test_positive_side_contains_a_true_interior_local_minimum(self):
        # Continuity and boundary values above an interior feasible value ensure
        # an interior minimum on [0,6], hence a local minimum of the full problem.
        upper=F(self.right['candidate_upper'])
        for a in [F(0),F(6)]:
            record=next(x for x in self.right['points'] if x['at']==[str(a)])
            lower=F(record['spectral_certificate']['lower'])+g.penalty(self.right['problem'],[a])
            self.assertGreater(lower,upper)

    def test_global_declaration_cannot_outlive_its_cover(self):
        bad=deepcopy(self.full);bad['partition_tree'].pop('right');self.assertFalse(g.verify(bad))
        bad=deepcopy(self.full);bad['partition_tree']['split']='7';self.assertFalse(g.verify(bad))
        bad=deepcopy(self.right);bad['problem']['box']=[['-6','6']];self.assertFalse(g.verify(bad))

    def test_forged_bound_objective_and_trial_rejected(self):
        for key,value in [('global_lower','0'),('candidate_upper','-100'),('global_gap','0'),
                          ('optimization_scope','all geometries'),('candidate_parameters',['2'])]:
            bad=deepcopy(self.full);bad[key]=value;self.assertFalse(g.verify(bad),key)
        bad=deepcopy(self.full);bad['problem']['penalty']['linear']=['0'];self.assertFalse(g.verify(bad))
        bad=deepcopy(self.full);bad['points'][0]['trial_legendre']=['1'];self.assertFalse(g.verify(bad))

    def test_second_eigenvalue_cannot_use_ground_concavity(self):
        bad=deepcopy(self.full);bad['points'][0]['spectral_certificate']['eigenvalue_index']=2
        self.assertFalse(g.verify(bad))

    def test_corner_directions_and_cross_term(self):
        p=g.problem([0],[[1],[2]],[[-2,2],[-1,1]],linear=[-1,-2],hessian=[[2,1],[1,2]])
        c=g.search(p,F(1,10**6),modes=4,max_leaves=1)['certificate']
        self.assertEqual(c['status'],'epsilon_global')
        self.assertLessEqual(F(c['global_lower']),0);self.assertGreaterEqual(F(c['candidate_upper']),0)
        self.assertLess(F(c['global_gap']),F(1,10**6))

    def test_two_nonconstant_directions(self):
        p=g.problem([0],[[0,1],[0,0,1]],[[-4,4],[-1,1]],linear=['1/50','-1/6'],
                    hessian=[['1/5','1/50'],['1/50',1]])
        c=g.search(p,F(1,100),modes=10,max_leaves=128)['certificate']
        self.assertTrue(g.verify(c));self.assertEqual(c['status'],'epsilon_global')

    def test_budget_failure_remains_a_valid_global_bracket(self):
        c=g.search(tilted(),F(1,1000),modes=8,max_leaves=1)['certificate']
        self.assertEqual(c['status'],'global_gap_open');self.assertTrue(g.verify(c))
        bad=deepcopy(c);bad['status']='epsilon_global';self.assertFalse(g.verify(bad))

    def test_concavity_bound_encloses_interior_values(self):
        p=tilted();oracle=g.Oracle(p,10,F(1,10**8));box=[[F(-6),F(6)]]
        ids=[oracle.point([x]) for x in [-6,6]]
        lower,_=g.relaxation(p,box,oracle.records,ids,'concavity',[(F(-1),F(1))])
        for a in [-5,-3,0,2,4]:self.assertLessEqual(lower,oracle.values[oracle.point([a])][0])

    def test_invalid_global_inputs(self):
        with self.assertRaises(ValueError):g.problem([0],[[0,1]],[[2,1]])
        with self.assertRaises(ValueError):g.search(tilted(),seeds=[[7]])
        with self.assertRaises(ValueError):g.problem([0],[[0,1]]*3,[[-1,1]]*3)
        with self.assertRaises(ValueError):g.problem([0],[[1],[2]],[[-1,1]]*2,hessian=[[1,2],[0,1]])

    def test_constant_objective_keeps_all_optimizer_boxes(self):
        p=g.problem([0],[[1]],[[-2,2]],linear=[-1])
        c=g.search(p,F(1,10**6),modes=4,max_leaves=1)['certificate']
        self.assertEqual(c['excluded_leaf_count'],0)
        self.assertEqual(c['possible_optimizer_boxes'],[[['-2','2']]])


if __name__=='__main__':unittest.main()
