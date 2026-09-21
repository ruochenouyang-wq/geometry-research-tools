from fractions import Fraction as F
from copy import deepcopy
import unittest
import moment_v12 as v


class MomentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result=v.search([0],[[0,1]],['1/50'],[['1/5']],modes=10,tolerance=F(1,1000))

    def test_confinement_is_exact_and_global(self):
        c=self.result['certificate'];self.assertTrue(v.verify(c))
        self.assertEqual(c['moment_confined_box'],[['-51/10','49/10']])
        self.assertEqual(c['status'],'epsilon_global')
        self.assertLessEqual(F(c['candidate_upper']),F(c['boxed_certificate']['candidate_upper']))

    def test_parameter_elimination_matches_joint_energy(self):
        p=v.g.problem([0],[[0,1]],[[0,1]],linear=['1/50'],hessian=[['1/5']])
        coeff=[F(1),F(1,3)];inverse,_,_=v.confinement(p,[v.g.v5.make_range([0,1])])
        energy,parameters,m=v.eliminated_value(p,coeff,inverse)
        mu0,_=v.moments(p,coeff)
        self.assertEqual(energy,mu0+sum(a*b for a,b in zip(parameters,m))+v.g.penalty(p,parameters))
        for delta in [F(-1),F(1),F(1,3)]:
            x=[parameters[0]+delta]
            self.assertGreater(mu0+x[0]*m[0]+v.g.penalty(p,x),energy)

    def test_wrong_inverse_moment_or_domain_rejected(self):
        for key,value in [('inverse_penalty_hessian',[['4']]),('candidate_moments',['0']),
                          ('moment_confined_box',[['-1','1']]),('candidate_upper','-100')]:
            c=deepcopy(self.result['certificate']);c[key]=value;self.assertFalse(v.verify(c))

    def test_positive_definite_not_just_gershgorin(self):
        # SPD, but its first row is not diagonally dominant.
        p=v.g.problem([0],[[0,1],[0,0,1]],[[-10,10]]*2,hessian=[[1,2],[2,5]])
        inverse,_,_=v.confinement(p,[v.g.v5.make_range([0,1]),v.g.v5.make_range([0,0,1])])
        self.assertEqual(inverse,[[5,-2],[-2,1]])
        with self.assertRaises(ValueError):v.search([0],[[0,1]],[0],[[0]])

    def test_constant_direction_does_not_remove_the_optimum(self):
        r=v.search([0],[[1]],[0],[[2]],modes=4,tolerance=F(1,10**6))
        self.assertEqual(r['certificate']['moment_confined_box'],[['-1/2','-1/2']])
        self.assertTrue(v.verify(r['certificate']))
        self.assertEqual(r['certificate']['candidate_upper'],'-1/4')

    def test_two_parameter_global_function_certificate(self):
        r=v.search([0],[[0,1],[0,0,1]],['1/50','-1/6'],[['1/5','1/50'],['1/50',1]],
                   modes=10,tolerance=F(1,100))
        c=r['certificate'];self.assertTrue(v.verify(c));self.assertEqual(c['status'],'epsilon_global')
        p=c['boxed_certificate']['problem'];z=[F(a)+F(b) for a,b in zip(p['penalty']['linear'],c['candidate_moments'])]
        a=list(map(F,c['candidate_parameters']));h=[list(map(F,row)) for row in p['penalty']['hessian']]
        for i in range(2):self.assertEqual(sum(h[i][j]*a[j] for j in range(2))+z[i],0)


if __name__=='__main__':unittest.main()
