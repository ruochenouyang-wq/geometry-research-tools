from fractions import Fraction as F
from copy import deepcopy
import unittest
import v08_cluster as v

class ClusterTests(unittest.TestCase):
    def test_exact_sphere_subspace(self):
        c=v.cluster([0],2,5)['certificate']
        self.assertTrue(v.verify(c));self.assertEqual(c['projector_hilbert_schmidt_squared_upper'],'0')

    def test_basis_change_invariance(self):
        c=v.cluster([0,1],2,8)['certificate'];cols=[[F(x) for x in col] for col in c['trial_columns']]
        transformed=[[x+2*y for x,y in zip(*cols)],[x-y for x,y in zip(*cols)]]
        d=v.certificate([0,1],transformed,c['ritz_upper'],c['gap_certificate'])
        self.assertEqual(c['full_residual_frobenius_squared'],d['full_residual_frobenius_squared'])
        self.assertEqual(c['projector_hilbert_schmidt_squared_upper'],d['projector_hilbert_schmidt_squared_upper'])

    def test_rank_deficient_basis_rejected(self):
        with self.assertRaises(ValueError):v.statistics([F(0)],[[F(1),F(0)],[F(2),F(0)]])

    def test_cluster_avoids_internal_small_gap(self):
        q=[0,'1/1000000',-100]
        one=v.cluster(q,1,20)['certificate'];two=v.cluster(q,2,20)['certificate']
        self.assertGreater(F(two['external_lower'])-F(two['ritz_upper']),1)
        self.assertLess(F(two['projector_hilbert_schmidt_squared_upper']),
                        F(one['projector_hilbert_schmidt_squared_upper'])/10**6)

    def test_forged_gap_and_residual_rejected(self):
        original=v.cluster([0,1],2,6)['certificate']
        for key,value in [('cluster_indices',[2,3]),('full_residual_frobenius_squared','0'),
                          ('ritz_upper','-100'),('trial_columns',[['1'],['1']])]:
            c=deepcopy(original);c[key]=value;self.assertFalse(v.verify(c))
        c=deepcopy(original);c['gap_certificate']['eigenvalue_index']=4;self.assertFalse(v.verify(c))

    def test_tail_is_not_omitted(self):
        cols=[[F(1),F(0)],[F(0),F(1)]]
        _,_,r=v.statistics([F(0),F(1)],cols)
        self.assertEqual(r,F(4,15))

if __name__=='__main__':unittest.main()
