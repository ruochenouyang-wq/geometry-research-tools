from fractions import Fraction as F
import unittest
import v06_parity as v

class ParityTests(unittest.TestCase):
    def test_counts_equal_global_counts(self):
        for q in [[0],[0,0,-100],[0,0,10,0,-3]]:
            proof=v.v5.make_range(q)
            a=v.ParityKernel(q,8,range_proof=proof)
            b=v.v3.Kernel(q,8,range_proof=proof)
            for x in [F(-110),F(-80),F(0),F(1)]:
                for kind in ['lower','upper_ritz','upper_tail']:
                    if kind!='upper_ritz' and x>=a.beta:continue
                    self.assertEqual(a.count(x,kind),b.count(x,kind,dense=True))

    def test_all_indices_include_odd_functions(self):
        for k in [1,2,3,4]:
            r=v.search([0],k,max_modes=8)
            exact=k*(k-1);c=r['certificate']
            self.assertLessEqual(F(c['lower']),exact);self.assertGreaterEqual(F(c['upper']),exact)
        # The second eigenvalue is 2, not the second even eigenvalue 6.

    def test_nearly_symmetric_is_rejected(self):
        with self.assertRaises(ValueError):v.search([0,'1/1000000000',-100])

    def test_global_verifier_accepts_without_parity_assumption(self):
        r=v.search([0,0,-100],2,max_modes=24)
        self.assertEqual(r['status'],'target_met')
        self.assertTrue(v.v3.verify(r['certificate'],independent=True))

if __name__=='__main__':unittest.main()
