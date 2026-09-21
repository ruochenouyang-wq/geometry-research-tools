from fractions import Fraction as F
from unittest.mock import patch
import unittest
import v10_adaptive as v

class AdaptiveTests(unittest.TestCase):
    def test_global_index_and_independent_verification(self):
        for k in [1,2,3]:
            r=v.adaptive([0,1],k,max_modes=16)
            self.assertEqual(r['status'],'target_met');self.assertTrue(v.v3.verify(r['certificate'],independent=True))

    def test_malicious_size_predictor_cannot_fake_accuracy(self):
        for proposal in [-100,10**9,None]:
            with patch.object(v,'next_size',return_value=proposal):r=v.adaptive([0,20,-100,0,1],max_modes=8)
            self.assertEqual(r['status'],'target_not_met')
            self.assertTrue(v.v3.verify(r['certificate']))
            self.assertTrue(all(x['modes']<=8 for x in r['attempts']))

    def test_function_precision_routed_by_residual(self):
        r=v.function([0,1],start_modes=14,max_modes=14)
        self.assertEqual(r['status'],'target_met')
        self.assertIn('coefficient_refinement',[x['action'] for x in r['attempts']])
        self.assertTrue(v.v7.verify(r['certificate']))

    def test_spatial_error_gets_larger_space(self):
        r=v.function([0,1],F(1,10**12),4,12)
        self.assertEqual(r['status'],'target_met')
        self.assertEqual(r['attempts'][0]['modes'],4)
        self.assertGreater(r['attempts'][-1]['modes'],4)
        self.assertNotIn('coefficient_refinement',[x['action'] for x in r['attempts']])

    def test_parity_is_exact(self):
        self.assertTrue(v.adaptive([0,0,1],max_modes=12)['parity_blocks'])
        self.assertFalse(v.adaptive([0,'1/1000000000',1],max_modes=12)['parity_blocks'])

    def test_function_budget_failure(self):
        r=v.function([0,1],start_modes=4,max_modes=4)
        self.assertEqual(r['status'],'target_not_met');self.assertTrue(v.v7.verify(r['certificate']))

if __name__=='__main__':unittest.main()
