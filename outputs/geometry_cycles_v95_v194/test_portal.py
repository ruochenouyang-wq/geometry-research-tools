"""End-to-end checks of scope, proof reuse, content identity and continuation."""
from copy import deepcopy
from fractions import Fraction as F
import tempfile
import unittest
from unittest.mock import patch
import portal
from common import projected


class PortalTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.service = portal.Service(self.folder.name)

    def test_ground_two_reused_decisions(self):
        result = self.service.call({'op':'ground','q':[0,0,2],'tolerance':'1/1000000'})
        self.assertTrue(result['verified'])
        with patch.object(portal.precision, 'full_ground', side_effect=AssertionError('Repeated solve')), \
             patch.object(projected, 'full_ground', side_effect=AssertionError('Repeated solve')):
            for t, status in [('119/50','proved'), ('12/5','refuted_by_spectral_existence')]:
                answer = self.service.call({'op':'threshold','certificate_id':result['certificate_id'],
                                            'q':[0,0,2],'threshold':t})
                self.assertEqual(answer['status'],status)
                self.assertEqual(answer['new_spectral_solves'],0)
        with self.assertRaises(ValueError):
            self.service.call({'op':'threshold','certificate_id':result['certificate_id'],
                               'q':[0,0,3],'threshold':'12/5'})

    def test_invalid_scope_rejected_for_reuse(self):
        for source in (projected.certify_sector([0],m=0,modes=2,bits=8),
                       projected.full_ground([0],mean_zero=False,modes=2,max_modes=2,bits=8)):
            self.assertTrue(portal.verify(source))
            with self.assertRaises(ValueError):
                portal.threshold(source,[0],1)

    def test_hash_rejects_modified_disk_proof(self):
        result = self.service.call({'op':'ground','q':[0],'max_steps':0})
        path = self.service.store.path(result['certificate_id'])
        path.write_text('{}')
        with self.assertRaises(ValueError):
            self.service.call({'op':'verify','certificate_id':result['certificate_id']})

    def test_no_arbitrary_paths_or_unknown_options(self):
        for request in ({'op':'fetch','certificate_id':'../../notes'},
                        {'op':'ground','q':[0],'directory':'/tmp/example'},
                        {'op':'eval','expression':'1+1'},
                        {'op':'ground','q':[0.0]}):
            with self.assertRaises((ValueError,TypeError)):
                self.service.call(request)

    def test_resume_keeps_bounds_and_problem_binding(self):
        old = self.service.call({'op':'ground','q':[0,1,1],'max_steps':0})
        options = {'op':'resume','certificate_id':old['certificate_id'],
                   'q':[0,1,1],'max_steps':3,'tolerance':'1/1000000'}
        result = self.service.call(options)
        self.assertGreaterEqual(F(result['lower']),F(old['lower']))
        self.assertLessEqual(F(result['upper']),F(old['upper']))
        with self.assertRaises(ValueError):
            self.service.call({**options,'q':[1]})

    def test_counterexample_is_separately_verified(self):
        result = self.service.call({'op':'counterexample','q':[0,0,2],'threshold':'12/5'})
        c = self.service.store.get(result['certificate_id'])
        self.assertTrue(portal.witness.verify(c,expected_q=[0,0,2],expected_threshold='12/5'))

    def test_all_main_backend_entry_points(self):
        requests = [
            {'op':'wide_ground','q':[0]*8+[1],'modes':2,'max_modes':2,'bits':16},
            {'op':'exponential','a':1,'degree':8,'modes':2,'max_modes':2,'bits':16},
            {'op':'family','q0':[0],'directions':[[0,0,1]],'box':[[0,1]],'threshold':2,'modes':2,'bits':16},
            {'op':'global_minimum','objective':{'q0':['0'],'direction':['0'],
                'domain':['-1','1'],'penalty':['0','0','1']},'max_leaves':2,'bits':16,'modes':2},
            {'op':'constraints','q':[0],'threshold':6,'cutoff':1,'modes':2,'bits':16},
            {'op':'ordered_spectrum','q':[0],'k':9,'modes':2,'max_modes':2,'max_radial':2,'bits':16},
            {'op':'gn_ratio','raw':[0,1]},
            {'op':'gn_refine','raw':[0,1],'iterations':0},
            {'op':'gn_diagnostic','raw':[0,1],'modes':2,'max_modes':2,'bits':16},
            {'op':'anisotropic','q':{'1,1,0':'1'},'L':1,'bits':16},
            {'op':'radius','q':[0],'radius':2,'threshold':'1/2','modes':2,'max_modes':2,'bits':16},
        ]
        for request in requests:
            with self.subTest(op=request['op']):
                answer = self.service.call(request)
                c = self.service.store.get(answer['certificate_id'])
                self.assertTrue(portal.verify(c))

    def test_radius_and_gn_summaries_retain_meaning(self):
        result = self.service.call({'op':'radius','q':[0],'radius':2,'threshold':'1/2'})
        self.assertEqual(result['status'],'proved')
        self.assertEqual(result['problem']['coordinates'],'normalized_t')
        self.assertEqual(result['radius_request'],{'radius':'2'})
        result = self.service.call({'op':'gn_ratio','raw':[0,1]})
        self.assertEqual(F(result['pi_K_lower']),F(9,40))
        self.assertFalse(result['universal_GN_sharp_constant_proved'])
        result = self.service.call({'op':'gn_refine','raw':[0,1],'iterations':0})
        self.assertEqual(F(result['final_pi_K_lower']),F(9,40))
        self.assertEqual(result['accepted_steps'],0)

    def test_threshold_tampering_does_not_verify(self):
        source = portal.precision.full_ground([0],max_steps=0)
        c = portal.threshold(source,[0],1)
        for key, value in [('scope','single_sector'),('new_spectral_solves',1),
                           ('status','refuted_by_spectral_existence'),('threshold','3')]:
            bad = deepcopy(c); bad[key] = value
            self.assertFalse(portal.verify(bad),key)


if __name__ == '__main__':
    unittest.main()
