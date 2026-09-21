"""Original-problem binding across the four new mathematical interfaces."""
from copy import deepcopy
from fractions import Fraction as F
from math import comb
import tempfile
import unittest
from unittest.mock import patch
import research_tool as tool

C12=[{'0,0,0':1},{'0,0,1':1}]
C18={'2,0,0':1,'0,2,0':2,'0,0,2':3}


class ResearchToolTests(unittest.TestCase):
    def setUp(self):
        self.directory=tempfile.TemporaryDirectory();self.addCleanup(self.directory.cleanup)
        self.service=tool.Service(self.directory.name)

    def test_c12_actual_full_space(self):
        out=self.service.call({'op':'constraints','q':{},'constraints':C12})
        self.assertEqual((F(out['lower']),F(out['upper'])),(2,2))
        self.assertIn('all_real_H1',out['scope'])
        out=self.service.call({'op':'constrained_inequality','q':{},'constraints':C12,'threshold':'21/10'})
        self.assertEqual(out['status'],'refuted_with_trial')

    def test_c14_both_original_spaces(self):
        out=self.service.call({'op':'spectrum','q':[-7]})
        self.assertTrue(out['all_requested_targets_met'])
        for row in out['spaces']:
            n,total=(8,20) if row['space']=='full_mean_zero' else (9,27)
            self.assertEqual((row['count']['count_lower'],row['count']['count_upper']),(n,n))
            self.assertEqual((F(row['negative_trace']['lower']),F(row['negative_trace']['upper'])),(total,total))

    def test_c16_original_monomial_polynomial(self):
        n=64;u=[F(comb(n,k)) for k in range(n+1)];u[0]-=F(2**n,n+1)
        out=self.service.call({'op':'cap_original','n':n,'raw_t_coefficients':list(map(str,u))})
        self.assertEqual(F(out['probability_ratio']),F(402777216,209564225))
        self.assertEqual(4*F(out['pi_times_area_K']),F(out['probability_ratio']))
        self.assertFalse(out['universal_GN_sharp_constant_proved'])

    def test_c18_target_preserved(self):
        out=self.service.call({'op':'parity_ground','q':C18,'target':'1/100000000'})
        self.assertLessEqual(F(out['exact_width']),F(1,10**8))
        self.assertEqual(out['requested_width'],'1/100000000')

    def test_record_does_not_accept_different_potential_or_constraints(self):
        out=self.service.call({'op':'constraints','q':{},'constraints':C12})
        saved=self.service.store.get(out['certificate_id'])
        for name,value in [('q',{'0,0,0':1}),('constraints',[]),('L',3)]:
            bad=deepcopy(saved);bad['arguments'][name]=value
            self.assertFalse(tool.verify(bad),name)

    def test_spectrum_scope_and_quantity_binding(self):
        out=self.service.call({'op':'spectrum','q':[-7],'spaces':['full_unprojected'],'quantities':['trace']})
        saved=self.service.store.get(out['certificate_id'])
        for name,value in [('spaces',['full_mean_zero']),('quantities',['ordered']),('k',10),('q',[-6])]:
            bad=deepcopy(saved);bad['arguments'][name]=value
            self.assertFalse(tool.verify(bad),name)

    def test_cap_original_scale_binding(self):
        out=self.service.call({'op':'cap','n':64,'original':True})
        saved=self.service.store.get(out['certificate_id'])
        for name,value in [('n',63),('original',False)]:
            bad=deepcopy(saved);bad['arguments'][name]=value
            self.assertFalse(tool.verify(bad),name)

    def test_replay_performs_no_new_search(self):
        out=self.service.call({'op':'parity_ground','q':C18})
        with patch.object(tool.parity,'adaptive_ground',side_effect=AssertionError('new search')):
            result=self.service.call({'op':'verify','certificate_id':out['certificate_id']})
            self.assertTrue(result['mathematical_proof_verified'])

    def test_hash_or_request_corruption_rejected(self):
        out=self.service.call({'op':'cap','n':4})
        self.service.store.path(out['certificate_id']).write_text('{}')
        with self.assertRaises((ValueError,KeyError)):
            self.service.call({'op':'fetch','certificate_id':out['certificate_id']})
        for request in ({'op':'fetch','certificate_id':'../../x'},
                        {'op':'cap','n':64,'directory':'/tmp'},
                        {'op':'evaluate','expression':'1+1'}):
            with self.assertRaises((ValueError,TypeError)):
                self.service.call(request)

    def test_bad_producer_cannot_return_another_valid_proof(self):
        wrong=tool.gn.evaluate_cap(4)
        with patch.dict(tool.OPERATIONS,{'cap':(lambda n,original=True:wrong,{'n','original'})}):
            with self.assertRaises(ValueError):self.service.call({'op':'cap','n':64})


if __name__=='__main__':unittest.main()
