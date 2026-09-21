import unittest
from unittest.mock import patch
from copy import deepcopy
import precision_tool as tool


class ToolTests(unittest.TestCase):
    def test_precise_singular_end_to_end_and_original_binding(self):
        q={'kind':'axis_profile','profile':'abs_power','exponent':'-1/4','amplitude':'-1'}
        response=tool.handle({'op':'precise_singular','function':q})
        cert=response['certificate']
        self.assertEqual(cert['status'],'target_met')
        self.assertEqual(len(cert['trial']['powers']),10)
        with patch.object(tool.enriched,'enrich',side_effect=AssertionError('no search during replay')):
            self.assertTrue(tool.handle({'op':'verify','certificate':cert,'function':q,'mean_zero':True})['verified'])
        self.assertFalse(tool.handle({'op':'verify','certificate':cert,'function':{**q,'amplitude':'1'}})['verified'])
        self.assertFalse(tool.handle({'op':'verify','certificate':cert,'mean_zero':False})['verified'])
        with self.assertRaises(ValueError):tool.handle({'op':'precise_singular','function':{**q,'axis':0}})

    def test_certificate_roundtrip_without_search(self):
        q={'kind':'axis_profile','profile':'step'}
        c=tool.handle({'op':'spectrum','function':q,'modes':2,'bits':16})['certificate']
        with patch.object(tool.direct,'full_ground',side_effect=AssertionError('no solving during replay')):
            self.assertTrue(tool.handle({'op':'verify','certificate':c,'function':q,'mean_zero':True})['verified'])
        self.assertFalse(tool.handle({'op':'verify','certificate':c,'function':{**q,'amplitude':'-1'}})['verified'])

    def test_approximation_is_not_a_spectrum(self):
        response=tool.handle({'op':'approximate','function':{'kind':'axis_profile','profile':'step'}})
        self.assertEqual(response['kind'],'function_L2_approximation_only')
        self.assertTrue(tool.handle({'op':'verify','certificate':response['certificate']})['verified'])
        with self.assertRaises(ValueError):
            tool.handle({'op':'verify','certificate':response['certificate'],'mean_zero':True})

    def test_unknown_fields_and_missing_function(self):
        for request in ({'op':'spectrum','function':{},'error_upper':'0'},
                        {'op':'approximate'},{'op':'verify','certificate':{},'sampled':True},
                        {'op':'solve'},[],None):
            with self.subTest(request=request),self.assertRaises(ValueError):tool.handle(request)

    def test_bad_certificate_never_becomes_verified(self):
        for cert in ({},None,[],{'format':'untrusted'}):
            self.assertFalse(tool.handle({'op':'verify','certificate':cert})['verified'])


if __name__=='__main__':unittest.main()
