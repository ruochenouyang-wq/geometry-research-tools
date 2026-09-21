import copy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch
from dependencies import ROOT
import research_tool as tool
import spectral_transfer as s


class ResearchToolTests(unittest.TestCase):
    def setUp(self):
        self.temporary=tempfile.TemporaryDirectory()
        self.service=tool.Service(self.temporary.name)
        self.request={'op':'solve','method':'l2','degree':0,
                      'function':{'kind':'axis_profile','profile':'step'},
                      'source_options':{'modes':2,'max_modes':2,'bits':12}}
    def tearDown(self):self.temporary.cleanup()

    def test_complete_original_function_jsonl_protocol(self):
        result=self.service.call(self.request)
        self.assertEqual((result['lower'],result['upper']),('1','4'))
        self.assertEqual(result['status'],'certified_open')
        record=self.service.get(result['certificate_id'])
        self.assertTrue(tool.verify_record(json.loads(json.dumps(record))))
        with patch.object(s,'solve',side_effect=AssertionError('no search')):
            replay=self.service.call({'op':'verify','certificate_id':result['certificate_id']})
            self.assertEqual(replay,result)

    def test_original_inputs_scope_order_and_target_cannot_change(self):
        result=self.service.call(self.request);record=self.service.get(result['certificate_id'])
        variants=[('degree',1),('degree',True),('mean_zero',False),('target_width','1'),
                  ('method','linf'),('source_options',{'ignored':1}),
                  ('function',{'kind':'axis_profile','profile':'step','axis':0})]
        for key,value in variants:
            bad=copy.deepcopy(record);bad['arguments'][key]=value
            self.assertFalse(tool.verify_record(bad),key)

    def test_integrity_identifier_and_unknown_fields_fail_closed(self):
        result=self.service.call(self.request);identifier=result['certificate_id']
        for bad in ['../x','a'*63,7]:
            with self.assertRaises(ValueError):self.service.get(bad)
        with self.assertRaises(ValueError):self.service.call({'op':'fetch','certificate_id':identifier,'path':'anything'})
        p=self.service.path(identifier);record=json.loads(p.read_text());record['certificate']['upper']='99'
        p.write_text(json.dumps(record))
        with self.assertRaises(ValueError):self.service.get(identifier)

    def test_invalid_solver_result_cannot_enter_store(self):
        result=self.service.call(self.request);record=self.service.get(result['certificate_id'])
        wrong=copy.deepcopy(record['certificate']);wrong['mean_zero']=False
        with patch.object(s,'solve',return_value=wrong):
            with self.assertRaises(ValueError):self.service.call(self.request)

    def test_executable_jsonl_rejection_and_success(self):
        # Exact constant avoids introducing changing timing fields into the proof.
        request={'op':'solve','function':{'kind':'analytic_sum','polynomial':{'0,0,0':'0'}},
                 'order':0,'source_options':{'modes':2,'max_modes':2,'bits':12}}
        p=subprocess.run([sys.executable,'-B',str(ROOT/'research_tool.py')],
                         input=json.dumps(request)+'\n'+json.dumps({'op':'not_supported'})+'\n',
                         text=True,capture_output=True,check=True,timeout=30)
        rows=[json.loads(line) for line in p.stdout.splitlines()]
        self.assertEqual(len(rows),2)
        self.assertEqual((rows[0]['lower'],rows[0]['upper']),('2','2'))
        self.assertEqual(rows[1]['status'],'not_certified')


if __name__=='__main__':unittest.main()
