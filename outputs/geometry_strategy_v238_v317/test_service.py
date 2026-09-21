"""Cross-component tests, including the public process boundary."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

import backend
from research_service import ResearchService, parse_request
from verification import verify_any

Q={'kind':'axis_profile','profile':'step','amplitude':'0','offset':'-5','axis':2}
R={'op':'research','target':'spectrum','function':Q,'tolerance':'1/1000'}


class ServiceTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.service=ResearchService(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_public_jsonl(self):
        data='\n'.join(['{',json.dumps(R|{'view':'full'}),
                        json.dumps({'op':'batch','requests':[{'op':'unknown'},R]}),
                        '{"op":"research","op":"fetch"}',
                        '{"op":"fetch","reference":"bad"}',json.dumps(R)])+'\n'
        run=subprocess.run([sys.executable,'-B',str(backend.ROOT/'research_service.py')],
                           input=data,text=True,capture_output=True,timeout=30,
                           env={**os.environ,'GEOMETRY_EVIDENCE_ROOT':self.tmp.name})
        self.assertEqual(run.returncode,0,run.stderr)
        rows=list(map(json.loads,run.stdout.splitlines()))
        self.assertEqual(len(rows),6)
        self.assertFalse(rows[0]['ok'])
        self.assertTrue(rows[1]['target_met'])
        self.assertFalse(rows[2]['responses'][0]['ok'])
        self.assertTrue(rows[2]['responses'][1]['certificate_valid'])
        self.assertFalse(rows[3]['ok'])
        self.assertFalse(rows[4]['ok'])
        self.assertTrue(rows[5]['target_met'])

    def test_exact_constant_fetch_verify_binding(self):
        summary=self.service.handle(R)
        before=self.service.runner.ledger.snapshot()['counts']['executor_calls']
        full=self.service.handle({'op':'fetch','reference':summary['evidence_ref']['sha256']})
        self.assertEqual(full['usage']['counts']['executor_calls'],before+1)
        self.assertIn('evidence_read',full['usage']['stages'])
        self.assertTrue(full['certificate_valid'])
        self.assertIsNone(full['target_met'])
        cert=full['certificate']
        self.assertEqual(cert['lower'],'-3')
        request={'op':'verify','certificate':cert,'function':Q,'expected_kind':'spectrum',
                 'mean_zero':True,'tolerance':'1/1000'}
        checked=self.service.handle(request)
        self.assertTrue(checked['verified'])
        self.assertTrue(checked['binding_complete'])
        self.assertTrue(checked['target_met'])
        wrong=self.service.handle(request|{'mean_zero':False})
        self.assertTrue(wrong['certificate_valid'])
        self.assertFalse(wrong['verified'])
        self.assertEqual(wrong['status'],'binding_mismatch')
        self.assertEqual(self.service.handle({'op':'verify','certificate':cert})['status'],'verified_unbound')

    def test_constant_checks_controls_before_shortcut(self):
        for field,value in [('budget',{'max_attempts':-1}),('budget',{'unknown':1}),
                            ('policy',{'min_relative_gain':'-1'}),('policy',{'unknown':1}),
                            ('budget',{'wall_seconds':True})]:
            with self.subTest(field=field,value=value),self.assertRaises(ValueError):
                self.service.handle(R|{field:value})

    def test_budget_state_survives_mathematical_success(self):
        raw=self.service.handle(R|{'view':'full'})
        fake={**raw,'status':'budget_exceeded','within_budget':False,
              'stop_reason':'wall_budget','cost':{'solver_calls':1},'certificate':raw['certificate']}
        # Simulate an over-budget computation at the execution/summary boundary.
        with patch.object(self.service,'_research',return_value=fake):
            result=self.service.handle(R|{'request_id':'over-budget'})
        self.assertTrue(result['target_met'])
        self.assertFalse(result['within_budget'])
        self.assertEqual(result['status'],'budget_exceeded')
        self.assertEqual(result['stop_reason'],'wall_budget')
        self.assertEqual(result['cost']['solver_calls'],1)
        zero=self.service.handle(R|{'budget':{'wall_seconds':0},'view':'full'})
        self.assertEqual(zero['status'],'budget_exceeded')

    def test_replay_time_is_inside_total_budget(self):
        from verification import assess
        def delayed(*args,**kw):
            time.sleep(0.02)
            return assess(*args,**kw)
        with patch('research_service.assess',side_effect=delayed):
            result=self.service.handle(R|{'budget':{'wall_seconds':0.01}})
        self.assertTrue(result['target_met'])
        self.assertEqual(result['status'],'budget_exceeded')
        self.assertFalse(result['within_budget'])

    def test_cache_does_not_cross_goals(self):
        first=self.service.handle(R)
        again=self.service.handle(R)
        changed=self.service.handle(R|{'mean_zero':False})
        self.assertFalse(first['interaction']['cache_hit'])
        self.assertTrue(again['interaction']['cache_hit'])
        self.assertFalse(changed['interaction']['cache_hit'])
        self.assertEqual(changed['bounds']['lower'],'-5')
        bad=self.service.handle({'op':'batch','requests':[{'op':'research'},R]})
        self.assertFalse(bad['responses'][0]['ok'])
        self.assertTrue(bad['responses'][1]['target_met'])

    def test_approximation_is_not_spectral_result(self):
        request={'op':'research','target':'approximation','function':Q,'tolerance':'1/1000'}
        result=self.service.handle(request|{'view':'full'})
        self.assertTrue(result['target_met'])
        self.assertFalse(verify_any(result['certificate'],expected_kind='spectrum'))
        with self.assertRaises(ValueError):
            self.service.handle(request|{'mean_zero':True})

    def test_json_boundary(self):
        for text in ('[]','{"x":NaN}','{"x":1e999}','{"x":1,"x":2}'):
            with self.subTest(text=text),self.assertRaises(ValueError):
                parse_request(text)

    def test_loose_spectral_goal_uses_cheap_probe_with_original_binding(self):
        q={'kind':'axis_profile','profile':'abs_power','amplitude':'-1/3',
           'offset':'0','axis':2,'exponent':'-1/5'}
        result=self.service.evaluate_task({'kind':'spectrum','function':q,'tolerance':'1/1000'})
        self.assertEqual(result['attempts'][0]['action']['route'],'spectrum')
        self.assertTrue(result['target_met'])
        self.assertTrue(verify_any(result['certificate'],expected_function=q,
                        expected_kind='spectrum',expected_mean_zero=True,expected_tolerance='1/1000'))


if __name__=='__main__':unittest.main()
