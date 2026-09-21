"""Regression and attack checks for the acceptance revision."""
from copy import deepcopy
from contextvars import copy_context
import json
import os
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from release_backend import ROOT,PREVIOUS,backend,verify_any,previous_service
from release_service import ResearchService
from request_runtime import ReplaySession,RequestRuntime

Q=deepcopy(backend.enriched.FUNCTION)
CONSTANT={'kind':'axis_profile','profile':'step','amplitude':'0','offset':'-5','axis':2}


class ScopedReplayTests(unittest.TestCase):
    def setUp(self):
        self.runtime=RequestRuntime('test-checker-version')
        self.source=json.loads((backend.FROZEN/'certificates/singular_n4.json').read_text())

    def test_content_binding_mutation_and_boolean_scope(self):
        r=self.runtime
        with r.session() as session:
            self.assertTrue(r.verify(self.source,expected_function=Q,expected_mean_zero=True))
            self.assertTrue(r.verify(deepcopy(self.source),expected_function=Q,expected_mean_zero=True))
            self.assertEqual(session.counts['hits'],1)
            for kw in ({'expected_function':dict(Q,amplitude='-1/2')},
                       {'expected_mean_zero':False},{'expected_mean_zero':1},
                       {'expected_tolerance':'1/2'},{'expected_kind':'approximation'}):
                with self.subTest(kw=kw):self.assertFalse(r.verify(self.source,**kw))
            modified=deepcopy(self.source)
            modified['lower']=modified['upper'];modified['exact_width']='0'
            self.assertFalse(r.verify(modified))
            modified=deepcopy(self.source)
            modified['verified']=True
            self.assertFalse(r.verify(modified))

    def test_request_cache_lifetime_and_exception_cleanup(self):
        r=self.runtime
        with r.session() as first:
            self.assertTrue(r.verify(self.source))
        with r.session() as second:
            self.assertTrue(r.verify(self.source))
            self.assertEqual(second.counts['hits'],0)
            self.assertIsNot(first,second)
        with self.assertRaises(RuntimeError):
            with r.session():raise RuntimeError('interrupted request')
        self.assertIsNone(r.current.get())

    def test_copied_context_cannot_extend_request_cache_lifetime(self):
        r=self.runtime
        with r.session() as first:
            self.assertTrue(r.verify(self.source))
            inherited=copy_context()
        self.assertTrue(first.closed)
        self.assertEqual(first.entries,{})
        with patch.object(r,'uncached_verify',wraps=r.uncached_verify) as checker:
            self.assertTrue(inherited.run(r.verify,self.source))
            self.assertEqual(checker.call_count,1)

    def test_cache_is_bounded_and_never_caches_failure(self):
        session=ReplaySession('v',max_entries=1,max_bytes=10000)
        calls=[]
        def checker(certificate,expected_function=None):
            calls.append(certificate)
            return certificate.get('valid') is True
        self.assertTrue(session.call('check',checker,{'valid':True,'id':1},{}))
        self.assertTrue(session.call('check',checker,{'valid':True,'id':2},{}))
        self.assertTrue(session.call('check',checker,{'valid':True,'id':2},{}))
        for _ in range(2):self.assertFalse(session.call('check',checker,{'valid':False},{}))
        self.assertEqual(len(calls),5)
        self.assertEqual(session.counts['rejected'],2)
        self.assertEqual(len(session.entries),1)

    def test_checker_cannot_mutate_the_owned_cache_key_or_callers_evidence(self):
        session=ReplaySession('v')
        evidence={'valid':True}
        def checker(certificate):
            certificate['valid']=False
            return True
        self.assertTrue(session.call('mutating-test-only',checker,evidence,{}))
        self.assertEqual(evidence,{'valid':True})
        self.assertTrue(session.call('mutating-test-only',checker,evidence,{}))
        self.assertEqual(session.counts['hits'],1)

    def test_no_shared_module_mutation(self):
        import adaptive_singular
        import verification
        original=(adaptive_singular.direct,verification.verify_any,
                  previous_service.verify_any,previous_service.assess,backend.direct.verify_full)
        r=RequestRuntime('isolated-other')
        with r.session():self.assertTrue(r.verify(self.source))
        current=(adaptive_singular.direct,verification.verify_any,
                 previous_service.verify_any,previous_service.assess,backend.direct.verify_full)
        self.assertEqual(original,current)
        self.assertIsNot(r.adaptive,adaptive_singular)


class ReleaseServiceTests(unittest.TestCase):
    def setUp(self):self.tmp=tempfile.TemporaryDirectory()
    def tearDown(self):self.tmp.cleanup()

    def test_singular_candidate_sequence_and_full_evidence_unchanged(self):
        task={'kind':'spectrum','function':Q,'tolerance':'1/100000000'}
        old=previous_service.ResearchService(self.tmp.name,cache=False).evaluate_task(task)
        new=ResearchService(self.tmp.name,cache=False).evaluate_task(task)
        self.assertEqual(old['certificate'],new['certificate'])
        self.assertEqual([a['action'] for a in old['attempts']],
                         [a['action'] for a in new['attempts']])
        self.assertTrue(new['target_met'])
        self.assertGreater(new['replay_reuse']['hits'],0)
        self.assertTrue(verify_any(new['certificate'],expected_function=Q,
                        expected_mean_zero=True,expected_tolerance=task['tolerance']))

    def test_disabled_reuse_is_exact_ablation(self):
        task={'kind':'spectrum','function':dict(Q,exponent='-1/8'),'tolerance':'1/100000000'}
        enabled=ResearchService(self.tmp.name,cache=False).evaluate_task(task)
        disabled=ResearchService(self.tmp.name,cache=False,replay_reuse=False).evaluate_task(task)
        self.assertEqual(enabled['certificate'],disabled['certificate'])
        self.assertEqual(disabled['replay_reuse']['hits'],0)
        self.assertFalse(disabled['replay_reuse']['enabled'])

    def test_invalid_and_exhausted_budgets_stay_distinct_from_proof(self):
        service=ResearchService(self.tmp.name,cache=False)
        request={'op':'research','target':'spectrum','function':CONSTANT}
        for extra in ({'budget':{'max_attempts':0}},{'mean_zero':1},{'policy':None}):
            with self.subTest(extra=extra),self.assertRaises((ValueError,TypeError)):
                service.handle(request|extra)
        late=service.handle(request|{'budget':{'wall_seconds':0}})
        self.assertTrue(late['certificate_valid']);self.assertTrue(late['target_met'])
        self.assertFalse(late['within_budget']);self.assertEqual(late['status'],'budget_exceeded')

    def test_summary_fetch_and_wrong_original_goal(self):
        service=ResearchService(self.tmp.name,cache=False)
        result=service.handle({'op':'research','target':'spectrum','function':CONSTANT})
        retrieved=service.handle({'op':'fetch','reference':result['evidence_ref']['sha256']})
        self.assertTrue(retrieved['certificate_valid']);self.assertIsNone(retrieved['target_met'])
        checked=service.handle({'op':'verify','certificate':retrieved['certificate'],
                                'function':CONSTANT,'mean_zero':False,'expected_kind':'spectrum',
                                'tolerance':'1/1000'})
        self.assertTrue(checked['certificate_valid']);self.assertFalse(checked['verified'])

    def test_actual_process_boundary_continues_after_bad_request(self):
        requests=['{',json.dumps({'op':'research','target':'approximation','function':CONSTANT}),
                  json.dumps({'op':'batch','requests':[{'op':'unknown'},
                    {'op':'research','target':'spectrum','function':CONSTANT}]})]
        result=subprocess.run([sys.executable,'-B',str(ROOT/'release_service.py')],
                   input='\n'.join(requests)+'\n',text=True,capture_output=True,timeout=20,
                   env={**os.environ,'GEOMETRY_EVIDENCE_ROOT':self.tmp.name})
        self.assertEqual(result.returncode,0,result.stderr)
        rows=list(map(json.loads,result.stdout.splitlines()))
        self.assertEqual(len(rows),3);self.assertFalse(rows[0]['ok']);self.assertTrue(rows[1]['target_met'])
        self.assertFalse(rows[2]['responses'][0]['ok']);self.assertTrue(rows[2]['responses'][1]['target_met'])


if __name__=='__main__':unittest.main()
