"""Post-reveal guard controls; never modifies the frozen C3 implementation."""
from contextlib import nullcontext
from copy import deepcopy
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
import subprocess
import signal
import sys
import time
from types import SimpleNamespace
import unittest
from unittest.mock import Mock, patch

ROOT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('c31_guard_under_test',ROOT/'controller.py')
guard=importlib.util.module_from_spec(spec);spec.loader.exec_module(guard)
TASK={'kind':'spectrum','function':{'kind':'axis_profile','profile':'abs_power','axis':0,
    'amplitude':'4/5','offset':'0','exponent':'-3/8'},'tolerance':'1/1000000',
    'mean_zero':True,'budget':{'wall_seconds':10}}


def certificate(width='2',task=None):
    t=TASK if task is None else task
    return {'exact_width':width,'function':deepcopy(t['function']),
            'mean_zero':t['mean_zero'],'tolerance':t['tolerance']}


def original(width='2',met=False):
    return {'version':'frozen_c3','route':'native_c3','certificate':certificate(width),
            'certificate_valid':True,'target_met':met,'status':'target_met' if met else 'certified_open',
            'attempts':[],'counters':{},'total_elapsed_seconds':0.01,'within_budget':True}


def assess(cert,task):
    if (type(cert) is not dict or cert.get('function')!=task['function'] or
        cert.get('mean_zero') is not task['mean_zero'] or cert.get('tolerance')!=task['tolerance'] or cert.get('damaged')):
        return {'certificate_valid':False,'target_met':False,'status':'verification_failed'}
    width=F(cert['exact_width']);met=width<=F(task['tolerance'])
    return {'certificate_valid':True,'target_met':met,'status':'target_met' if met else 'certified_open',
            'metric':str(width),'tolerance':task['tolerance']}


def child_response(cert):
    return {'outcome':'returned','response':{'certificate':cert,'certificate_valid':True,'target_met':True,
            'status':'untrusted_child_status'},'allocated_seconds':8,'elapsed_seconds':0.01}


class GuardControls(unittest.TestCase):
    def run_mock(self,base=None,child=None,task=None,child_error=None,replay=None):
        base=original() if base is None else base;task=deepcopy(TASK if task is None else task)
        baseline=SimpleNamespace(normalize_task=lambda x:deepcopy(x),assess=Mock(side_effect=replay or assess))
        c3=SimpleNamespace(solve=Mock(return_value=deepcopy(base)))
        worker=Mock(return_value=child if child is not None else child_response(certificate('1')),side_effect=child_error)
        events=[]
        with patch.object(guard,'_frozen_scope',return_value=nullcontext()), \
             patch.object(guard,'_modules',return_value=(c3,baseline)),patch.object(guard,'_baseline_child',worker):
            result=guard.solve(task,progress=events.append)
        return result,worker,baseline,events

    def test_narrower_verified_open_replaces(self):
        result,worker,baseline,events=self.run_mock()
        self.assertTrue(result['guard']['replaced']);self.assertEqual(result['certificate']['exact_width'],'1')
        self.assertTrue(result['certificate_valid']);self.assertFalse(result['target_met'])
        worker.assert_called_once();baseline.assess.assert_called_once()
        passed,seconds=worker.call_args.args
        self.assertEqual(passed['function'],TASK['function']);self.assertEqual(passed['tolerance'],TASK['tolerance'])
        self.assertAlmostEqual(passed['budget']['wall_seconds'],seconds)
        self.assertLessEqual(seconds,10-guard.FINAL_RESERVE)
        self.assertTrue(any(e['action']=='guard_parent_frozen_replay' for e in events))

    def test_equal_and_wider_preserve_original(self):
        for width in ('2','3'):
            result,worker,_,_=self.run_mock(child=child_response(certificate(width)))
            self.assertEqual(result['certificate'],original()['certificate']);self.assertFalse(result['guard']['replaced'])
            self.assertEqual(result['guard']['reason'],'fallback_no_strict_width_improvement');worker.assert_called_once()

    def test_damaged_or_wrong_task_evidence_preserves_original(self):
        for cert in (dict(certificate('0'),damaged=True),dict(certificate('0'),tolerance='1'),
                     dict(certificate('0'),function=dict(TASK['function'],axis=1))):
            result,_,baseline,_=self.run_mock(child=child_response(cert))
            self.assertEqual(result['certificate'],original()['certificate']);baseline.assess.assert_called_once()
            self.assertEqual(result['guard']['reason'],'fallback_original_verifier_rejected')

    def test_timeout_preserves_original_without_assessing_partial_response(self):
        result,worker,baseline,_=self.run_mock(child={'outcome':'timeout','response':{'certificate':certificate('0')}})
        self.assertEqual(result['certificate'],original()['certificate']);self.assertTrue(result['certificate_valid'])
        self.assertEqual(result['guard']['reason'],'fallback_timeout');worker.assert_called_once();baseline.assess.assert_not_called()

    def test_worker_exception_preserves_original(self):
        result,worker,_,_=self.run_mock(child_error=RuntimeError('mock failed launch'))
        self.assertEqual(result['certificate'],original()['certificate'])
        self.assertEqual(result['guard']['reason'],'fallback_exception');worker.assert_called_once()

    def test_parent_replay_exception_preserves_original(self):
        def failed(*args):raise ArithmeticError('mock damaged replay')
        result,_,_,_=self.run_mock(replay=failed)
        self.assertEqual(result['certificate'],original()['certificate']);self.assertEqual(result['guard']['reason'],'fallback_exception')

    def test_slow_parent_replay_times_out_and_preserves_original(self):
        def slow(cert,task):time.sleep(0.3);return assess(cert,task)
        task=deepcopy(TASK);task['budget']['wall_seconds']=0.1
        previous=signal.getsignal(signal.SIGALRM)
        with patch.object(guard,'FINAL_RESERVE',0.03),patch.object(guard,'OUTPUT_RESERVE',0.02), \
             patch.object(guard,'MIN_FALLBACK',0.01):
            result,worker,_,_=self.run_mock(task=task,replay=slow)
        self.assertEqual(result['guard']['reason'],'fallback_parent_replay_timeout')
        self.assertEqual(result['certificate'],original()['certificate'])
        self.assertTrue(result['certificate_valid']);self.assertFalse(result['guard']['replaced'])
        self.assertIs(signal.getsignal(signal.SIGALRM),previous)
        self.assertEqual(signal.getitimer(signal.ITIMER_REAL),(0.0,0.0));worker.assert_called_once()

    def test_missing_certificate_preserves_original(self):
        result,_,baseline,_=self.run_mock(child=child_response(None))
        self.assertEqual(result['certificate'],original()['certificate']);baseline.assess.assert_called_once()

    def test_insufficient_budget_does_not_start(self):
        task=deepcopy(TASK);task['budget']['wall_seconds']=1.2
        result,worker,baseline,_=self.run_mock(task=task)
        worker.assert_not_called();baseline.assess.assert_not_called()
        self.assertEqual(result['guard']['reason'],'insufficient_remaining_budget_for_protected_fallback')

    def test_existing_success_does_not_start(self):
        result,worker,_,_=self.run_mock(base=original('0',True))
        worker.assert_not_called();self.assertTrue(result['target_met'])

    def test_scope_and_validity_gate(self):
        negative=deepcopy(TASK);negative['function']['amplitude']='-1'
        unconstrained=deepcopy(TASK);unconstrained['mean_zero']=False
        for task in (negative,unconstrained):
            result,worker,_,_=self.run_mock(task=task);worker.assert_not_called();self.assertFalse(result['guard']['eligible'])
        result,worker,_,_=self.run_mock(base=dict(original(),certificate_valid=False))
        worker.assert_not_called();self.assertFalse(result['guard']['eligible'])

    def test_final_cost_includes_calls_and_serialization(self):
        def replay(cert,task):time.sleep(0.004);return assess(cert,task)
        result,_,_,_=self.run_mock(replay=replay)
        self.assertGreaterEqual(result['total_elapsed_seconds'],0.004)
        self.assertGreaterEqual(result['total_cpu_seconds'],result['parent_cpu_seconds'])
        self.assertTrue(any(x['action']=='response_json_boundary' for x in result['attempts']))
        self.assertIsNone(result['actual_model_tokens'])

    def test_direct_child_timeout_inherits_group_and_is_reaped(self):
        real_popen=subprocess.Popen;created=[]
        def launch(command,**kwargs):
            self.assertIs(kwargs.get('start_new_session'),False)
            proc=real_popen([sys.executable,'-B','-c','import time; time.sleep(2)'],**kwargs)
            created.append(proc);return proc
        with patch.object(guard.subprocess,'Popen',side_effect=launch):
            result=guard._baseline_child(TASK,0.05)
        self.assertEqual(result['outcome'],'timeout');self.assertIsNotNone(created[0].poll())
        self.assertLess(result['elapsed_seconds'],0.8)

    def test_wrong_runtime_name_restored_after_private_loading(self):
        marker=SimpleNamespace(__file__='/not/the/frozen/runtime.py')
        before=sys.modules.get('runtime')
        try:
            sys.modules['runtime']=marker
            with guard._frozen_scope():
                c3,baseline=guard._modules()
                self.assertEqual(Path(sys.modules['runtime'].__file__).resolve(),guard.C3_ROOT/'runtime.py')
                self.assertEqual(Path(baseline.__file__).resolve(),guard.BASELINE_ROOT/'solver.py')
            self.assertIs(sys.modules['runtime'],marker)
        finally:
            if before is None:sys.modules.pop('runtime',None)
            else:sys.modules['runtime']=before


class PublicPostRevealRegression(unittest.TestCase):
    def test_n29_and_e18_original_tolerance_ten_seconds_once(self):
        e18=deepcopy(TASK);e18['function'].update(amplitude='2/3',exponent='-2/7')
        diagnostics=[]
        for name,task in (('N29',deepcopy(TASK)),('E18',e18)):
            result=guard.solve(task)
            self.assertTrue(result['certificate_valid'],result.get('status'))
            self.assertTrue(result['within_budget'],result['total_elapsed_seconds'])
            if name=='N29':
                self.assertTrue(result['guard']['attempted']);self.assertTrue(result['guard']['replaced'],result['guard']['reason'])
                self.assertLess(F(result['certificate']['exact_width']),F(result['guard']['original_exact_width']))
            else:
                self.assertTrue(result['target_met']);self.assertFalse(result['guard']['attempted'])
            diagnostics.append({'case':name,'label':'post-reveal functional regression, not blind comparison',
                'status':result['status'],'certificate_valid':result['certificate_valid'],'target_met':result['target_met'],
                'exact_width_float':float(F(result['certificate']['exact_width'])),
                'guard_reason':result['guard']['reason'],'guard_attempted':result['guard']['attempted'],
                'guard_replaced':result['guard']['replaced'],'elapsed_seconds':result['total_elapsed_seconds'],
                'total_cpu_seconds':result['total_cpu_seconds'],'parent_cpu_seconds':result['parent_cpu_seconds'],
                'child_cpu_seconds':result['child_cpu_seconds']})
        print('POST_REVEAL_FUNCTIONAL='+json.dumps(diagnostics,allow_nan=False))


if __name__=='__main__':unittest.main()
