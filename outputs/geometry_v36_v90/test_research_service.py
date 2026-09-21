"""V90 cross-format acceptance, mathematical scope and actual JSONL integration."""
import copy,json,subprocess,sys,tempfile,unittest
from pathlib import Path
from fractions import Fraction as F
import research_registry as registry
import protocol,run,token_meter
import exact_extrema,uniform_refine,positive_search,precision_bridge,matrix_refine,original_goal

ROOT=Path(__file__).resolve().parent
GOAL={'profile':'sphere_ground','parameters':['a'],'potential':'a*t','penalty':'a*a/5',
      'domain':{'kind':'all_real'},'threshold':'0'}


class ResearchServiceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        lo=positive_search.lower_certificate([0,0,'20/9'],[[0,0,1]],['-1/3'],range_backend=precision_bridge.ExactRangeBackend())
        hi=positive_search.ansatz_upper_certificate([0,0,'20/9'],[[0,0,1]],[{'t':'1/2','weight':'1'}])
        cls.positive=positive_search.optimization_certificate(lo,hi,range_backend=precision_bridge.ExactRangeBackend())
        cls.matrix=precision_bridge.precise_matrix(matrix_refine.refine([[0,1]],budget=0))
        cls.extrema=exact_extrema.maximize([0,1,0,-1],max_nodes=8)
        cls.uniform=uniform_refine.synthesize_log_sobolev()
        cls.goal=original_goal.solve(protocol.expand_request(GOAL))

    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.s=protocol.Service(self.tmp.name)
    def tearDown(self):self.tmp.cleanup()

    def test_registered_new_formats_verify_without_search(self):
        for cert in (self.extrema,self.positive,self.matrix,self.uniform,self.goal):
            self.assertTrue(registry.verify(cert))
            self.assertTrue(protocol.checked(cert))

    def test_malformed_or_unknown_proofs_rejected(self):
        for value in (None,[],{'format':'made_up'},{'method':uniform_refine.LOG_SOBOLEV_METHOD}):
            self.assertFalse(registry.verify(value))

    def test_unknown_backend_not_loaded_from_json(self):
        c=copy.deepcopy(self.positive);c['lower_certificate']['range_backend']='some_external.module'
        self.assertFalse(registry.verify(c))
        with self.assertRaises(ValueError):registry.backend('some_external.module')

    def test_positive_ansatz_upper_never_becomes_spectral_upper(self):
        r=self.s.receipt(self.positive)
        self.assertIsNone(r['bounds']['spectral_upper'])
        self.assertFalse(r['scope']['ansatz_upper_is_spectral_upper'])
        self.assertLessEqual(F(r['bounds']['spectral_lower']),F(23,36))
        self.assertGreaterEqual(F(r['bounds']['ansatz_upper']),F(23,36))
        c=copy.deepcopy(self.positive);c['spectral_upper']=c['ansatz_upper']
        self.assertFalse(registry.verify(c))

    def test_positive_actual_rayleigh_bound_remains_separate(self):
        c=positive_search.optimization_certificate(self.positive['lower_certificate'],self.positive['ansatz_upper_certificate'],
                rayleigh=positive_search.rayleigh_certificate([0,0,'20/9'],[1]),range_backend=precision_bridge.ExactRangeBackend())
        r=self.s.receipt(c)
        self.assertGreater(F(r['bounds']['spectral_upper']),F(r['bounds']['ansatz_upper']))

    def test_new_extrema_is_polynomial_not_spectral(self):
        r=self.s.receipt(self.extrema)
        self.assertEqual(r['kind'],'polynomial_maximum')
        self.assertNotIn('geometry',r['scope'])
        self.assertEqual(self.s.inspect(r['task'])['interval'],['-1','1'])

    def test_original_goal_exact_infimum_and_trust_disclosed(self):
        r=self.s.receipt(self.goal)
        self.assertEqual(r['bounds'],{'lower':'0','upper':'0'})
        self.assertEqual(r['target'],'0')
        self.assertIn('log-Sobolev',r['scope']['analytic_dependency'])
        self.assertEqual(self.s.inspect(r['task']),protocol.expand_request(GOAL))

    def test_uniform_counterexample_is_actual_disproof(self):
        c=uniform_refine.synthesize_log_sobolev(coefficient='1/7');r=self.s.receipt(c)
        self.assertEqual(r['status'],'disproved')
        self.assertLess(F(r['bounds']['counterexample_objective_upper']),0)
        self.assertEqual(c['counterexample']['spectral_objective_upper'],'-1/91')

    def test_ansatz_obstruction_never_marked_original_disproof(self):
        c=uniform_refine.synthesize(coefficient='1/5',order=1);r=self.s.receipt(c)
        self.assertEqual(r['status'],'ansatz_obstruction')
        self.assertEqual(r['bounds'],{})

    def test_matrix_globality_scope_preserved(self):
        r=self.s.receipt(self.matrix)
        self.assertEqual(r['scope']['ansatz'],'poisson_exponential')
        self.assertEqual(r['scope']['objective'],'trace(C G)')

    def test_each_new_receipt_roundtrip_and_full_wire_budget(self):
        for cert in (self.extrema,self.positive,self.matrix,self.uniform,self.goal):
            r=self.s.receipt(cert);b=self.s.budget(r,500)
            self.assertTrue(b['fits']);self.assertLessEqual(b['wire_tokens'],500)
            self.assertEqual(b['wire_tokens'],token_meter.count(token_meter.wire(b)))
            self.assertEqual(json.loads(b['text']),r)
            self.assertEqual(self.s.inspect(r['evidence'])['certificate'],cert)

    def test_tiny_budget_has_no_invented_verdict(self):
        b=self.s.budget(self.s.receipt(self.uniform),16)
        self.assertFalse(b['fits']);self.assertNotIn('status',b)

    def test_all_new_task_references_restore_in_fresh_service(self):
        receipts=[self.s.receipt(c) for c in (self.extrema,self.positive,self.matrix,self.uniform,self.goal)]
        goals=[r['task'] for r in receipts];cp=self.s.checkpoint(goals,receipts)
        restored=protocol.Service(self.tmp.name).restore(cp)
        self.assertEqual(restored['results'],receipts);self.assertEqual(restored['goals'],goals)

    def test_new_evidence_cannot_be_rebound_to_another_goal(self):
        wrong=self.s.register(dict(GOAL,potential='2*a*t'))
        with self.assertRaises(ValueError):self.s.receipt(self.goal,wrong)

    def test_rounding_direction_for_every_bound_kind(self):
        c=positive_search.optimization_certificate(self.positive['lower_certificate'],self.positive['ansatz_upper_certificate'],
                rayleigh=positive_search.rayleigh_certificate([0,0,'20/9'],[1]),range_backend=precision_bridge.ExactRangeBackend())
        exact=protocol.summary(c);display=protocol.summary(c,compact_numbers=True)
        for key,value in exact['bounds'].items():
            if key.endswith('_lower'):self.assertLessEqual(F(display['bounds'][key]),F(value))
            else:self.assertGreaterEqual(F(display['bounds'][key]),F(value))

    def test_service_research_original_goal(self):
        out=run.call(self.s,{'op':'research','goal':GOAL})
        self.assertEqual(json.loads(out['text'])['status'],'proved')

    def test_service_barta_no_log_sobolev_dependency(self):
        out=run.call(self.s,{'op':'research','goal':GOAL,'rule':'barta'})
        packet=json.loads(out['text']);self.assertEqual(packet['status'],'proved')
        self.assertNotIn('analytic_dependency',packet['scope'])

    def test_service_refine_uses_exact_range_and_old_dual(self):
        out=run.call(self.s,{'op':'refine','spec':{'directions':[[0,1]]},'budget':0})
        c=self.s.inspect(json.loads(out['text'])['evidence'])['certificate']
        self.assertEqual(c['method'],'precise_global_matrix_v88');self.assertEqual(c['gap'],'0')

    def test_legacy_optimize_and_lemma_request_operations(self):
        for request in ({'op':'optimize','spec':{'directions':[[0,1]]},'max_rounds':1},
                        {'op':'lemma','spec':{'directions':[[0,1]]}}):
            out=run.call(self.s,request);self.assertTrue(out['fits'])
            self.assertTrue(protocol.checked(self.s.inspect(json.loads(out['text'])['evidence'])))

    def test_actual_jsonl_new_operations_and_invalid_recovery(self):
        requests=[{'op':'uniform','spec':{'coefficient':'1/7'}},
                  {'op':'uniform','spec':{'direction':[0,0,1]}},
                  {'op':'extrema','spec':{'polynomial':[0,1,0,-1]}},
                  {'op':'research','goal':GOAL},
                  {'op':'positive','spec':{'potential':[0,0,'20/9']},
                   'options':{'degree':2,'grid_size':9,'max_exchanges':1,'max_newton':1,'range_backend':'exact_sturm_minimum_v1'}}]
        p=subprocess.run([sys.executable,str(ROOT/'run.py'),'serve','--store',self.tmp.name],
                         input='\n'.join(map(json.dumps,requests))+'\n',text=True,capture_output=True)
        self.assertEqual(p.returncode,0,p.stderr);rows=list(map(json.loads,p.stdout.splitlines()))
        self.assertEqual(len(rows),len(requests));self.assertIsNone(rows[1]['mathematical_verdict'])
        self.assertEqual(json.loads(rows[0]['text'])['status'],'disproved')
        for i in (0,2,3,4):self.assertTrue(rows[i]['fits'])

    def test_cli_verifies_new_and_rejects_tampered_formats(self):
        file=Path(self.tmp.name)/'cert.json';file.write_text(json.dumps(self.uniform))
        p=subprocess.run([sys.executable,str(ROOT/'run.py'),'verify',str(file)],text=True,capture_output=True)
        self.assertEqual(p.returncode,0,p.stdout+p.stderr)
        bad=copy.deepcopy(self.uniform);bad['coefficient']='1/7';file.write_text(json.dumps(bad))
        p=subprocess.run([sys.executable,str(ROOT/'run.py'),'verify',str(file)],text=True,capture_output=True)
        self.assertEqual(p.returncode,1)

    def test_unsupported_geometry_and_hidden_defaults_rejected(self):
        with self.assertRaises(ValueError):run.call(self.s,{'op':'research','goal':dict(GOAL,geometry='plane')})
        with self.assertRaises(ValueError):run.call(self.s,{'op':'research','goal':{k:v for k,v in GOAL.items() if k!='profile'}})


if __name__=='__main__':unittest.main()
