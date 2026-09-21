import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import evidence_store as e

ROOT = Path(__file__).resolve().parent


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='.evidence-test-', dir=ROOT)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.store = e.EvidenceStore(self.root)

    def test_298_full_content_dedup(self):
        a = {'format':'candidate', 'nested':{'upper':'2', 'lower':'1'}}
        b = {'nested':{'lower':'1', 'upper':'2'}, 'format':'candidate'}
        self.assertEqual(self.store.put(a), self.store.put(b))
        self.assertEqual(len(list(self.root.glob('*.json'))), 1)
        c = {'format':'candidate', 'nested':{'upper':'2', 'lower':'0'}}
        self.assertNotEqual(self.store.put(c), self.store.put(a))
        self.assertEqual(self.store.get(self.store.put(a)), a)

    def test_299_strict_json_and_complete_summary(self):
        for bad in ({1:'numeric key'}, {'x':float('nan')}, {'x':float('inf')}, ('tuple',)):
            with self.subTest(bad=str(bad)):
                with self.assertRaises(ValueError): self.store.put(bad)
        a={'format':'proposal','large_evidence':list(range(1000))}
        b={'format':'proposal','large_evidence':list(range(999))+[-1]}
        sa=e.summary(a); sb=e.summary(b)
        self.assertNotEqual(sa['sha256'],sb['sha256'])
        self.assertFalse(sa['summary_is_certificate'])
        self.assertTrue(sa['complete_content_hashed'])
        self.assertEqual(sa,self.store.summary(self.store.put(a)))
        self.assertNotEqual(e.digest({'n':True}),e.digest({'n':1}))
        self.assertNotEqual(e.digest({'n':1.0}),e.digest({'n':1}))

    def test_300_tamper_and_failed_atomic_publish(self):
        value={'format':'proposal','bound':'7/8'}
        ref=self.store.put(value)
        (self.root/(ref+'.json')).write_text('{"bound":"1","format":"proposal"}')
        with self.assertRaises(ValueError): self.store.get(ref)
        with self.assertRaises(ValueError): self.store.put(value)
        new={'format':'proposal','bound':'8/9'}
        with patch.object(e.os,'link',side_effect=OSError('injected publish failure')):
            with self.assertRaises(OSError): self.store.put(new)
        self.assertFalse((self.root/(e.digest(new)+'.json')).exists())
        self.assertFalse(list(self.root.glob('.pending-*')))
        self.assertEqual(self.store.get(self.store.put(new)),new)

    def test_301_paths_sizes_and_unknown_format(self):
        for ref in ('../escape','A'*64,'0'*63,'/tmp/example',True):
            with self.subTest(ref=ref):
                with self.assertRaises(ValueError): self.store.get(ref)
        ref='0'*64
        target=self.root/'ordinary.json'; target.write_text('{}')
        (self.root/(ref+'.json')).symlink_to(target)
        with self.assertRaises((ValueError,OSError)): self.store.get(ref)
        small=e.EvidenceStore(self.root/'small',max_bytes=64)
        with self.assertRaises(ValueError): small.put({'large':'x'*100})
        oversized=small.root/(ref+'.json'); oversized.write_bytes(b'x'*65)
        with self.assertRaises(ValueError): small.get(ref)
        unknown=self.store.put({'format':'future_math_certificate','verified':True})
        self.assertEqual(self.store.format_status(unknown)['status'],'unknown_format_unverified')
        self.assertFalse(self.store.format_status(unknown)['verified'])
        chain=[]
        for _ in range(100): chain=[chain]
        with self.assertRaises(ValueError): self.store.put(chain)

    def test_302_default_math_replay_and_callback_boundary(self):
        import backend
        valid=backend.pieces.piecewise_model({'kind':'axis_profile','profile':'step'})
        ref=self.store.put(valid)
        with patch.object(backend.direct,'full_ground',side_effect=AssertionError('spectral search disabled')), \
             patch.object(backend.enriched,'enrich',side_effect=AssertionError('trial search disabled')):
            self.assertTrue(self.store.verify(ref)['verified'])
        invalid=json.loads(json.dumps(valid)); invalid['error_upper']='-1'
        self.assertFalse(self.store.verify(self.store.put(invalid))['verified'])
        unknown=self.store.put({'format':'not-known','verified':True})
        self.assertEqual(self.store.verify(unknown)['status'],'unknown_format')
        for callback in (lambda c:1, lambda c:{'verified':True}):
            custom=e.EvidenceStore(self.root,callback)
            self.assertFalse(custom.verify(ref)['verified'])
        custom=e.EvidenceStore(self.root,lambda c:c==valid)
        self.assertTrue(custom.verify(ref)['verified'])
        def broken(c): raise RuntimeError('callback fault')
        self.assertEqual(e.EvidenceStore(self.root,broken).verify(ref)['status'],'verifier_error')

    def test_303_cache_binds_full_content_version_and_integrity(self):
        calls=[]
        def check(c): calls.append(c); return c.get('value')==1
        cached=e.EvidenceStore(self.root,check,verifier_version='rule-v1')
        ref=cached.put({'format':'test-only','value':1})
        self.assertFalse(cached.verify(ref)['cache_hit'])
        self.assertTrue(cached.verify(ref)['cache_hit']); self.assertEqual(len(calls),1)
        cached.verifier_version='rule-v2'
        self.assertFalse(cached.verify(ref)['cache_hit']); self.assertEqual(len(calls),2)
        changed=cached.put({'format':'test-only','value':2})
        self.assertFalse(cached.verify(changed)['verified']); self.assertEqual(len(calls),3)
        (self.root/(ref+'.json')).write_bytes(e.canonical({'format':'test-only','value':0}))
        with self.assertRaises(ValueError): cached.verify(ref)
        self.assertEqual(len(calls),3)
        uncached=e.EvidenceStore(self.root,check)
        clean=uncached.put({'format':'test-only','value':1,'new':True})
        uncached.verify(clean); uncached.verify(clean)
        self.assertEqual(len(calls),5)
        cached.verifier=lambda c:False
        self.assertFalse(cached.verify(clean)['verified'])

    def test_304_problem_binding_including_cached_results(self):
        import backend
        value=json.loads((backend.FROZEN/'certificates'/'step_n4.json').read_text())
        ref=self.store.put(value); problem=e.problem_claim(value)
        self.assertTrue(self.store.verify(ref,problem)['verified'])

        self.assertTrue(self.store.verify(ref,problem)['cache_hit'])
        for key,new in (('geometry','unit_circle'),('scope','only_m0'),
                        ('mean_zero',False),('eigenvalue_index',2),
                        ('tolerance','1/1000000000000'),('measure','unnormalized_area')):
            bad=json.loads(json.dumps(problem)); bad[key]=new
            with self.subTest(key=key):
                self.assertEqual(self.store.verify(ref,bad)['status'],'problem_mismatch')
        bad=json.loads(json.dumps(problem)); bad['function']['amplitude']='2'
        self.assertEqual(self.store.verify(ref,bad)['status'],'problem_mismatch')
        partial={'function':problem['function']}
        self.assertFalse(self.store.verify(ref,partial)['verified'])
        bad=json.loads(json.dumps(problem)); bad['eigenvalue_index']=True
        self.assertFalse(self.store.verify(ref,bad)['verified'])
        self.assertEqual(value['status'],'certified_open')
        self.assertTrue(self.store.verify(ref,problem)['verified'])

    def test_305_typed_graph_does_not_promote_claims(self):
        graph=e.ResearchGraph(self.store)
        source=graph.add_node('problem',{'name':'original'})
        target=graph.add_node('candidate',{'name':'trial'})
        for kind in e.RELATION_TYPES:
            edge=graph.add_relation(kind,[source],target,conditions={'assumed':True})
            result=graph.relation_status(edge)
            self.assertFalse(result['verified'])
            self.assertEqual(result['status'],'heuristic' if kind=='heuristic' else 'unverified_no_rule')
        with self.assertRaises(ValueError): graph.add_relation('proved',[source],target)
        with self.assertRaises(ValueError): graph.add_relation('implies',['0'*64],target)
        with self.assertRaises(ValueError): graph.add_node('theorem',{})
        self.assertEqual(len(graph.export()['relations']),6)


    def test_306_cycles_cannot_justify_the_original_claim(self):
        graph=e.ResearchGraph(self.store)
        a,b,c=[graph.add_node('claim',{'name':name}) for name in ('A','B','C')]
        graph.add_relation('implies',[a],b)
        graph.add_relation('upper',[a,b],c)
        before=graph.export()
        for sources,target in (([c],a),([b],a),([a],a),([a,c],b)):
            with self.subTest(sources=sources,target=target):
                with self.assertRaisesRegex(ValueError,'Circular'):
                    graph.add_relation('equivalent',sources,target)
        self.assertEqual(before,graph.export())
        order=graph.topological_order()
        self.assertLess(order.index(a),order.index(b)); self.assertLess(order.index(b),order.index(c))
        with self.assertRaises(ValueError): graph.add_relation('implies',[[a]],b)
        with self.assertRaises(ValueError): graph.add_relation('implies',[a],{})


    def _bound_graph(self):
        import backend
        cert=json.loads((backend.FROZEN/'certificates'/'step_n4.json').read_text())
        ref=self.store.put(cert); problem=e.problem_claim(cert)
        graph=e.ResearchGraph(self.store)
        source=graph.add_node('certificate',problem,evidence=ref)
        target=graph.add_node('claim',{'problem':problem,'value':cert['upper']})
        return graph,cert,problem,source,target

    def test_307_bound_rules_require_conditions_scope_and_source(self):
        graph,cert,p,source,target=self._bound_graph()
        edge=graph.add_relation('upper',[source],target,scope=p,
                                conditions={'mean_zero':True},rule='certified_bound_v1')
        self.assertTrue(graph.relation_status(edge)['verified'])
        bad=graph.add_relation('upper',[source],target,scope=p,
                               conditions={'mean_zero':False},rule='certified_bound_v1')
        self.assertEqual(graph.relation_status(bad)['status'],'condition_not_established')
        bad=graph.add_relation('upper',[source],target,scope={'geometry':'any_manifold'},rule='certified_bound_v1')
        self.assertEqual(graph.relation_status(bad)['status'],'scope_mismatch')
        smaller=graph.add_node('claim',{'problem':p,'value':'0'})
        bad=graph.add_relation('upper',[source],smaller,scope=p,rule='certified_bound_v1')
        self.assertFalse(graph.relation_status(bad)['verified'])
        low=graph.add_node('claim',{'problem':p,'value':'0'})
        good=graph.add_relation('lower',[source],low,scope=p,rule='certified_bound_v1')
        self.assertTrue(graph.relation_status(good)['verified'])
        bad=graph.add_relation('equivalent',[source],target,scope=p,rule='certified_bound_v1')
        self.assertEqual(graph.relation_status(bad)['status'],'rule_kind_mismatch')
        missing=graph.add_relation('implies',[source],target,scope=p,rule='unregistered_theorem')
        self.assertEqual(graph.relation_status(missing)['status'],'unverified_unknown_rule')

    def test_307_approximation_never_becomes_spectral_bound(self):
        import backend
        cert=backend.pieces.piecewise_model({'kind':'axis_profile','profile':'step'})
        p=e.problem_claim(cert); graph=e.ResearchGraph(self.store)
        source=graph.add_node('certificate',p,evidence=self.store.put(cert))
        target=graph.add_node('claim',{'problem':p,'error_upper':'0','approximant_digest':e.digest(cert['cells'])})
        edge=graph.add_relation('approximation',[source],target,scope=p,rule='certified_approximation_v1')
        self.assertTrue(graph.relation_status(edge)['verified'])
        false=graph.add_node('claim',{'problem':p,'value':'0'})
        bad=graph.add_relation('upper',[source],false,scope=p,rule='certified_bound_v1')
        self.assertEqual(graph.relation_status(bad)['status'],'bound_claim_mismatch')
        other=graph.add_node('claim',{'problem':p,'error_upper':'0','approximant_digest':'0'*64})
        bad=graph.add_relation('approximation',[source],other,scope=p,rule='certified_approximation_v1')
        self.assertFalse(graph.relation_status(bad)['verified'])

    def test_307_read_only_replay_tamper_cycles_and_explicit_custom_rules(self):
        graph,cert,p,source,target=self._bound_graph()
        good=graph.add_relation('upper',[source],target,scope=p,rule='certified_bound_v1')
        original=graph.export()
        with patch.object(self.store,'put',side_effect=AssertionError('replay cannot write')):
            recovered=e.ResearchGraph.from_export(self.store,original)
            self.assertTrue(recovered.replay()['all_relations_verified'])
        unknown=graph.add_relation('equivalent',[source],target,scope=p,rule='my_equivalence')
        self.assertFalse(graph.relation_status(unknown)['verified'])
        # Even a trusted extension checker cannot promote a heuristic edge.
        rules={'my_equivalence':{'version':'v1','kinds':['equivalent','heuristic'],'check':lambda *args:1}}
        recovered=e.ResearchGraph.from_export(self.store,graph.export(),rules)
        self.assertFalse(recovered.relation_status(unknown)['verified'])
        heuristic=recovered.add_relation('heuristic',[source],target,scope=p,rule='my_equivalence')
        self.assertEqual(recovered.relation_status(heuristic)['status'],'heuristic')
        forged=e.canonical({'format':e.RELATION_FORMAT,'kind':'implies','sources':[target],
                            'target':source,'conditions':{},'scope':p,'rule':None})
        cycle=self.store.put(json.loads(forged)); bad=dict(original); bad['relations']=original['relations']+[cycle]
        with self.assertRaisesRegex(ValueError,'Circular'): e.ResearchGraph.from_export(self.store,bad)
        bad=dict(original); bad['verified']=True
        with self.assertRaises(ValueError): e.ResearchGraph.from_export(self.store,bad)
        (self.root/(good+'.json')).write_bytes(b'{}')
        with self.assertRaises(ValueError): e.ResearchGraph.from_export(self.store,original)

    def test_307_custom_rule_needs_verified_premise_and_explicit_boolean(self):
        import backend
        graph,cert,p,source,target=self._bound_graph()
        calls=[]
        def check(relation,sources,target,context):
            calls.append(relation['kind'])
            return (backend.pieces.rational(target['statement']['value']) >=
                    backend.pieces.rational(context['certificates'][0]['upper']))
        rules={'upper_implies_v1':{'version':'1','kinds':['implies'],'check':check}}
        graph=e.ResearchGraph.from_export(self.store,graph.export(),rules)
        edge=graph.add_relation('implies',[source],target,scope=p,rule='upper_implies_v1')
        self.assertTrue(graph.relation_status(edge)['verified'])
        fake_source=graph.add_node('candidate',p)
        bad=graph.add_relation('implies',[fake_source],target,scope=p,rule='upper_implies_v1')
        self.assertEqual(graph.relation_status(bad)['status'],'source_not_certified')
        bad=graph.add_relation('implies',[source],target,scope=p,
                               conditions={'unknown_assumption':True},rule='upper_implies_v1')
        self.assertEqual(graph.relation_status(bad)['status'],'condition_not_established')
        self.assertEqual(calls,['implies'])
        mismatch=dict(p); mismatch['mean_zero']=False
        invalid=graph.add_node('certificate',mismatch,evidence=graph.nodes[source]['evidence_reference'])
        bad=graph.add_relation('implies',[invalid],target,scope=p,rule='upper_implies_v1')
        self.assertEqual(graph.relation_status(bad)['status'],'source_problem_mismatch')


if __name__ == '__main__':
    unittest.main()
