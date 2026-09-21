"""Independent bounded audit of the final research registry and JSONL facade."""
import copy
from fractions import Fraction as F
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import exact_extrema
import matrix_refine
import original_goal
import positive_search
import precision_bridge
import protocol
import research_registry as registry
import run
import token_meter as meter
import uniform_refine


ROOT = Path(__file__).resolve().parent
RAW = {'profile': 'sphere_ground', 'parameters': ['a'], 'potential': 'a*t',
       'penalty': 'a*a/5', 'domain': {'kind': 'all_real'}, 'threshold': '0'}


class IntegrationAudit(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        backend = precision_bridge.ExactRangeBackend()
        lower = positive_search.lower_certificate([0, 0, '20/9'], [[0, 0, 1]], ['-1/3'], range_backend=backend)
        upper = positive_search.ansatz_upper_certificate([0, 0, '20/9'], [[0, 0, 1]], [{'t': '1/2', 'weight': '1'}])
        cls.positive = positive_search.optimization_certificate(lower, upper,
                          rayleigh=positive_search.rayleigh_certificate([0, 0, '20/9'], [1]), range_backend=backend)
        cls.extrema = exact_extrema.maximize([0, 1, 0, -1], interval=[-1, 2], max_nodes=4)
        cls.family = precision_bridge.precise_family([[0, 1]], [[1]], max_nodes=0)
        cls.matrix = precision_bridge.precise_matrix(matrix_refine.refine([[0, 1]], budget=0), max_nodes=0)
        cls.classical = uniform_refine.synthesize_log_sobolev()
        cls.counterexample = uniform_refine.synthesize_log_sobolev(coefficient='1/7')
        cls.obstruction = uniform_refine.synthesize(coefficient='1/5', order=1)
        cls.proved = original_goal.solve(protocol.expand_request(RAW))
        cls.refuted = original_goal.solve(protocol.expand_request(dict(RAW, threshold='1')))
        cls.unresolved = original_goal.solve(protocol.expand_request(dict(RAW, penalty='a*a/10', threshold='-1/10')))
        cls.certificates = [cls.positive, cls.extrema, cls.family, cls.matrix, cls.classical,
                            cls.counterexample, cls.obstruction, cls.proved, cls.refuted, cls.unresolved]

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='geometry-integration-audit-')
        self.addCleanup(self.tmp.cleanup)
        self.service = protocol.Service(self.tmp.name)

    def serve(self, requests):
        result = subprocess.run([sys.executable, str(ROOT/'run.py'), 'serve', '--store', self.tmp.name],
                                input='\n'.join(meter.wire(r) for r in requests)+'\n',
                                text=True, capture_output=True, timeout=20, check=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        lines = result.stdout.splitlines(keepends=True)
        self.assertEqual(len(lines), len(requests))
        return [(line, json.loads(line)) for line in lines]

    def test_every_new_decision_rounds_each_bound_in_the_safe_direction(self):
        for cert in self.certificates:
            exact = protocol.summary(cert)
            shown = protocol.summary(cert, compact_numbers=True)
            self.assertEqual(exact['status'], shown['status'])
            for key, value in exact['bounds'].items():
                with self.subTest(method=cert.get('method', cert.get('format')), key=key):
                    if value is None:
                        self.assertIsNone(shown['bounds'][key])
                    elif key == 'lower' or key.endswith('_lower'):
                        self.assertLessEqual(F(shown['bounds'][key]), F(value))
                    else:
                        self.assertGreaterEqual(F(shown['bounds'][key]), F(value))

    def test_positive_summary_separates_closed_ansatz_and_open_spectral_gap(self):
        shown = protocol.summary(self.positive)
        self.assertEqual(shown['status'], 'ansatz_gap_closed')
        self.assertEqual(F(shown['bounds']['ansatz_gap']), 0)
        self.assertGreater(F(shown['bounds']['spectral_gap']), 0)
        self.assertFalse(shown['scope']['ansatz_upper_is_spectral_upper'])
        self.assertEqual(shown['scope']['space'], 'full_sphere')

    def test_extrema_keeps_original_interval_without_spectral_scope(self):
        shown = protocol.summary(self.extrema)
        self.assertEqual(shown['scope']['interval'], ['-1', '2'])
        self.assertEqual(shown['task']['interval'], ['-1', '2'])
        self.assertEqual(shown['kind'], 'polynomial_maximum')
        self.assertNotIn('geometry', shown['scope'])

    def test_matrix_and_family_scopes_retain_poisson_limitation(self):
        for cert in (self.matrix, self.family):
            shown = protocol.summary(cert)
            self.assertEqual(shown['scope']['ansatz'], 'poisson_exponential')
            self.assertEqual(shown['scope']['amplitudes'], 'all_real')
            self.assertFalse(shown['formal_kernel'])
        self.assertEqual(protocol.summary(self.matrix)['scope']['objective'], 'trace(C G)')

    def test_classical_analytic_dependency_survives_budget_and_restore(self):
        for cert in (self.classical, self.proved, self.refuted, self.unresolved):
            receipt = self.service.receipt(cert)
            self.assertIn('not proved by this program', receipt['scope']['analytic_dependency'])
            self.assertFalse(receipt['formal_kernel'])
            packed = self.service.budget(receipt, 100000)
            self.assertEqual(json.loads(packed['text'])['scope'], receipt['scope'])
            restored = self.service.restore(self.service.checkpoint([], [receipt]))
            self.assertEqual(restored['results'][0]['scope'], receipt['scope'])

    def test_classical_disclosure_cannot_be_removed_from_verified_receipt(self):
        receipt = self.service.receipt(self.classical)
        receipt['scope'].pop('analytic_dependency')
        with self.assertRaises(ValueError): self.service.budget(receipt)
        with self.assertRaises(ValueError): self.service.poll(receipt, receipt['evidence'])

    def test_obstruction_is_not_a_spectral_counterexample(self):
        obstruction = protocol.summary(self.obstruction)
        refutation = protocol.summary(self.counterexample)
        self.assertEqual(obstruction['status'], 'ansatz_obstruction')
        self.assertEqual(obstruction['bounds'], {})
        self.assertNotIn('analytic_dependency', obstruction['scope'])
        self.assertEqual(refutation['status'], 'disproved')
        self.assertLess(F(refutation['bounds']['counterexample_objective_upper']), 0)

    def test_original_goal_all_three_statuses_keep_true_objective_and_source(self):
        for status, cert in (('proved', self.proved), ('refuted', self.refuted), ('unresolved', self.unresolved)):
            receipt = self.service.receipt(cert)
            self.assertEqual(receipt['status'], status)
            self.assertEqual(receipt['scope']['objective'], 'inf_parameters(lambda_1(-Delta+q)+penalty)')
            self.assertEqual(self.service.inspect(receipt['task']), cert['compiled_goal']['source'])
            self.assertTrue(registry.verify(self.service.inspect(receipt['evidence'])['certificate']))
        self.assertIsNone(protocol.summary(self.unresolved)['bounds']['lower'])

    def test_backend_registry_is_allowlisted_and_does_not_trust_json_names(self):
        for name in ('evil.module', '__import__', '../precision_bridge', {'name': 'exact_sturm_minimum_v1'}):
            with self.subTest(name=name), self.assertRaises(ValueError): registry.backend(name)
        bad = copy.deepcopy(self.positive)
        bad['lower_certificate']['range_backend'] = 'evil.module'
        self.assertFalse(registry.verify(bad))

    def test_exact_backend_rejects_valid_range_for_wrong_domain(self):
        bad = copy.deepcopy(self.positive)
        lower = bad['lower_certificate']
        lower['range_proof'] = exact_extrema.maximize(lower['range_proof']['polynomial'], interval=[0, 1], max_nodes=4)
        self.assertTrue(exact_extrema.verify(lower['range_proof']))
        self.assertFalse(registry.verify(bad))

    def test_precise_family_cannot_enter_legacy_library_by_format_confusion(self):
        self.assertTrue(registry.verify(self.family))
        with self.assertRaises(ValueError): self.service.add_lemma(self.family)
        self.assertEqual(self.service.lemma_index(), [])

    def test_all_new_tasks_restore_with_exact_bound_evidence(self):
        receipts = [self.service.receipt(c) for c in self.certificates]
        goals = [r['task'] for r in receipts]
        checkpoint = self.service.checkpoint(goals, receipts)
        fresh = protocol.Service(self.tmp.name)
        restored = fresh.restore(checkpoint)
        self.assertEqual(restored['goals'], goals)
        self.assertEqual(restored['results'], receipts)
        self.assertEqual(fresh.calls, 0)

    def test_context_uses_per_task_scope_for_extrema_and_mixed_research(self):
        for certificates in ([self.extrema], [self.extrema, self.proved, self.positive]):
            receipts = [self.service.receipt(c) for c in certificates]
            context = self.service.context(receipts, [])
            self.assertEqual(context['scope'], 'per_task')
            self.assertEqual(context['available_input_profile'], protocol.PROFILE)
            self.assertNotIn('profile', context)
            restored = self.service.restore(context['checkpoint'])
            self.assertEqual(restored['results'], receipts)
            self.assertNotIn('geometry', receipts[0]['scope'])

    def test_each_research_task_rejects_different_exact_content_binding(self):
        for cert in self.certificates:
            original = protocol.summary(cert)['task']
            wrong = copy.deepcopy(original)
            wrong['audit_wrong_task'] = True
            ref = self.service.store.put('goal', wrong)
            with self.assertRaises(ValueError): self.service.receipt(cert, ref)

    def test_new_evidence_tamper_is_rejected_before_projection(self):
        receipt = self.service.receipt(self.positive)
        path = Path(self.tmp.name)/(receipt['evidence'].replace(':', '_')+'.json')
        path.write_text('{}')
        with self.assertRaises(ValueError): self.service.inspect(receipt['evidence'], '/certificate/status')

    def test_complete_budget_envelope_counts_both_encodings_and_exact_boundary(self):
        for encoding in meter.ENCODINGS:
            receipt = self.service.receipt(self.classical)
            large = self.service.budget(receipt, 100000, encoding)
            self.assertEqual(large['wire_tokens'], meter.count(meter.wire(large), encoding))
            self.assertEqual(large['text_tokens'], meter.count(large['text'], encoding))
            boundary = self.service.budget(receipt, large['wire_tokens'], encoding)
            self.assertTrue(boundary['fits'])
            self.assertEqual(json.loads(boundary['text']), receipt)
            self.assertGreater(large['wire_tokens'], large['text_tokens'])

    def test_insufficient_budget_reports_smallest_complete_available_envelope(self):
        receipt = self.service.receipt(self.classical)
        concise = copy.deepcopy(receipt); concise['diagnostics'] = 'available_in_evidence'
        for encoding in meter.ENCODINGS:
            full = self.service.budget(receipt, 100000, encoding)
            short = self.service.budget(concise, 100000, encoding)
            minimum = min(full['wire_tokens'], short['wire_tokens'])
            result = self.service.budget(receipt, minimum-1, encoding)
            self.assertFalse(result['fits'])
            self.assertEqual(result['required_wire_tokens'], minimum)
            self.assertNotIn('text', result)
            self.assertNotIn('status', result)

    def test_jsonl_count_includes_complete_emitted_result_and_preserves_disclosure(self):
        receipt = self.service.receipt(self.classical)
        requests = [{'op': 'budget', 'receipt': receipt, 'limit': 100000, 'encoding': name}
                    for name in meter.ENCODINGS]
        for encoding, (line, result) in zip(meter.ENCODINGS, self.serve(requests)):
            self.assertEqual(meter.count(line, encoding), result['wire_tokens'])
            self.assertIn('analytic_dependency', json.loads(result['text'])['scope'])

    def test_jsonl_null_reference_errors_do_not_abort_later_requests(self):
        requests = [{'op': 'solve', 'goal_ref': None}, {'op': 'restore', 'ref': None},
                    {'op': 'solve', 'goal_ref': 123}, {'op': 'extrema', 'spec': {'polynomial': [0, 1]}}]
        rows = self.serve(requests)
        for _, row in rows[:3]:
            self.assertIn('error', row); self.assertIsNone(row['mathematical_verdict'])
        self.assertTrue(rows[-1][1]['fits'])

    def test_five_new_operations_run_in_jsonl_and_each_complete_packet_fits(self):
        requests = [
            {'op': 'extrema', 'spec': {'polynomial': [0, 1, 0, -1]}, 'max_nodes': 4},
            {'op': 'positive', 'spec': {'potential': [0, 0, '20/9']},
             'options': {'degree': 2, 'grid_size': 9, 'max_exchanges': 1, 'max_newton': 1}},
            {'op': 'uniform', 'spec': {'coefficient': '1/7'}},
            {'op': 'refine', 'spec': {'directions': [[0, 1]]}, 'budget': 0, 'max_nodes': 0},
            {'op': 'research', 'goal': RAW}]
        for line, response in self.serve(requests):
            self.assertTrue(response['fits'])
            self.assertEqual(meter.count(line), response['wire_tokens'])
            self.assertLessEqual(response['wire_tokens'], 500)
            receipt = json.loads(response['text'])
            self.assertEqual(self.service.poll(receipt, receipt['evidence']), {'unchanged': receipt['evidence']})

    def test_jsonl_rejected_new_inputs_have_no_verdict_and_recover(self):
        requests = [
            {'op': 'positive', 'spec': {'potential': [0]}, 'options': {'range_backend': 'evil.module'}},
            {'op': 'uniform', 'spec': {'direction': [0, 0, 1]}},
            {'op': 'refine', 'spec': {'directions': [[0, 1]], 'cost': [[0]]}, 'budget': 0},
            {'op': 'research', 'goal': dict(RAW, geometry='plane')},
            {'op': 'research', 'goal': RAW}]
        rows = self.serve(requests)
        for _, response in rows[:-1]: self.assertIsNone(response['mathematical_verdict'])
        self.assertTrue(rows[-1][1]['fits'])


if __name__ == '__main__':
    unittest.main()
