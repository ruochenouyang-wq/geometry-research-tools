import json
from pathlib import Path
import unittest
from fractions import Fraction

import interaction as subject


def example_certificate():
    # Synthetic contract fixture: callbacks below test transport, not mathematics.
    return {'format': 'fixture_spectrum', 'function': {'q': '0'},
            'scope': 'fixture_space', 'mean_zero': True,
            'tolerance': '1/10', 'lower': '1', 'upper': '2',
            'exact_width': '1', 'status': 'certified_open'}


class EvidenceTests(unittest.TestCase):
    def test_retrieval_checks_hash_before_verifier(self):
        cert = example_certificate()
        summary = subject.present_result({'certificate': cert}, store=lambda k,v:'saved')
        calls = []
        def check(c):
            calls.append(c)
            return True
        recovered = subject.retrieve_evidence(summary['evidence_ref'], lambda p:cert,
                                               check, require_verified=True)
        self.assertEqual(recovered['certificate'], cert)
        self.assertEqual(recovered['evidence_state'], 'verified')
        changed = dict(cert, upper='3')
        with self.assertRaises(ValueError):
            subject.retrieve_evidence(summary['evidence_ref'], lambda p:changed, check)
        self.assertEqual(len(calls), 1)

    def test_hash_match_is_not_a_proof(self):
        cert = example_certificate()
        ref = {'sha256':subject.digest(cert), 'locator':'saved'}
        self.assertEqual(subject.retrieve_evidence(ref,lambda p:cert)['evidence_state'], 'unchecked')
        with self.assertRaises(ValueError):
            subject.retrieve_evidence(ref,lambda p:cert,require_verified=True)
        with self.assertRaises(ValueError):
            subject.retrieve_evidence(ref,lambda p:cert,lambda c:False,True)
        ref['bytes'] = 0
        with self.assertRaises(ValueError):
            subject.retrieve_evidence(ref,lambda p:cert)

    def test_store_reference_preserves_evidence_and_reduces_transport(self):
        path = Path(__file__).resolve().parent.parent/'geometry_precision_followup'/'approximations'/'same25slots_half.json'
        certificate = json.loads(path.read_text())
        memory = {}
        def store(key, value):
            memory[key] = value
            return 'memory:'+key
        summary = subject.present_result({'certificate': certificate}, store=store)
        self.assertNotIn('certificate', summary)
        ref = summary['evidence_ref']
        self.assertEqual(memory[ref['sha256']], certificate)
        self.assertLess(len(subject.canonical(summary)), len(subject.canonical(certificate)))

    def test_no_store_cannot_silently_drop_evidence(self):
        certificate = example_certificate()
        response = subject.present_result({'certificate': certificate})
        self.assertEqual(response['certificate'], certificate)
        with self.assertRaises(ValueError):
            subject.present_result({'certificate': certificate}, store=lambda k,v: None)


class RequestRunnerTests(unittest.TestCase):
    def test_duplicate_execution_and_copy_isolation(self):
        calls = []
        def execute(request):
            calls.append(request)
            return {'ok': True, 'values':[request['value']]}
        runner = subject.RequestRunner(execute, 'fixture-code-v1')
        first = runner.execute({'value': 3})
        first['result']['values'][0] = 999
        second = runner.execute({'value': 3})
        self.assertEqual(len(calls), 1)
        self.assertTrue(second['cache_hit'])
        self.assertEqual(second['result']['values'], [3])

    def test_full_semantics_and_namespace_bound(self):
        request = {'op':'spectrum','function':{'q':'0'},'mean_zero':True,'tolerance':'1/10'}
        original = subject.request_key(request, 'v1')
        for change in [{'mean_zero':False},{'tolerance':'1/100'},{'function':{'q':'1'}}]:
            self.assertNotEqual(subject.request_key(dict(request,**change),'v1'), original)
        self.assertNotEqual(subject.request_key(request, 'v2'), original)
        self.assertEqual(subject.request_key(dict(reversed(list(request.items()))),'v1'), original)
        with self.assertRaises(ValueError):
            subject.request_key({1:'ambiguous'},'v1')

    def test_failure_not_cached_and_disable_supported(self):
        calls = []
        runner = subject.RequestRunner(lambda r:(calls.append(r) or {'ok':False}), 'v1')
        runner.execute({});runner.execute({})
        self.assertEqual(len(calls),2)
        runner = subject.RequestRunner(lambda r:{'ok':True}, 'v1', cache=False)
        self.assertFalse(runner.execute({})['cache_hit'])
        self.assertFalse(runner.execute({})['cache_hit'])



class SummaryTests(unittest.TestCase):
    def test_directed_rounding_negative_and_tiny_values(self):
        for value in ['1/3', '-1/3', '1/1000000000000000000', '-7/2', '0']:
            for places in [0, 3, 24]:
                low = subject.outward_decimal(value, places)
                high = subject.outward_decimal(value, places, upper=True)
                self.assertLessEqual(Fraction(low), Fraction(value))
                self.assertGreaterEqual(Fraction(high), Fraction(value))
                self.assertLessEqual(Fraction(high)-Fraction(low), Fraction(1, 10**places))
        with self.assertRaises(ValueError):
            subject.outward_decimal(0.1)

    def test_claimed_met_is_recomputed_from_exact_width(self):
        cert = example_certificate()
        cert['status'] = 'certified_met'
        summary = subject.summarize_result({'certificate': cert}, verifier=lambda c: True)
        self.assertEqual(summary['status'], 'certified_open')
        self.assertEqual(summary['bounds']['lower'], '1')

    def test_false_width_and_reversed_bounds_rejected(self):
        cert = example_certificate()
        cert['exact_width'] = '0'
        with self.assertRaises(ValueError):
            subject.summarize_result({'certificate': cert})
        cert.update(lower='3', exact_width='-1')
        with self.assertRaises(ValueError):
            subject.summarize_result({'certificate': cert})

    def test_scope_target_and_unchecked_state_survive(self):
        summary = subject.summarize_result({'certificate': example_certificate()})
        self.assertEqual(summary['scope']['mean_zero'], True)
        self.assertEqual(summary['target']['tolerance'], '1/10')
        self.assertEqual(summary['status'], 'unchecked')

    def test_rejection_cannot_be_reported_as_success(self):
        summary = subject.summarize_result({'certificate': example_certificate()},
                                           verifier=lambda c: False)
        self.assertFalse(summary['ok'])
        self.assertEqual(summary['status'], 'verification_failed')

    def test_literal_verifier_acceptance_required(self):
        summary = subject.summarize_result({'certificate': example_certificate()},
                                           verifier=lambda c: {'verified': False})
        self.assertFalse(summary['ok'])

    def test_invalid_and_failed_result(self):
        with self.assertRaises(ValueError):
            subject.summarize_result({'certificate': {}})
        result = subject.summarize_result({'ok': False, 'error': 'ValueError', 'reason': 'bad'})
        self.assertEqual(result['status'], 'execution_failed')

    def test_existing_piece_certificate_replay_and_no_spectral_claim(self):
        import backend
        path = Path(__file__).resolve().parent.parent/'geometry_precision_followup'/'approximations'/'step_exact.json'
        certificate = json.loads(path.read_text())
        summary = subject.summarize_result({'certificate': certificate}, verifier=backend.pieces.verify_piecewise)
        self.assertEqual(summary['evidence_state'], 'verified')
        self.assertEqual(summary['target']['quantity'], 'function_L2_error')
        self.assertFalse(summary['scope']['spectral_transfer_claimed'])


if __name__ == '__main__':
    unittest.main()
