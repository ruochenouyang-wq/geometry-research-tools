import json
import unittest
from unittest.mock import Mock, patch

import request_contract as contract

STEP = {'kind': 'axis_profile', 'profile': 'step'}
SINGULAR = {'kind': 'axis_profile', 'profile': 'abs_power',
            'exponent': '-1/4', 'amplitude': '-1'}


class ContractTests(unittest.TestCase):
    def test_duplicate_key_cannot_replace_expected_function(self):
        line = '{"op":"verify","function":{"profile":"step"},"function":null}'
        self.assertIsNone(json.loads(line)['function'])
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            contract.parse_request(line)

    def test_nested_duplicate_and_valid_roundtrip(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            contract.parse_request('{"op":"spectrum","function":{"axis":0,"axis":2}}')
        request = {'op': 'spectrum', 'function': STEP}
        self.assertEqual(contract.parse_request(json.dumps(request)), request)

    def test_nonfinite_json_and_overflow_are_rejected(self):
        for value in ('NaN', 'Infinity', '-Infinity', '1e999'):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'Non-finite'):
                contract.parse_request('{"op":"verify","certificate":{"value":' + value + '}}')
        with self.assertRaisesRegex(ValueError, 'Non-finite'):
            contract.validate_request({'op': 'verify', 'certificate': {'x': float('nan')}})

    def test_non_object_roots_are_rejected_before_execution(self):
        for value in ('[]', 'null', '1', '"spectrum"', 'true'):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, 'object'):
                contract.parse_request(value)

    def test_null_cannot_silently_disable_certificate_binding(self):
        cert = contract.backend.legacy.handle({'op': 'spectrum', 'function': STEP,
                                              'modes': 2, 'bits': 16})['certificate']
        for field in ('function', 'mean_zero', 'tolerance'):
            request = {'op': 'verify', 'certificate': cert, field: None}
            self.assertTrue(contract.backend.legacy.handle(request)['verified'])
            with self.subTest(field=field), self.assertRaisesRegex(ValueError, 'null binding'):
                contract.validate_request(request)
        self.assertTrue(contract.execute({'op': 'verify', 'certificate': cert})['verified'])

    def test_valid_approximation_is_not_a_spectral_claim(self):
        cert = contract.backend.legacy.handle({'op': 'approximate', 'function': STEP})['certificate']
        old = {'op': 'verify', 'certificate': cert, 'function': STEP}
        self.assertTrue(contract.backend.legacy.handle(old)['verified'])
        response = contract.execute({**old, 'expected_kind': 'spectrum'})
        self.assertFalse(response['verified'])
        self.assertFalse(response['kind_matches'])
        self.assertTrue(contract.execute({**old, 'expected_kind': 'approximation'})['verified'])

    def test_kind_metadata_is_removed_and_wrong_generation_not_started(self):
        spy = Mock(return_value={'ok': True, 'certificate': {'format': contract.backend.pieces.FORMAT}})
        contract.execute({'op': 'approximate', 'function': STEP, 'expected_kind': 'approximation'}, executor=spy)
        self.assertNotIn('expected_kind', spy.call_args.args[0])
        spy.reset_mock()
        with self.assertRaises(ValueError):
            contract.execute({'op': 'approximate', 'function': STEP, 'expected_kind': 'spectrum'}, executor=spy)
        spy.assert_not_called()

    def test_invalid_singular_budgets_do_not_compute_source(self):
        for field, value in (('max_terms', 0), ('tolerance', '0'),
                             ('precision_bits', 31), ('iterations', 65)):
            request = {'op': 'precise_singular', 'function': SINGULAR, field: value}
            with patch.object(contract.backend.direct, 'full_ground', wraps=contract.backend.direct.full_ground) as old_source:
                with self.assertRaises(ValueError):
                    contract.backend.legacy.handle(request)
                self.assertEqual(old_source.call_count, 1)
            with patch.object(contract.backend.direct, 'full_ground') as new_source:
                with self.subTest(field=field), self.assertRaises(ValueError):
                    contract.execute(request)
                new_source.assert_not_called()

    def test_zero_iteration_is_a_valid_budget_not_a_false_error(self):
        request = {'op': 'precise_singular', 'function': SINGULAR,
                   'max_terms': 1, 'iterations': 0}
        self.assertEqual(contract.validate_request(request), request)
        cert = contract.execute(request)['certificate']
        self.assertTrue(contract.backend.enriched.verify(cert))

    def test_all_route_budgets_and_types_fail_before_executor(self):
        cases = [('spectrum', 'modes', True), ('spectrum', 'bits', 97),
                 ('spectrum', 'max_m', -1), ('spectrum', 'sqrt_bits', 7),
                 ('spectrum', 'near_tail', 33), ('spectrum', 'mean_zero', 1),
                 ('spectrum', 'tolerance', 0.01), ('approximate', 'levels', 0),
                 ('approximate', 'degree', 33), ('approximate', 'sqrt_bits', 257),
                 ('approximate', 'root_ratio', '1'), ('verify', 'expected_kind', [])]
        for op, field, value in cases:
            spy = Mock()
            with self.subTest(op=op, field=field), self.assertRaises(ValueError):
                contract.execute({'op': op, 'function': STEP, field: value}, executor=spy)
            spy.assert_not_called()
        for request in ({'op': []}, {'op': 'spectrum'}, {'op': 'spectrum', 'function': STEP, 'typo': 3}):
            with self.assertRaises(ValueError):
                contract.validate_request(request)

    def test_bad_spectral_tolerance_prevents_all_sector_work(self):
        request = {'op': 'spectrum', 'function': STEP, 'modes': 1,
                   'bits': 8, 'max_m': 0, 'tolerance': '0'}
        with patch.object(contract.backend.direct, 'certify_sector', wraps=contract.backend.direct.certify_sector) as old:
            with self.assertRaises(ValueError):
                contract.backend.legacy.handle(request)
            self.assertEqual(old.call_count, 1)
        with patch.object(contract.backend.direct, 'certify_sector') as new:
            with self.assertRaises(ValueError):
                contract.execute(request)
            new.assert_not_called()

    def test_valid_boundary_budgets_remain_usable(self):
        cert = contract.execute({'op': 'approximate', 'function': STEP,
                                 'levels': 1, 'degree': 0, 'sqrt_bits': 0})['certificate']
        self.assertTrue(contract.backend.pieces.verify_piecewise(cert))

    def test_semantic_digest_normalizes_defaults_and_exact_rationals(self):
        a = {'op': 'spectrum', 'function': STEP, 'request_id': 'a'}
        b = {'op': 'spectrum', 'function': {**STEP, 'axis': 2, 'amplitude': '2/2', 'offset': 0},
             **contract.DEFAULTS['spectrum'], 'request_id': 'b'}
        self.assertEqual(contract.request_digest(a), contract.request_digest(b))
        self.assertNotEqual(contract.request_digest(a), contract.request_digest({**a, 'tolerance': '1/100'}))

    def test_request_identity_and_expected_digest_are_bound_before_dispatch(self):
        request = {'op': 'approximate', 'function': STEP, 'request_id': 'run-2'}
        digest = contract.request_digest(request)
        spy = Mock(wraps=contract.backend.legacy.handle)
        response = contract.execute({**request, 'expected_request_digest': digest}, executor=spy)
        self.assertEqual(response['request_id'], 'run-2')
        self.assertEqual(response['request_digest'], digest)
        self.assertFalse(set(spy.call_args.args[0]) & contract.METADATA)
        spy.reset_mock()
        with self.assertRaisesRegex(ValueError, 'digest mismatch'):
            contract.execute({**request, 'expected_request_digest': '0' * 64}, executor=spy)
        spy.assert_not_called()

    def test_executor_cannot_mutate_the_callers_original_request(self):
        request = {'op': 'approximate', 'function': dict(STEP)}
        def executor(payload):
            payload['function']['amplitude'] = '-1'
            return {'ok': True}
        contract.execute(request, executor=executor)
        self.assertEqual(request['function'], STEP)

    def test_valid_but_unfinished_spectrum_is_not_target_success(self):
        request = {'op': 'spectrum', 'function': STEP, 'modes': 2, 'bits': 16,
                   'tolerance': '1/100000000'}
        old = contract.backend.legacy.handle(request)
        self.assertTrue(old['ok'])
        self.assertEqual(old['certificate']['status'], 'certified_open')
        response = contract.execute(request)
        self.assertTrue(response['execution_ok'])
        self.assertTrue(response['certificate_valid'])
        self.assertFalse(response['target_met'])
        self.assertEqual(response['status'], 'certified_open')

    def test_solver_claim_cannot_bypass_actual_certificate_replay(self):
        cert = contract.backend.legacy.handle({'op': 'approximate', 'function': STEP})['certificate']
        cert['error_squared'] = '999'
        response = contract.execute({'op': 'approximate', 'function': STEP},
                                    executor=lambda _: {'ok': True, 'verified': True, 'status': 'target_met', 'certificate': cert})
        self.assertFalse(response['certificate_valid'])
        self.assertFalse(response['verified'])
        self.assertIsNone(response['target_met'])

    def test_replay_marks_missing_external_bindings(self):
        cert = contract.backend.legacy.handle({'op': 'spectrum', 'function': STEP, 'modes': 2, 'bits': 16})['certificate']
        forbidden = Mock(side_effect=AssertionError('verification must not solve'))
        response = contract.execute({'op': 'verify', 'certificate': cert}, executor=forbidden)
        forbidden.assert_not_called()
        self.assertTrue(response['certificate_valid'])
        self.assertTrue(response['verified'])
        self.assertFalse(response['binding_complete'])
        self.assertIsNone(response['target_met'])
        self.assertEqual(response['status'], 'verified_unbound')

    def test_wrong_original_keeps_validity_separate_from_binding(self):
        cert = contract.backend.legacy.handle({'op': 'approximate', 'function': STEP})['certificate']
        response = contract.execute({'op': 'verify', 'certificate': cert,
                                     'function': {**STEP, 'amplitude': '-1'}, 'expected_kind': 'approximation'})
        self.assertTrue(response['certificate_valid'])
        self.assertFalse(response['bindings_match'])
        self.assertFalse(response['verified'])

    def test_custom_verifier_receives_original_bindings(self):
        request = {'op': 'approximate', 'function': STEP}
        verifier = Mock(wraps=contract.verify_certificate)
        response = contract.execute(request, verifier=verifier)
        self.assertTrue(response['certificate_valid'])
        self.assertEqual(verifier.call_count, 1)
        self.assertEqual(verifier.call_args.kwargs['expected_function'], STEP)

    def test_stable_rejection_codes_do_not_depend_on_error_messages(self):
        cases = [('[]', 'invalid_root'), ('{', 'malformed_json'),
                 ('{"op":"verify","op":"spectrum"}', 'duplicate_key'),
                 ('{"op":"verify","mean_zero":null}', 'null_binding'),
                 ('{"op":"unknown"}', 'unknown_operation'),
                 ('{"op":"verify","tolerance":NaN}', 'nonfinite_number')]
        for line, code in cases:
            with self.subTest(code=code), self.assertRaises(contract.ContractError) as caught:
                contract.parse_request(line)
            self.assertEqual(caught.exception.code, code)

    def test_invalid_evidence_and_wrong_goal_have_distinct_reasons(self):
        cert = contract.backend.legacy.handle({'op': 'approximate', 'function': STEP})['certificate']
        wrong_goal = contract.execute({'op': 'verify', 'certificate': cert, 'expected_kind': 'spectrum'})
        self.assertEqual(wrong_goal['reason_code'], 'conclusion_kind_mismatch')
        self.assertTrue(wrong_goal['certificate_valid'])
        bad = {**cert, 'error_squared': '1'}
        self.assertEqual(contract.execute({'op': 'verify', 'certificate': bad})['reason_code'], 'invalid_evidence')
        self.assertEqual(contract.execute({'op': 'verify'})['reason_code'], 'missing_certificate')
        self.assertEqual(contract.execute({'op': 'verify', 'certificate': {'format': 'future'}})['reason_code'], 'unsupported_certificate_format')

    def test_request_stream_survives_parse_and_executor_failure(self):
        calls = []
        def executor(request):
            calls.append(request)
            if len(calls) == 1:
                raise RuntimeError('Injected worker failure')
            return contract.backend.legacy.handle(request)
        lines = ['{', json.dumps({'op': 'approximate', 'function': STEP, 'request_id': 'failed'}),
                 '\n', json.dumps({'op': 'approximate', 'function': STEP, 'request_id': 'succeeded'})]
        responses = list(contract.process_lines(lines, executor=executor))
        self.assertEqual([r['line_number'] for r in responses], [1, 2, 4])
        self.assertEqual(responses[0]['reason_code'], 'malformed_json')
        self.assertEqual(responses[1]['reason_code'], 'internal_error')
        self.assertEqual(responses[1]['stage'], 'executor')
        self.assertEqual(responses[1]['request_id'], 'failed')
        self.assertIsNone(responses[1]['certificate_valid'])
        self.assertIsNone(responses[1]['target_met'])
        self.assertTrue(responses[2]['certificate_valid'])

    def test_verifier_failure_is_not_a_bad_mathematical_certificate(self):
        cert = contract.backend.legacy.handle({'op': 'approximate', 'function': STEP})['certificate']
        calls = []
        def verifier(*args, **kwargs):
            calls.append(1)
            if len(calls) == 1:
                raise RuntimeError('Injected checker failure')
            return contract.verify_certificate(*args, **kwargs)
        line = json.dumps({'op': 'verify', 'certificate': cert})
        responses = list(contract.process_lines([line, line], verifier=verifier))
        self.assertEqual(responses[0]['stage'], 'verifier')
        self.assertIsNone(responses[0]['certificate_valid'])
        self.assertTrue(responses[1]['verified'])

    def test_explicit_backend_failure_and_resource_failure_are_separate(self):
        request = {'op': 'approximate', 'function': STEP}
        response = contract.execute(request, executor=lambda _: {'ok': False, 'reason': 'worker unavailable'})
        self.assertEqual(response['reason_code'], 'backend_failure')
        self.assertIsNone(response['certificate_valid'])
        response = contract.execute(request, executor=Mock(side_effect=TimeoutError('budget ended')))
        self.assertEqual(response['reason_code'], 'resource_limit')
        self.assertIsNone(response['target_met'])

    def test_user_interrupt_is_not_swallowed(self):
        with self.assertRaises(KeyboardInterrupt):
            list(contract.process_lines([json.dumps({'op': 'approximate', 'function': STEP})],
                                        executor=Mock(side_effect=KeyboardInterrupt())))


if __name__ == '__main__':
    unittest.main()
