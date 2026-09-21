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


if __name__ == '__main__':
    unittest.main()
