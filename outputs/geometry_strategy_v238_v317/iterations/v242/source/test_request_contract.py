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


if __name__ == '__main__':
    unittest.main()
