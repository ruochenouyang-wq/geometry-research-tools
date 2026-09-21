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


if __name__ == '__main__':
    unittest.main()
