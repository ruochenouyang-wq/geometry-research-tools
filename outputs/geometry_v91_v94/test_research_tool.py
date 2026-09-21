"""Small real-kernel JSONL tests with isolated content-addressed evidence."""
from fractions import Fraction as F
import copy
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import projected_spectrum as spectrum
import research_tool as tool


SMALL = {'modes': 2, 'max_modes': 2, 'bits': 12, 'max_m': 1}


class ResearchToolTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='research-tool-audit-')
        self.addCleanup(self.tmp.cleanup)
        self.service = tool.Service(self.tmp.name)

    def request_lines(self, requests):
        output = io.StringIO()
        tool.serve(io.StringIO('\n'.join(tool.wire(r) for r in requests)+'\n'), output, self.service)
        return [json.loads(line) for line in output.getvalue().splitlines()]

    def test_sector_default_projection_and_exact_zero_potential(self):
        result = self.service.call({'op': 'sector', 'q': [0], 'modes': 2, 'bits': 12})
        self.assertEqual(result['status'], 'certified_sector_bound')
        self.assertEqual(result['lower'], '2'); self.assertEqual(result['upper'], '2')
        self.assertEqual(result['scope']['space'], 'single_azimuth_sector_only')
        self.assertTrue(result['scope']['mean_zero'])
        self.assertEqual(result['scope']['projection'], 'remove_l0')

    def test_sector_keeps_m_index_and_unprojected_scope(self):
        result = self.service.call({'op': 'sector', 'q': [3], 'm': 2, 'mean_zero': False,
                                    'k': 2, 'modes': 2, 'bits': 12})
        self.assertEqual(result['lower'], '15')
        self.assertEqual(result['scope']['m'], 2)
        self.assertEqual(result['scope']['k'], 2)
        self.assertFalse(result['scope']['mean_zero'])

    def test_mean_zero_full_sphere_is_not_a_sector(self):
        result = self.service.call({'op': 'mean_zero_ground', 'q': [0], **SMALL})
        self.assertEqual(result['lower'], '2'); self.assertEqual(result['width'], '0')
        self.assertEqual(result['scope']['space'], 'all_mean_zero_H1_functions_on_unit_S2')
        self.assertNotIn('m', result['scope'])

    def test_weighted_poincare_proved_and_existence_refutation_are_distinct(self):
        for threshold, status in ((2, 'proved'), (3, 'refuted_by_spectral_existence')):
            result = self.service.call({'op': 'weighted_poincare', 'q': [0], 'threshold': threshold, **SMALL})
            self.assertEqual(result['status'], status)
            self.assertEqual(result['scope']['threshold'], str(threshold))
            self.assertEqual(result['refutation_does_not_include_an_explicit_function'], threshold == 3)

    def test_weighted_poincare_preserves_undetermined(self):
        options = {'modes': 1, 'max_modes': 1, 'bits': 8, 'max_m': 0}
        result = self.service.call({'op': 'mean_zero_ground', 'q': [0, 1], **options})
        self.assertGreater(F(result['width']), 0)
        threshold = (F(result['lower'])+F(result['upper']))/2
        inequality = self.service.call({'op': 'weighted_poincare', 'q': [0, 1], 'threshold': str(threshold), **options})
        self.assertEqual(inequality['status'], 'undetermined')
        self.assertFalse(inequality['refutation_does_not_include_an_explicit_function'])

    def test_quadratic_potential_threshold_examples_cover_all_mean_zero_functions(self):
        requests = [{'op': 'weighted_poincare', 'q': [0, 0, 2], 'threshold': threshold}
                    for threshold in ('238265540615/100000000000', '12/5')]
        results = self.request_lines(requests)
        self.assertEqual([r['status'] for r in results], ['proved', 'refuted_by_spectral_existence'])
        for result in results:
            self.assertEqual(result['scope']['space'], 'all_mean_zero_H1_functions_on_unit_S2')
            self.assertEqual(result['scope']['quantifier'], 'every_real_H1_function_u_with_integral_u_zero')
        self.assertTrue(results[1]['refutation_does_not_include_an_explicit_function'])

    def test_two_jsonl_requests_can_verify_and_fetch_previous_evidence(self):
        first = self.request_lines([{'op': 'sector', 'q': [0], 'modes': 2}])[0]
        identifier = first['certificate_id']
        second = self.request_lines([{'op': 'verify', 'certificate_id': identifier},
                                     {'op': 'fetch', 'certificate_id': identifier}])
        self.assertTrue(second[0]['verified']); self.assertTrue(second[1]['verified'])
        self.assertEqual(second[0]['scope'], first['scope'])
        self.assertTrue(spectrum.verify_sector(second[1]['certificate']))
        fresh = tool.Service(self.tmp.name)
        self.assertEqual(fresh.call({'op': 'verify', 'certificate_id': identifier}), second[0])

    def test_duplicate_identical_certificates_have_full_stable_content_hash(self):
        request = {'op': 'sector', 'q': [0], 'modes': 2}
        left = self.service.call(request); right = self.service.call(request)
        self.assertEqual(left, right)
        self.assertRegex(left['certificate_id'], r'^[0-9a-f]{64}$')
        data = (Path(self.tmp.name)/(left['certificate_id']+'.json')).read_bytes()
        self.assertEqual(hashlib.sha256(data).hexdigest(), left['certificate_id'])

    def test_content_hash_tamper_rejected_for_fetch_and_verify(self):
        result = self.service.call({'op': 'sector', 'q': [0], 'modes': 2})
        identifier = result['certificate_id']
        (Path(self.tmp.name)/(identifier+'.json')).write_text('{}')
        for op in ('fetch', 'verify'):
            with self.assertRaisesRegex(ValueError, 'hash mismatch'):
                self.service.call({'op': op, 'certificate_id': identifier})

    def test_rehashed_but_invalid_certificate_still_rejected(self):
        certificate = spectrum.certify_sector([0], modes=2)
        certificate['scope'] = 'all_mean_zero_H1_functions_on_unit_S2'
        data = tool.wire(certificate).encode(); identifier = hashlib.sha256(data).hexdigest()
        (Path(self.tmp.name)/(identifier+'.json')).write_bytes(data)
        with self.assertRaisesRegex(ValueError, 'proof verification failed'):
            self.service.call({'op': 'fetch', 'certificate_id': identifier})

    def test_invalid_proof_cannot_be_inserted_in_store(self):
        certificate = spectrum.certify_sector([0], modes=2)
        certificate['lower'] = '1000'
        with self.assertRaises(ValueError): self.service.store.put(certificate)

    def test_certificate_id_cannot_read_arbitrary_paths_or_symlinks(self):
        for identifier in ('../projected_spectrum.py', '/etc/passwd', 'a'*63, 'A'*64, None, 3):
            with self.subTest(identifier=identifier), self.assertRaises(ValueError):
                self.service.call({'op': 'fetch', 'certificate_id': identifier})
        link = Path(self.tmp.name)/('0'*64+'.json'); link.symlink_to(Path(__file__))
        with self.assertRaises(ValueError): self.service.call({'op': 'fetch', 'certificate_id': '0'*64})

    def test_unknown_format_is_not_verified(self):
        self.assertFalse(tool.verify_certificate({'format': 'pretend_full_scope'}))

    def test_full_operations_cannot_override_mean_zero_even_with_true(self):
        for operation in ('mean_zero_ground', 'weighted_poincare'):
            for value in (False, True, 0, 'false'):
                request = {'op': operation, 'q': [0], 'mean_zero': value}
                if operation == 'weighted_poincare': request['threshold'] = 2
                with self.assertRaisesRegex(ValueError, 'Unknown fields'): self.service.call(request)

    def test_client_cannot_inject_scope_or_certificate(self):
        for field in ('scope', 'certificate', 'range_proof', 'path', 'options'):
            with self.assertRaises(ValueError): self.service.call({'op': 'sector', 'q': [0], field: {}})

    def test_generated_sector_must_bind_requested_m_q_mean_zero_and_k(self):
        wrong = spectrum.certify_sector([0], m=1, modes=2)
        with patch.object(tool.spectrum, 'certify_sector', return_value=wrong), self.assertRaises(ValueError):
            self.service.call({'op': 'sector', 'q': [0], 'm': 0, 'modes': 2})
        for key, value in (('q', [1]), ('mean_zero', False), ('k', 2)):
            actual = spectrum.certify_sector([0], modes=2)
            request = {'op': 'sector', 'q': [0], 'modes': 2, key: value}
            with patch.object(tool.spectrum, 'certify_sector', return_value=actual), self.assertRaises(ValueError):
                self.service.call(request)

    def test_generated_full_operation_cannot_accept_sector_certificate(self):
        wrong = spectrum.certify_sector([0], modes=2)
        with patch.object(tool.spectrum, 'full_ground', return_value=wrong), self.assertRaises(ValueError):
            self.service.call({'op': 'mean_zero_ground', 'q': [0], **SMALL})

    def test_generated_inequality_must_bind_threshold(self):
        wrong = spectrum.weighted_poincare([0], 2, **SMALL)
        with patch.object(tool.spectrum, 'weighted_poincare', return_value=wrong), self.assertRaises(ValueError):
            self.service.call({'op': 'weighted_poincare', 'q': [0], 'threshold': 3, **SMALL})

    def test_generated_sector_must_bind_requested_mode_count(self):
        other_size = spectrum.certify_sector([0], modes=3)
        self.assertTrue(spectrum.verify_sector(other_size))
        with patch.object(tool.spectrum, 'certify_sector', return_value=other_size), self.assertRaises(ValueError):
            self.service.call({'op': 'sector', 'q': [0], 'modes': 2})

    def test_full_tolerance_cannot_be_replaced_by_a_wider_valid_certificate(self):
        loose = spectrum.full_ground([0, 0, 2], tolerance=F(1, 100), **SMALL)
        self.assertTrue(spectrum.verify_full(loose))
        self.assertEqual(loose['status'], 'certified_target_met')
        self.assertGreater(F(loose['exact_width']), F(1, 10**10))
        with patch.object(tool.spectrum, 'full_ground', return_value=loose), self.assertRaises(ValueError):
            self.service.call({'op': 'mean_zero_ground', 'q': [0, 0, 2],
                               'tolerance': '1/10000000000', **SMALL})

    def test_inequality_tolerance_must_bind_nested_full_certificate(self):
        loose = spectrum.weighted_poincare([0, 0, 2], 2, tolerance=F(1, 100), **SMALL)
        self.assertTrue(spectrum.verify_inequality(loose))
        with patch.object(tool.spectrum, 'weighted_poincare', return_value=loose), self.assertRaises(ValueError):
            self.service.call({'op': 'weighted_poincare', 'q': [0, 0, 2], 'threshold': 2,
                               'tolerance': '1/10000000000', **SMALL})

    def test_full_radial_mode_counts_stay_within_requested_interval(self):
        for size in (1, 3):
            full = spectrum.full_ground([0], modes=size, max_modes=size, bits=12, max_m=1)
            self.assertTrue(spectrum.verify_full(full))
            with patch.object(tool.spectrum, 'full_ground', return_value=full), self.assertRaises(ValueError):
                self.service.call({'op': 'mean_zero_ground', 'q': [0], **SMALL})
            inequality = spectrum.inequality_certificate(full, 2)
            with patch.object(tool.spectrum, 'weighted_poincare', return_value=inequality), self.assertRaises(ValueError):
                self.service.call({'op': 'weighted_poincare', 'q': [0], 'threshold': 2, **SMALL})

    def test_full_azimuth_sectors_cannot_exceed_requested_max_m(self):
        full = spectrum.full_ground([0, 0, 2], **SMALL)
        self.assertGreater(max(c['azimuth_m'] for c in full['sectors']), 0)
        request = {'op': 'mean_zero_ground', 'q': [0, 0, 2], **dict(SMALL, max_m=0)}
        with patch.object(tool.spectrum, 'full_ground', return_value=full), self.assertRaises(ValueError):
            self.service.call(request)
        inequality = spectrum.inequality_certificate(full, 2)
        with patch.object(tool.spectrum, 'weighted_poincare', return_value=inequality), self.assertRaises(ValueError):
            self.service.call(dict(request, op='weighted_poincare', threshold=2))

    def test_arithmetic_failure_is_separate_from_invalid_request_and_recovers(self):
        valid = spectrum.certify_sector([0], modes=2)
        request = {'op': 'sector', 'q': [0], 'modes': 2}
        with patch.object(tool.spectrum, 'certify_sector', side_effect=[ArithmeticError('injected internal failure'), valid]):
            results = self.request_lines([request, request])
        self.assertEqual(results[0]['status'], 'computation_failed')
        self.assertIsNone(results[0]['mathematical_verdict'])
        self.assertEqual(results[1]['status'], 'certified_sector_bound')
        invalid = self.request_lines([dict(request, bits=0)])[0]
        self.assertEqual(invalid['status'], 'invalid_request')

    def test_floats_booleans_and_decimal_strings_are_not_rational_coefficients(self):
        for value in (0.0, True, False, '0.5', 'NaN', '1/0'):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.service.call({'op': 'sector', 'q': [value]})

    def test_invalid_boolean_and_integer_sizes_are_rejected(self):
        for key, value in (('mean_zero', 1), ('m', True), ('m', -1), ('m', 65),
                           ('modes', 0), ('modes', 65), ('k', 0), ('k', 13), ('bits', 7), ('bits', 161), ('bits', 44.0)):
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                self.service.call({'op': 'sector', 'q': [0], key: value})

    def test_full_search_budget_and_tolerance_validation(self):
        for changes in ({'max_modes': 7}, {'max_m': True}, {'max_m': -1}, {'max_m': 65},
                        {'tolerance': 1e-10}, {'tolerance': True}, {'tolerance': '0'},
                        {'tolerance': '2'}, {'tolerance': '1/'+str(10**31)}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.service.call({'op': 'mean_zero_ground', 'q': [0], **changes})

    def test_threshold_float_and_boolean_rejected(self):
        for value in (2.0, True, '2.0'):
            with self.assertRaises(ValueError):
                self.service.call({'op': 'weighted_poincare', 'q': [0], 'threshold': value, **SMALL})

    def test_bad_json_unknown_operation_and_unknown_fields_recover_without_logs(self):
        output = io.StringIO()
        requests = '{broken}\n'+tool.wire({'op': 'unknown'})+'\n'+tool.wire({'op': 'sector', 'q': [0], 'typo': 1})+'\n'+tool.wire({'op': 'sector', 'q': [0], 'modes': 2})+'\n'
        tool.serve(io.StringIO(requests), output, self.service)
        rows = [json.loads(line) for line in output.getvalue().splitlines()]
        self.assertEqual(len(rows), 4)
        for row in rows[:3]: self.assertIsNone(row['mathematical_verdict'])
        self.assertEqual(rows[-1]['status'], 'certified_sector_bound')

    def test_quoted_command_payload_is_literal_and_never_executed(self):
        marker = Path(self.tmp.name)/'must-not-exist'
        literal = '__import__("pathlib").Path("'+str(marker)+'").touch(); `touch '+str(marker)+'`; $(touch '+str(marker)+')'
        rows = self.request_lines([{'op': 'sector', 'q': [literal]}, {'op': literal},
                                   {'op': 'sector', 'q': [0], 'modes': 2}])
        self.assertFalse(marker.exists())
        self.assertIsNone(rows[0]['mathematical_verdict']); self.assertIsNone(rows[1]['mathematical_verdict'])
        self.assertEqual(rows[2]['lower'], '2')


if __name__ == '__main__':
    unittest.main()
