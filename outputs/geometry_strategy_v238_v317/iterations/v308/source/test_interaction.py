import json
from pathlib import Path
import unittest

import interaction as subject


def example_certificate():
    # Synthetic contract fixture: callbacks below test transport, not mathematics.
    return {'format': 'fixture_spectrum', 'function': {'q': '0'},
            'scope': 'fixture_space', 'mean_zero': True,
            'tolerance': '1/10', 'lower': '1', 'upper': '2',
            'exact_width': '1', 'status': 'certified_open'}


class SummaryTests(unittest.TestCase):
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
