import hashlib
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import compat_controller as repaired


def public_rows():
    path = ROOT.parent / 'geometry_strategy_v317_runtime/evaluation/holdout_20260920T042408501655Z/release_r0_output.json'
    return json.loads(path.read_text())['rows']


class ControllerRepairs(unittest.TestCase):
    def test_singular_fallback_refinement_stays_legal(self):
        for row in public_rows():
            if row['case_id'] not in ('H07', 'H09'):
                continue
            request = {k: v for k, v in row['task'].items() if k != 'kind'}
            request['target'] = 'spectrum'
            action = repaired._controller.first_action(request, lambda _: True)
            updated = repaired.next_action(action, {'kind': 'spectral_lower_budget',
                                                   'bottleneck': 'radial_sector'})
            legal, rejected = repaired.legal_actions(request, [updated])
            self.assertEqual(rejected, [])
            self.assertEqual(legal[0]['max_terms'], 6)
            self.assertNotIn('near_tail', legal[0])

    def test_new_certificate_residual_schema_and_frozen_isolation(self):
        import research_controller as frozen
        frozen_file = repaired.OLD / 'research_controller.py'
        before = hashlib.sha256(frozen_file.read_bytes()).hexdigest()
        row = next(r for r in public_rows() if r['case_id'] == 'H05')
        certificate = row['attempts'][0]['response']['certificate']
        self.assertEqual(certificate['format'], 'adaptive_singular_temple_v1')
        request = {'target': 'spectrum', **row['task']}
        diagnosis = repaired.diagnose(certificate, request)
        self.assertEqual(diagnosis['evidence'], 'exact_orthogonal_decomposition')
        self.assertEqual(repaired.backend.F(diagnosis['projected_squared']) +
                         repaired.backend.F(diagnosis['outside_trial_squared']),
                         repaired.backend.F(diagnosis['complete_squared']))
        self.assertIsNot(frozen.diagnose, repaired.diagnose)
        self.assertEqual(hashlib.sha256(frozen_file.read_bytes()).hexdigest(), before)


if __name__ == '__main__':
    unittest.main()
