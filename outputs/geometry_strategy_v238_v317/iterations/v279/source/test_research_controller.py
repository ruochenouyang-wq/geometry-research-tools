import unittest
from copy import deepcopy
import research_controller as rc

STEP = {'kind': 'axis_profile', 'profile': 'step', 'amplitude': '-1'}
SINGULAR = deepcopy(rc.backend.enriched.FUNCTION)


class ControllerTests(unittest.TestCase):
    def test_real_step_approximation_is_exact(self):
        result = rc.solve({'target': 'approximation', 'function': STEP})
        self.assertEqual(result['status'], 'target_met')
        self.assertEqual(result['certificate']['error_squared'], '0')

    def test_real_spectral_goal_is_distinct(self):
        result = rc.solve({'target': 'spectrum', 'function': STEP, 'tolerance': '1'})
        self.assertEqual(result['status'], 'target_met')
        self.assertIn('lower', result['certificate'])
        self.assertNotIn('error_upper', result['certificate'])

    def test_forged_verified_flag_cannot_accept(self):
        bad = lambda action: {'verified': True, 'certificate': {'error_upper': '0',
                              'spectral_transfer_claimed': False}}
        with self.assertRaises(ValueError):
            rc.solve({'target': 'approximation', 'function': STEP}, {'approximate': bad})

    def test_unsupported_approximation_never_calls_solver(self):
        q = dict(SINGULAR, exponent='-1/3')
        called = []
        result = rc.solve({'target': 'approximation', 'function': q},
                          {'approximate': lambda action: called.append(action)})
        self.assertEqual(result['status'], 'unsupported')
        self.assertEqual(called, [])

    def test_actions_cannot_change_scope_or_domain(self):
        req = rc.normalize({'target': 'spectrum', 'function': STEP})
        actions = [dict(rc.first_action(req), mean_zero=False),
                   dict(rc.first_action(req), modes=0),
                   dict(rc.first_action(req), function=SINGULAR)]
        good, bad = rc.legal_actions(req, actions)
        self.assertEqual(good, [])
        self.assertEqual(len(bad), 3)


if __name__ == '__main__':
    unittest.main()
