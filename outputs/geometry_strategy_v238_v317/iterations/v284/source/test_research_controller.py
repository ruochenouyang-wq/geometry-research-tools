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

    def test_real_approximation_partition_changes_the_right_parameter(self):
        req = rc.normalize({'target': 'approximation', 'function': SINGULAR})
        for levels, degree, expected in ((1, 8, 'core'), (20, 0, 'annuli')):
            cert = rc.backend.pieces.piecewise_model(SINGULAR, levels=levels, degree=degree)
            d = rc.diagnose(cert, req)
            self.assertEqual(d['bottleneck'], expected)
            self.assertEqual(rc.F(d['core_squared'])+rc.F(d['annuli_squared']),
                             rc.F(cert['error_squared']))
            action = dict(rc.first_action(req), levels=levels, degree=degree)
            changed = rc.next_action(action, d)
            self.assertEqual(changed['degree' if expected == 'core' else 'levels'],
                             action['degree' if expected == 'core' else 'levels'])

    def test_real_full_residual_separates_coefficients_and_representation(self):
        req = rc.normalize({'target': 'spectrum', 'function': SINGULAR})
        source = rc.backend.direct.full_ground(SINGULAR, modes=2, bits=12)
        kinds = []
        for iterations in (0, 14):
            proposal = rc.backend.enriched.propose_trial(
                rc.backend.enriched.DEFAULT_EXPONENTS[:3], iterations=iterations)
            cert = rc.backend.enriched.trial_certificate(proposal['powers'],
                proposal['coefficients'], source)
            d = rc.diagnose(cert, req)
            self.assertEqual(rc.F(d['projected_squared'])+rc.F(d['outside_trial_squared']),
                             rc.F(d['complete_squared']))
            kinds.append(d['bottleneck'])
        self.assertEqual(kinds, ['coefficients', 'representation'])

    def test_real_spectral_bottleneck_selects_angular_or_radial_work(self):
        req = rc.normalize({'target': 'spectrum', 'function': SINGULAR, 'route': 'spectrum'})
        for max_m, expected in ((0, 'angular_tail'), (2, 'radial_sector')):
            cert = rc.backend.direct.full_ground(SINGULAR, modes=2, bits=12, max_m=max_m)
            d = rc.diagnose(cert, req)
            self.assertEqual(d['bottleneck'], expected)
            action = dict(rc.first_action(req), max_m=max_m)
            changed = rc.next_action(action, d)
            self.assertEqual(changed['max_m'] > max_m, expected == 'angular_tail')
            self.assertIn('not strictly separated', d['hypothesis'])

    def test_automatic_routes_obey_structure_and_goal(self):
        cases = [('approximation', STEP, 'approximate'), ('spectrum', STEP, 'spectrum'),
                 ('spectrum', SINGULAR, 'singular'),
                 ('spectrum', dict(SINGULAR, exponent='-1/3'), 'spectrum')]
        for target, q, expected in cases:
            a = rc.first_action(rc.normalize({'target': target, 'function': q}))
            self.assertEqual(a['route'], expected)
        result = rc.solve({'target': 'spectrum', 'function': SINGULAR, 'tolerance': '1'})
        self.assertEqual(result['attempts'][0]['action']['route'], 'singular')
        self.assertEqual(result['status'], 'target_met')

    def test_real_incremental_budget_reaches_approximation_goal(self):
        result = rc.solve({'target': 'approximation', 'function': SINGULAR,
                           'tolerance': '1/10', 'budget': {'max_attempts': 3}})
        self.assertEqual(result['status'], 'target_met')
        self.assertEqual(len(result['attempts']), 2)
        self.assertEqual(result['attempts'][1]['action']['levels'], 4)

    def test_work_budget_prevents_starting_unaffordable_action(self):
        result = rc.solve({'target': 'spectrum', 'function': SINGULAR,
                           'budget': {'work_units': 1}})
        self.assertEqual(result['cost']['solver_calls'], 0)
        self.assertEqual(result['stop_reason'], 'work_budget')

    def test_capability_callback_extends_dispatch_without_self_certification(self):
        req = rc.normalize({'target': 'spectrum', 'function': dict(SINGULAR, exponent='-1/3')})
        supports = lambda r: r['function']['exponent'] == '-1/3'
        action = rc.first_action(req, supports)
        self.assertEqual(action['route'], 'singular')
        self.assertEqual(len(rc.legal_actions(req, [action], supports)[0]), 1)
        self.assertEqual(rc.legal_actions(req, [action])[0], [])


if __name__ == '__main__':
    unittest.main()
