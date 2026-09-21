import unittest
from time import sleep
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
        result = rc.solve({'target': 'approximation', 'function': STEP}, {'approximate': bad})
        self.assertIsNone(result['certificate'])
        self.assertEqual(result['status'], 'verification_failure')
        self.assertEqual(result['cost']['failed_attempts'], 1)

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

    def test_no_gain_stops_and_preserves_the_best_replayed_result(self):
        good = rc.backend.pieces.piecewise_model(SINGULAR, levels=6, degree=4)
        worse = rc.backend.pieces.piecewise_model(SINGULAR, levels=2, degree=2)
        remaining = [good, worse, worse]
        result = rc.solve({'target': 'approximation', 'function': SINGULAR,
                           'budget': {'max_attempts': 8}},
                          {'approximate': lambda action: remaining.pop(0)})
        self.assertEqual(result['stop_reason'], 'no_meaningful_gain')
        self.assertEqual(result['certificate'], good)
        self.assertEqual(result['cost']['solver_calls'], 3)
        self.assertEqual(result['attempts'][-1]['best_metric'], good['error_upper'])

    def test_real_route_comparison_uses_two_probes_then_better_branch(self):
        result = rc.solve({'target': 'spectrum', 'function': SINGULAR,
                           'budget': {'max_attempts': 3}})
        routes = [a['action']['route'] for a in result['attempts']]
        self.assertEqual(routes, ['singular', 'spectrum', 'singular'])
        self.assertEqual(set(result['route_comparison']), {'singular', 'spectrum'})
        self.assertLess(rc.F(result['route_comparison']['singular']['best_metric']),
                        rc.F(result['route_comparison']['spectrum']['best_metric']))

    def test_failed_specialized_solver_falls_back_and_is_costed(self):
        def broken(action):
            raise RuntimeError('Candidate search failed')
        result = rc.solve({'target': 'spectrum', 'function': SINGULAR, 'tolerance': '1',
                           'policy': {'compare_routes': False}},
                          {'singular': broken, 'spectrum': rc.backend.legacy.handle})
        self.assertEqual(result['status'], 'target_met')
        self.assertEqual([a['action']['route'] for a in result['attempts']], ['singular', 'spectrum'])
        self.assertEqual(result['cost']['solver_calls'], 2)
        self.assertEqual(result['cost']['failed_attempts'], 1)
        self.assertGreater(result['cost']['work_units_charged'], result['attempts'][1]['work_units_charged'])
        self.assertGreater(result['cost']['verification_seconds'], 0)

    def test_invalid_certificate_falls_back_without_trusting_self_report(self):
        cert = rc.backend.direct.full_ground(SINGULAR, modes=2, bits=12, tolerance='1')
        cert['lower'] = cert['upper']
        cert['exact_width'] = '0'
        result = rc.solve({'target': 'spectrum', 'function': SINGULAR, 'tolerance': '1'},
                          {'singular': lambda action: {'verified': True, 'certificate': cert},
                           'spectrum': rc.backend.legacy.handle})
        self.assertEqual(result['status'], 'target_met')
        self.assertFalse(result['attempts'][0]['verified'])
        self.assertEqual(result['attempts'][0]['failed_phase'], 'verify')
        self.assertNotEqual(result['certificate']['exact_width'], '0')

    def test_late_valid_certificate_is_not_budget_success(self):
        certificate = rc.backend.pieces.piecewise_model(STEP)
        def slow(action):
            sleep(0.01)
            return certificate
        result = rc.solve({'target': 'approximation', 'function': STEP,
                           'budget': {'wall_seconds': 0.005}}, {'approximate': slow})
        self.assertEqual(result['status'], 'budget_exceeded')
        self.assertTrue(result['mathematical_target_met'])
        self.assertFalse(result['within_budget'])
        self.assertEqual(result['certificate'], certificate)

    def test_zero_wall_budget_starts_no_backend(self):
        result = rc.solve({'target': 'approximation', 'function': STEP,
                           'budget': {'wall_seconds': 0}})
        self.assertEqual(result['cost']['solver_calls'], 0)
        self.assertEqual(result['status'], 'budget_exceeded')

    def test_injected_verifier_gets_exact_original_bindings(self):
        bindings = []
        q = dict(SINGULAR, exponent='-1/3')
        def verifier(cert, **kwargs):
            bindings.append(kwargs)
            return cert == {'format': 'test_only', 'lower': '0', 'upper': '1/10', 'exact_width': '1/10'}
        solvers = {'singular_supports': lambda r: True,
                   'singular': lambda action: {'format': 'test_only', 'lower': '0',
                                               'upper': '1/10', 'exact_width': '1/10'}}
        result = rc.solve({'target': 'spectrum', 'function': q, 'tolerance': '1/5'}, solvers, verifier)
        self.assertEqual(result['status'], 'target_met')
        self.assertEqual(bindings, [{'expected_function': q, 'expected_mean_zero': True,
                                     'expected_tolerance': '1/5'}])

    def test_adaptive_same_shape_diagnostic_fixture_and_refinement(self):
        # Controller compatibility fixture, not a new mathematical verification claim.
        source = rc.backend.direct.full_ground(SINGULAR, modes=2, bits=12)
        proposal = rc.backend.enriched.propose_trial(rc.backend.enriched.DEFAULT_EXPONENTS[:3])
        cert = rc.backend.enriched.trial_certificate(proposal['powers'], proposal['coefficients'], source)
        cert['format'] = 'adaptive_singular_temple_v1'
        req = rc.normalize({'target': 'spectrum', 'function': SINGULAR})
        d = rc.diagnose(cert, req)
        self.assertEqual(d['kind'], 'trial_residual_budget')
        self.assertEqual(rc.next_action(rc.first_action(req), d)['max_terms'], 6)

    def test_malformed_budgets_are_not_ignored(self):
        for extra in ({'budget': {'wall_seconds': float('inf')}},
                      {'budget': {'wall_seconds': True}},
                      {'budget': {'invented_budget': 10}},
                      {'policy': {'compare_routes': 'false'}}):
            with self.assertRaises((ValueError, TypeError)):
                rc.solve(dict(target='approximation', function=STEP, **extra))


if __name__ == '__main__':
    unittest.main()
