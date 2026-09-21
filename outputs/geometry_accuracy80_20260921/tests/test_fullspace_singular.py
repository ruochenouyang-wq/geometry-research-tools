"""Independent binding, mathematics and budget checks for the full-space route."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import signal
import sys
import time
import unittest

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import fullspace_singular as fs
import spectral_gap
import power_forms


def function(alpha='-2/5', amplitude='-1/3', offset='0', axis=2):
    return {'kind': 'axis_profile', 'axis': axis, 'profile': 'abs_power',
            'exponent': alpha, 'amplitude': amplitude, 'offset': offset}


def bounded_solve(q, tolerance, **options):
    started = time.perf_counter()
    previous = signal.getsignal(signal.SIGALRM)
    def timeout(signum, frame):
        raise TimeoutError('Ten-second solve plus independent replay cap')
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, 10)
    try:
        result = fs.solve(q, tolerance, **options)
        cert = result['certificate']
        valid = fs.verify(cert, q, False, tolerance)
        return result, valid, time.perf_counter()-started
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
        signal.signal(signal.SIGALRM, previous)


class FullspaceSingularTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.inputs = {
            'H06': (function(), '1/1000000'),
            'H08': (function('-3/10'), '1/100000000'),
            'H10': (function('-3/7', '-3/5'), '1/10000000000'),
        }
        cls.results = {name: bounded_solve(q, tol)
                       for name, (q, tol) in cls.inputs.items()}

    def test_all_three_original_fullspace_tasks_meet_bound_and_replay(self):
        for name, (q, tol) in self.inputs.items():
            with self.subTest(name=name):
                result, valid, elapsed = self.results[name]
                self.assertTrue(valid)
                self.assertTrue(result['target_met'])
                self.assertLess(elapsed, 10)
                cert = result['certificate']
                self.assertIs(cert['mean_zero'], False)
                self.assertLessEqual(F(cert['exact_width']), F(tol))
                self.assertEqual(cert['gap']['azimuth_m'], 0)
                self.assertIs(cert['gap']['mean_zero'], False)
                self.assertEqual(cert['gap']['eigenvalue_index'], 2)
                self.assertLess(F(cert['upper']), F(cert['gap']['beta']))
                self.assertTrue(cert['full_infinite_space_covered'])

    def test_failed_precision_stages_are_preserved(self):
        result = self.results['H10'][0]
        self.assertEqual([r['terms'] for r in result['attempts']], [3, 6, 10])
        self.assertEqual([r['status'] for r in result['attempts']],
                         ['certified_open', 'certified_open', 'target_met'])
        self.assertTrue(all(r['verified'] for r in result['attempts']))
        self.assertTrue(all(r['elapsed_seconds'] > 0 for r in result['attempts']))

    def test_constant_moments_are_known_exactly(self):
        M, H, R = fs.trial_matrices(function(), [0])
        self.assertEqual(M, [[F(2)]])
        self.assertEqual(H, [[F(-10, 9)]])
        self.assertEqual(R, [[F(10, 9)]])
        stats = fs.statistics(function(), [0], [1])
        self.assertEqual(F(stats['rayleigh']), F(-5, 9))
        self.assertEqual(F(stats['residual_squared']), F(20, 81))

    def test_operator_uses_m0_coefficient_and_cancels_leading_singularity(self):
        q = function()
        delta = F(8, 5)
        self.assertEqual(fs.operator_terms(q, delta)[delta], delta*(delta+1))
        ss = fs.generated_powers(q, 3)
        C = fs.cancellation_transform(q, ss)
        v = fs.matvec(C, [F(1), F(0)])
        combined = {}
        for s, coeff in zip(ss, v):
            for r, value in fs.operator_terms(q, s).items():
                combined[r] = combined.get(r, F(0))+coeff*value
        self.assertEqual(combined[F(q['exponent'])], 0)

    def test_strong_domain_rejects_nonsquare_integrable_residual(self):
        for ss in ([0, '3/2'], [0, '1'], [0, '7/5']):
            with self.subTest(powers=ss), self.assertRaises(ValueError):
                fs.trial_matrices(function(), ss)
        with self.assertRaises(ValueError):
            fs.normalize(function('-1/2'))

    def test_clipped_form_h10_exact_integral_and_positive_coercivity(self):
        q = self.inputs['H10'][0]
        cert = fs.form_certificate(q)
        self.assertTrue(fs.verify_form(cert, q))
        self.assertEqual(cert['method'], 'clipped_negative_L2')
        self.assertEqual(F(cert['proof']['delta']), F(1, 128))
        self.assertEqual(F(cert['proof']['cap']), F(24, 5))
        self.assertEqual(F(cert['proof']['residual_L2_squared']), F(81, 100))
        # The outward dyadic square root is slightly above exact 9/10.
        self.assertGreaterEqual(F(cert['proof']['eta_upper']), F(9, 10))
        self.assertGreater(F(cert['form']['energy_factor']), 0)
        self.assertLessEqual(F(cert['form']['energy_factor']), F(1, 10))
        with self.assertRaises(ValueError):
            fs.form_certificate(q, method='centered_L2')

    def test_form_evidence_and_function_binding_reject_tampering(self):
        q = self.inputs['H10'][0]
        original = fs.form_certificate(q)
        for section, key, value in [
            ('proof', 'residual_L2_squared', '0'),
            ('proof', 'cap', '0'), ('proof', 'eta_upper', '0'),
            ('form', 'energy_factor', '1'), ('form', 'mass_offset', '0'),
        ]:
            changed = deepcopy(original)
            changed[section][key] = value
            with self.subTest(key=key):
                self.assertFalse(fs.verify_form(changed, q))
        self.assertFalse(fs.verify_form(original, function('-3/7', '-2/3')))

    def test_full_certificate_binds_original_quantity_and_space(self):
        q, tol = self.inputs['H10']
        cert = self.results['H10'][0]['certificate']
        self.assertFalse(fs.verify(cert, q, True, tol))
        self.assertFalse(fs.verify(cert, function('-3/7', '-2/3'), False, tol))
        self.assertFalse(fs.verify(cert, q, False, '1/1000000000'))
        self.assertFalse(fs.verify(cert, dict(q, axis=1), False, tol))
        changed = deepcopy(cert); changed['mean_zero'] = True
        self.assertFalse(fs.verify(changed, q, False, tol))
        with self.assertRaises(ValueError):
            fs.solve(q, tol, mean_zero=True)

    def test_numeric_proof_fields_cannot_be_edited_into_success(self):
        q, tol = self.inputs['H06']
        cert = self.results['H06'][0]['certificate']
        for field, value in [('lower', cert['upper']), ('exact_width', '0'),
                             ('matrix_digest', '0'*64), ('gap_digest', '0'*64),
                             ('scope', 'finite_trial_only')]:
            changed = deepcopy(cert); changed[field] = value
            with self.subTest(field=field):
                self.assertFalse(fs.verify(changed, q, False, tol))
        changed = deepcopy(cert)
        changed['statistics']['residual_squared'] = '0'
        self.assertFalse(fs.verify(changed, q, False, tol))
        changed = deepcopy(cert)
        changed['trial']['coefficients'][0] = '2'
        self.assertFalse(fs.verify(changed, q, False, tol))
        changed = deepcopy(cert)
        changed['full_space_domination']['additional_energy_nonnegative'] = False
        self.assertFalse(fs.verify(changed, q, False, tol))

    def test_gap_kernel_comparison_and_custom_form_are_replayed(self):
        q, tol = self.inputs['H10']
        original = self.results['H10'][0]['certificate']
        for target in ('comparison', 'moments', 'form', 'm', 'index', 'barrier'):
            changed = deepcopy(original)
            gap = changed['gap']
            if target == 'comparison': gap['comparison_matrix'][0][0] = '0'
            elif target == 'moments': gap['kernel_evidence']['A'][0][0] = '0'
            elif target == 'form': gap['form_certificate']['form']['mass_offset'] = '0'
            elif target == 'm': gap['azimuth_m'] = 1
            elif target == 'index': gap['eigenvalue_index'] = 1
            else: gap['lower'] = gap['radial_tail_lower']; gap['beta'] = gap['lower']
            changed['gap_digest'] = fs.digest(gap)
            with self.subTest(target=target):
                self.assertFalse(fs.verify(changed, q, False, tol))
        self.assertFalse(spectral_gap.verify_gap(original['gap']))

    def test_exact_offset_covariance_of_all_operator_matrices(self):
        q = function('-1/3', '-2/5')
        ss = fs.generated_powers(q, 6)
        M, H, R = fs.trial_matrices(q, ss)
        b = F(7, 3)
        shifted = fs.trial_matrices(dict(q, offset=str(b)), ss)
        self.assertEqual(shifted[0], M)
        self.assertEqual(shifted[1], [[h+b*M[i][j] for j, h in enumerate(row)]
                                     for i, row in enumerate(H)])
        self.assertEqual(shifted[2], [[r+2*b*H[i][j]+b*b*M[i][j] for j, r in enumerate(row)]
                                     for i, row in enumerate(R)])

    def test_additional_public_parameter_and_positive_offset(self):
        q = function('-1/3', '-2/5', '7/3', axis=0)
        result, valid, elapsed = bounded_solve(q, '1/10000000')
        self.assertTrue(valid)
        self.assertTrue(result['target_met'])
        self.assertLess(elapsed, 10)
        self.assertGreater(F(result['certificate']['lower']), 0)

    def test_new_h06_interval_sits_in_independent_frozen_coarse_interval(self):
        q, _ = self.inputs['H06']
        coarse = fs.direct.full_ground(q, mean_zero=False, modes=4, bits=24, near_tail=4)
        cert = self.results['H06'][0]['certificate']
        self.assertLessEqual(F(coarse['lower']), F(cert['lower']))
        self.assertLessEqual(F(cert['upper']), F(coarse['upper']))

    def test_invalid_budget_does_not_start_a_search(self):
        for budget in (True, -1, float('nan'), float('inf'), 61):
            with self.subTest(budget=budget), self.assertRaises(ValueError):
                fs.solve(function(), wall_seconds=budget)

    def test_outward_cube_root_has_exact_direction_and_one_unit_resolution(self):
        for value in (F(0), F(1), F(8), F(27, 125), F(81, 100), F(1, 10**30), F(10**30)):
            for bits in (8, 40, 80):
                with self.subTest(value=value, bits=bits):
                    root = power_forms.cube_root_upper(value, bits)
                    self.assertGreaterEqual(root**3, value)
                    if value:
                        self.assertLess((root-F(1, 2**bits))**3, value)
                    else:
                        self.assertEqual(root, 0)
        with self.assertRaises(ValueError):
            power_forms.cube_root_upper(F(-1))

    def test_sqrt_tangent_is_an_upper_bound_on_the_declared_core(self):
        for t in (F(0), F(1, 100), F(1, 3), F(3, 4), F(1)):
            tangent = 1-t/2
            self.assertGreaterEqual(tangent, 0)
            self.assertGreaterEqual(tangent*tangent, 1-t)
        # An independently simplified integral coefficient checks signs and
        # the probability-measure normalization without using production D.
        q = function('-13/27', '-22/23')
        cert = power_forms.certificate(q, '67/50')
        p = F(13, 27)
        independent_D = 2*p*(2+5*p)/((2-3*p)*(4-p*p))
        self.assertEqual(F(cert['proof']['integral_coefficient_D']), independent_D)
        base, A = F(67, 50), F(22, 23)
        self.assertEqual(F(cert['proof']['residual_norm_cube_upper']),
                         A**3*base**(-15)*independent_D**2)

    def test_selected_holder_form_replays_and_opens_the_strong_initial_tail(self):
        q = function('-13/27', '-22/23')
        cert = power_forms.select_form(q, 12)
        self.assertTrue(power_forms.verify_form(cert, q))
        self.assertTrue(fs.verify_form(cert, q))
        self.assertGreater(fs.direct.tail_bound(cert['form'], 12), 2)
        candidates = cert['selection']['attempts']
        selected = next(x for x in candidates if x['base'] == cert['selection']['selected_base'])
        self.assertEqual(F(selected['initial_tail_lower']),
                         max(F(x['initial_tail_lower']) for x in candidates if x['admissible']))
        proof = cert['selected_form']['proof']
        self.assertGreaterEqual(F(proof['eta_upper'])**3, F(proof['residual_norm_cube_upper']))
        self.assertLess(F(proof['eta_upper']), F(1, 2))

    def test_holder_form_rejects_direction_measure_and_selection_tampering(self):
        q = function('-13/27', '-22/23')
        original = power_forms.select_form(q, 12)
        for field, value in [('eta_upper', '0'), ('eta_upper_cube', '0'),
                             ('residual_norm_cube_upper', '0'),
                             ('integral_coefficient_D', '1'),
                             ('integration_region', 'whole_sphere'),
                             ('dual_function_exponent', '4')]:
            cert = deepcopy(original)
            cert['selected_form']['proof'][field] = value
            with self.subTest(field=field):
                self.assertFalse(power_forms.verify_form(cert, q))
        cert = deepcopy(original); cert['form']['energy_factor'] = '1'
        self.assertFalse(power_forms.verify_form(cert, q))
        cert = deepcopy(original); cert['measure'] = 'd_sigma'
        self.assertFalse(power_forms.verify_form(cert, q))
        cert = deepcopy(original); cert['selection']['attempts'][0]['initial_tail_lower'] = '999'
        self.assertFalse(power_forms.verify_form(cert, q))
        self.assertFalse(power_forms.verify_form(original, dict(q, offset='1')))

    def test_predeclared_strong_public_fullspace_stress_cases(self):
        for alpha in ('-5/12', '-13/27'):
            with self.subTest(alpha=alpha):
                q = function(alpha, '-22/23')
                result, valid, elapsed = bounded_solve(q, '1/10000000000')
                self.assertTrue(valid)
                self.assertTrue(result['target_met'])
                self.assertLess(elapsed, 10)
                self.assertEqual(result['gap_attempts'][0]['modes'], 12)
                self.assertEqual(result['form_certificate']['format'], power_forms.SELECTED_FORMAT)

    def test_second_gap_stage_retains_first_open_attempt_without_restarting_history(self):
        q = function('-5/12', '-22/23')
        result, valid, elapsed = bounded_solve(q, '1/10000000000', max_terms=3)
        self.assertTrue(valid)
        self.assertFalse(result['target_met'])
        self.assertLess(elapsed, 10)
        self.assertEqual([x['modes'] for x in result['gap_attempts']], [12, 16])
        self.assertEqual([x['gap_modes'] for x in result['attempts']], [12, 16])
        self.assertTrue(all(x['status'] == 'certified_open' for x in result['attempts']))

    def test_centered_weak_parameter_route_keeps_the_original_gap_stage(self):
        for name in ('H06', 'H08'):
            result = self.results[name][0]
            self.assertEqual([x['modes'] for x in result['gap_attempts']], [8])
            self.assertEqual(result['form_certificate']['method'], 'centered_L2')


if __name__ == '__main__':
    unittest.main()
