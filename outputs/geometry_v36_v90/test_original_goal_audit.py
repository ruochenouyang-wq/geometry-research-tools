"""Independent bounded V89 binding/optimization audit; no production edits."""
from fractions import Fraction as F
import copy
import unittest

import original_goal as original
import uniform_refine as uniform


def goal(potential, penalty, parameters, box=None, threshold='0', space='full_sphere'):
    return {'geometry': 'unit_sphere', 'function_space': space, 'eigenvalue_index': 1,
            'parameters': parameters, 'potential': potential, 'penalty': penalty,
            'domain': {'kind': 'all_real'} if box is None else {'kind': 'box', 'box': box},
            'threshold': threshold}


class OriginalGoalIndependentAudit(unittest.TestCase):
    def test_four_parameter_offset_constants_slopes_and_cross_terms(self):
        potential = '7/5+3*t/2+a*(1+2*t)+b*(-2-t)+c*(3+3*t/2)+d*(-5-4*t)'
        penalty = ('(3/2+2*a-b+3*c/2-4*d)**2/6-a+2*b-3*c+5*d-7/5'
                   '+2*a*a+3*b*b+4*c*c+5*d*d+a*b+2*b*c-3*c*d'
                   '+(a-2*b+3*c-4*d)/7+11/13')
        cert = original.solve(goal(potential, penalty, ['a', 'b', 'c', 'd'], threshold='-100'))
        self.assertTrue(original.verify(cert))
        expected = {'constant': '11/13', 'linear': ['1/7', '-2/7', '3/7', '-4/7'],
                    'hessian': [['4', '1', '0', '0'], ['1', '6', '2', '0'],
                                ['0', '2', '8', '-3'], ['0', '0', '-3', '10']]}
        self.assertEqual(cert['derived_quadratic'], expected)
        self.assertEqual(cert['amplitude_offset'], '3/2')
        self.assertEqual(cert['amplitude_slopes'], ['2', '-1', '3/2', '-4'])
        self.assertEqual(cert['parameter_constants'], ['1', '-2', '3', '-5'])
        for point in [(F(1, 3), F(-2, 5), F(3, 7), F(-4, 9)), (F(0),)*4]:
            a, b, c, d = point
            expanded = 2*a*a + 3*b*b + 4*c*c + 5*d*d + a*b + 2*b*c - 3*c*d
            expanded += (a-2*b+3*c-4*d)/7 + F(11, 13)
            h = [[F(x) for x in row] for row in expected['hessian']]
            stored = F(expected['constant']) + sum(F(x)*v for x, v in zip(expected['linear'], point))
            stored += sum(point[i]*h[i][j]*point[j]/2 for i in range(4) for j in range(4))
            self.assertEqual(stored, expanded)

    def test_nonconvex_box_minimum_on_face_interior(self):
        raw = goal('4+2*t+a*(3+2*t)+b*(-5+t)',
                   '(2+2*a+b)**2/6-3*a+5*b-4-a*a+(b-a/2)**2',
                   ['a', 'b'], [['-1', '1'], ['-1', '1']], '-1')
        cert = original.solve(raw)
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert['lower_bound'], '-1')
        a, b = map(F, cert['relaxation_minimizer'])
        self.assertEqual(abs(a), 1)
        self.assertEqual(b, a/2)
        self.assertEqual(cert['status'], 'proved')

    def test_singular_face_does_not_hide_nonconvex_global_minimum(self):
        raw = goal('(2*a-b+3*c)*t', '(2*a-b+3*c)**2/6+(a+b)**2-c*c',
                   ['a', 'b', 'c'], [['-1', '1'], ['-1', '1'], ['-2', '3']], '-9')
        cert = original.solve(raw)
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert['lower_bound'], '-9')
        a, b, c = map(F, cert['relaxation_minimizer'])
        self.assertEqual(a+b, 0)
        self.assertEqual(c, 3)
        self.assertEqual(cert['status'], 'proved')

    def test_degenerate_parameter_box_is_rejected_by_existing_compiler(self):
        # The shared compiler requires strictly positive side widths. This is
        # an explicit input limitation, not a supported singleton-box case.
        with self.assertRaises(ValueError):
            original.solve(goal('a*t', 'a*a/6+(a-b)**2', ['a', 'b'],
                                [['1', '1'], ['2', '3']], '1'))

    def test_degree_six_noneven_remainder_binding_and_reintegration(self):
        potential = '9/7+5*t/3+(t**3-2*t*t/3)**2+a*(4+2*t)'
        penalty = '(5/3+2*a)**2/6-9/7-4*a+a*a'
        cert = original.solve(goal(potential, penalty, ['a']))
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert['remainder_polynomial'], ['9/7', '0', '0', '0', '4/9', '-4/3', '1'])
        self.assertEqual(cert['remainder_range']['coordinate'], 't')
        self.assertEqual(cert['remainder_range']['polynomial'], ['-9/7', '0', '0', '0', '-4/9', '4/3', '-1'])
        self.assertLessEqual(F(cert['remainder_lower']), F(9, 7))
        self.assertEqual(cert['derived_quadratic']['linear'], ['0'])
        self.assertEqual(cert['derived_quadratic']['hessian'], [['2']])
        for witness in cert['constant_trial_candidates']:
            p = list(map(F, witness['original_potential']))
            mean = sum((p[k]/(k+1) for k in range(0, len(p), 2)), F(0))
            self.assertEqual(mean, F(witness['spectral_rayleigh_upper']))
            self.assertEqual(mean + F(witness['original_penalty']), F(witness['objective_upper']))

    def test_constant_energy_shift_changes_bounds_and_threshold_together(self):
        raw = goal('t**3+t*t+2*t+a*(3-t)', '(a-1)**2', ['a'], [['-2', '2']], '-1')
        first = original.solve(raw)
        shifted = copy.deepcopy(raw)
        shifted['potential'] = '('+raw['potential']+')+7/11'
        shifted['threshold'] = str(F(raw['threshold'])+F(7, 11))
        second = original.solve(shifted)
        self.assertTrue(original.verify(second))
        for bound in ('lower_bound', 'upper_bound'):
            self.assertEqual(F(second[bound])-F(first[bound]), F(7, 11))
        self.assertEqual(first['status'], second['status'])
        self.assertEqual(first['infimum_gap'], second['infimum_gap'])

    def test_true_goal_with_unbounded_lower_relaxation_stays_unresolved(self):
        # Independently, lambda(-Delta+a*t) >= -abs(a), hence the objective
        # is >= 10-abs(a)+a*a/7 >= 33/4. V89's chosen quadratic is weaker.
        cert = original.solve(goal('a*t', 'a*a/7+10', ['a']))
        self.assertTrue(original.verify(cert))
        self.assertIsNone(cert['lower_bound'])
        self.assertEqual(cert['upper_bound'], '10')
        self.assertEqual(cert['status'], 'unresolved')

    def test_actual_rayleigh_refutation_for_nonlinear_fixed_potential(self):
        cert = original.solve(goal('t*t-3', '0', [], threshold='-5/2'))
        self.assertTrue(original.verify(cert))
        self.assertEqual(cert['constant_trial']['probability_mass'], '1')
        self.assertEqual(cert['constant_trial']['dirichlet_energy'], '0')
        self.assertEqual(cert['upper_bound'], '-8/3')
        self.assertEqual(cert['status'], 'refuted')

    def test_spaces_and_strict_threshold_boundary(self):
        for space in ('full_sphere', 'axisymmetric'):
            exact = original.solve(goal('0', '0', [], space=space))
            self.assertEqual(exact['status'], 'proved')
            self.assertEqual(exact['operator_binding']['function_space'], space)
            above = original.solve(goal('0', '0', [], threshold='1/999999937', space=space))
            self.assertEqual(above['status'], 'refuted')
            touching = original.solve(goal('t', '0', [], space=space))
            self.assertEqual(touching['upper_bound'], touching['threshold'])
            self.assertEqual(touching['status'], 'unresolved')
        with self.assertRaises(ValueError):
            original.solve(goal('a*t', 'a*a/6', ['a'], space='mean_zero'))

    def test_genuine_failed_or_finite_lemma_is_not_success(self):
        cert = original.solve(goal('a*t', 'a*a/5', ['a']))
        failed = uniform.synthesize(coefficient='1/5', order=1, all_real=True)
        self.assertTrue(uniform.verify(failed))
        self.assertFalse(uniform.is_proved(failed))
        with self.assertRaises(ValueError):
            original._build(cert['compiled_goal'], F(1, 5), 'barta', failed, cert['remainder_range'])
        finite = uniform.synthesize(coefficient='1/5', order=3, amplitude_box=['-2', '2'])
        self.assertTrue(uniform.is_proved(finite))
        with self.assertRaises(ValueError):
            original._build(cert['compiled_goal'], F(1, 5), 'barta', finite, cert['remainder_range'])


if __name__ == '__main__':
    unittest.main()
