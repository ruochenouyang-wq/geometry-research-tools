"""Independent original-function, norm, source and transfer acceptance checks.

Rodrigues coefficients and analytic moments below do not call the new model
generator or its verification helpers. Frozen imports disable bytecode writes.
"""
from copy import deepcopy
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json
import sys
import unittest

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from dependencies import wide
import function_models as fm
import spectral_transfer as st


def rodrigues(l):
    return {l - 2*k: F((-1)**k * factorial(2*l-2*k),
                       2**l * factorial(k) * factorial(l-k) * factorial(l-2*k))
            for k in range(l//2+1)}


def profile_projection(degree, exponent=F(-1, 4), amplitude=F(-1), offset=F(0)):
    """Probability-L2 projection of offset+amplitude*|t|**exponent."""
    p = [F(0)] * (degree+1)
    coeffs = []
    for l in range(degree+1):
        leg = rodrigues(l)
        inner = sum((c*(amplitude/(j+exponent+1)+offset/F(j+1))
                     for j, c in leg.items() if j % 2 == 0), F(0))
        a = (2*l+1)*inner
        coeffs.append(a)
        for j, c in leg.items():
            p[j] += a*c
    norm = offset**2 + 2*offset*amplitude/(exponent+1) + amplitude**2/(2*exponent+1)
    error = norm - sum((a*a/F(2*l+1) for l, a in enumerate(coeffs)), F(0))
    return p, coeffs, error


def radial_inner(p, q):
    return sum((a*b/F(i+j+1) for i, a in enumerate(p) for j, b in enumerate(q)
                if (i+j) % 2 == 0), F(0))


class IndependentSourceAndMath(unittest.TestCase):
    def test_rodrigues_orthogonality_and_probability_normalization(self):
        polys = []
        for l in range(9):
            p = [F(0)]*(l+1)
            for j, c in rodrigues(l).items():
                p[j] = c
            polys.append(p)
        for i, p in enumerate(polys):
            for j, q in enumerate(polys):
                self.assertEqual(radial_inner(p, q), F(1, 2*i+1) if i == j else 0)

    def test_singular_exact_error_not_sampling(self):
        expected = {4: F(5202, 43681), 8: F(14450, 159201),
                    12: F(4118450, 54066609),
                    24: F(1335263295935522, 24084192784943929)}
        for d, error in expected.items():
            p, cs, actual = profile_projection(d)
            self.assertEqual(actual, error)
            self.assertEqual(cs[0], F(-4, 3))
            self.assertTrue(all(c == 0 for c in cs[1::2]))
            self.assertEqual(radial_inner(p, p), 2-error)

    def test_fixed_constant_reference_uniform_theta(self):
        for d in (0, 4, 8, 12, 24):
            p, _, error = profile_projection(d)
            p[0] += F(4, 3)
            theta2 = radial_inner(p, p)
            self.assertEqual(theta2, F(2, 9)-error)
            self.assertLessEqual(theta2, F(2, 9))
            self.assertLess(theta2, F(1, 4))

    def test_nested_projection_reference_error_identity(self):
        for r, n in ((4, 8), (4, 12), (8, 12), (8, 24)):
            pr, _, er = profile_projection(r)
            pn, _, en = profile_projection(n)
            difference = [a-(pr[j] if j < len(pr) else 0) for j, a in enumerate(pn)]
            self.assertEqual(radial_inner(difference, difference), er-en)
            self.assertLessEqual(er-en, er)

    def test_amplitude_offset_moments(self):
        p, _, error = profile_projection(8, F(1, 3), F(-3, 2), F(5, 7))
        self.assertEqual(radial_inner(p, [1]), F(5, 7)-F(3, 2)*F(3, 4))
        _, _, base_error = profile_projection(8, F(1, 3), F(1), F(0))
        self.assertEqual(error, F(9, 4)*base_error)

    def test_original_source_certificates_replay_and_bind(self):
        for degree in (4, 8, 12):
            cert = json.loads((HERE/f'SOURCE_PROBE_degree{degree}.json').read_text())
            q, _, _ = profile_projection(degree)
            self.assertTrue(wide.verify_full(cert, expected_q=q, expected_mean_zero=True))
            q[0] += 1
            self.assertFalse(wide.verify_full(cert, expected_q=q, expected_mean_zero=True))
            self.assertFalse(wide.verify_full(cert, expected_mean_zero=False))

    def test_source_budget_24_is_real_rejection(self):
        q, _, _ = profile_projection(24)
        self.assertGreater(max(map(abs, q)), 10**6)
        with self.assertRaisesRegex(ValueError, 'Coefficient magnitude'):
            wide.full_ground(q, modes=4, max_modes=4, bits=12, max_m=2)

    def test_source_open_status_and_angular_tail(self):
        cert = json.loads((HERE/'SOURCE_PROBE_degree12.json').read_text())
        self.assertEqual(cert['status'], 'certified_bound_open_gap')
        self.assertGreater(F(cert['exact_width']), F(cert['tolerance']))
        self.assertTrue(cert['angular_tail_cannot_improve_best_upper'])
        bad = deepcopy(cert)
        bad['angular_tail_lower'] = str(F(cert['angular_tail_lower'])+1)
        self.assertFalse(wide.verify_full(bad))


def singular(axis=2, **changes):
    value = {'kind': 'axis_profile', 'axis': axis, 'profile': 'abs_power',
             'exponent': '-1/4', 'amplitude': '-1', 'offset': '0'}
    value.update(changes)
    return value


def xyz_coefficients(q, axis=2):
    result = {}
    for j, c in enumerate(q):
        if c:
            e = [0, 0, 0]
            e[axis] = j
            result[','.join(map(str, e))] = str(c)
    return result


class IndependentFunctionModels(unittest.TestCase):
    def test_all_axes_full_rodrigues_coefficients_and_residual(self):
        for axis in (0, 1, 2):
            for degree in (4, 8, 12, 24):
                original = singular(axis)
                model = fm.l2_model(original, degree=degree)
                p, cs, error = profile_projection(degree)
                self.assertEqual(model['polynomial'], xyz_coefficients(p, axis))
                self.assertEqual(F(model['error_squared']), error)
                self.assertEqual(F(model['target_norm_squared']), 2)
                self.assertEqual([F(e['projection_coefficient']) for e in model['legendre_projection']], cs)
                self.assertTrue(fm.verify_model(model, expected_function=original,
                                               expected_polynomial=xyz_coefficients(p, axis)))

    def test_step_projection_is_probability_not_area_norm(self):
        original = {'kind': 'axis_profile', 'axis': 1, 'profile': 'step'}
        model = fm.l2_model(original, degree=1)
        self.assertEqual(model['polynomial'], {'0,0,0': '1/2', '0,1,0': '3/4'})
        self.assertEqual(F(model['target_norm_squared']), F(1, 2))
        self.assertEqual(F(model['error_squared']), F(1, 16))
        self.assertEqual(F(model['error_upper']), F(1, 4))

    def test_original_function_axis_exponent_amplitude_offset_binding(self):
        model = fm.l2_model(singular(), degree=4)
        for changed in (singular(axis=0), singular(exponent='-1/5'),
                        singular(amplitude='1'), singular(offset='1')):
            self.assertFalse(fm.verify_model(model, expected_function=changed))
            bad = deepcopy(model)
            bad['function'] = fm.normalize_function(changed)
            self.assertFalse(fm.verify_model(bad))

    def test_projection_and_norm_payload_tampering(self):
        model = fm.l2_model(singular(), degree=8)
        for field in ('error_squared', 'error_upper', 'target_norm_squared', 'projection_norm_squared'):
            bad = deepcopy(model)
            bad[field] = str(F(bad[field])+F(1, 1024))
            self.assertFalse(fm.verify_model(bad))
        bad = deepcopy(model)
        bad['legendre_projection'][-1]['basis_mass'] = '1'
        self.assertFalse(fm.verify_model(bad))
        bad = deepcopy(model)
        bad['norm'] = 'Linf'
        self.assertFalse(fm.verify_model(bad))

    def test_l2_integrability_boundary_and_ae_representative(self):
        for exponent in ('-1/2', '-2'):
            with self.assertRaises(ValueError):
                fm.l2_model(singular(exponent=exponent), degree=4)
        model = fm.l2_model(singular(exponent='-499/1000'), degree=0)
        self.assertTrue(fm.verify_model(model))
        self.assertGreater(F(model['error_upper']), 1)
        self.assertEqual(model['point_value_at_axis_zero'], '0')
        zero_power = fm.l2_model(singular(exponent='0', amplitude='2', offset='3'), degree=4)
        self.assertEqual(zero_power['polynomial'], {'0,0,0': '5'})
        self.assertEqual(zero_power['point_value_at_axis_zero'], '5')
        self.assertEqual(F(zero_power['error_squared']), 0)

    def test_exact_profile_amplitude_offset_and_non_singular_power(self):
        original = singular(axis=0, exponent='1/3', amplitude='-3/2', offset='5/7')
        model = fm.l2_model(original, degree=8)
        p, _, error = profile_projection(8, F(1, 3), F(-3, 2), F(5, 7))
        self.assertEqual(model['polynomial'], xyz_coefficients(p, 0))
        self.assertEqual(F(model['error_squared']), error)

    def test_sqrt_rounding_is_outward_and_tight(self):
        for z in (F(0), F(1), F(2, 9), F(5202, 43681), F(1, 10**30)):
            for bits in (0, 8, 40):
                upper = fm.sqrt_upper(z, bits)
                self.assertGreaterEqual(upper**2, z)
                if upper:
                    self.assertLess((upper-F(1, 1 << bits))**2, z)

    def test_analytic_signed_combination_exact_taylor_coefficients(self):
        original = {'kind': 'analytic_sum', 'terms': [
            {'function': 'exp', 'argument': {'0,0,1': '1/4'}, 'coefficient': '-2'},
            {'function': 'sin', 'argument': {'1,1,0': '1/4'}},
            {'function': 'log1p', 'argument': {'0,0,1': '1/4'}}]}
        model = fm.analytic_model(original, order=2)
        self.assertEqual(model['polynomial'], {'0,0,0': '-2', '0,0,1': '-1/4',
                                               '0,0,2': '-3/32', '1,1,0': '1/4'})
        exp_error = F(1, 4)**3/F(factorial(3))/(1-F(1, 16))
        sin_error = F(1, 4)**3/F(factorial(3))
        log_error = F(1, 4)**3/3/(1-F(1, 4))
        self.assertEqual(F(model['error_upper']), 2*exp_error+sin_error+log_error)
        self.assertTrue(fm.verify_model(model, expected_function=original))

    def test_log_domain_and_large_exp_fallback(self):
        with self.assertRaises(ValueError):
            fm.analytic_model({'kind': 'analytic_sum', 'terms': [
                {'function': 'log1p', 'argument': {'0,0,1': '1'}}]}, order=4)
        model = fm.analytic_model({'kind': 'analytic_sum', 'terms': [
            {'function': 'exp', 'argument': {'0,0,1': '2'}}]}, order=0)
        self.assertEqual(model['polynomial'], {'0,0,0': '1'})
        self.assertGreaterEqual(F(model['error_upper']), 2*3**2)
        self.assertTrue(fm.verify_model(model))

    def test_original_input_unknown_field_float_and_bool_rejection(self):
        for changed in (singular(axis=True), singular(amplitude=0.5),
                        {**singular(), 'unproved_error': '0'}):
            with self.assertRaises(ValueError):
                fm.l2_model(changed)

    def test_raw_log_domain_survives_cancellation_and_zero_coefficient(self):
        zero_model = fm.analytic_model({'kind': 'analytic_sum'})
        for coefficients in (('1', '-1'), ('0',)):
            for argument in ({'0,0,1': '2'}, {'0,0,0': '-1'}):
                original = {'kind': 'analytic_sum', 'terms': [
                    {'function': 'log1p', 'argument': argument, 'coefficient': coefficient}
                    for coefficient in coefficients]}
                with self.subTest(coefficients=coefficients, argument=argument):
                    with self.assertRaises(ValueError):
                        fm.normalize_function(original)
                    with self.assertRaises(ValueError):
                        fm.analytic_model(original)
                    self.assertFalse(fm.verify_model(zero_model, expected_function=original))
        valid_cancellation = {'kind': 'analytic_sum', 'terms': [
            {'function': 'log1p', 'argument': {'0,0,1': '1/2'}, 'coefficient': coefficient}
            for coefficient in ('1', '-1')]}
        model = fm.analytic_model(valid_cancellation)
        self.assertEqual(model['polynomial'], {})
        self.assertEqual(F(model['error_upper']), 0)
        self.assertTrue(fm.verify_model(model, expected_function=valid_cancellation))


class IndependentTransfers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.models = {d: fm.l2_model(singular(), degree=d) for d in (4, 8, 12)}
        cls.sources = {}
        for d in cls.models:
            proof = json.loads((HERE/f'SOURCE_PROBE_degree{d}.json').read_text())
            cls.sources[d] = st._source_record(cls.models[d]['polynomial'], True,
                                              'axis_isometry_wide', proof)
        cls.step = fm.l2_model({'kind': 'axis_profile', 'profile': 'step'}, degree=0)
        cls.step_source = st.polynomial_source(cls.step['polynomial'], mean_zero=True,
                                               modes=2, max_modes=2, bits=8)

    def test_additive_linf_exact_minmax(self):
        original = {'kind': 'analytic_sum', 'terms': [
            {'function': 'exp', 'argument': {'0,0,1': '1/4'}}]}
        model = fm.analytic_model(original, order=0)
        source = st.polynomial_source(model['polynomial'], modes=2, max_modes=2, bits=8)
        cert = st.transfer_linf(model, source, target_width='1/1000')
        self.assertEqual(F(model['error_upper']), F(2, 7))
        self.assertEqual([F(cert[k]) for k in ('lower', 'upper')], [F(19, 7), F(23, 7)])
        self.assertEqual(cert['status'], 'certified_open')
        self.assertTrue(st.verify(cert, expected_function=original, expected_mean_zero=True))

    def test_l2_step_exact_formula(self):
        cert = st.transfer_l2(self.step, self.step_source)
        self.assertEqual([F(cert[k]) for k in ('lower', 'upper')], [F(1), F(4)])
        self.assertTrue(st.verify(cert))

    def test_l2_full_space_and_mean_zero_are_distinct(self):
        source = st.polynomial_source(self.step['polynomial'], mean_zero=False,
                                      modes=2, max_modes=2, bits=8)
        cert = st.transfer_l2(self.step, source)
        self.assertEqual([F(cert[k]) for k in ('lower', 'upper')], [F(0), F(1)])
        self.assertTrue(st.verify(cert, expected_mean_zero=False))
        self.assertFalse(st.verify(cert, expected_mean_zero=True))

    def test_fixed_reference_sphere_l2_cross_terms(self):
        cert = st.fixed_reference({'1,1,0': '1/2', '0,0,2': '1/3'}, {})
        self.assertEqual(F(cert['distance_squared']), F(7, 180))
        self.assertEqual(F(cert['shift']), 1)
        self.assertTrue(st.verify_fixed_reference(cert))

    def test_fixed_reference_original_projection_identity(self):
        for d in (4, 8, 12):
            model, source = self.models[d], self.sources[d]
            cert = st.transfer_fixed(model, source, {'0,0,0': '-4/3'})
            coercivity = cert['coercivity']
            self.assertEqual(F(coercivity['distance_squared']), F(2, 9)-F(model['error_squared']))
            self.assertEqual(F(coercivity['shift']), F(7, 3))
            rho = F(model['error_upper'])/(1-F(coercivity['distance_upper']))
            self.assertEqual(F(cert['lower']), (1-rho)*(F(source['lower'])+F(7, 3))-F(7, 3))
            self.assertEqual(F(cert['upper']), (1+rho)*(F(source['upper'])+F(7, 3))-F(7, 3))
            self.assertTrue(st.verify(cert, expected_function=singular()))

    def test_relative_error_equal_one_is_rejected(self):
        model = fm.l2_model({'kind': 'axis_profile', 'profile': 'step', 'amplitude': '2'}, degree=0)
        source = st.polynomial_source(model['polynomial'], modes=2, max_modes=2, bits=8)
        self.assertEqual(F(model['error_upper']), 1)
        with self.assertRaisesRegex(ValueError, 'strictly below one'):
            st.transfer_l2(model, source)
        with self.assertRaisesRegex(ValueError, 'strictly below one'):
            st.transfer_fixed(self.step, self.step_source, {})
        with self.assertRaisesRegex(ValueError, 'theta < 1'):
            st.fixed_reference({'0,0,0': '1'}, {})

    def test_rotation_uses_same_full_sphere_but_binds_axis(self):
        original_proof = self.sources[4]['source_proof']
        for axis in (0, 1, 2):
            model = fm.l2_model(singular(axis), degree=4)
            source = st._source_record(model['polynomial'], True, 'axis_isometry_wide', original_proof)
            cert = st.transfer_l2(model, source)
            self.assertEqual(source['conversion']['axis'], axis)
            self.assertTrue(st.verify(cert, expected_function=singular(axis)))
            self.assertFalse(st.verify(cert, expected_function=singular((axis+1) % 3)))
            bad = deepcopy(source)
            bad['conversion']['axis'] = (axis+1) % 3
            self.assertFalse(st.verify_source(bad))

    def test_single_sector_cannot_replace_full_spectrum(self):
        with self.assertRaises(ValueError):
            st._source_record(self.models[4]['polynomial'], True, 'axis_isometry_wide',
                              self.sources[4]['source_proof']['sectors'][0])

    def test_cross_polynomial_model_source_pair_is_rejected(self):
        with self.assertRaises(ValueError):
            st.transfer_l2(self.models[4], self.sources[8])
        bad = deepcopy(self.sources[4])
        bad['polynomial'] = self.models[8]['polynomial']
        self.assertFalse(st.verify_source(bad))

    def test_source_model_coercivity_scope_and_status_tampering(self):
        cert = st.transfer_fixed(self.models[4], self.sources[4], {'0,0,0': '-4/3'})
        for path, value in ((('mean_zero',), False), (('scope',), 'finite_Ritz_only'),
                            (('measure',), 'd_sigma'), (('eigenvalue_index',), 2),
                            (('status',), 'target_met'), (('lower',), '100'),
                            (('coercivity', 'shift'), '0'),
                            (('error_evidence', 'embedding', 'energy_coefficient'), '1/2'),
                            (('model', 'function', 'amplitude'), '1')):
            bad = deepcopy(cert)
            cur = bad
            for key in path[:-1]:
                cur = cur[key]
            cur[path[-1]] = value
            self.assertFalse(st.verify(bad), path)
        self.assertFalse(st.verify(cert, expected_width='1/10'))

    def test_open_source_remains_usable_without_claiming_target(self):
        cert = st.transfer_l2(self.models[12], self.sources[12], target_width='1/1000')
        self.assertEqual(cert['source']['source_proof']['status'], 'certified_bound_open_gap')
        self.assertEqual(cert['status'], 'certified_open')
        self.assertGreater(F(cert['exact_width']), F(1, 1000))
        self.assertTrue(st.verify(cert))

    def test_zero_remainder_and_norm_route_separation(self):
        model = fm.analytic_model({'kind': 'analytic_sum', 'polynomial': {'0,0,0': '-3'}})
        source = st.polynomial_source(model['polynomial'], modes=2, max_modes=2, bits=8)
        cert = st.transfer_linf(model, source)
        self.assertEqual((cert['lower'], cert['upper'], cert['status']), ('-1', '-1', 'target_met'))
        with self.assertRaises(ValueError):
            st.transfer_l2(model, source)
        with self.assertRaises(ValueError):
            st.transfer_linf(self.step, self.step_source)

    def test_malformed_verifier_inputs_fail_closed(self):
        for value in (None, [], 'not a certificate', 1, True):
            self.assertFalse(st.verify_source(value))
            self.assertFalse(st.verify(value))
            self.assertFalse(st.verify_range(value))
            self.assertFalse(st.verify_fixed_reference(value))
            self.assertFalse(fm.verify_model(value))
        cert = st.transfer_l2(self.step, self.step_source)
        for field in ('model', 'source', 'coercivity'):
            bad = deepcopy(cert)
            bad[field] = []
            self.assertFalse(st.verify(bad))

    def test_coarse_valid_source_lower_uses_certified_coercive_floor(self):
        model = self.models[4]
        original = self.sources[4]['source_proof']
        sectors = []
        for old in original['sectors']:
            kernel = wide.Kernel(old['q_coefficients'], old['azimuth_m'], True,
                                 old['modes'], original['range_proof'])
            sectors.append(wide.sector_certificate(kernel, 1, F(-100), F(old['upper']),
                                                   old['upper_kind'], independent=True))
        proof = wide.full_certificate(original['q_coefficients'], True, sectors,
                                      original['range_proof'], F(original['tolerance']))
        self.assertEqual(F(proof['lower']), -100)
        source = st._source_record(model['polynomial'], True, 'axis_isometry_wide', proof)
        direct = st.transfer_l2(model, source)
        self.assertEqual(F(direct['error_evidence']['source_shifted_lower_used']), 1)
        self.assertEqual(F(direct['lower']), 1-F(model['error_upper'])-F(direct['error_evidence']['shift']))
        fixed = st.transfer_fixed(model, source, {'0,0,0': '-4/3'})
        factor = F(fixed['coercivity']['coercivity_factor'])
        rho = F(model['error_upper'])/factor
        self.assertEqual(F(fixed['error_evidence']['source_shifted_lower_used']), factor)
        self.assertEqual(F(fixed['lower']), (1-rho)*factor-F(7, 3))
        self.assertTrue(st.verify(direct))
        self.assertTrue(st.verify(fixed))

    def test_non_axis_source_covers_cartesian_full_and_mean_zero_spaces(self):
        p = {'1,0,0': '1/8', '0,1,0': '1/16'}
        for mean_zero in (False, True):
            source = st.polynomial_source(p, mean_zero=mean_zero, L=2, bits=8)
            self.assertEqual(source['backend'], 'cartesian_constraints')
            self.assertTrue(st.verify_source(source, expected_polynomial=p,
                                             expected_mean_zero=mean_zero))
            self.assertFalse(st.verify_source(source, expected_mean_zero=not mean_zero))
            bad = deepcopy(source)
            bad['source_proof']['constraints'] = [{'0,0,1': '1'}]
            self.assertFalse(st.verify_source(bad))

    def test_saved_top_level_certificates_replay(self):
        verifiers = {wide.FULL_FORMAT: wide.verify_full,
                     fm.ANALYTIC_FORMAT: fm.verify_model, fm.L2_FORMAT: fm.verify_model,
                     st.SOURCE: st.verify_source, st.COERCIVITY: st.verify_fixed_reference,
                     st.LINF: st.verify, st.L2: st.verify, st.FIXED: st.verify}
        count = 0
        for path in HERE.rglob('*.json'):
            value = json.loads(path.read_text())
            if isinstance(value, dict) and value.get('format') in verifiers:
                self.assertTrue(verifiers[value['format']](value), str(path))
                count += 1
        self.assertGreaterEqual(count, 3)


if __name__ == '__main__':
    unittest.main()
