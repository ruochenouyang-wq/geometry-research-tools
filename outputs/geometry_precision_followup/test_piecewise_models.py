import copy
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import subprocess
import sys
import unittest
import piecewise_models as pm


SOURCE = {'kind': 'axis_profile', 'axis': 2, 'profile': 'abs_power',
          'exponent': '-1/4', 'amplitude': '-1', 'offset': '0'}


def independent_t_coefficients(cell):
    a, b = F(cell['left']), F(cell['right'])
    local = list(map(F, cell['polynomial_in_s']))
    out = [F(0)]*len(local)
    for k, coefficient in enumerate(local):
        for j in range(k+1):
            out[j] += coefficient*comb(k, j)*(-a)**(k-j)/(b-a)**k
    return out


class PiecewiseModelTests(unittest.TestCase):
    def test_defaults_and_canonical_exact_input(self):
        self.assertEqual(pm.normalize_function({'kind': 'axis_profile', 'profile': 'abs_power', 'exponent': '-1/4'}),
                         dict(SOURCE, amplitude='1'))
        self.assertTrue(pm.verify_piecewise(pm.piecewise_model(SOURCE), SOURCE))

    def test_omitted_exponent_keeps_old_default_and_is_unsupported(self):
        source = {'kind': 'axis_profile', 'profile': 'abs_power'}
        with self.assertRaises(ValueError): pm.normalize_function(source)
        with self.assertRaises(ValueError): pm.piecewise_model(source)
        self.assertFalse(pm.verify_piecewise(pm.piecewise_model(SOURCE), source))

    def test_unknown_and_wrong_source_rejected(self):
        for bad in ({'kind': 'analytic_sum'}, dict(SOURCE, exponent='-1/3'),
                    dict(SOURCE, extra=0), dict(SOURCE, axis=3), dict(SOURCE, axis=True),
                    dict(SOURCE, amplitude=1.0), dict(SOURCE, offset=False),
                    dict(SOURCE, profile='step')):
            with self.subTest(bad=bad), self.assertRaises(ValueError): pm.normalize_function(bad)

    def test_budget_and_grid_boundary_rejected(self):
        for kwargs in ({'levels': 0}, {'levels': 65}, {'levels': True}, {'degree': 33},
                       {'degree': -1}, {'root_ratio': 0}, {'root_ratio': 1},
                       {'root_ratio': '-1/2'}, {'root_ratio': '1/4097'},
                       {'root_ratio': 0.5}, {'sqrt_bits': 257}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError): pm.piecewise_model(SOURCE, **kwargs)

    def test_step_exact_source_and_probability_mass(self):
        source = {'kind': 'axis_profile', 'profile': 'step', 'amplitude': '-3', 'offset': '2'}
        cert = pm.piecewise_model(source)
        self.assertEqual(cert['error_squared'], '0')
        self.assertEqual(cert['error_upper'], '0')
        self.assertEqual(cert['source_norm_squared'], '5/2')
        self.assertEqual([pm.evaluate(cert, x) for x in (-1, 0, 1)], [2, 2, -1])
        self.assertEqual(cert['expanded_polynomial_coefficient_slots'], 2)
        self.assertEqual(cert['shape_coefficient_count'], 2)

    def test_step_axis_binding(self):
        for axis in range(3):
            source = {'kind': 'axis_profile', 'profile': 'step', 'axis': axis}
            cert = pm.piecewise_model(source)
            self.assertTrue(pm.verify_piecewise(cert, source))
            self.assertFalse(pm.verify_piecewise(cert, dict(source, axis=(axis+1)%3)))

    def test_geometric_fourth_power_cover(self):
        cert = pm.piecewise_model(SOURCE, levels=6)
        cells = cert['cells']
        self.assertEqual(cells[0]['left'], '0')
        self.assertEqual(cells[-1]['right'], '1')
        for i, cell in enumerate(cells):
            self.assertEqual(F(cell['left_fourth_root'])**4, F(cell['left']))
            self.assertEqual(F(cell['right_fourth_root'])**4, F(cell['right']))
            if i: self.assertEqual(cells[i-1]['right'], cell['left'])
            self.assertLess(F(cell['left']), F(cell['right']))

    def test_core_exact_nonzero_error_not_omitted(self):
        cert = pm.piecewise_model(SOURCE, levels=6, degree=3)
        cell = cert['cells'][0]
        root = F(cell['right_fourth_root'])
        self.assertEqual(cell['degree'], 0)
        self.assertEqual(F(cell['polynomial_in_s'][0]), -F(4, 3)/root)
        self.assertEqual(F(cell['error_squared_dt']), F(2, 9)*root**2)
        self.assertGreater(F(cell['error_squared_dt']), 0)
        self.assertGreater(F(cert['error_squared']), F(cell['error_squared_dt']))

    def test_total_source_norm_not_half_or_double(self):
        cert = pm.piecewise_model(SOURCE)
        self.assertEqual(cert['source_norm_squared'], '2')
        self.assertEqual(sum(F(c['source_norm_squared_dt']) for c in cert['cells']), 2)

    def test_exact_rational_source_moments_independent(self):
        cert = pm.piecewise_model(SOURCE, levels=3, degree=4)
        for cell in cert['cells']:
            A, B = F(cell['left_fourth_root']), F(cell['right_fourth_root'])
            for k, value in enumerate(cell['source_monomial_moments_dt']):
                # Substitute t=x^4: -t^k*t^-1/4 dt = -4*x^(4k+2) dx.
                expected = -F(4, 4*k+3)*(B**(4*k+3)-A**(4*k+3))
                self.assertEqual(F(value), expected)

    def test_shape_transport_matches_separate_full_projection(self):
        for amplitude, offset in ((F(-1), F(0)), (F(3, 2), F(7, 5))):
            source = dict(SOURCE, amplitude=str(amplitude), offset=str(offset))
            cert = pm.piecewise_model(source, levels=3, degree=4, root_ratio='2/3')
            for cell in cert['cells'][1:]:
                direct = pm._cell(F(cell['left_fourth_root']), F(cell['right_fourth_root']), 4, amplitude, offset)
                self.assertEqual(cell, direct)

    def test_complete_error_geometric_sum_and_saturation_floor(self):
        for ratio in (F(1, 2), F(3, 4)):
            cert = pm.piecewise_model(SOURCE, levels=6, degree=3, root_ratio=ratio)
            shape_error = F(cert['shape_proof']['error_squared_dt'])
            expected = F(2, 9)*ratio**12 + shape_error*(1-ratio**12)/(1-ratio**2)
            self.assertEqual(F(cert['error_squared']), expected)
            self.assertGreater(shape_error, 0)

    def test_shifted_legendre_rodrigues_independent(self):
        cert = pm.piecewise_model(SOURCE, levels=1, degree=8)
        for row in cert['cells'][-1]['projection']:
            n = row['degree']
            expected = [F((-1)**(n-k)*comb(n, k)*comb(n+k, k)) for k in range(n+1)]
            self.assertEqual(list(map(F, row['shifted_legendre_in_s'])), expected)

    def test_direct_full_squared_residual_independent(self):
        source = dict(SOURCE, amplitude='-3/2', offset='7/5')
        cert = pm.piecewise_model(source, levels=3, degree=3)
        total = F(0)
        for cell in cert['cells']:
            a, b = F(cell['left']), F(cell['right'])
            A, B = F(cell['left_fourth_root']), F(cell['right_fourth_root'])
            p = independent_t_coefficients(cell)
            norm = sum((c*d*(b**(i+j+1)-a**(i+j+1))/F(i+j+1)
                        for i, c in enumerate(p) for j, d in enumerate(p)), F(0))
            cross = sum((c*(F(7, 5)*(b**(k+1)-a**(k+1))/F(k+1)
                            - F(3, 2)*F(4, 4*k+3)*(B**(4*k+3)-A**(4*k+3)))
                         for k, c in enumerate(p)), F(0))
            target = F(49, 25)*(b-a)-F(21, 5)*F(4, 3)*(B**3-A**3)+F(9, 2)*(B**2-A**2)
            error = target-2*cross+norm
            self.assertEqual(error, F(cell['error_squared_dt']))
            total += error
        self.assertEqual(total, F(cert['error_squared']))

    def test_offset_cross_terms_and_error_translation_invariance(self):
        basic = pm.piecewise_model(SOURCE)
        shifted = pm.piecewise_model(dict(SOURCE, offset='9/7'))
        self.assertEqual(basic['error_squared'], shifted['error_squared'])
        self.assertEqual(F(shifted['source_norm_squared']), F(81, 49)-F(24, 7)+2)
        self.assertEqual(pm.evaluate(shifted, F(1, 3))-pm.evaluate(basic, F(1, 3)), F(9, 7))

    def test_signed_amplitude_error_scaling(self):
        basic = pm.piecewise_model(SOURCE)
        changed = pm.piecewise_model(dict(SOURCE, amplitude='3/2'))
        self.assertEqual(F(changed['error_squared']), F(9, 4)*F(basic['error_squared']))

    def test_zero_amplitude_constant_exact(self):
        cert = pm.piecewise_model(dict(SOURCE, amplitude='0', offset='7/3'))
        self.assertEqual(cert['error_squared'], '0')
        self.assertEqual(pm.evaluate(cert, 0), F(7, 3))
        self.assertEqual(pm.evaluate(cert, 1), F(7, 3))

    def test_even_representation_and_core_value_convention(self):
        cert = pm.piecewise_model(SOURCE)
        for t in (F(1), F(1, 16), F(1, 10000000)):
            self.assertEqual(pm.evaluate(cert, t), pm.evaluate(cert, -t))
        self.assertEqual(cert['source_value_at_axis_zero'], '0')
        self.assertNotEqual(pm.evaluate(cert, 0), 0)  # source and approximant may differ on a null set.

    def test_approximation_defined_at_all_seams_and_endpoints(self):
        cert = pm.piecewise_model(SOURCE)
        for cell in cert['cells']:
            t = F(cell['left'])
            self.assertEqual(pm.evaluate(cert, t), F(cell['polynomial_in_s'][0]))
        self.assertIsInstance(pm.evaluate(cert, 1), F)
        for bad in (2, '-3/2', 0.5):
            with self.assertRaises(ValueError): pm.evaluate(cert, bad)

    def test_degree_improvement_for_fixed_cover(self):
        errors = [F(pm.piecewise_model(SOURCE, levels=3, degree=d)['error_squared']) for d in (0, 1, 2, 3, 4)]
        self.assertTrue(all(b < a for a, b in zip(errors, errors[1:])))

    def test_refining_core_improves_whole_domain(self):
        errors = [F(pm.piecewise_model(SOURCE, levels=j, degree=3)['error_squared']) for j in (1, 2, 4, 6, 8)]
        self.assertTrue(all(b < a for a, b in zip(errors, errors[1:])))

    def test_sqrt_outward_minimal_dyadic(self):
        for value in (F(0), F(1, 16), F(2), F(1, 10**50)):
            bound = pm.sqrt_upper(value)
            self.assertGreaterEqual(bound**2, value)
            if bound: self.assertLess((bound-F(1, 2**40))**2, value)

    def test_same_25_slot_error_beats_frozen_global_projection(self):
        cert = pm.piecewise_model(SOURCE, levels=6, degree=3)
        old_error = F(1335263295935522, 24084192784943929)
        self.assertEqual(cert['expanded_polynomial_coefficient_slots'], 25)
        self.assertEqual(cert['shape_coefficient_count'], 4)
        self.assertLess(F(cert['error_squared']), old_error/64)
        self.assertGreater(F(cert['error_squared']), F(1, 10**12))  # no 1e-6 error target attained.

    def test_predeclared_high_precision_first_attempt_not_target(self):
        cert = pm.piecewise_model(SOURCE, levels=32, degree=24, root_ratio='1/2')
        self.assertTrue(pm.verify_piecewise(cert, SOURCE))
        self.assertEqual(F(cert['error_upper']), F(75073, 549755813888))
        self.assertGreater(F(cert['error_upper']), F(1, 10**8))
        self.assertEqual(cert['expanded_polynomial_coefficient_slots'], 801)

    def test_predeclared_high_precision_fallback_meets_target(self):
        cert = pm.piecewise_model(SOURCE, levels=48, degree=20, root_ratio='2/3')
        self.assertTrue(pm.verify_piecewise(cert, SOURCE))
        self.assertEqual(F(cert['error_upper']), F(1833, 1099511627776))
        self.assertLess(F(cert['error_upper']), F(1, 10**8))
        self.assertGreater(F(cert['error_squared']), 0)
        self.assertEqual(cert['expanded_polynomial_coefficient_slots'], 1009)
        self.assertFalse(cert['spectral_transfer_claimed'])

    def test_13_active_slot_comparison(self):
        cert = pm.piecewise_model(SOURCE, levels=4, degree=2)
        self.assertEqual(cert['expanded_polynomial_coefficient_slots'], 13)
        self.assertLess(F(cert['error_squared']), F(1335263295935522, 24084192784943929)/9)

    def test_actual_frozen_degree24_baseline_replay(self):
        previous = Path(__file__).resolve().parent.parent/'geometry_general_v235_v237'
        program = "import json;import function_models as f;c=f.l2_model("+repr(SOURCE)+",degree=24);print(json.dumps({'error':c['error_squared'],'nonzero':len(c['polynomial']),'verified':f.verify_model(c)}))"
        completed = subprocess.run([sys.executable, '-B', '-c', program], cwd=previous,
                                   text=True, capture_output=True, check=True)
        actual = json.loads(completed.stdout)
        self.assertEqual(actual, {'error': '1335263295935522/24084192784943929', 'nonzero': 13, 'verified': True})

    def test_expected_source_binding_all_parameters(self):
        cert = pm.piecewise_model(SOURCE)
        for field, value in (('axis', 0), ('amplitude', '1'), ('offset', '1'), ('exponent', '-1/3')):
            self.assertFalse(pm.verify_piecewise(cert, dict(SOURCE, **{field: value})))

    def test_tampered_error_tail_cell_cover_and_unknown_rejected(self):
        cert = pm.piecewise_model(SOURCE)
        changes = [('error_squared', '0'), ('error_upper', '0'), ('norm', 'Linf'),
                   ('scope', 'sphere_area'), ('levels', 5), ('root_ratio', '3/4'),
                   ('global_polynomial', True), ('unknown', 0)]
        for field, value in changes:
            bad = copy.deepcopy(cert); bad[field] = value
            self.assertFalse(pm.verify_piecewise(bad), field)
        for field, value in (('left', '1/8'), ('right_fourth_root', '1/7'), ('error_squared_dt', '0')):
            bad = copy.deepcopy(cert); bad['cells'][0][field] = value
            self.assertFalse(pm.verify_piecewise(bad), field)
        bad = copy.deepcopy(cert); bad['cells'].pop(0); self.assertFalse(pm.verify_piecewise(bad))
        bad = copy.deepcopy(cert); bad['cells'][-1]['projection'].pop(); self.assertFalse(pm.verify_piecewise(bad))
        bad = copy.deepcopy(cert); bad['shape_proof']['error_squared_dt'] = '0'; self.assertFalse(pm.verify_piecewise(bad))

    def test_json_roundtrip_replay_and_grid_variants(self):
        for ratio in ('1/2', '3/4', '7/8'):
            cert = pm.piecewise_model(SOURCE, levels=6, degree=3, root_ratio=ratio)
            self.assertTrue(pm.verify_piecewise(json.loads(json.dumps(cert)), SOURCE))
            self.assertGreater(F(cert['error_squared']), 0)

    def test_reproducible_budget_table_no_spectral_claim(self):
        rows = pm.budget_table(budgets=((2, 2), (6, 3)))
        self.assertEqual([r['coefficient_slots'] for r in rows], [7, 25])
        self.assertEqual(rows, pm.budget_table(budgets=((2, 2), (6, 3))))
        cert = pm.piecewise_model(SOURCE)
        self.assertFalse(cert['spectral_transfer_claimed'])
        self.assertFalse(cert['global_polynomial'])

    def test_saved_note_recipes_replay_exact_values_and_hashes(self):
        import hashlib
        note = (Path(__file__).resolve().parent/'PIECEWISE_NOTES.md').read_text()
        block = note.split('<!-- SAVED_RECIPES_BEGIN -->\n```json\n', 1)[1].split('\n```\n<!-- SAVED_RECIPES_END -->', 1)[0]
        recipes = json.loads(block)
        for row in recipes:
            cert = pm.piecewise_model(row['function'], **row['parameters'])
            self.assertTrue(pm.verify_piecewise(cert, row['function']))
            self.assertEqual(cert['error_squared'], row['error_squared'])
            self.assertEqual(cert['error_upper'], row['error_upper'])
            serialized = pm._canonical(cert).encode()
            self.assertEqual(len(serialized), row['certificate_bytes'])
            self.assertEqual(hashlib.sha256(serialized).hexdigest(), row['certificate_sha256'])


if __name__ == '__main__': unittest.main()
