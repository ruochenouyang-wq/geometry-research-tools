"""Independent original-moment tests; recurrence versus implementation Rodrigues.

All expected matrix entries are computed from the three-term Legendre recurrence,
repeated differentiation, repeated multiplication by (1-t^2), and direct moments.
No expected value calls direct_moments' basis, product, moment or assembly helpers.
"""
from copy import deepcopy
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import sys
import unittest

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import direct_moments as d
import enriched_trial as enriched
from bridge import base


def profile(name='step', axis=2, amplitude='1', offset='0', exponent='-1/4'):
    q = {'kind': 'axis_profile', 'axis': axis, 'profile': name,
         'amplitude': amplitude, 'offset': offset}
    if name == 'abs_power':
        q['exponent'] = exponent
    return q


STEP = profile()
SINGULAR = profile('abs_power', amplitude='-1')


def plus(p, q):
    return [(p[k] if k < len(p) else F(0))+(q[k] if k < len(q) else F(0))
            for k in range(max(len(p), len(q)))]


def scale(p, c):
    return [c*x for x in p]


def times(p, q):
    out = [F(0)]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return out


def recurrence(l):
    previous, current = [F(1)], [F(0), F(1)]
    if l == 0:
        return previous
    for k in range(1, l):
        previous, current = current, scale(plus(scale([F(0)]+current, 2*k+1),
                                                    scale(previous, -k)), F(1, k+1))
    return current


def derivative(p, m):
    for _ in range(m):
        p = [j*p[j] for j in range(1, len(p))]
    return p


def independent_moment(q, k, power):
    amplitude, offset = F(q['amplitude']), F(q['offset'])
    if q['profile'] == 'step':
        return ((offset+amplitude)**power+(-1)**k*offset**power)/F(k+1)
    if k % 2:
        return F(0)
    exponent = F(q['exponent'])
    return sum((F(2*comb(power, r))*amplitude**r*offset**(power-r)/
                (k+r*exponent+1) for r in range(power+1)), F(0))


def independent_product(q, m, l, k, power):
    product = times(derivative(recurrence(l), m), derivative(recurrence(k), m))
    for _ in range(m):
        product = times(product, [F(1), F(0), F(-1)])
    return sum((c*independent_moment(q, j, power) for j, c in enumerate(product)), F(0))


def independent_assembly(q, m, mean_zero, modes, near_tail=0):
    start = 1 if m == 0 and mean_zero else m
    degrees = list(range(start, start+modes))
    mass = [independent_product(q, m, l, l, 0) for l in degrees]
    v = [[independent_product(q, m, l, k, 1) for k in degrees] for l in degrees]
    t = [[independent_product(q, m, l, k, 2) for k in degrees] for l in degrees]
    a = [[v[i][j]+(degrees[i]*(degrees[i]+1)*mass[i] if i == j else 0)
          for j in range(modes)] for i in range(modes)]
    projected = list(range(start, start+modes+near_tail))
    if m == 0 and mean_zero:
        projected = [0]+projected
    c = [[t[i][j]-sum((independent_product(q, m, l, k, 1)*
                       independent_product(q, m, h, k, 1)/
                       independent_product(q, m, k, k, 0) for k in projected), F(0))
          for j, h in enumerate(degrees)] for i, l in enumerate(degrees)]
    return {'degrees': degrees, 'mass': mass, 'A': a, 'V': v, 'T': t, 'C': c}


class IndependentMatrices(unittest.TestCase):
    def test_raw_dt_moments_with_amplitude_offset_cross_terms(self):
        cases = [STEP, SINGULAR, profile(amplitude='-3/2', offset='5/7'),
                 profile('abs_power', amplitude='2/3', offset='-1/5', exponent='1/3')]
        for q in cases:
            for power in (0, 1, 2):
                for k in range(9):
                    self.assertEqual(d.moment(q, k, power), independent_moment(q, k, power))

    def test_mass_and_matrices_from_independent_recurrence(self):
        cases = [STEP, SINGULAR,
                 profile('abs_power', amplitude='2/3', offset='-1/5', exponent='1/3')]
        for q in cases:
            for m, mean_zero in ((0, False), (0, True), (1, True), (2, True)):
                expected = independent_assembly(q, m, mean_zero, 3)
                actual = d.matrix_assembly(q, m=m, mean_zero=mean_zero, modes=3)
                for field in expected:
                    self.assertEqual(actual[field], expected[field], (q, m, mean_zero, field))

    def test_five_exact_scalar_calibrations(self):
        rows = [(STEP, 0, True, F(2, 3), F(1, 3), F(1, 3), F(1, 24), F(5, 2)),
                (STEP, 0, False, F(2), F(1), F(1), F(1, 2), F(1, 2)),
                (STEP, 1, True, F(4, 3), F(2, 3), F(2, 3), F(1, 3), F(5, 2)),
                (SINGULAR, 0, True, F(2, 3), F(-8, 11), F(4, 5), F(4, 605), F(10, 11)),
                (SINGULAR, 1, True, F(4, 3), F(-64, 33), F(16, 5), F(688, 1815), F(6, 11))]
        for q, m, zero, mass, v, t, c, quotient in rows:
            actual = d.matrix_assembly(q, m=m, mean_zero=zero, modes=1)
            self.assertEqual((actual['mass'][0], actual['V'][0][0], actual['T'][0][0],
                              actual['C'][0][0], actual['A'][0][0]/actual['mass'][0]),
                             (mass, v, t, c, quotient))

    def test_mean_zero_removes_constant_exactly_once(self):
        compressed = d.matrix_assembly(STEP, m=0, mean_zero=True, modes=1)
        self.assertEqual(compressed['removed_constant'], {'mass': F(2), 'column': [F(1, 2)]})
        self.assertEqual(compressed['C'][0][0]+F(1, 8), F(1, 6))
        self.assertIsNone(d.matrix_assembly(STEP, m=0, mean_zero=False, modes=1)['removed_constant'])
        self.assertIsNone(d.matrix_assembly(STEP, m=1, mean_zero=True, modes=1)['removed_constant'])

    def test_constant_potential_has_zero_head_tail_gram(self):
        q = profile(amplitude='0', offset='-3/7')
        for m in (0, 1, 2):
            actual = d.matrix_assembly(q, m=m, modes=3, near_tail=3)
            self.assertEqual(actual['C'], [[F(0)]*3 for _ in range(3)])
            self.assertTrue(all(not any(t['column']) for t in actual['near_tail']))
            for i, l in enumerate(actual['degrees']):
                self.assertEqual(actual['A'][i][i]/actual['mass'][i], l*(l+1)-F(3, 7))

    def test_near_tail_remainder_matches_complete_projection(self):
        for q in (STEP, SINGULAR):
            for m in (0, 1, 2):
                expected = independent_assembly(q, m, True, 2, near_tail=3)
                actual = d.matrix_assembly(q, m=m, modes=2, near_tail=3)
                self.assertEqual(actual['C'], expected['C'])
                self.assertEqual(base.inertia(actual['C'])[0], 0)
                for entry in actual['near_tail']:
                    k = entry['degree']
                    self.assertEqual(entry['mass'], independent_product(q, m, k, k, 0))
                    self.assertEqual(entry['column'], [independent_product(q, m, l, k, 1)
                                                      for l in actual['degrees']])

    def test_axis_permutation_preserves_radial_data(self):
        expected = d.matrix_assembly(SINGULAR, m=1, modes=3, near_tail=2)
        for axis in (0, 1, 2):
            self.assertEqual(d.matrix_assembly({**SINGULAR, 'axis': axis}, m=1, modes=3, near_tail=2), expected)

    def test_form_bound_uses_probability_variance(self):
        form = d.form_bound(SINGULAR)
        self.assertEqual(F(form['center']), F(-4, 3))
        self.assertEqual(F(form['variance']), F(2, 9))
        self.assertGreaterEqual(F(form['eta_upper'])**2, F(2, 9))
        self.assertEqual(F(form['mass_offset']), F(-4, 3)-F(form['eta_upper']))
        step = d.form_bound(profile(amplitude='-3', offset='1'))
        self.assertEqual((F(step['energy_factor']), F(step['mass_offset'])), (F(1), F(-2)))

    def test_form_integrability_and_relative_bound_failures(self):
        for exponent in ('-1/2', '-1'):
            with self.assertRaises(ValueError):
                d.matrix_assembly(profile('abs_power', exponent=exponent), modes=2)
        q = profile('abs_power', exponent='-499/1000', amplitude='-1')
        self.assertEqual(len(d.matrix_assembly(q, modes=1)['A']), 1)
        with self.assertRaisesRegex(ValueError, 'eta < 1'):
            d.Kernel(q, modes=1)

    def test_input_boundaries_and_exact_rational_only(self):
        for value in ({}, [], None, {'kind': 'analytic_sum'}):
            with self.assertRaises(ValueError):
                d.normalize(value)
        for kwargs in ({'m': True}, {'m': 13}, {'modes': 0}, {'modes': 33},
                       {'near_tail': -1}, {'mean_zero': 1}):
            with self.assertRaises(ValueError):
                d.matrix_assembly(STEP, **kwargs)
        for change in ({'amplitude': 0.5}, {'axis': True}, {'unexpected': 'value'}):
            with self.assertRaises(ValueError):
                d.normalize({**STEP, **change})


class IndependentSpectrum(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.step = d.full_ground(STEP, modes=3, bits=14, max_m=2, tolerance='1/100000000')
        cls.singular = d.full_ground(SINGULAR, modes=3, bits=14, max_m=2,
                                     tolerance='1/100000000', near_tail=2)

    def test_strict_tail_denominator_and_exact_zero_inertia(self):
        kernel = d.Kernel(STEP, m=1, modes=2, near_tail=2)
        with self.assertRaises(ValueError):
            kernel.matrix(kernel.beta)
        with self.assertRaises(ValueError):
            kernel.matrix(kernel.beta+1)
        self.assertEqual(len(kernel.matrix(kernel.beta, 'upper')), 2)
        q = profile(amplitude='0', offset='-3/7')
        kernel = d.Kernel(q, m=1, modes=2)
        eigenvalue = F(11, 7)
        cert = d.sector_certificate(kernel, eigenvalue, eigenvalue)
        self.assertTrue(d.verify_sector(cert))
        self.assertGreater(cert['lower_inertia'][1], 0)
        self.assertGreater(cert['upper_inertia'][1], 0)

    def test_schur_matrix_against_independent_formula(self):
        q, m, modes, near, x = SINGULAR, 1, 2, 3, F(0)
        expected = independent_assembly(q, m, True, modes, near)
        kernel = d.Kernel(q, m=m, modes=modes, near_tail=near)
        energy = F(kernel.form['energy_factor'])
        offset = F(kernel.form['mass_offset'])
        d_k = lambda k: energy*k*(k+1)+offset
        start = expected['degrees'][-1]+1
        matrix = [[expected['A'][i][j]-(x*expected['mass'][i] if i == j else 0)
                   -expected['C'][i][j]/(d_k(start+near)-x)
                   -sum((independent_product(q, m, l, k, 1)*independent_product(q, m, h, k, 1)/
                         (independent_product(q, m, k, k, 0)*(d_k(k)-x))
                         for k in range(start, start+near)), F(0))
                   for j, h in enumerate(expected['degrees'])] for i, l in enumerate(expected['degrees'])]
        self.assertEqual(kernel.matrix(x), matrix)

    def test_near_tail_comparison_is_monotone_in_psd_order(self):
        for q in (STEP, SINGULAR):
            for m in (0, 1):
                coarse = d.Kernel(q, m=m, modes=2)
                refined = d.Kernel(q, m=m, modes=2, near_tail=3)
                x = F(-1)
                a, b = coarse.matrix(x), refined.matrix(x)
                difference = [[b[i][j]-a[i][j] for j in range(2)] for i in range(2)]
                self.assertEqual(base.inertia(difference)[0], 0)

    def test_complete_sphere_includes_m1_and_remaining_angular_tail(self):
        for cert, original, trial_upper in ((self.step, STEP, F(5, 2)),
                                             (self.singular, SINGULAR, F(6, 11))):
            self.assertTrue(d.verify_full(cert, expected_function=original, expected_mean_zero=True))
            self.assertEqual([s['azimuth_m'] for s in cert['sectors']], [0, 1])
            self.assertTrue(cert['angular_tail_cannot_improve_best_upper'])
            self.assertLessEqual(F(cert['upper']), trial_upper)
            self.assertGreaterEqual(F(cert['angular_tail_lower']), F(cert['upper']))
            self.assertEqual(cert['status'], 'certified_open')

    def test_uncomputed_m1_is_retained_in_global_open_interval(self):
        cert = d.full_ground(STEP, modes=2, bits=10, max_m=0, tolerance='1/100000000')
        self.assertEqual(len(cert['sectors']), 1)
        self.assertEqual(F(cert['angular_tail_lower']), 2)
        self.assertEqual(F(cert['lower']), 2)
        self.assertFalse(cert['angular_tail_cannot_improve_best_upper'])
        self.assertEqual(cert['status'], 'certified_open')
        self.assertTrue(d.verify_full(cert))

    def test_full_space_constant_mode_differs_from_mean_zero(self):
        q = profile(amplitude='0', offset='-3/7')
        full = d.full_ground(q, mean_zero=False, modes=1, bits=10, max_m=1)
        zero = d.full_ground(q, mean_zero=True, modes=1, bits=10, max_m=1)
        self.assertEqual(F(full['upper']), F(-3, 7))
        self.assertEqual(F(zero['upper']), F(11, 7))
        self.assertTrue(d.verify_full(full, expected_mean_zero=False))
        self.assertFalse(d.verify_full(full, expected_mean_zero=True))

    def test_original_function_and_norm_scope_binding(self):
        for changed in ({**SINGULAR, 'axis': 0}, {**SINGULAR, 'amplitude': '1'},
                        {**SINGULAR, 'offset': '1'}, {**SINGULAR, 'exponent': '-1/5'}):
            self.assertFalse(d.verify_full(self.singular, expected_function=changed))
            bad = deepcopy(self.singular)
            bad['function'] = changed
            self.assertFalse(d.verify_full(bad))
        self.assertFalse(d.verify_full(self.singular, expected_tolerance='1/10'))

    def test_sector_constant_far_gram_and_tail_payload_tampering(self):
        sector = self.singular['sectors'][0]
        for path, value in ((('kernel_evidence', 'C', 0, 0), '0'),
                            (('kernel_evidence', 'removed_constant', 'mass'), '1'),
                            (('kernel_evidence', 'near_tail', 0, 'degree'), 50),
                            (('kernel_evidence', 'mass', 0), '1'),
                            (('radial_tail_lower',), '100000'),
                            (('remainder_tail_lower',), '100000'),
                            (('function_approximation_error',), '1/100'),
                            (('projection',), 'none')):
            bad = deepcopy(sector)
            cursor = bad
            for key in path[:-1]:
                cursor = cursor[key]
            cursor[path[-1]] = value
            self.assertFalse(d.verify_sector(bad), path)

    def test_missing_duplicated_sectors_and_forged_global_status(self):
        for sectors in ([], [self.step['sectors'][1]], [self.step['sectors'][0]]*2):
            with self.assertRaises(ValueError):
                d.full_certificate(STEP, True, sectors)
        for field, value in (('angular_tail_lower', '1000'), ('omitted_azimuth_m_start', 3),
                             ('status', 'target_met'), ('mean_zero', False),
                             ('measure', 'd_sigma'), ('eigenvalue_index', 2),
                             ('scope', 'finite_Ritz_only'), ('upper', '-100')):
            bad = deepcopy(self.step)
            bad[field] = value
            self.assertFalse(d.verify_full(bad), field)

    def test_malformed_verifiers_and_sector_not_global(self):
        for bad in (None, [], {}, 1, True, 'certificate'):
            self.assertFalse(d.verify_sector(bad))
            self.assertFalse(d.verify_full(bad))
        self.assertFalse(d.verify_full(self.step['sectors'][0]))
        self.assertFalse(d.verify_sector(self.step))

    def test_budget_exhaustion_keeps_valid_original_function_interval(self):
        cert = d.full_ground(SINGULAR, modes=1, bits=8, max_m=0, tolerance='1/1000000000000')
        self.assertTrue(d.verify_full(cert, expected_function=SINGULAR))
        self.assertEqual(cert['status'], 'certified_open')
        self.assertEqual(cert['function_approximation_error'], '0')
        self.assertGreater(F(cert['exact_width']), F(cert['tolerance']))

    def test_saved_original_function_certificates_replay(self):
        def visit(value, path):
            if isinstance(value, dict):
                if value.get('format') == d.FULL:
                    self.assertTrue(d.verify_full(value), str(path))
                elif value.get('format') == d.SECTOR:
                    self.assertTrue(d.verify_sector(value), str(path))
                for child in value.values():
                    visit(child, path)
            elif isinstance(value, list):
                for child in value:
                    visit(child, path)
        for path in HERE.rglob('*.json'):
            visit(json.loads(path.read_text()), path)

    def test_enriched_weak_energy_and_complete_operator_gram(self):
        powers = [F(0), F(8, 5), F(7, 4), F(2), F(13, 5)]
        moment = lambda r: F(4)/((r+1)*(r+3))
        def weak_energy(s, t):
            r = s+t
            mixed = (2*s*t*(1/(r-1)-2/(r+1)+1/(r+3))) if s*t else F(0)
            return mixed-r*moment(r)+2/(r+3)+2/(r+1)
        def strong(s):
            result = [(s, (s+1)*(s+2)), (s-F(1, 4), F(-1))]
            if s:
                result.append((s-2, -s*(s-1)))
            return result
        mass, energy, norm = enriched.trial_matrices(powers)
        for i, s in enumerate(powers):
            for j, t in enumerate(powers):
                self.assertEqual(mass[i][j], moment(s+t))
                self.assertEqual(energy[i][j], weak_energy(s, t)-moment(s+t-F(1, 4)))
                expected_norm = sum((a*b*moment(r+v) for r, a in strong(s)
                                     for v, b in strong(t)), F(0))
                self.assertEqual(norm[i][j], expected_norm)
        action = {}
        for s, c in ((F(0), F(1)), (F(7, 4), F(-16, 21))):
            for r, a in enriched.action_terms(s).items():
                action[r] = action.get(r, F(0))+c*a
        self.assertEqual({r: a for r, a in action.items() if a},
                         {F(0): F(2), F(3, 2): F(16, 21), F(7, 4): F(-55, 7)})
        self.assertNotIn(F(-2), enriched.action_terms(0))

    def test_enriched_domain_temple_gap_and_nonprojected_residual(self):
        for powers in ([], [F(3, 2)], [1], [0, 0], [2, 0]):
            with self.assertRaises(ValueError):
                enriched.trial_matrices(powers)
        with self.assertRaises(ValueError):
            enriched.weighted_moment(-1)
        with self.assertRaises(ValueError):
            enriched.trial_statistics([0], [0])
        stats = enriched.trial_statistics([0], [1])
        self.assertEqual(F(stats['mass']), F(4, 3))
        self.assertEqual(F(stats['rayleigh']), F(6, 11))
        self.assertEqual(F(stats['operator_norm_squared']), F(128, 165))
        self.assertEqual(F(stats['residual_squared']), F(172, 605))
        # The residual of a one-dimensional projected operator would be zero.
        self.assertGreater(F(stats['residual_squared']), 0)
        proof = enriched.second_eigenvalue_bound()
        eta = F(proof['form_bound']['eta_upper'])
        self.assertEqual(F(proof['lower']), (1-eta)*6-F(4, 3)-eta)
        self.assertEqual(proof['real_component'], 'cos_phi_only')
        self.assertEqual(proof['radial_parities'], 'both_even_and_odd')
        self.assertTrue(proof['not_the_second_eigenvalue_of_cos_plus_sin'])
        self.assertGreater(F(enriched.trial_statistics([2], [1])['rayleigh']), F(proof['lower']))
        with self.assertRaisesRegex(ValueError, 'Temple requires'):
            enriched.trial_certificate([2], [1], self.singular)

    def test_enriched_global_binding_rounding_and_saved_evidence(self):
        proposal = enriched.propose_trial([0, F(7, 4), 2], precision_bits=32, iterations=2)
        cert = enriched.trial_certificate(proposal['powers'], proposal['coefficients'], self.singular)
        self.assertTrue(enriched.verify(cert, expected_source=self.singular))
        self.assertFalse(enriched.verify(cert, expected_source=self.step))
        lower = min(F(cert['m1_lower']), F(cert['angular_tail_lower']),
                    *(F(s['lower']) for s in self.singular['sectors'] if s['azimuth_m'] != 1))
        self.assertEqual(F(cert['lower']), lower)
        self.assertEqual(F(cert['upper']), min(F(cert['statistics']['rayleigh']), F(self.singular['upper'])))
        for path, value in ((('statistics', 'residual_squared'), '0'),
                            (('matrices', 'R', 0, 0), '0'),
                            (('second_eigenvalue_proof', 'eigenvalue_index'), 1),
                            (('second_eigenvalue_proof', 'real_component'), 'cos_plus_sin'),
                            (('second_eigenvalue_proof', 'lower'), cert['statistics']['rayleigh']),
                            (('function', 'amplitude'), '1'), (('function', 'axis'), 0),
                            (('angular_tail_lower',), '1000'), (('mean_zero',), False),
                            (('source_digest',), 'unbound'), (('status',), 'target_met')):
            bad = deepcopy(cert)
            cursor = bad
            for key in path[:-1]:
                cursor = cursor[key]
            cursor[path[-1]] = value
            self.assertFalse(enriched.verify(bad), path)
        def visit(value, path):
            if isinstance(value, dict):
                if value.get('format') == enriched.FORMAT:
                    self.assertTrue(enriched.verify(value), str(path))
                    trial = value['trial']
                    self.assertEqual(enriched.trial_statistics(trial['powers'], trial['coefficients']),
                                     value['statistics'])
                for child in value.values():
                    visit(child, path)
            elif isinstance(value, list):
                for child in value:
                    visit(child, path)
        for path in HERE.rglob('*.json'):
            visit(json.loads(path.read_text()), path)


if __name__ == '__main__':
    unittest.main()
