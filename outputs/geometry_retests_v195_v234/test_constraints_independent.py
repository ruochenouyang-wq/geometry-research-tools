"""Independent full-sphere finite-moment review and analytic calibration."""
import copy
import json
from pathlib import Path
from fractions import Fraction as F
from unittest.mock import patch
import unittest
import global_constraints as g
from support import old_anisotropic as a, old_constraints

ROOT = Path(__file__).resolve().parent
ONE = {'0,0,0': '1'}
X = {'1,0,0': '1'}
Y = {'0,1,0': '1'}
Z = {'0,0,1': '1'}
XY = {'1,1,0': '1'}
C12 = [ONE, Z]


def _format_nodes(value):
    if isinstance(value, dict):
        if value.get('format') in (g.TRIAL, g.SPECTRAL, g.TRANSFER, g.ADAPTIVE, g.INEQUALITY):
            yield value
        for child in value.values():
            yield from _format_nodes(child)
    elif isinstance(value, list):
        for child in value:
            yield from _format_nodes(child)


def combine(*terms):
    return a.encode(a.add(*(a.polynomial(term) for term in terms)))


class GlobalConstraintIndependentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.original = g.certify({}, C12, L=1, bits=28)
        cls.nontrivial = g.certify(XY, [ONE, combine(X, Z)], L=2, bits=24)

    def test_original_full_sphere_two_and_axisymmetric_six(self):
        self.assertEqual((self.original['lower'], self.original['upper']), ('2', '2'))
        axisymmetric = old_constraints.certify_moment_sector([0], [[1], [0, 1]],
                                                             m=0, modes=3, bits=32)
        self.assertLessEqual(F(axisymmetric['lower']), 6)
        self.assertGreaterEqual(F(axisymmetric['upper']), 6)
        self.assertGreater(F(axisymmetric['lower']), 5)
        self.assertNotEqual(axisymmetric['scope'], self.original['scope'])
        for trial in (X, Y):
            c = g.trial_enclosure({}, C12, trial, support_L=1)
            self.assertEqual(c['upper'], '2')
            self.assertEqual(c['moment_inner_products'], ['0', '0'])

    def test_complete_real_basis_dimension_and_three_linear_modes(self):
        for L in range(5):
            data = g.full_harmonics(L)
            self.assertEqual(data['dimension'], (L+1)**2)
            for l in range(L+1):
                self.assertEqual(sum(b['degree'] == l for b in data['basis']), 2*l+1)
        basis = g.full_harmonics(1)['basis']
        self.assertEqual([b['polynomial'] for b in basis], [ONE, Z, X, Y])
        self.assertEqual(g.full_harmonics(1)['mass'], ['1', '1/3', '1/3', '1/3'])

    def test_single_constant_plus_x_constraint_has_exact_ground_three_halves(self):
        rows = [combine(ONE, X)]
        c = g.certify({}, rows, L=1, bits=32)
        self.assertFalse(c['forms']['constraint_system']['constant_in_row_span'])
        self.assertLessEqual(F(c['lower']), F(3, 2))
        self.assertGreaterEqual(F(c['upper']), F(3, 2))
        self.assertLess(F(c['upper']), 2)
        trial = g.trial_enclosure({}, rows, {'1,0,0': 1, '0,0,0': '-1/3'}, support_L=1)
        self.assertEqual((trial['mass'], trial['energy'], trial['upper']), ('4/9', '2/3', '3/2'))
        self.assertEqual(trial['laplace_lower'], 0)

    def test_cross_cos_sin_m_constraint_has_ground_six_seventeenths(self):
        row = {'0,0,0': 1, '1,0,0': 1, '0,1,0': 2, '0,0,1': 3}
        c = g.certify({}, [row], L=1, bits=32)
        self.assertLessEqual(F(c['lower']), F(6, 17))
        self.assertGreaterEqual(F(c['upper']), F(6, 17))
        self.assertFalse(c['forms']['constraint_system']['constant_in_row_span'])
        trial = g.trial_enclosure({}, [row],
                    {'0,0,0': '-14/3', '1,0,0': 1, '0,1,0': 2, '0,0,1': 3}, support_L=1)
        self.assertEqual(F(trial['upper']), F(6, 17))

    def test_two_rows_jointly_remove_constant(self):
        rows = [combine(ONE, X), {'0,0,0': 1, '1,0,0': -1}]
        system = g.constraint_nullspace(rows, L=1)
        self.assertTrue(system['constant_in_row_span'])
        self.assertEqual(list(map(F, system['constant_span_coefficients'])), [F(1, 2), F(1, 2)])
        self.assertEqual(g.certify({}, rows, L=1)['lower'], '2')

    def test_empty_redundant_and_complete_linear_exclusions(self):
        self.assertEqual(g.certify({}, [], L=1)['lower'], '0')
        self.assertEqual(g.certify({}, [ONE, X, Y, Z], L=2)['lower'], '6')
        self.assertEqual(g.certify({}, [ONE, X, Y], L=1)['lower'], '2')
        redundant = [ONE, Z, {'0,0,0': 2}, {}, {'0,0,1': -3}]
        self.assertEqual(g.certify({}, redundant, L=1)['lower'], '2')

    def test_dense_mass_closed_form_and_root_calibration(self):
        row = {'0,0,0': '1/2', '0,0,1': 1, '0,0,2': '3/2'}  # 1+P1+P2
        forms = g.reduced_forms({}, [row], L=2)
        free = forms['constraint_system']['free_columns']
        iz, ip2 = free.index(1), free.index(4)
        mass = [[F(x) for x in r] for r in forms['reduced_mass']]
        energy = [[F(x) for x in r] for r in forms['reduced_A']]
        self.assertEqual([[mass[i][j] for j in (iz, ip2)] for i in (iz, ip2)],
                         [[F(4, 9), F(1, 15)], [F(1, 15), F(6, 25)]])
        self.assertEqual([[energy[i][j] for j in (iz, ip2)] for i in (iz, ip2)],
                         [[F(2, 3), 0], [0, F(6, 5)]])
        self.assertEqual(g.Kernel({}, [row], L=2).matrix(F(1), 'upper_ritz')[iz][ip2], F(-1, 15))
        c = g.certify({}, [row], L=2, bits=32)
        lo, hi = F(c['lower']), F(c['upper'])
        self.assertTrue(1 < lo <= hi < 2)
        self.assertGreaterEqual(23*lo*lo-156*lo+180, 0)
        self.assertLessEqual(23*hi*hi-156*hi+180, 0)

    def test_complete_constraint_support_rejects_truncation(self):
        p2 = {'0,0,0': '-1/2', '0,0,2': '3/2'}
        with self.assertRaises(ValueError):
            g.constraint_nullspace([p2], L=1)
        with self.assertRaises(ValueError):
            g.certify({}, [p2], L=1)
        result = g.adaptive({}, [p2], start_L=0, max_L=1)
        self.assertEqual(result['status'], 'unsupported_within_budget')
        self.assertTrue(g.verify(result))
        # Raw degree two may represent a constant on S2 and must be accepted.
        system = g.constraint_nullspace([a.encode(a.R2)], L=0)
        self.assertEqual(system['rank'], 1)
        self.assertTrue(system['constant_in_row_span'])

    def test_entire_head_eliminated_requires_growth(self):
        with self.assertRaises(ValueError):
            g.certify({}, [ONE, X, Y, Z], L=1)
        c = g.adaptive({}, [ONE, X, Y, Z], start_L=1, max_L=2)
        self.assertEqual(c['attempts'][0]['status'], 'constraints_eliminate_head')
        self.assertEqual(c['final_certificate']['lower'], '6')
        self.assertTrue(g.verify(c))

    def test_transformed_complete_tail_columns_match_direct_integrals(self):
        rows = [{'0,0,0': 1, '1,0,0': 1, '0,1,0': 2, '0,0,1': 3}]
        forms = g.reduced_forms(XY, rows, L=1)
        system = forms['constraint_system']
        t = [[F(x) for x in row] for row in system['T']]
        basis = [a.polynomial(b['polynomial']) for b in forms['harmonics']['basis']]
        trials = [a.add(*(a.scale(b, t[i][j]) for i, b in enumerate(basis))) for j in range(system['dimension'])]
        tail = g.transformed_tail(XY, rows, L=1)
        self.assertEqual(len(tail['columns']), 5+7)
        for column in tail['columns']:
            b = a.solid_harmonic(column['degree'], column['m'], column['part'])
            self.assertEqual(F(column['mass']), a.inner(b, b))
            self.assertEqual(list(map(F, column['entries'])),
                             [a.inner(b, a.multiply(a.polynomial(XY), u)) for u in trials])
        for u in trials:
            self.assertEqual(a.inner(u, a.polynomial(rows[0])), 0)
            qu = a.multiply(a.polynomial(XY), u)
            omitted_norm = a.inner(qu, qu)-sum(a.inner(qu, b)**2/a.inner(b, b) for b in basis)
            column_index = trials.index(u)
            reconstructed = sum((F(c['entries'][column_index])**2/F(c['mass']) for c in tail['columns']), F(0))
            self.assertEqual(omitted_norm, reconstructed)

    def test_equivalent_rowspaces_transfer_without_research(self):
        target = [combine(ONE, Z), {'0,0,0': 1, '0,0,1': -1}, {'0,0,0': 3}]
        with patch.object(g, 'certify', side_effect=AssertionError('Should replay, not search')):
            c = g.equivalent_constraints(self.original, target)
            self.assertTrue(g.verify(c))
        with self.assertRaises(ValueError):
            g.equivalent_constraints(self.original, [ONE, X])

    def test_scope_inputs_mass_and_tail_tampering_rejected(self):
        for name in ['scope', 'q', 'constraints', 'constant_span', 'mass', 'tail']:
            bad = copy.deepcopy(self.nontrivial)
            if name == 'scope': bad['scope'] = 'axisymmetric_only'
            elif name == 'q': bad['potential'] = X
            elif name == 'constraints': bad['constraints'] = [ONE]
            elif name == 'constant_span': bad['forms']['constraint_system']['constant_in_row_span'] = False
            elif name == 'mass': bad['forms']['reduced_mass'][0][0] = '999'
            else: bad['tail']['columns'].pop()
            self.assertFalse(g.verify(bad), name)
        self.assertFalse(g.verify(self.original, expected_constraints=[ONE, X]))
        self.assertFalse(g.verify(self.original, expected_q=X))

    def test_original_inequality_has_explicit_transverse_refutation(self):
        proof = g.inequality({}, C12, 2, L=1)
        refutation = g.inequality({}, C12, 3, L=1)
        self.assertEqual(proof['status'], 'proved')
        self.assertEqual(refutation['status'], 'refuted_with_trial')
        self.assertEqual(refutation['trial_certificate']['moment_inner_products'], ['0', '0'])
        self.assertEqual(refutation['trial_certificate']['upper'], '2')
        self.assertTrue(g.verify(refutation))
        bad = copy.deepcopy(refutation)
        bad['trial_certificate']['trial'] = Z
        self.assertFalse(g.verify(bad))

    def test_saved_original_proof_and_derived_certificate_replay(self):
        found, original_found = 0, False
        for path in sorted((ROOT/'C12'/'certificates').rglob('*.json')):
            if path.name == 'INDEPENDENT_TESTS.json': continue
            for c in _format_nodes(json.loads(path.read_text())):
                self.assertTrue(g.verify(c), str(path)+':'+c['format'])
                found += 1
                if c['format'] == g.SPECTRAL and c['potential'] == {} and c['constraints'] == C12:
                    self.assertEqual((c['lower'], c['upper']), ('2', '2'))
                    original_found = True
        self.assertGreaterEqual(found, 3)
        self.assertTrue(original_found, 'Must replay the unchanged original C12 problem')


if __name__ == '__main__':
    unittest.main()
