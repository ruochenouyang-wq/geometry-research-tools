"""Independent reflection-block completeness, exact moments and scope review."""
import copy
import json
from pathlib import Path
from fractions import Fraction as F
from unittest.mock import patch
import unittest
import parity_spectrum as p
from support import old_anisotropic as a

ROOT = Path(__file__).resolve().parent
Q = {'2,0,0': '1', '0,2,0': '2', '0,0,2': '3'}
XY = {'1,1,0': '1'}


def _certificate_nodes(value):
    if isinstance(value, dict):
        if (value.get('format') in (p.BLOCK_FORMAT, p.FULL_FORMAT, p.DRIVER_FORMAT, a.FORMAT)
                and 'potential' in value
                and any(key in value for key in ('finite_form','characters','certificate'))):
            yield value
        for child in value.values():
            yield from _certificate_nodes(child)
    elif isinstance(value, list):
        for child in value:
            yield from _certificate_nodes(child)


class ParityIndependentTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.block = p.block_ground(Q, (1, 0, 0), L=3, bits=32)
        cls.original = p.adaptive_ground(Q, target=F(1, 10**8), bits=40)

    def test_simplex_range_and_original_independent_guard(self):
        c = p.potential_range(Q)
        self.assertEqual((c['lower'], c['upper']), ('1', '3'))
        self.assertTrue(p.verify_range(c, expected_q=Q))
        self.assertEqual(F(a.rayleigh(Q, {'1,0,0': 1})['upper']), F(18, 5))
        final = self.original['certificate']
        self.assertTrue(3 <= F(final['lower']) <= F(final['upper']) <= F(18, 5))
        bad = copy.deepcopy(c)
        bad['lower'] = '2'
        self.assertFalse(p.verify_range(bad))

    def test_eight_character_projectors_sum_to_identity(self):
        poly = a.polynomial({'0,0,0': 1, '1,0,0': 2, '0,1,1': -3,
                             '1,1,1': '1/7', '2,0,2': '2/5'})
        pieces = [p.parity_project(poly, character) for character in p.PARITIES]
        self.assertEqual(a.add(*pieces), poly)
        for i, character in enumerate(p.PARITIES):
            self.assertEqual(p.parity_project(pieces[i], character), pieces[i])
            for j in range(i):
                self.assertEqual(a.inner(pieces[i], pieces[j]), 0)

    def test_all_harmonic_dimensions_and_characters_through_degree_eight(self):
        for degree in range(9):
            total = 0
            for character in p.PARITIES:
                k = sum(character)
                expected = (degree-k)//2+1 if degree >= k and (degree-k) % 2 == 0 else 0
                basis = p.degree_basis(character, degree)
                self.assertEqual(len(basis), expected)
                total += len(basis)
                for b in basis:
                    self.assertTrue(all(tuple(v % 2 for v in e) == character for e in b['polynomial']))
                    self.assertTrue(all(sum(e) == degree for e in b['polynomial']))
                    self.assertGreater(b['mass'], 0)
            self.assertEqual(total, 2*degree+1)

    def test_zero_potential_calibrates_all_eight_block_floors(self):
        for character in p.PARITIES:
            first = sum(character) or 2
            c = p.block_ground({}, character, L=max(3, first), bits=16)
            expected = str(first*(first+1))
            self.assertEqual((c['lower'], c['upper']), (expected, expected))
            self.assertTrue(p.verify_block(c, expected_parity=character))
        self.assertEqual(p.parity_basis((0, 0, 0), 1), [])
        self.assertTrue(all(b['degree'] >= 2 for b in p.parity_basis((0, 0, 0), 4)))

    def test_same_parity_next_omitted_degree_and_floor(self):
        for character, L, expected in [((1,0,0),3,5), ((0,0,0),2,4),
                                        ((1,1,1),3,5), ((1,1,0),4,6)]:
            self.assertEqual(p.next_omitted_degree(character, L), expected)
            k = p.Kernel(Q, character, L=L)
            self.assertEqual(k.beta, expected*(expected+1)+1)
            with self.assertRaises(ValueError):
                k.matrix(k.beta, 'lower')

    def test_quadratic_tail_complete_norm_for_x(self):
        columns = p.tail_couplings(Q, (1,0,0), L=1)['columns']
        # E[q²x²]=32/35 and <qx,x>/<x,x>=8/5.
        # Norm of all omitted cubic harmonics =32/35-64/75=32/525.
        self.assertEqual(len(columns), 2)
        self.assertTrue(all(c['degree'] == 3 for c in columns))
        norm = sum((c['entries'][0]**2/c['mass'] for c in columns), F(0))
        self.assertEqual(norm, F(32, 525))
        k = p.Kernel(Q, (1,0,0), L=1)
        correction = k.matrix(0, 'upper_ritz')[0][0]-k.matrix(0, 'lower')[0][0]
        self.assertEqual(correction, F(32, 6825))

    def test_zero_character_projection_removes_constant_norm(self):
        form = p.block_form(Q, (0,0,0), L=2)
        u = form['basis'][0]['polynomial']
        qu = a.multiply(a.polynomial(Q), u)
        self.assertEqual(a.integrate(qu), F(1, 5))
        finite_norm = sum(a.inner(qu, b['polynomial'])**2/b['mass'] for b in form['basis'])
        omitted_norm = a.inner(qu, qu)-a.integrate(qu)**2-finite_norm
        tail = p.tail_couplings(Q, (0,0,0), L=2)['columns']
        self.assertEqual(sum((c['entries'][0]**2/c['mass'] for c in tail), F(0)), omitted_norm)
        self.assertTrue(all(c['degree'] == 4 for c in tail))

    def test_qxy_cannot_be_split_into_eight_invariant_blocks(self):
        self.assertFalse(p.reflection_symmetry(XY)['all_eight_reflections'])
        cross = p.cross_block_form(XY, (1,0,0), (0,1,0), L=1)
        self.assertEqual(cross, [[F(1, 15)]])
        for action in (lambda:p.block_ground(XY, (1,0,0)), lambda:p.full_ground(XY),
                       lambda:p.adaptive_ground(XY)):
            with self.assertRaises(ValueError): action()
        self.assertTrue(all(value == 0 for row in p.cross_block_form(Q, (1,0,0), (0,1,0)) for value in row))

    def test_symmetry_is_checked_modulo_sphere_relation(self):
        sphere_zero = a.multiply(a.polynomial({'1,0,0':1}), a.add(a.R2, {a.ZERO:-1}))
        equivalent = a.add(a.polynomial(Q), sphere_zero)
        self.assertTrue(p.reflection_symmetry(equivalent)['all_eight_reflections'])
        self.assertEqual((p.potential_range(equivalent)['lower'],p.potential_range(equivalent)['upper']), ('1','3'))
        c = p.block_ground(equivalent, (1,0,0), L=3, bits=32)
        self.assertEqual((c['lower'],c['upper']), (self.block['lower'],self.block['upper']))
        self.assertTrue(p.verify_block(c, expected_q=equivalent))
        self.assertFalse(p.verify_block(c, expected_q=Q))

    def test_cache_is_bound_to_q_character_and_degree_and_returns_copies(self):
        cache = p.FormCache()
        first, _ = cache.form(Q, (1,0,0), 3)
        first['A'][0][0] = F(999)
        replay, _ = cache.form(Q, (1,0,0), 3)
        self.assertNotEqual(replay['A'][0][0], 999)
        other, _ = cache.form(Q, (0,1,0), 3)
        self.assertNotEqual(replay['parity'], other['parity'])
        bigger, _ = cache.form(Q, (1,0,0), 5)
        self.assertGreater(len(bigger['mass']),len(replay['mass']))
        shifted = a.add(a.polynomial(Q),{a.ZERO:F(1)})
        changed, _ = cache.form(shifted,(1,0,0),3)
        self.assertEqual(changed['A'][0][0]-replay['A'][0][0],replay['mass'][0])

    def test_all_eight_other_block_lowers_remain_in_final_certificate(self):
        c = self.original['certificate']
        self.assertEqual(len(c['characters']),8)
        self.assertEqual([tuple(e['parity']) for e in c['characters']], list(p.PARITIES))
        self.assertEqual(F(c['lower']), min(F(e['lower']) for e in c['characters']))
        self.assertEqual(F(c['upper']), min(F(e['upper']) for e in c['characters'] if e['format']==p.BLOCK_FORMAT))
        bad = copy.deepcopy(c); bad['characters'].pop()
        self.assertFalse(p.verify_full(bad))
        bad = copy.deepcopy(c); bad['characters'][0]=copy.deepcopy(bad['characters'][1])
        self.assertFalse(p.verify_full(bad))
        bad = copy.deepcopy(c)
        index=next(i for i,e in enumerate(bad['characters']) if e['format']!=p.BLOCK_FORMAT)
        bad['characters'][index]['lower']='999'
        self.assertFalse(p.verify_full(bad))

    def test_original_target_reached_without_claiming_runtime_proof(self):
        c=self.original['certificate']
        self.assertEqual(self.original['status'],'target_met')
        self.assertLessEqual(F(c['exact_width']),F(1,10**8))
        self.assertTrue(p.verify_result(self.original,expected_q=Q,expected_width=F(1,10**8)))
        self.assertFalse(p.verify_result(self.original,expected_width=F(1,10**10)))
        bad=copy.deepcopy(self.original);bad['status']='certified_open_gap'
        self.assertFalse(p.verify_result(bad))

    def test_000_constant_tail_mass_q_and_scope_tampering_rejected(self):
        for mode in ('tail','mass','scope','q','projection','degree_floor'):
            bad=copy.deepcopy(self.block)
            if mode=='tail':bad['complete_tail']['columns'].pop()
            elif mode=='mass':bad['mass'][0]='1'
            elif mode=='scope':bad['scope']='all_mean_zero_H1'
            elif mode=='q':bad['potential']={'2,0,0':'1'}
            elif mode=='projection':bad['mean_zero']=False
            else:bad['complete_tail']['first_omitted_degree']=4
            self.assertFalse(p.verify_block(bad),mode)

    def test_max_degree_budget_is_respected_in_initial_character_coverage(self):
        c=p.adaptive_ground({'4,0,0':-100},start_L=1,max_L=1,bits=16,wall_budget_seconds=10)
        retained=[e['retained_degree'] for e in c['certificate']['characters'] if e['format']==p.BLOCK_FORMAT]
        self.assertTrue(all(L<=1 for L in retained),'Initial full coverage must not silently exceed max_L')
        self.assertEqual(c['status'],'certified_open_gap')
        self.assertTrue(p.verify_result(c))

    def test_saved_original_proofs_replay_without_search(self):
        found,original_found=0,False
        with patch.object(p,'block_ground',side_effect=AssertionError('Verifier must not search')):
            for path in sorted((ROOT/'C18').rglob('*.json')):
                if path.name=='INDEPENDENT_TESTS.json':continue
                for c in _certificate_nodes(json.loads(path.read_text())):
                    verify=(p.verify_block if c['format']==p.BLOCK_FORMAT else
                            p.verify_full if c['format']==p.FULL_FORMAT else
                            a.verify if c['format']==a.FORMAT else p.verify_result)
                    self.assertTrue(verify(c),str(path)+':'+c['format'])
                    found+=1
                    if c['format']==p.DRIVER_FORMAT and c['potential']==Q and c['requested_width']=='1/100000000':
                        self.assertEqual(c['status'],'target_met')
                        original_found=True
        self.assertGreaterEqual(found,4)
        self.assertTrue(original_found,'Replay the unchanged C18 potential and original 1e-8 goal')


if __name__=='__main__':
    unittest.main()
