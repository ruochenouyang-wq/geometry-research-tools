"""Public development examples; no holdout or private-seed access."""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from fractions import Fraction as F
from math import comb
from copy import deepcopy
from unittest.mock import patch
import json
import unittest

import general_approx as g


def source(alpha='-2/9',amplitude='-1/3',offset='0',axis=2):
    return {'kind':'axis_profile','profile':'abs_power','exponent':alpha,
            'amplitude':amplitude,'offset':offset,'axis':axis}


def independent_integral(alpha,k,power,left_root,right_root):
    """Substitute t=u^d directly; do not call the implementation's moments."""
    d = alpha.denominator
    primitive_exponent = F(k+1)+power*alpha
    integer_exponent = primitive_exponent*d
    assert integer_exponent.denominator == 1 and integer_exponent > 0
    return (right_root**int(integer_exponent)-left_root**int(integer_exponent))/primitive_exponent


def direct_cell_error(q,left_root,right_root,polynomial_in_s,scale):
    """Expand the actual returned local polynomial in global t and integrate."""
    alpha,A,B = (F(q[k]) for k in ('exponent','amplitude','offset'))
    d = alpha.denominator
    left,right = left_root**d,right_root**d
    width = right-left
    coefficients = [sum((F(c)*comb(i,k)*(-left)**(i-k)/width**i
                         for i,c in enumerate(polynomial_in_s) if i>=k),F(0))*scale
                    for k in range(len(polynomial_in_s))]
    coefficients[0] += B
    source_norm = (A*A*independent_integral(alpha,0,2,left_root,right_root)
                   +2*A*B*independent_integral(alpha,0,1,left_root,right_root)+B*B*width)
    cross = sum((c*(A*independent_integral(alpha,k,1,left_root,right_root)
                    +B*(right**(k+1)-left**(k+1))/F(k+1))
                 for k,c in enumerate(coefficients)),F(0))
    norm = sum((c*e*(right**(k+l+1)-left**(k+l+1))/F(k+l+1)
                for k,c in enumerate(coefficients) for l,e in enumerate(coefficients)),F(0))
    return source_norm-2*cross+norm


class GeneralApproxTests(unittest.TestCase):
    def test_public_H11_H12_H13_reach_their_original_targets(self):
        for name,alpha,A,tol in [('H11','-2/9','-1/3','1/1000000'),
                                 ('H12','-2/11','-3/5','1/100000000'),
                                 ('H13','-2/9','-3/5','1/10000000000')]:
            with self.subTest(name=name):
                q = source(alpha,A)
                c = g.solve(q,tol)
                self.assertTrue(g.verify(c,q,tol))
                self.assertLessEqual(g.rational(c['error_upper']),F(tol))
                self.assertEqual(c['status'],'target_met')
                self.assertGreater(g.rational(c['error_squared']),0)
                self.assertEqual(c['model']['kind'],'finite_even_piecewise_polynomial_scaled_templates')

    def test_new_public_development_parameters(self):
        for q,tol in [(source('-1/3','-2/5','7/3',1),'1/100000000'),
                      (source('-3/8','3/7','-1',0),'1/10000000000'),
                      (source('-5/13','-1/2'),'1/100000000'),
                      (source('-1/10','1'),'1/1000000')]:
            c = g.solve(q,tol)
            self.assertTrue(g.verify(c,q,tol))
            self.assertLessEqual(g.rational(c['error_upper']),F(tol))

    def test_general_moments_by_direct_root_substitution(self):
        for alpha in (F(-1,4),F(-2,9),F(-2,11),F(-5,13)):
            for k in (0,1,3):
                for power in (0,1,2):
                    self.assertEqual(g.power_moment(alpha,k,F(1,3),F(2,3),power),
                                     independent_integral(alpha,k,power,F(1,3),F(2,3)))

    def test_core_exact_calibrations(self):
        for alpha,error in [(F(-1,4),F(2,9)),(F(-2,9),F(36,245)),(F(-2,11),F(44,567))]:
            self.assertEqual(g._core_error(alpha,g._core_coefficient(alpha)),error)
            self.assertLess(error,1/(1+2*alpha))

    def test_local_polynomial_error_by_independent_global_expansion(self):
        q,ratio = source('-2/9','-2/3','5/7'),F(17,18)
        shape = list(g.shape_series(F(q['exponent']),ratio,3))[-1]
        alpha = F(q['exponent']);m,d = -alpha.numerator,alpha.denominator
        for j in (0,2):
            right,left = ratio**j,ratio**(j+1)
            actual = direct_cell_error(q,left,right,shape[0],F(q['amplitude'])*right**(-m))
            expected = F(q['amplitude'])**2*right**(d-2*m)*shape[1]
            self.assertEqual(actual,expected)

    def test_whole_domain_error_is_the_sum_of_actual_polynomial_cell_errors(self):
        q = source('-2/11','3/5','-2/3')
        c = g.solve(q,'1/100000000',max_degree=2,max_levels=3)
        r,L,d,m = F(c['root_ratio']),c['levels'],c['root_degree'],c['singular_numerator']
        total = direct_cell_error(q,F(0),r**L,[F(c['core_normal_coefficient'])],F(q['amplitude'])*r**(-m*L))
        for block in c['blocks']:
            coeffs = list(map(g.rational,c['templates'][str(block['degree'])]['polynomial_in_s']))
            for j in range(block['start'],block['start']+block['count']):
                total += direct_cell_error(q,r**(j+1),r**j,coeffs,F(q['amplitude'])*r**(-m*j))
        self.assertEqual(total,g.rational(c['error_squared']))

    def test_offset_and_amplitude_scaling_are_exact(self):
        q = source('-2/9','1','0')
        c = g.solve(q,'1/1000000',max_levels=8,max_degree=3)
        shifted = g.solve(dict(q,offset='-19/7'),'1/1000000',max_levels=8,max_degree=3)
        self.assertEqual(c['error_squared'],shifted['error_squared'])
        self.assertEqual(c['templates'],shifted['templates'])
        t = F(1,2)
        self.assertEqual(g.evaluate(shifted,t),g.evaluate(c,t)-F(19,7))
        constant = g.solve(source('-2/9','0','-19/7'),'1/10000000000')
        self.assertEqual(g.rational(constant['error_squared']),0)
        self.assertEqual(g.evaluate(constant,0),-F(19,7))

    def test_root_scale_exponents_change_with_the_source_exponent(self):
        for alpha,gamma,degree in [('-2/9',5,9),('-2/11',7,11),('-1/4',2,4)]:
            c = g.solve(source(alpha),'1/1000000')
            self.assertEqual(c['squared_error_scale_exponent'],gamma)
            self.assertEqual(c['root_degree'],degree)

    def test_true_finite_nonuniform_polynomial_representation(self):
        c = g.solve(source(), '1/10000000000')
        self.assertGreater(len(c['templates']),1)
        self.assertLess(c['stored_template_coefficient_slots'],c['expanded_polynomial_coefficient_slots'])
        cursor = 0
        for block in c['blocks']:
            self.assertEqual(block['start'],cursor)
            cursor += block['count']
            self.assertTrue(all(type(x) is str for x in c['templates'][str(block['degree'])]['polynomial_in_s']))
        self.assertEqual(cursor,c['levels'])
        self.assertEqual(g.evaluate(c,F(2,3)),g.evaluate(c,-F(2,3)))

    def test_core_boundary_convention_matches_evaluation(self):
        q = source()
        c = g.solve(q,'1/100',max_levels=4,max_degree=2)
        r,L,d,m = F(c['root_ratio']),c['levels'],c['root_degree'],c['singular_numerator']
        core_value = F(q['offset'])+F(q['amplitude'])*r**(-m*L)*F(c['core_normal_coefficient'])
        self.assertEqual(g.evaluate(c,0),core_value)
        self.assertEqual(g.evaluate(c,r**(d*L)),core_value)
        self.assertIn('core includes t_levels',c['model']['seam_rule'])

    def test_independent_verify_uses_neither_planner_nor_shape_cache(self):
        q,tol = source(), '1/1000000'
        c = json.loads(json.dumps(g.solve(q,tol)))
        with patch.object(g,'shape_series',side_effect=AssertionError('projection search')):
            with patch.object(g.ShapeCache,'get',side_effect=AssertionError('planning cache')):
                self.assertTrue(g.verify(c,q,tol))

    def test_cache_key_includes_exponent_and_cached_shapes_remain_checked(self):
        cache = g.ShapeCache()
        ratio = F(17,18)
        first = list(g.shape_series(F(-2,9),ratio,2,cache))
        second = list(g.shape_series(F(-2,11),ratio,2,cache))
        self.assertNotEqual(first[-1],second[-1])
        self.assertIsNotNone(cache.get(F(-2,9),ratio,2))
        self.assertIsNotNone(cache.get(F(-2,11),ratio,2))
        with patch.object(g,'_local_moments',side_effect=AssertionError('cached search needlessly recomputed')):
            self.assertEqual(list(g.shape_series(F(-2,9),ratio,2,cache)),first)
        q = source()
        cold = g.solve(q,'1/1000000',root_ratio=ratio)
        warm = g.solve(q,'1/1000000',root_ratio=ratio,shape_cache=cache)
        self.assertEqual(cold,warm)

    def test_cache_poisoning_cannot_create_a_false_certificate(self):
        cache = g.ShapeCache()
        cache.put(F(-2,9),F(17,18),0,[F(0)],F(0))
        with self.assertRaises(ArithmeticError):
            g.solve(source(),'1/1000000',root_ratio='17/18',shape_cache=cache)

    def test_explicit_budget_open_is_still_a_valid_bound(self):
        q,tol = source(), '1/10000000000'
        c = g.solve(q,tol,max_degree=0,max_levels=1)
        self.assertTrue(g.verify(c,q,tol))
        self.assertEqual(c['status'],'certified_open')
        self.assertGreater(g.rational(c['error_upper']),F(tol))

    def test_near_L2_boundary_needs_more_layers_and_large_fraction_support(self):
        q,tol = source('-49/100','-1/10'),'1/1000000'
        insufficient = g.solve(q,tol)
        enough = g.solve(q,tol,max_levels=8192)
        self.assertTrue(g.verify(insufficient,q,tol))
        self.assertTrue(g.verify(enough,q,tol))
        self.assertEqual(insufficient['status'],'certified_open')
        self.assertEqual(enough['status'],'target_met')
        self.assertGreater(len(enough['error_squared']),4300)

    def test_large_rational_io_does_not_change_global_python_limits(self):
        before = sys.get_int_max_str_digits()
        value = F(10**6000+7,10**5000+3)
        self.assertEqual(g.rational(g.rational_text(value)),value)
        self.assertEqual(sys.get_int_max_str_digits(),before)

    def test_original_function_and_requested_tolerance_are_bound(self):
        q,tol = source(), '1/1000000'
        c = g.solve(q,tol)
        for change in ({'axis':0},{'amplitude':'1/3'},{'offset':'1'},{'exponent':'-2/11'}):
            self.assertFalse(g.verify(c,dict(q,**change),tol))
        self.assertFalse(g.verify(c,q,'1/10000000'))

    def test_math_and_representation_tampering_is_rejected(self):
        c = g.solve(source(),'1/1000000')
        for key,value in [('error_squared','0'),('error_upper','0'),('root_degree',4),
                          ('singular_numerator',1),('squared_error_scale_exponent',2),
                          ('scope','one_sample_only'),('status','certified_open'),
                          ('formal_proof_assistant_checked',0)]:
            bad = deepcopy(c);bad[key]=value
            self.assertFalse(g.verify(bad),key)
        bad = deepcopy(c);bad['blocks'][0]['start']=1
        self.assertFalse(g.verify(bad))
        bad = deepcopy(c);bad['blocks'][-1]['count']-=1
        self.assertFalse(g.verify(bad))
        bad = deepcopy(c);name=next(iter(bad['templates']))
        bad['templates'][name]['polynomial_in_s'][0]='0'
        self.assertFalse(g.verify(bad))
        bad = deepcopy(c);bad['model']['kind']='original_singular_function'
        self.assertFalse(g.verify(bad))

    def test_invalid_domains_and_inexact_inputs_are_rejected(self):
        for alpha in ('-1/2','-3/4','0','1/4'):
            with self.assertRaises(ValueError):
                g.solve(source(alpha))
        for q in (dict(source(),axis=True),dict(source(),amplitude=0.1),dict(source(),extra='field')):
            with self.assertRaises(ValueError):
                g.solve(q)
        for ratio in ('0','1','-1/2'):
            with self.assertRaises(ValueError):
                g.solve(source(),root_ratio=ratio)


if __name__ == '__main__':
    unittest.main()
