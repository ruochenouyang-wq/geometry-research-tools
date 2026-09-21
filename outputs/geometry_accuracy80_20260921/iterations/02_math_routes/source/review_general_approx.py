"""Independent public-input checks for the rational-power polynomial model.

Uses global t-polynomial expansion and elementary antiderivatives, rather than
the implementation's scaled-error assembly or Legendre projection machinery.
No evaluation seed or holdout is read.
"""
from copy import deepcopy
from fractions import Fraction as F
from math import comb,isqrt
import json
import sys
import unittest
from unittest.mock import patch

import general_approx as g


PUBLIC_CASES = (
    ({'kind':'axis_profile','profile':'abs_power','axis':2,
      'exponent':'-1/3','amplitude':'-2/3','offset':'3/7'},'2/3',4,5),
    ({'kind':'axis_profile','profile':'abs_power','axis':0,
      'exponent':'-2/5','amplitude':'5/4','offset':'-7/6'},'3/4',3,4),
    ({'kind':'axis_profile','profile':'abs_power','axis':1,
      'exponent':'-1/4','amplitude':'0','offset':'-3'},'1/2',2,3),
)


def exact_nth_root(value,degree):
    if value in (0,1): return value
    low,high = 0,1 << ((value.bit_length()+degree-1)//degree)
    while low+1<high:
        mid = (low+high)//2
        if mid**degree<=value: low=mid
        else: high=mid
    root = high if high**degree==value else low
    if root**degree!=value: raise ValueError('not an exact power')
    return root


def exact_power(value,exponent):
    value,exponent = F(value),F(exponent)
    if value == 0:
        if exponent<=0: raise ValueError('divergent endpoint')
        return F(0)
    root = F(exact_nth_root(value.numerator,exponent.denominator),
             exact_nth_root(value.denominator,exponent.denominator))
    return root**exponent.numerator


def monomial_integral(exponent,left,right):
    antiderivative = F(exponent)+1
    if antiderivative<=0 and left==0: raise ValueError('divergent integral')
    return (exact_power(right,antiderivative)-exact_power(left,antiderivative))/antiderivative


def polynomial_value(coefficients,t):
    return sum((c*t**k for k,c in enumerate(coefficients)),F(0))


def polynomial_norm(coefficients,left,right):
    return sum((a*b*monomial_integral(i+j,left,right)
                for i,a in enumerate(coefficients) for j,b in enumerate(coefficients)),F(0))


def independent_cells(cert):
    """Expand every compact template to an ordinary polynomial in physical t."""
    q = cert['function']
    alpha,A,b = (F(q[k]) for k in ('exponent','amplitude','offset'))
    d,m = alpha.denominator,-alpha.numerator
    r,L = F(cert['root_ratio']),cert['levels']
    cells = [(F(0),r**(d*L),[b+A*r**(-m*L)*F(cert['core_normal_coefficient'])])]
    for block in cert['blocks']:
        template = list(map(F,cert['templates'][str(block['degree'])]['polynomial_in_s']))
        for j in range(block['start'],block['start']+block['count']):
            left,right = r**(d*(j+1)),r**(d*j)
            # s=(t-left)/(right-left); coefficients are integer powers of t.
            polynomial = [F(0)]*len(template)
            for k,c in enumerate(template):
                for ell in range(k+1):
                    polynomial[ell] += A*r**(-m*j)*c*comb(k,ell)*(-left)**(k-ell)/(right-left)**k
            polynomial[0] += b
            cells.append((left,right,polynomial))
    return sorted(cells,key=lambda x:x[0])


def independent_error_squared(function,cells):
    alpha,A,b = (F(function[k]) for k in ('exponent','amplitude','offset'))
    total = F(0)
    for left,right,polynomial in cells:
        source = A*A*monomial_integral(2*alpha,left,right)
        source += 2*A*b*monomial_integral(alpha,left,right)+b*b*(right-left)
        cross = sum((c*(A*monomial_integral(alpha+k,left,right)
                        +b*monomial_integral(k,left,right))
                     for k,c in enumerate(polynomial)),F(0))
        total += source-2*cross+polynomial_norm(polynomial,left,right)
    return total


class GeneralApproxIndependentReview(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.certificates = [g.solve(q,tolerance='1/1000',root_ratio=r,
                                   max_degree=degree,max_levels=levels)
                            for q,r,degree,levels in PUBLIC_CASES]

    def test_whole_domain_error_by_expanding_every_physical_polynomial(self):
        for c in self.certificates:
            with self.subTest(function=c['function']):
                cells = independent_cells(c)
                self.assertEqual(cells[0][0],0)
                self.assertEqual(cells[-1][1],1)
                self.assertTrue(all(a[1]==b[0] for a,b in zip(cells,cells[1:])))
                error = independent_error_squared(c['function'],cells)
                self.assertEqual(error,F(c['error_squared']))
                self.assertEqual(error,F(c['core_error_squared'])+F(c['annuli_error_squared']))

    def test_evaluation_is_the_claimed_finite_piecewise_polynomial(self):
        for c in self.certificates:
            cells = independent_cells(c)
            self.assertEqual(g.evaluate(c,0),cells[0][2][0])
            self.assertEqual(g.evaluate(c,1),polynomial_value(cells[-1][2],F(1)))
            for left,right,p in cells:
                t = (left+right)/2
                expected = polynomial_value(p,t)
                self.assertEqual(g.evaluate(c,t),expected)
                self.assertEqual(g.evaluate(c,-t),expected)
            core = cells[0]
            self.assertEqual(g.evaluate(c,core[1]/3),g.evaluate(c,core[1]/2))
            if F(c['function']['amplitude']):
                # Nonzero-alpha source is nonconstant in this core; the output
                # is a finite constant and cannot be a copy of the source.
                alpha = F(c['function']['exponent'])
                d = alpha.denominator
                t1 = core[1]/(F(2)**d)
                t2 = core[1]/(F(3)**d)
                self.assertNotEqual(exact_power(t1,alpha),exact_power(t2,alpha))
                self.assertEqual(g.evaluate(c,t1),g.evaluate(c,t2))

    def test_core_constant_is_the_exact_least_squares_constant(self):
        for c in self.certificates:
            q = c['function'];alpha,A,b = (F(q[k]) for k in ('exponent','amplitude','offset'))
            left,right,p = independent_cells(c)[0]
            mean = (A*monomial_integral(alpha,0,right)+b*right)/right
            self.assertEqual(p[0],mean)
            baseline = independent_error_squared(q,[(left,right,p)])
            perturbed = independent_error_squared(q,[(left,right,[p[0]+F(2,7)])])
            self.assertEqual(perturbed-baseline,right*F(4,49))

    def test_nonoptimal_core_with_honest_reintegrated_error_is_valid(self):
        c = deepcopy(self.certificates[0])
        c['core_normal_coefficient'] = str(F(c['core_normal_coefficient'])+F(1,10))
        cells = independent_cells(c)
        error = independent_error_squared(c['function'],cells)
        core_error = independent_error_squared(c['function'],cells[:1])
        den = 2**c['sqrt_bits']
        n = isqrt(error.numerator*den*den//error.denominator)
        if F(n*n,den*den)<error: n+=1
        upper = F(n,den)
        c.update(error_squared=str(error),core_error_squared=str(core_error),
                 error_upper=str(upper),sqrt_upper_squared=str(upper*upper),
                 status='target_met' if upper<=F(c['tolerance']) else 'certified_open')
        self.assertTrue(g.verify(c,c['function'],c['tolerance']))

    def test_seam_values_match_declared_core_exception_and_outer_annuli(self):
        for c in self.certificates:
            self.assertEqual(c['model']['seam_rule'],
                             'outer cell at internal annulus seams; core includes t_levels and zero')
            cells = independent_cells(c)
            self.assertEqual(g.evaluate(c,cells[0][1]),cells[0][2][0])
            for inner,outer in zip(cells[1:],cells[2:]):
                seam = inner[1]
                self.assertEqual(seam,outer[0])
                self.assertEqual(g.evaluate(c,seam),polynomial_value(outer[2],seam))

    def test_moments_by_physical_endpoint_antiderivatives(self):
        for q,r,_,_ in PUBLIC_CASES:
            alpha,root = F(q['exponent']),F(r)
            for k in (0,1,4):
                for power in (0,1,2):
                    expected = monomial_integral(k+power*alpha,root**alpha.denominator,F(1))
                    self.assertEqual(g.power_moment(alpha,k,root,F(1),power),expected)

    def test_verifier_does_not_run_planner_projection_or_cache(self):
        c = self.certificates[0]
        with (patch.object(g,'solve',side_effect=AssertionError('planner')),
             patch.object(g,'shape_series',side_effect=AssertionError('projection')),
             patch.object(g.ShapeCache,'get',side_effect=AssertionError('cache')),
             patch.object(g.ShapeCache,'put',side_effect=AssertionError('cache'))):
            self.assertTrue(g.verify(c,c['function'],'1/1000'))

    def test_poisoned_planning_cache_cannot_create_a_false_certificate(self):
        q,r,degree,levels = PUBLIC_CASES[0]
        cache = g.ShapeCache()
        cache.put(F(q['exponent']),F(r),0,[F(0)],F(0))
        with self.assertRaises(ArithmeticError):
            g.solve(q,tolerance='1/1000',root_ratio=r,max_degree=degree,
                    max_levels=levels,shape_cache=cache)

    def test_tampering_cannot_hide_polynomial_or_source_error(self):
        original = self.certificates[0]
        mutations = [
            lambda c:c.update(error_squared='0'),
            lambda c:c.update(error_upper='0'),
            lambda c:c.update(core_error_squared='0'),
            lambda c:c.update(core_normal_coefficient='0'),
            lambda c:c.update(root_ratio='3/4'),
            lambda c:c.update(root_degree=c['root_degree']+1),
            lambda c:c['function'].update(exponent='-2/5'),
            lambda c:c['function'].update(amplitude='-1/3'),
            lambda c:c['function'].update(offset='0'),
            lambda c:c['templates'][next(iter(c['templates']))]['polynomial_in_s'].__setitem__(0,'0'),
            lambda c:c['templates'][next(iter(c['templates']))].update(error_squared='0'),
            lambda c:c['model'].update(polynomial_rule='source_function'),
            lambda c:c.update(spectral_transfer_claimed=True),
        ]
        for i,mutation in enumerate(mutations):
            c = deepcopy(original);mutation(c)
            with self.subTest(mutation=i): self.assertFalse(g.verify(c,original['function'],'1/1000'))

    def test_gaps_overlaps_and_missing_annuli_are_rejected(self):
        original = self.certificates[0]
        for mutation in ('gap','overlap','missing','extra','unused'):
            c = deepcopy(original)
            if mutation=='gap': c['blocks'][0]['start']=1
            elif mutation=='overlap': c['blocks'].append(deepcopy(c['blocks'][0]))
            elif mutation=='missing': c['blocks'][-1]['count']-=1
            elif mutation=='extra': c['blocks'][-1]['count']+=1
            else: c['templates']['32']={'polynomial_in_s':['0']*33,'error_squared':'1'}
            with self.subTest(mutation=mutation): self.assertFalse(g.verify(c))

    def test_external_request_function_and_tolerance_are_bound(self):
        c = self.certificates[0]
        self.assertFalse(g.verify(c,dict(c['function'],axis=1)))
        self.assertFalse(g.verify(c,tolerance='1/100000000'))
        self.assertFalse(g.verify(None))
        self.assertFalse(g.verify([]))

    def test_metadata_type_changes_and_hidden_extra_fields_are_rejected(self):
        original = self.certificates[0]
        for field,value in [('spectral_transfer_claimed',0),('global_polynomial',0),
                            ('formal_proof_assistant_checked',0),('hidden_function','original')]:
            bad = deepcopy(original);bad[field] = value
            with self.subTest(field=field): self.assertFalse(g.verify(bad))
        bad = deepcopy(original);bad['blocks'][0]['count'] = True
        self.assertFalse(g.verify(bad))

    def test_large_rational_text_roundtrip_preserves_global_integer_limit(self):
        limit_query = getattr(sys,'get_int_max_str_digits',None)
        limit = limit_query() if limit_query else None
        value = -F(10**5000+17,10**5001+19)
        encoded = g.rational_text(value)
        self.assertGreater(len(encoded),10000)
        self.assertEqual(g.rational(json.loads(json.dumps(encoded))),value)
        self.assertEqual(limit_query() if limit_query else None,limit)

    def test_dyadic_error_upper_is_outward_and_status_is_truthful(self):
        for c in self.certificates:
            error,upper,tolerance = F(c['error_squared']),F(c['error_upper']),F(c['tolerance'])
            step = F(1,2**c['sqrt_bits'])
            self.assertGreaterEqual(upper*upper,error)
            if upper: self.assertLess((upper-step)**2,error)
            self.assertEqual(c['status'],'target_met' if upper<=tolerance else 'certified_open')
        self.assertEqual(self.certificates[2]['error_squared'],'0')
        self.assertEqual(g.evaluate(self.certificates[2],F(17,29)),-3)


if __name__ == '__main__': unittest.main(verbosity=2)
