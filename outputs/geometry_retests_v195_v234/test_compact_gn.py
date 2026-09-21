import copy
from fractions import Fraction as F
from math import comb
import unittest
import compact_gn as g
from support import old_gn


def cap(n,original=False):
    return g.centered_function([[n,1]],2**n if original else 1)['terms']


def dense_mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return c


def dense_mean(p): return sum((x/F(i+1) for i,x in enumerate(p) if i%2==0),F(0))


class CompactGNTests(unittest.TestCase):
    def test_sparse_canonical_merge_and_mean(self):
        c=g.mean_certificate([[4,2],[0,3],[4,-1],[2,0]])
        self.assertEqual(c['terms'],[[0,'3'],[4,'1']])
        self.assertEqual(F(c['mean']),F(16,5))
        self.assertTrue(g.verify(c))

    def test_invalid_sparse_exponents_and_inexact_coefficients(self):
        for p in [[[True,1]],[[-1,1]],[[1,1.0]],[[1,True]],[[g.MAX_POWER+1,1]]]:
            with self.assertRaises(ValueError):g.mean_certificate(p)

    def test_centering_and_scale_representative(self):
        c=g.centered_function([[4,1],[0,2]],-3)
        self.assertEqual(c['terms'],[[0,'3/5'],[4,'-3']])
        self.assertEqual(c['mean'],'0')
        self.assertTrue(g.verify(c))
        with self.assertRaises(ValueError):g.centered_function([[4,1]],0)

    def test_mass_pair_kernel_direct_covariance(self):
        for a,b in [(0,64),(1,1),(4,16),(16,64),(64,64)]:
            self.assertEqual(g.mass_kernel(a,b),F(1,a+b+1)-F(1,(a+1)*(b+1)))
        c=g.mass_form([[4,2],[16,-3]])
        direct=4*g.mass_kernel(4,4)-12*g.mass_kernel(4,16)+9*g.mass_kernel(16,16)
        self.assertEqual(F(c['value']),direct)
        self.assertTrue(g.verify(c))

    def test_energy_pair_kernel_sphere_factor(self):
        self.assertEqual(g.energy_kernel(0,0),0)
        self.assertEqual(g.energy_kernel(64,0),0)
        self.assertEqual(g.energy_kernel(1,1),F(1,6))
        self.assertEqual(g.energy_kernel(4,16),F(16,105))
        c=g.energy_form([[4,2],[16,-3]])
        function=c['function']['terms']
        self.assertEqual(F(c['value']),F(g.moments(function)['energy']))
        self.assertTrue(g.verify(c))

    def test_quartic_kernel_support_and_operation_count(self):
        a=g.quartic_moment(cap(4));b=g.quartic_moment(cap(64))
        self.assertEqual(a['operation_counts'],b['operation_counts'])
        self.assertEqual(b['operation_counts']['coefficient_products'],13)
        self.assertEqual(b['support_sizes'],{'input':2,'square':3,'fourth':5})
        self.assertEqual(b['fourth_power'][-1][0],256)
        self.assertTrue(g.verify(b))

    def test_quartic_independent_binomial_formula(self):
        for n in (4,16,64):
            expected=sum(F(comb(4,j))*(-F(1,n+1))**(4-j)/F(n*j+1) for j in range(5))
            self.assertEqual(F(g.quartic_moment(cap(n))['quartic']),expected)

    def test_all_original_trials_exact_scale(self):
        for n in (4,16,64):
            c=g.evaluate_cap(n,True);m=c['ratio']['moments'];law=c['family_value']
            for field,factor in [('mass',2**(2*n)),('energy',2**(2*n)),('quartic',2**(4*n))]:
                self.assertEqual(F(m[field]),F(law['normalized_'+field])*factor)
            self.assertTrue(g.verify(c,expected_n=n))
            self.assertFalse(g.verify(c,expected_n=n+1))

    def test_n64_exact_ratio_reference(self):
        r=g.ratio(cap(64,True))
        self.assertEqual(F(r['probability_ratio']),F(402777216,209564225))
        self.assertEqual(F(r['pi_times_area_K']),F(100694304,209564225))
        self.assertEqual(F(r['probability_ratio']),4*F(r['pi_times_area_K']))

    def test_signed_scale_ratio_invariance(self):
        p=cap(16);q=[[n,-F(7,3)*F(c)] for n,c in p]
        a,b=g.ratio(p),g.ratio(q)
        self.assertEqual(a['probability_ratio'],b['probability_ratio'])
        self.assertFalse(g.verify(a,expected_terms=q))

    def test_old_dense_supported_low_n_agreement(self):
        for n in (4,16):
            old=old_gn.moments(g.original_cap(n));new=g.evaluate_cap(n)['ratio']['moments']
            for field in ('mean','mass','energy','quartic'):self.assertEqual(F(old[field]),F(new[field]))

    def test_old_degree64_rejection_is_real(self):
        original=g.original_cap(64)
        self.assertEqual(len(original),65)
        with self.assertRaises(ValueError):old_gn.moments(original)
        self.assertTrue(g.verify(g.bridge_original(original,64),expected_original_coefficients=original))

    def test_n64_dense_original_independent_integration(self):
        n=64;p=[F(comb(n,j)) for j in range(n+1)];p[0]-=F(2**n,n+1)
        square=dense_mul(p,p);fourth=dense_mul(square,square)
        dp=[j*p[j] for j in range(1,len(p))]
        energy=dense_mean(dense_mul([1,0,-1],dense_mul(dp,dp)))
        c=g.evaluate_cap(n)['ratio']['moments']
        self.assertEqual(F(c['mass']),dense_mean(square))
        self.assertEqual(F(c['quartic']),dense_mean(fourth))
        self.assertEqual(F(c['energy']),energy)

    def test_constant_and_zero_energy_reject(self):
        for p in ([],[[0,1]],[[4,1]]):
            with self.assertRaises(ValueError):g.ratio(p)

    def test_full_residual_keeps_degree192(self):
        r=g.residual(cap(64))
        self.assertEqual(r['full_degree'],192)
        self.assertEqual(r['residual'][-1][0],192)
        self.assertNotEqual(F(r['residual'][-1][1]),0)
        self.assertEqual(r['residual_mean'],'0')
        self.assertEqual(r['inner_trial'],'0')
        self.assertFalse(r['degree_truncation_used'])
        self.assertTrue(g.verify(r))

    def test_laplacian_matches_old_dense_low_n(self):
        n=4;r=g.residual(cap(n));old=old_gn.residual(old_gn.centered_cap(n))
        # Expand only in the independent test; implementation residual stays sparse.
        out=[F(0)]*(3*n+1)
        for power,c in r['residual']:
            for j in range(power+1):out[j]+=F(c)*F(comb(power,j),2**power)
        while len(out)>1 and out[-1]==0:out.pop()
        self.assertEqual(out,list(map(F,old['residual'])))

    def test_sparse_direction_sign_and_tamper(self):
        d=g.direction(cap(64))
        self.assertGreater(F(d['derivative_pi_times_K']),0)
        self.assertEqual(d['derivative_pi_times_K'],d['identity_derivative'])
        v=[[n,-F(c)] for n,c in d['direction']]
        self.assertLess(F(g.direction(cap(64),v)['derivative_pi_times_K']),0)
        self.assertTrue(g.verify(d))
        d['derivative_pi_times_K']='0';self.assertFalse(g.verify(d))

    def test_plane_n64_complete_charts(self):
        basis=[cap(16),cap(64)];c=g.plane(basis,F(1,10**7))
        self.assertTrue(g.verify_plane(c,basis))
        for x,y in [(1,0),(0,1),(1,3),(-4,1),(7,-2)]:
            u=g._add(g._scale(g.sparse(basis[0]),F(x)),g._scale(g.sparse(basis[1]),F(y)))
            self.assertLessEqual(F(g.ratio(u)['pi_times_area_K']),F(c['upper_pi_times_K']))
        self.assertFalse(c['universal_GN_upper_bound'])
        self.assertTrue(c['target_met'])
        c['charts'].pop();self.assertFalse(g.verify_plane(c))

    def test_plane_original_scaling_and_binding(self):
        basis=[cap(4,True),cap(64,True)]
        c=g.plane(basis,F(1,10**6))
        self.assertTrue(g.verify_plane(c,basis))
        self.assertFalse(g.verify_plane(c,[cap(4),cap(64)]))
        self.assertGreaterEqual(F(c['lower_pi_times_K']),F(g.ratio(basis[1])['pi_times_area_K']))

    def test_dependent_plane_rejected(self):
        with self.assertRaises(ValueError):g.plane([cap(64),cap(64)])

    def test_joint_plane_support_budget_rejects_before_return(self):
        f=g.centered_function([[i,1 if i==1 else F(1,1024)] for i in range(1,33)])['terms']
        h=g.centered_function([[i,1 if i==33 else F(1,1024)] for i in range(33,65)])['terms']
        self.assertEqual(len(f),33);self.assertEqual(len(h),33)
        with self.assertRaisesRegex(ValueError,'union of plane supports'):
            g.plane([f,h],F(1,10**5))

    def test_centered_function_respects_complete_support_budget(self):
        with self.assertRaisesRegex(ValueError,'Centering exceeds'):
            g.centered_function([[i,1] for i in range(1,65)])

    def test_family_symbolic_identity_and_integer_monotonicity(self):
        law=g.family_law();self.assertTrue(g.verify(law))
        self.assertEqual(law['n_domain'],'all_integers_n>=1')
        self.assertTrue(all(F(x)>0 for x in law['difference_after_n_equals_x_plus_1']))
        for n in (1,2,4,16,64,1024):
            a,b=g.family_value(n),g.family_value(n+1)
            self.assertLess(F(a['pi_times_area_K']),F(b['pi_times_area_K']))
            self.assertLess(F(a['pi_times_area_K']),F(1,2))
        law['limit_pi_times_area_K']='1';self.assertFalse(g.verify(law))

    def test_no_false_continuous_monotonicity(self):
        def j(n):return F(3)*n*(2*n+1)*(2*n*n-n+3)/(2*(4*n+1)*(3*n+1)*(n+1)**2)
        self.assertLess(j(F(11,10)),j(F(1)))
        with self.assertRaises(ValueError):g.family_value(F(11,10))

    def test_n1024_and_million_closed_values(self):
        for n in (1024,4096,10**6):
            c=g.family_value(n)
            self.assertTrue(g.verify(c,expected_n=n))
            self.assertEqual(F(c['pi_times_area_K']),F(g.ratio(cap(n))['pi_times_area_K']))
        self.assertTrue(g.verify(g.evaluate_cap(1024,False)))

    def test_original_bridge_all_coefficients_bound(self):
        p=g.original_cap(64);c=g.bridge_original(p,64)
        self.assertEqual(c['verified_coefficient_count'],65)
        self.assertTrue(g.verify(c,expected_n=64,expected_original_coefficients=p))
        for i in (0,31,64):
            altered=list(p);altered[i]+=1
            with self.assertRaises(ValueError):g.bridge_original(altered,64)
            self.assertFalse(g.verify(c,expected_original_coefficients=altered))

    def test_original_bridge_does_not_accept_normalized_cap(self):
        p=[F(x,2**16) for x in g.original_cap(16)]
        with self.assertRaises(ValueError):g.bridge_original(p,16)


if __name__=='__main__':unittest.main()
