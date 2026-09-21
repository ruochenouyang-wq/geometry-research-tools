from copy import deepcopy
from fractions import Fraction as F
import random
import unittest
import error_bounds as e
import run


class Algebra(unittest.TestCase):
    def test_bernstein_interior_extremum(self):
        p = [F(1,9),F(-2,3),F(1)]  # (t-1/3)^2
        lo,hi = e.polynomial_range(p,3)
        self.assertLessEqual(lo,0)
        self.assertGreaterEqual(hi,F(16,9))
        self.assertGreater(e.evaluate(p,F(-1)),0)
        self.assertGreater(e.evaluate(p,F(1)),0)

    def test_subdivision_encloses_and_tightens(self):
        p = list(map(F,[1,-2,3,4,-1,5]))
        bounds = [e.polynomial_range(p,d) for d in range(5)]
        for (lo,hi),(a,b) in zip(bounds,bounds[1:]):
            self.assertLessEqual(lo,a)
            self.assertGreaterEqual(hi,b)
        for j in range(101):
            value = e.evaluate(p,F(j,50)-1)
            self.assertLessEqual(bounds[-1][0],value)
            self.assertGreaterEqual(bounds[-1][1],value)

    def test_bernstein_reconstruction_at_rational_points(self):
        from math import comb
        p = list(map(F,[1,2,0,-4,2]))
        b = e.bernstein_coefficients(p)
        for x in [F(0),F(1,7),F(1,2),F(1)]:
            value = sum(b[i]*comb(4,i)*x**i*(1-x)**(4-i) for i in range(5))
            self.assertEqual(value,e.evaluate(p,2*x-1))

    def test_banded_inertia_matches_dense(self):
        rng = random.Random(508)
        for n in [2,5,8]:
            for width in [0,1,2,4]:
                a = [[F(0)]*n for _ in range(n)]
                for i in range(n):
                    for j in range(i,min(n,i+width+1)):
                        a[i][j] = a[j][i] = F(rng.randint(-5,5),rng.randint(1,4))
                self.assertEqual(e.banded_inertia(a),e.base.inertia(a))

    def test_zero_pivot_exact_fallback(self):
        for a in [[[0,1],[1,0]],[[1,1],[1,1]],[[0,0],[0,0]]]:
            self.assertEqual(e.banded_inertia(a),e.base.inertia(a))

    def test_basis_change_preserves_ritz_counts(self):
        d = e.base.assemble([F(0),F(1)],3)
        t = [[F(1),0,F(-1,3)],[0,F(1),0],[0,0,F(4,3)]]  # T_2=(4P_2-P_0)/3
        for x in [-2,0,2,7]:
            b = e.base.shifted(d,F(x))
            transformed = [[sum(t[k][i]*b[k][l]*t[l][j] for k in range(3) for l in range(3)) for j in range(3)] for i in range(3)]
            self.assertEqual(e.base.inertia(b),e.base.inertia(transformed))

    def test_invalid_parameters(self):
        for d in [-1,8,True]:
            with self.assertRaises(ValueError): e.polynomial_range([F(1)],d)
        for degree in [0,17,True]:
            with self.assertRaises(ValueError): e.exponential_trial([F(0)],degree)
        with self.assertRaises(ValueError): e.banded_inertia([[1,2],[0,1]])


class Certificates(unittest.TestCase):
    def test_dense_and_banded_certificates_identical(self):
        for q,k,n in [([0,1],1,6),([1,-2,3],2,6),([0,1000],1,2),([0],2,4)]:
            self.assertEqual(e.base.certify(q,k,n,20),e.certify_banded(q,k,n,20))

    def test_exact_exponential_family(self):
        for a in [F(-3,2),F(0),F(1),F(10)]:
            cert = e.barta_certificate([0,2*a,-a*a],[0,-a])
            self.assertEqual(F(cert['lower']),-a*a)
            self.assertEqual(cert['lower'],cert['upper'])
            self.assertTrue(e.verify(cert))

    def test_exponential_generator_rediscovers_matched_family(self):
        c = e.certify_exponential([0,20,-100],1)
        self.assertEqual(c['exact_width'],'0')
        self.assertEqual(c['lower'],'-100')

    def test_gaussian_function_family(self):
        c = e.certify_exponential([0,0,28,0,-16],2)
        self.assertEqual(c['exact_width'],'0')
        self.assertEqual(c['lower'],'4')

    def test_even_sector_sphere_second_eigenvalue(self):
        q = [F(0)]
        gap = e.even_separation(q,6,30)
        lower = e.gap_value(q,gap)
        self.assertLessEqual(lower,6)
        self.assertGreater(lower,F(599,100))

    def test_even_gap_cannot_apply_to_odd_trial(self):
        gap = {'type':'even_potential_minmax','subdivision_depth':3,'potential_lower':'0'}
        with self.assertRaises(ValueError): e.temple_certificate([0],[0,1],gap)
        with self.assertRaises(ValueError): e.temple_certificate([0,1],[1],gap)

    def test_even_ground_certificate_overlaps_full_schur(self):
        q = [0,0,-100]
        c = e.certify_temple(q,20,symmetry=True)
        d = e.certify_banded(q,1,20,40)
        self.assertLessEqual(max(F(c['lower']),F(d['lower'])),min(F(c['upper']),F(d['upper'])))
        self.assertGreater(F(c['gap_lower'])-F(c['rayleigh_quotient']),1)
        self.assertTrue(e.verify(c))

    def test_forged_even_sector_gap_rejected(self):
        c = e.certify_temple([0,0,-100],20,symmetry=True)
        c['gap_proof']['tail_lower'] = '1000000'
        self.assertFalse(e.verify(c))

    def test_normalization_does_not_change_local_energy(self):
        self.assertEqual(e.barta_certificate([1,2],[0,1]),e.barta_certificate([1,2],[100,1]))

    def test_barta_constant_trial_encloses_potential(self):
        c = e.barta_certificate([0,1],[0])
        self.assertEqual((c['lower'],c['upper']),('-1','1'))

    def test_local_energy_by_hand(self):
        # s=a*t+b*t^2; independent pointwise derivative formula.
        a,b = F(2,3),F(-1,5)
        p = e.local_energy([F(1),F(-2),F(3)],[0,a,b])
        for t in [F(-1),F(-1,3),F(0),F(2,5),F(1)]:
            expected = -(1-t*t)*(2*b+(a+2*b*t)**2)+2*t*(a+2*b*t)+1-2*t+3*t*t
            self.assertEqual(e.evaluate(p,t),expected)

    def test_full_residual_by_hand_includes_tail(self):
        self.assertEqual(e.residual_statistics([F(0),F(1)],[F(1)]),(F(0),F(0),F(1,3)))

    def test_sphere_and_constant_shift_temple(self):
        for shift in [F(0),F(-7,3),F(11,5)]:
            c = e.certify_temple([shift],4)
            self.assertEqual(F(c['lower']),shift)
            self.assertEqual(F(c['upper']),shift)
            self.assertEqual(c['excited_weight_upper'],'0')

    def test_residual_square_improves_weak_case(self):
        old = e.base.certify([0,1],1,6,36)
        c = e.certify_temple([0,1],6)
        self.assertLess(F(c['exact_width']),F(old['exact_width'])/100)
        self.assertTrue(e.verify(c))

    def test_independent_schur_and_temple_enclosures_overlap(self):
        for q in [[0,1],[1,-2,3],[0,20]]:
            a,b = e.certify_banded(q,1,10,40),e.certify_temple(q,10)
            self.assertLessEqual(max(F(a['lower']),F(b['lower'])),min(F(a['upper']),F(b['upper'])))

    def test_weak_trial_does_not_fabricate_gap(self):
        q = [F(0),F(1000)]
        gap = {'type':'potential_minmax','subdivision_depth':3,'potential_lower':'-1000'}
        with self.assertRaises(ValueError): e.temple_certificate(q,[1],gap)
        with self.assertRaises(ValueError): e.residual_statistics([F(0)],[F(0)])

    def test_temple_tampering_rejected(self):
        original = e.certify_temple([0,1],6)
        for key,value in [('tail_residual_squared','0'),('gap_lower','100'),('eigenvalue_index',2),
                          ('lower','0'),('excited_weight_upper','0'),('scope','all manifolds')]:
            c = deepcopy(original);c[key] = value
            self.assertFalse(e.verify(c),key)

    def test_forged_gap_proof_rejected(self):
        c = e.certify_temple([0,1],6)
        c['gap_proof']['potential_lower'] = '100'
        self.assertFalse(e.verify(c))

    def test_nested_gap_for_wrong_operator_rejected(self):
        c = e.certify_temple([0,20],8)
        c['gap_proof']['certificate'] = e.base.certify([100],2,4,16)
        self.assertFalse(e.verify(c))

    def test_barta_tampering_rejected(self):
        original = e.barta_certificate([0,20,-100],[0,-10])
        for key,value in [('s_power_coefficients',['0','0']),('local_energy_coefficients',['0']),
                          ('eigenvalue_index',2),('upper','-101'),('scope','all manifolds')]:
            c = deepcopy(original);c[key] = value
            self.assertFalse(e.verify(c),key)

    def test_reflection_symmetry(self):
        a,b = e.certify_temple([1,2,3],6),e.certify_temple([1,-2,3],6)
        self.assertEqual((a['lower'],a['upper']),(b['lower'],b['upper']))

    def test_auto_recognizes_exact_function(self):
        r = run.search([0,20,-100],tolerance=F(1,10**10))
        self.assertEqual(r['status'],'target_met')
        self.assertEqual(len(r['attempts']),1)
        self.assertEqual(r['certificate']['lower'],'-100')

    def test_unmet_target_is_explicit(self):
        r = run.search([0,1],'temple',F(1,10**30),max_modes=4)
        self.assertEqual(r['status'],'target_not_met')
        self.assertTrue(e.verify(r['certificate']))


if __name__ == '__main__': unittest.main()
