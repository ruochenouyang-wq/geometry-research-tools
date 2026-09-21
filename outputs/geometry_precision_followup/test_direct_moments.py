import unittest
from copy import deepcopy
from unittest.mock import patch
import direct_moments as d

STEP={'kind':'axis_profile','profile':'step'}
SING={'kind':'axis_profile','profile':'abs_power','exponent':'-1/4','amplitude':'-1'}


class DirectMomentsTests(unittest.TestCase):
    def test_exact_step_moments(self):
        for j in range(7):
            self.assertEqual(d.moment(STEP,j),d.F(1,j+1))
            self.assertEqual(d.moment(STEP,j,2),d.F(1,j+1))

    def test_exact_singular_moments(self):
        for j in range(7):
            self.assertEqual(d.moment(SING,j),d.F(-8,4*j+3) if j%2==0 else 0)
            self.assertEqual(d.moment(SING,j,2),d.F(4,2*j+1) if j%2==0 else 0)
        self.assertEqual(d.form_bound(SING)['variance'],'2/9')

    def test_negative_step_amplitude_pointwise_bound(self):
        q={**STEP,'amplitude':'-3','offset':'7/2'}
        self.assertEqual(d.form_bound(q)['mass_offset'],'1/2')

    def test_associated_basis_orthogonality(self):
        for m in range(3):
            for l in range(m,m+5):
                for k in range(m,m+5):
                    got=d.integral_product(STEP,m,l,k,0)
                    self.assertEqual(got,d.wide.basis_mass(m,l) if l==k else 0)

    def test_removed_constant_is_required(self):
        data=d.matrix_assembly(STEP,modes=1)
        self.assertEqual(data['removed_constant']['column'],[d.F(1,2)])
        self.assertEqual(data['C'],[[d.F(1,24)]])
        self.assertIsNone(d.matrix_assembly(STEP,m=1,modes=1)['removed_constant'])

    def test_constant_potential_has_zero_gram(self):
        for mean_zero in (True,False):
            k=d.Kernel({**STEP,'amplitude':'0','offset':'-7/3'},mean_zero=mean_zero,modes=3,near_tail=3)
            self.assertTrue(all(v==0 for row in k.data['C'] for v in row))
            self.assertTrue(all(v==0 for row in k.data['near_tail'] for v in row['column']))
            lam=d.F(2 if mean_zero else 0)-d.F(7,3)
            cert=d.sector_certificate(k,lam,lam)
            self.assertTrue(d.verify_sector(cert))

    def test_tail_refinement_is_safe_and_tighter(self):
        for q in (STEP,SING):
            a=d.Kernel(q,modes=4,near_tail=0)
            b=d.Kernel(q,modes=4,near_tail=5)
            ma,mb=a.matrix('0'),b.matrix('0')
            difference=[[mb[i][j]-ma[i][j] for j in range(4)] for i in range(4)]
            self.assertEqual(d.base.inertia(difference)[0],0)

    def test_strict_tail_boundary(self):
        k=d.Kernel(STEP,modes=2)
        with self.assertRaises(ValueError):k.matrix(k.beta)

    def test_input_rejections(self):
        for q in ({**STEP,'extra':1},{**SING,'exponent':'-1/2'},
                  {**SING,'amplitude':'-4'},{**STEP,'amplitude':0.1},
                  {**STEP,'axis':True},{'kind':'analytic_sum'},None):
            with self.subTest(q=q),self.assertRaises((ValueError,TypeError)):
                d.Kernel(q,modes=2)

    def test_integer_and_precision_budgets(self):
        for arguments in ({'m':True},{'modes':0},{'near_tail':33},{'mean_zero':1},{'sqrt_bits':1}):
            with self.subTest(arguments=arguments),self.assertRaises(ValueError):d.Kernel(STEP,**arguments)
        with self.assertRaises(ValueError):d.certify_sector(STEP,bits=1)
        with self.assertRaises(ValueError):d.full_ground(STEP,max_m=13)

    def test_endpoint_tampering(self):
        c=d.certify_sector(STEP,modes=3,bits=16)
        for key,value in [('lower',c['upper']),('lower_inertia',[0,0,99]),('function_approximation_error','1')]:
            bad=deepcopy(c);bad[key]=value
            self.assertFalse(d.verify_sector(bad))

    def test_moment_gram_and_form_tampering(self):
        c=d.certify_sector(SING,modes=3,bits=16,near_tail=2)
        for key in ('A','V','T','C'):
            bad=deepcopy(c);bad['kernel_evidence'][key][0][0]='0'
            self.assertFalse(d.verify_sector(bad))
        bad=deepcopy(c);bad['form_bound']['eta_upper']='0'
        self.assertFalse(d.verify_sector(bad))

    def test_full_space_and_binding(self):
        c=d.full_ground(STEP,modes=4,bits=16)
        self.assertTrue(d.verify_full(c,STEP,True,'1/100000000'))
        self.assertEqual(len(c['sectors']),2)
        self.assertFalse(d.verify_full(c,SING))
        self.assertFalse(d.verify_full(c,expected_mean_zero=False))
        self.assertFalse(d.verify_full(c,expected_tolerance='1/100'))

    def test_omitted_m_tail_cannot_be_ignored(self):
        c=d.full_ground(SING,modes=3,bits=16,max_m=0)
        self.assertFalse(c['angular_tail_cannot_improve_best_upper'])
        self.assertEqual(c['lower'],c['angular_tail_lower'])
        bad=deepcopy(c);bad['lower']=c['sectors'][0]['lower']
        self.assertFalse(d.verify_full(bad))

    def test_reordered_or_missing_sector_rejected(self):
        c=d.full_ground(STEP,modes=3,bits=16)
        bad=deepcopy(c);bad['sectors'].reverse()
        self.assertFalse(d.verify_full(bad))
        bad=deepcopy(c);bad['sectors']=bad['sectors'][1:]
        self.assertFalse(d.verify_full(bad))

    def test_axis_isometry(self):
        bounds=[]
        for axis in range(3):
            c=d.full_ground({**STEP,'axis':axis},modes=3,bits=16)
            bounds.append((c['lower'],c['upper']))
        self.assertEqual(bounds,[bounds[0]]*3)

    def test_additive_constant_shift(self):
        c=d.full_ground(STEP,modes=3,bits=16)
        b=d.full_ground({**STEP,'offset':'7/3'},modes=3,bits=16)
        for name in ('lower','upper'):self.assertEqual(d.F(b[name])-d.F(c[name]),d.F(7,3))

    def test_verifier_never_calls_search(self):
        c=d.full_ground(SING,modes=3,bits=16)
        with patch.object(d,'full_ground',side_effect=AssertionError('search forbidden')),patch.object(d,'certify_sector',side_effect=AssertionError('search forbidden')):
            self.assertTrue(d.verify_full(c))

    def test_wrong_type_certificates_return_false(self):
        for bad in (None,[],True,1,'proof',{}):
            self.assertFalse(d.verify_sector(bad))
            self.assertFalse(d.verify_full(bad))


if __name__=='__main__':unittest.main()
