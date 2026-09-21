import copy
from fractions import Fraction as F
from math import factorial
import unittest
import function_models as fm


def profile(name='abs_power',axis=2,**options):
    return dict(kind='axis_profile',axis=axis,profile=name,**options)


def analytic(name='exp',argument=None,coefficient=1):
    return {'kind':'analytic_sum','terms':[{'function':name,'argument':argument or {},'coefficient':coefficient}]}


class FunctionModelTests(unittest.TestCase):
    def test_normalization_defaults(self):
        self.assertEqual(fm.normalize_function({'kind':'analytic_sum'}),{'kind':'analytic_sum','polynomial':{},'terms':[]})
        self.assertEqual(fm.normalize_function({'kind':'axis_profile','profile':'abs_power'}),
            {'kind':'axis_profile','profile':'abs_power','axis':2,'amplitude':'1','offset':'0','exponent':'1'})

    def test_canonical_term_combination_and_zero_removal(self):
        c={'kind':'analytic_sum','terms':[{'function':'sin','argument':{'1,0,0':1},'coefficient':2},
            {'function':'sin','argument':{'1,0,0':'1'},'coefficient':-2},{'function':'exp'}]}
        result=fm.analytic_model(c)
        self.assertEqual(result['polynomial'],{'0,0,0':'1'})
        self.assertEqual(result['error_upper'],'0')
        self.assertEqual(len(result['function']['terms']),1)
        self.assertTrue(fm.verify_model(result,c))

    def test_unknown_and_inexact_inputs_rejected(self):
        bad=[{'kind':'analytic_sum','extra':1},analytic(coefficient=1.0),analytic(argument={'0,0,1':True}),
            profile(exponent=-0.25),profile(axis=True),profile(amplitude=False),
            {'kind':'axis_profile','profile':'step','exponent':'1'},
            {'kind':'analytic_sum','terms':[{'function':'sqrt'}]}]
        for value in bad:
            with self.assertRaises(ValueError): fm.normalize_function(value)

    def test_degree_and_duplicate_xyz_rejection(self):
        for p in [{'0,0,25':1},{'1,0,0':1,'01,0,0':1},{(1,0,0):1}]:
            with self.assertRaises(ValueError): fm.polynomial(p)

    def test_sqrt_outward_rounding(self):
        for q in [F(0),F(9,16),F(2),F(1,10**100),F(10**60,3)]:
            upper=fm.sqrt_upper(q,40)
            self.assertGreaterEqual(upper*upper,q)
            if upper>0:self.assertLess((upper-F(1,2**40))**2,q)
        self.assertEqual(fm.sqrt_upper(F(9,16),2),F(3,4))

    def test_sqrt_and_degree_budget_reject(self):
        for q,bits in [(-1,40),(1,-1),(1,257),(1,True),(1.0,40)]:
            with self.assertRaises(ValueError): fm.sqrt_upper(q,bits)
        with self.assertRaises(ValueError): fm.l2_model(profile(),degree=25)
        with self.assertRaises(ValueError): fm.analytic_model(analytic(),order=True)

    def test_step_degree_zero_signed_amplitude_offset(self):
        c=fm.l2_model(profile('step',amplitude=-3,offset=2),degree=0)
        self.assertEqual(c['polynomial'],{'0,0,0':'1/2'})
        self.assertEqual(F(c['target_norm_squared']),F(5,2))
        self.assertEqual(F(c['error_squared']),F(9,4))
        self.assertEqual(F(c['error_upper']),F(3,2))
        self.assertEqual(c['point_value_at_axis_zero'],'2')
        self.assertTrue(fm.verify_model(c))

    def test_step_linear_projection_probability_constants(self):
        c=fm.l2_model(profile('step'),degree=1)
        self.assertEqual(c['polynomial'],{'0,0,0':'1/2','0,0,1':'3/4'})
        self.assertEqual(F(c['projection_norm_squared']),F(7,16))
        self.assertEqual(F(c['error_squared']),F(1,16))
        self.assertEqual(F(c['error_upper']),F(1,4))
        even=fm.l2_model(profile('step'),degree=2)
        self.assertEqual(c['polynomial'],even['polynomial'])

    def test_absolute_value_quadratic_projection(self):
        c=fm.l2_model(profile(exponent=1),degree=2)
        self.assertEqual(c['polynomial'],{'0,0,0':'3/16','0,0,2':'15/16'})
        self.assertEqual(F(c['error_squared']),F(1,192))

    def test_polynomial_profile_zero_error(self):
        c=fm.l2_model(profile(exponent=2,amplitude=-2,offset=3),degree=2)
        self.assertEqual(c['polynomial'],{'0,0,0':'3','0,0,2':'-2'})
        self.assertEqual(c['error_squared'],'0');self.assertEqual(c['error_upper'],'0')
        zero=fm.l2_model(profile(exponent=0,amplitude=-3,offset=3),degree=0)
        self.assertEqual(zero['polynomial'],{})
        self.assertEqual(zero['error_squared'],'0')

    def test_negative_quarter_profile_exact_reference_errors(self):
        expected={4:F(5202,43681),8:F(14450,159201),12:F(4118450,54066609),24:F(1335263295935522,24084192784943929)}
        for degree,error in expected.items():
            c=fm.l2_model(profile(exponent='-1/4',amplitude=-1),degree=degree)
            self.assertEqual(F(c['target_norm_squared']),2)
            self.assertEqual(F(c['error_squared']),error)
            self.assertEqual(F(c['monomial_source_moments'][0]),F(-4,3))
            self.assertTrue(fm.verify_model(c))

    def test_degree24_frontend_has_no_old_coefficient_cap(self):
        c=fm.l2_model(profile(exponent='-1/4',amplitude=-1),degree=24)
        self.assertGreater(max(abs(F(x)) for x in c['polynomial'].values()),10**6)
        self.assertTrue(fm.verify_model(c))

    def test_power_integrability_boundary(self):
        for alpha in ('-1/2','-3/4'):
            with self.assertRaises(ValueError):fm.l2_model(profile(exponent=alpha))
        c=fm.l2_model(profile(exponent='-49/100'),degree=4)
        self.assertEqual(F(c['target_norm_squared']),50)
        self.assertEqual(c['point_value_at_axis_zero'],'0')

    def test_offset_cross_term_is_not_omitted(self):
        c=fm.l2_model(profile(exponent='-1/4',amplitude=-2,offset=3),degree=0)
        self.assertEqual(F(c['target_norm_squared']),1)
        self.assertEqual(F(c['error_squared']),F(8,9))
        self.assertEqual(c['polynomial'],{'0,0,0':'1/3'})

    def test_each_axis_embedding_and_source_binding(self):
        for axis in (0,1,2):
            function=profile('step',axis=axis);c=fm.l2_model(function,degree=1)
            key=','.join(str(int(j==axis)) for j in range(3))
            self.assertEqual(c['polynomial'],{'0,0,0':'1/2',key:'3/4'})
            self.assertTrue(fm.verify_model(c,function))
            self.assertFalse(fm.verify_model(c,profile('step',axis=(axis+1)%3)))

    def test_error_monotone_under_exact_projection(self):
        errors=[F(fm.l2_model(profile(exponent='-1/4'),degree=d)['error_squared']) for d in (0,2,4,8,12,24)]
        self.assertEqual(errors,sorted(errors,reverse=True))

    def test_rodrigues_coefficients_independent_of_recurrence(self):
        c=fm.l2_model(profile(exponent='-1/4'),degree=12)
        for row in c['legendre_projection']:
            n=row['degree'];expected=[F(0)]*(n+1)
            for k in range(n//2+1):
                expected[n-2*k]=F((-1)**k*factorial(2*n-2*k),2**n*factorial(k)*factorial(n-k)*factorial(n-2*k))
            self.assertEqual(list(map(F,row['legendre_coefficients'])),expected)

    def test_projection_norm_by_independent_monomial_product(self):
        c=fm.l2_model(profile(exponent='-1/4',amplitude=-1,offset=2),degree=8)
        p=list(map(F,c['axis_power_coefficients']))
        expected=sum((a*b/F(i+j+1) for i,a in enumerate(p) for j,b in enumerate(p) if (i+j)%2==0),F(0))
        self.assertEqual(F(c['projection_norm_squared']),expected)

    def test_exact_expected_polynomial_binding(self):
        c=fm.l2_model(profile('step'),degree=1)
        self.assertTrue(fm.verify_model(c,expected_polynomial={'0,0,0':F(1,2),'0,0,1':F(3,4)}))
        self.assertFalse(fm.verify_model(c,expected_polynomial={'0,0,0':F(1,2)}))

    def test_l2_tampering_rejected(self):
        source=fm.l2_model(profile(exponent='-1/4'),degree=4)
        for field,value in [('error_squared','0'),('error_upper','0'),('norm','Linf'),('scope','area'),('target_norm_squared','1'),('degree',3),('sqrt_bits',20)]:
            c=copy.deepcopy(source);c[field]=value;self.assertFalse(fm.verify_model(c),field)
        c=copy.deepcopy(source);c['legendre_projection'].pop();self.assertFalse(fm.verify_model(c))

    def test_exp_geometric_full_tail(self):
        c=fm.analytic_model(analytic('exp',{'0,0,1':'1/4'}),order=2)
        self.assertEqual(c['polynomial'],{'0,0,0':'1','0,0,1':'1/4','0,0,2':'1/32'})
        self.assertEqual(F(c['error_upper']),F(1,360))
        self.assertEqual(c['term_proofs'][0]['tail']['method'],'absolute_exponential_series_geometric_tail')
        self.assertTrue(fm.verify_model(c))

    def test_exp_large_radius_lagrange_fallback(self):
        c=fm.analytic_model(analytic('exp',{'0,0,1':4}),order=0)
        self.assertEqual(F(c['error_upper']),324)
        self.assertEqual(c['term_proofs'][0]['tail']['exponential_envelope'],'81')
        self.assertTrue(fm.verify_model(c))

    def test_sine_multivariate_full_taylor(self):
        c=fm.analytic_model(analytic('sin',{'1,1,0':'1/4'}),order=4)
        self.assertEqual(c['polynomial'],{'1,1,0':'1/4','3,3,0':'-1/384'})
        self.assertEqual(F(c['error_upper']),F(1,122880))

    def test_log1p_series_and_boundary(self):
        c=fm.analytic_model(analytic('log1p',{'0,0,1':'1/4'}),order=3)
        self.assertEqual(c['polynomial'],{'0,0,1':'1/4','0,0,2':'-1/32','0,0,3':'1/192'})
        self.assertEqual(F(c['error_upper']),F(1,768))
        for radius in (1,-1,2):
            with self.assertRaises(ValueError):fm.analytic_model(analytic('log1p',{'0,0,1':radius}))

    def test_log_domain_checked_before_zero_or_cancellation(self):
        zero=fm.analytic_model({'kind':'analytic_sum'})
        for argument in ({'0,0,1':2},{'0,0,0':-1}):
            for coefficients in ((1,-1),(0,)):
                source={'kind':'analytic_sum','terms':[
                    {'function':'log1p','argument':argument,'coefficient':c}
                    for c in coefficients]}
                with self.subTest(argument=argument,coefficients=coefficients):
                    with self.assertRaises(ValueError):fm.normalize_function(source)
                    with self.assertRaises(ValueError):fm.analytic_model(source)
                    self.assertFalse(fm.verify_model(zero,expected_function=source))
        safe={'kind':'analytic_sum','terms':[
            {'function':'log1p','argument':{'0,0,1':'1/2'},'coefficient':c}
            for c in (1,-1)]}
        self.assertEqual(fm.normalize_function(safe),fm.normalize_function({'kind':'analytic_sum'}))
        self.assertTrue(fm.verify_model(zero,expected_function=safe))

    def test_signed_sum_triangle_error(self):
        function={'kind':'analytic_sum','polynomial':{'0,0,0':5},'terms':[
            {'function':'exp','argument':{'0,0,1':'1/4'},'coefficient':-2},
            {'function':'sin','argument':{'1,1,0':'1/4'},'coefficient':3},
            {'function':'log1p','argument':{'0,0,1':'1/4'},'coefficient':-1}]}
        c=fm.analytic_model(function,order=4)
        self.assertEqual(F(c['error_upper']),F(2,117760)+F(3,122880)+F(1,3840))
        self.assertTrue(fm.verify_model(c,function))

    def test_expanded_degree_counts_actual_nonzero_sine_terms(self):
        c=fm.analytic_model(analytic('sin',{'8,0,0':1}),order=4)
        self.assertEqual(c['polynomial'],{'8,0,0':'1','24,0,0':'-1/6'})
        with self.assertRaises(ValueError):fm.analytic_model(analytic('exp',{'0,0,2':1}),order=13)

    def test_analytic_missing_tail_and_source_tamper(self):
        source=fm.analytic_model(analytic('exp',{'0,0,1':'1/4'}),order=4)
        c=copy.deepcopy(source);del c['term_proofs'][0]['tail'];self.assertFalse(fm.verify_model(c))
        c=copy.deepcopy(source);c['error_upper']='0';self.assertFalse(fm.verify_model(c))
        c=copy.deepcopy(source);c['function']['terms'][0]['coefficient']='2';self.assertFalse(fm.verify_model(c))
        c=copy.deepcopy(source);c['unknown']=1;self.assertFalse(fm.verify_model(c))

    def test_wrong_model_family_is_rejected(self):
        with self.assertRaises(ValueError):fm.analytic_model(profile('step'))
        with self.assertRaises(ValueError):fm.l2_model(analytic())


if __name__=='__main__':unittest.main()
