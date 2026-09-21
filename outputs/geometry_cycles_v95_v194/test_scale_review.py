"""Independent analytic calibration and adversarial replay of radius transfers."""
import copy
import json
from pathlib import Path
import unittest
from fractions import Fraction as F
import sphere_scale as s
from common import projected
import ordered_spectrum as ordered


SMALL=dict(modes=2,max_modes=2,max_m=2,bits=18)


def loose_zero_point(unit_lower,unit_upper,radius=2):
    """Legitimate deliberately non-nested enclosures, not forged certificates."""
    kernel=projected.Kernel([0],0,True,2)
    sector=projected.sector_certificate(kernel,1,F(unit_lower),F(unit_upper),'upper_ritz')
    full=projected.full_certificate([0],True,[sector],kernel.range_proof,F(1,10**6))
    return s.scale_ground([0],radius,full)


class ScaleIndependentReview(unittest.TestCase):
    def test_physical_coordinate_pullback_powers(self):
        c=s.physical_pullback([1,-2,3],F(3,2))
        self.assertEqual(list(map(F,c['pulled_back_potential'])),[F(1),F(-3),F(27,4)])
        self.assertEqual(list(map(F,c['unit_operator_potential'])),[F(9,4),F(-27,4),F(243,16)])
        self.assertTrue(s.verify(c))

    def test_pullback_external_problem_binding(self):
        c=s.physical_pullback([0,1],2)
        self.assertTrue(s.verify(c,expected_q=[0,1],expected_coordinates='ambient_x3',expected_radius=2))
        self.assertFalse(s.verify(c,expected_q=[0,2]))
        self.assertFalse(s.verify(c,expected_coordinates='normalized_t'))

    def test_constant_potential_ground_radius_formula(self):
        for r,c in [(F(2),F(-7)),(F(3,2),F(5,3)),(F(1,2),F(-9))]:
            proof=s.backend_ground([c],r,**SMALL)
            exact=2/r**2+c
            self.assertEqual(F(proof['lower']),exact)
            self.assertEqual(F(proof['upper']),exact)
            self.assertTrue(s.verify(proof,expected_q=[c],expected_radius=r))

    def test_unit_potential_binding_rejects_missing_radius_power(self):
        wrong=projected.full_ground([0,0,4],modes=2,max_modes=2,max_m=2,bits=18)
        with self.assertRaises(ValueError): s.scale_ground([0,0,1],2,wrong,'ambient_x3')
        correct=s.scale_ground([0,0,1],2,wrong,'normalized_t')
        self.assertTrue(s.verify(correct,expected_coordinates='normalized_t'))
        self.assertFalse(s.verify(correct,expected_coordinates='ambient_x3'))

    def test_natural_area_mass_and_energy_scale(self):
        c=s.scale_trial([3],2,[1],m=0,degrees=[1])
        # P1(t)=t; integrals in t are 2/3 and 4/3. Azimuth 2*pi cancels.
        self.assertEqual(F(c['mass_divided_by_azimuth_factor']),F(8,3))
        self.assertEqual(F(c['dirichlet_energy_divided_by_azimuth_factor']),F(4,3))
        self.assertEqual(F(c['potential_energy_divided_by_azimuth_factor']),F(8))
        self.assertEqual(F(c['rayleigh_quotient']),F(7,2))
        self.assertTrue(s.verify(c))

    def test_physical_vs_normalized_quadratic_trial(self):
        ambient=s.scale_trial([0,0,1],2,[1],m=0,degrees=[1],coordinates='ambient_x3')
        normal=s.scale_trial([0,0,1],2,[1],m=0,degrees=[1],coordinates='normalized_t')
        # Potential average with weight t² is 3/5; physical x3² adds R².
        self.assertEqual(F(ambient['rayleigh_quotient']),F(29,10))
        self.assertEqual(F(normal['rayleigh_quotient']),F(11,10))
        self.assertEqual(F(ambient['potential_energy_divided_by_azimuth_factor']),F(32,5))
        self.assertTrue(s.verify(ambient))

    def test_nonaxisymmetric_trial_azimuth_factor(self):
        c=s.scale_trial([0,0,1],2,[1],m=1,degrees=[1],coordinates='ambient_x3')
        # P1^1=-sqrt(1-t²), potential weighted mean t²=1/5.
        self.assertEqual(F(c['mass_divided_by_azimuth_factor']),F(16,3))
        self.assertEqual(F(c['dirichlet_energy_divided_by_azimuth_factor']),F(8,3))
        self.assertEqual(F(c['rayleigh_quotient']),F(13,10))
        self.assertTrue(s.verify(c))

    def test_mixed_legendre_trial_cross_terms_by_direct_integration(self):
        # P1+(2/5)P3=(2/5)t+t³. Independent power convolution includes every
        # cross term, instead of reusing the implementation's Legendre action.
        def product(a,b):
            out=[F(0)]*(len(a)+len(b)-1)
            for i,x in enumerate(a):
                for j,y in enumerate(b): out[i+j]+=x*y
            return out
        def integral(p):
            return sum((2*x/F(i+1) for i,x in enumerate(p) if i%2==0),F(0))
        r=F(3,2);u=[F(0),F(2,5),F(0),F(1)]
        square=product(u,u);du=[u[1],0,3*u[3]]
        mass=r*r*integral(square)
        kinetic=integral(product([1,0,-1],product(du,du)))
        pulled_operator=[r*r,-2*r**3,3*r**4]
        potential=integral(product(pulled_operator,square))
        c=s.scale_trial([1,-2,3],r,[1,F(2,5)],m=0,degrees=[1,3],coordinates='ambient_x3')
        self.assertEqual(F(c['mass_divided_by_azimuth_factor']),mass)
        self.assertEqual(F(c['dirichlet_energy_divided_by_azimuth_factor']),kinetic)
        self.assertEqual(F(c['potential_energy_divided_by_azimuth_factor']),potential)
        self.assertEqual(F(c['rayleigh_quotient']),(kinetic+potential)/mass)

    def test_ordered_negative_trace_scales_but_count_does_not(self):
        req=ordered.normalize_request([-28],k=8,modes=2,max_modes=2,max_m=0,max_radial=1,bits=18)
        original=ordered.ordered_bounds(req,ordered.assemble_sectors([-28],[]))
        c=s.scale_ordered([-7],2,original)
        self.assertEqual([F(x['lower']) for x in c['ordered_intervals']],[F(-13,2)]*3+[F(-11,2)]*5)
        # l=1..4 are negative: multiplicities 3+5+7+9=24, trace=93.
        self.assertEqual(c['negative_count_lower'],24)
        self.assertEqual(c['negative_count_upper'],24)
        self.assertEqual(F(c['negative_trace_lower']),93)
        self.assertEqual(F(c['negative_trace_upper']),93)
        self.assertTrue(s.verify(c))

    def test_zero_eigenvalue_is_not_negative_after_scaling(self):
        req=ordered.normalize_request([-2],k=3,modes=2,max_modes=2,max_m=0,max_radial=1,bits=18)
        c=s.scale_ordered([F(-1,2)],2,ordered.ordered_bounds(req,ordered.assemble_sectors([-2],[])))
        self.assertEqual(c['negative_count_upper'],0)
        self.assertEqual(F(c['negative_trace_upper']),0)
        self.assertEqual([F(x['upper']) for x in c['ordered_intervals']],[0,0,0])

    def test_negative_eigenvalue_four_corner_interval(self):
        anchor=s.backend_ground([-10,1],2,**SMALL)
        c=s.radius_cell([-10,1],[1,2],anchor)
        L,U=map(F,c['unit_eigenvalue_interval'])
        self.assertLess(U,0)
        corners=[L/4,L,U/4,U]
        self.assertEqual(list(map(F,c['four_signed_products'])),corners)
        self.assertEqual(F(c['lower']),min(corners))
        self.assertEqual(F(c['upper']),max(corners))
        self.assertNotEqual(F(c['upper']),U)  # The smaller positive multiplier wins.
        self.assertTrue(s.verify(c))

    def test_sign_crossing_interval_four_corners(self):
        anchor=s.backend_ground([-2,1],F(3,2),**SMALL)
        c=s.radius_cell([-2,1],[1,2],anchor)
        L,U=map(F,c['unit_eigenvalue_interval'])
        self.assertLess(L,0);self.assertGreater(U,0)
        self.assertEqual(F(c['lower']),L)
        self.assertEqual(F(c['upper']),U)

    def test_constant_uniform_cell_matches_analytic_radius_range(self):
        c=s.radius_cell([-7],[1,2],s.backend_ground([-7],F(3,2),**SMALL))
        self.assertEqual((F(c['lower']),F(c['upper'])),(F(-13,2),F(-5)))
        self.assertTrue(s.verify(c))

    def test_minimum_incumbent_is_point_not_maximum_upper(self):
        c=s.adaptive_minimum([-7],[1,2],max_splits=0,**SMALL)['certificate']
        self.assertEqual(c['best_feasible_radius'],'2')
        self.assertEqual(F(c['upper']),F(-13,2))
        self.assertEqual(F(c['maximum_upper']),F(-5))
        self.assertEqual(F(c['lower']),F(c['upper']))
        self.assertTrue(s.verify(c))

    def test_angular_budget_and_radius_budget_return_open_bound(self):
        c=s.adaptive_minimum([0,1],[1,2],max_splits=0,tolerance=F(1,10**6),
            modes=2,max_modes=2,max_m=0,bits=18)['certificate']
        self.assertEqual(c['status'],'certified_global_bound_open_gap')
        self.assertGreater(F(c['width']),F(c['tolerance']))
        self.assertTrue(s.verify(c))
        # Honest wide bounds still enclose independently evaluated interior points.
        for r in (F(1),F(4,3),F(2)):
            p=s.backend_ground([0,1],r,**SMALL)
            self.assertLessEqual(F(c['lower']),F(p['upper']))
            self.assertLessEqual(F(p['lower']),F(c['maximum_upper']))

    def test_complete_domain_and_checkpoint_binding(self):
        c=s.adaptive_minimum([0,1],[1,2],max_splits=1,**SMALL)['certificate']
        self.assertTrue(s.verify(c))
        missing=copy.deepcopy(c);missing['cells'].pop()
        self.assertFalse(s.verify(missing))
        changed=copy.deepcopy(c);changed['domain']=['1','3']
        self.assertFalse(s.verify(changed))
        saved=s.checkpoint(c)
        self.assertTrue(s.verify_checkpoint(saved,expected_q=[0,1],expected_domain=[1,2]))
        self.assertFalse(s.verify_checkpoint(saved,expected_domain=[1,3]))
        saved['certificate_digest']='0'*64
        self.assertFalse(s.verify_checkpoint(saved))

    def test_cache_problem_binding_and_copy_isolation(self):
        point=s.backend_ground([0,1],2,**SMALL)
        cache=s.AnchorCache([0,1],entries=[point],**SMALL)
        a=cache.get(2);a['upper']='999'
        self.assertTrue(s.verify(cache.get(2)))
        self.assertEqual(cache.computed,0)
        self.assertEqual(cache.reused,2)
        with self.assertRaises(ValueError): s.AnchorCache([0,2],entries=[point],**SMALL)
        with self.assertRaises(ValueError): s.AnchorCache([0,1],coordinates='ambient_x3',entries=[point],**SMALL)

    def test_nonnested_cache_cannot_discard_old_incumbent(self):
        old=loose_zero_point(0,2);narrow=loose_zero_point(F(19,10),F(21,10))
        old_cache=s.AnchorCache([0],entries=[old],**SMALL)
        before=s.adaptive_minimum([0],[1,2],max_splits=0,cache=old_cache)['certificate']
        self.assertEqual(F(before['upper']),F(1,2))
        cache=s.AnchorCache([0],entries=[narrow],**SMALL)
        after=s.adaptive_minimum([0],[1,2],max_splits=0,cache=cache,saved=s.checkpoint(before))['certificate']
        self.assertTrue(s.verify(after))
        self.assertLessEqual(F(after['upper']),F(before['upper']))
        self.assertGreaterEqual(F(after['lower']),F(before['lower']))

    def test_resume_keeps_all_old_bounds(self):
        first=s.adaptive_minimum([0,1],[1,2],max_splits=0,**SMALL)['certificate']
        later=s.resume_minimum(s.checkpoint(first),max_splits=1,**SMALL)['certificate']
        self.assertTrue(s.verify(later))
        self.assertGreaterEqual(F(later['lower']),F(first['lower']))
        self.assertLessEqual(F(later['upper']),F(first['upper']))

    def test_threshold_witness_strictness_and_radius(self):
        trial=s.scale_trial([0],2,[1],m=0,degrees=[1])
        equal=s.transfer_counterexample(trial,F(1,2))
        strict=s.transfer_counterexample(trial,F(3,4))
        self.assertEqual(equal['status'],'trial_does_not_refute')
        self.assertEqual(strict['status'],'explicit_counterexample')
        self.assertEqual(F(strict['strict_margin']),F(1,4))
        self.assertTrue(s.verify(strict,expected_q=[0],expected_radius=2,expected_threshold=F(3,4)))
        self.assertFalse(s.verify(strict,expected_threshold=F(1,2)))
        wrong=copy.deepcopy(strict);wrong['radius']='1'
        self.assertFalse(s.verify(wrong))

    def test_whole_radius_assertion_refuted_at_specific_radius(self):
        c=s.research_job([0],F(3,4),domain=[1,2],max_splits=0,
            trial={'coefficients':[1],'m':0,'degrees':[1],'radius':2},**SMALL)['certificate']
        self.assertEqual(c['status'],'refuted_by_explicit_function')
        self.assertTrue(s.verify(c,expected_q=[0],expected_threshold=F(3,4)))
        with self.assertRaises(ValueError):
            s.research_job([0],F(3,4),domain=[1,2],max_splits=0,
                trial={'coefficients':[1],'m':0,'degrees':[1],'radius':3},**SMALL)

    def test_point_job_external_radius_binding(self):
        c=s.research_job([0],F(1,2),radius=2,**SMALL)['certificate']
        self.assertTrue(s.verify(c,expected_radius=2))
        self.assertFalse(s.verify(c,expected_radius=1))
        d=s.research_job([0],F(1,2),domain=[1,2],max_splits=0,**SMALL)['certificate']
        self.assertFalse(s.verify(d,expected_radius=2))

    def test_nonpositive_radius_zero_trial_and_mean_not_allowed(self):
        for r in (0,-1):
            with self.assertRaises(ValueError): s.scale_trial([0],r,[1],m=0,degrees=[1])
        with self.assertRaises(ValueError): s.scale_trial([0],1,[0],m=0,degrees=[1])
        with self.assertRaises(ValueError): s.scale_trial([0],1,[1],m=0,degrees=[0])

    def test_saved_certificates_replay(self):
        directory=Path(__file__).resolve().parent/'round10'
        kinds={s.POINT,s.PULLBACK,s.TRIAL,s.ORDERED,s.CELL,s.GLOBAL,s.CHECKPOINT,s.VIOLATION,s.JOB,'radius_cell_intersection_v191'}
        seen=set();count=0
        def walk(obj):
            nonlocal count
            if isinstance(obj,dict):
                if obj.get('format') in kinds:
                    fingerprint=json.dumps(obj,sort_keys=True)
                    if fingerprint not in seen:
                        self.assertTrue(s.verify(obj),obj.get('format'))
                        seen.add(fingerprint);count+=1
                for value in obj.values(): walk(value)
            elif isinstance(obj,list):
                for value in obj: walk(value)
        # RESULTS.json contains abbreviated {format,path,verified} index rows,
        # not certificates. Replay the actual saved proof documents it indexes.
        paths=list((directory/'certificates').glob('*.json'))+list(directory.glob('*candidate.json'))
        for path in paths: walk(json.loads(path.read_text()))
        self.assertGreater(count,0,'Awaiting saved round10 certificates from the implementation author')


if __name__=='__main__': unittest.main()
