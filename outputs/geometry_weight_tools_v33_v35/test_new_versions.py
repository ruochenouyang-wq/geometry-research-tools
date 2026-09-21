from fractions import Fraction as F
from itertools import product
from copy import deepcopy
import unittest
import global_core as core
import point_oracle as point
import simplex_v15 as simplex
import bernstein_v17 as bernstein
import equality_v19 as equality
import robust_v20 as robust
import constrained_v21 as constrained


def tilt():return core.problem([0],[[0,1]],[[-6,6]],[F(1,50)],[[F(1,5)]])

def leaf(tree):
    while 'axis' in tree:tree=tree['left']
    return tree


class V13AdaptiveTests(unittest.TestCase):
    def test_adaptive_removes_spectral_floor(self):
        p=core.problem([0],[[0,1]],[[-20,20]])
        fixed=point.Oracle(p,F(1,100000),max_modes=4);adaptive=point.Oracle(p,F(1,100000),max_modes=20)
        a=fixed.records[fixed.point([20])];b=adaptive.records[adaptive.point([20])]
        self.assertGreater(point.width(a),F(1,100000));self.assertLessEqual(point.width(b),F(1,100000))
        self.assertTrue(core.v3.verify(b['spectral_certificate'],True))
    def test_small_budget_remains_open(self):
        r=core.search(tilt(),F(1,100000000),version=13,max_leaves=1,max_modes=4)
        self.assertEqual(r['certificate']['status'],'global_gap_open');self.assertTrue(core.verify(r['certificate']))
    def test_point_upper_is_actual_trial(self):
        o=point.Oracle(tilt(),F(1,1000));r=o.records[o.point([3])]
        self.assertEqual(F(r['rayleigh_quotient']),core.v3.old.residual_statistics([F(0),F(3)],list(map(F,r['trial_legendre'])))[0])
    def test_odd_maximum_is_attempted(self):
        o=point.Oracle(tilt(),F(1,10**15),max_modes=5);o.point([6]);self.assertEqual(o.attempts[-1]['modes'],5)
    def test_invalid_budget_rejected(self):
        with self.assertRaises(ValueError):point.Oracle(tilt(),F(1,100),start_modes=2)
    def test_false_global_claim_rejected(self):
        c=core.search(tilt(),version=13,max_leaves=1)['certificate'];c['status']='epsilon_global';self.assertFalse(core.verify(c))


class V14OrbitTests(unittest.TestCase):
    def test_reflection_reuses_search_and_checks_dense(self):
        o=point.Oracle(tilt(),F(1,10000),orbits=True);a=o.records[o.point([-3])];calls=o.solves;b=o.records[o.point([3])]
        self.assertEqual(calls,o.solves);self.assertEqual(a['rayleigh_quotient'],b['rayleigh_quotient'])
        self.assertTrue(core.v3.verify(b['spectral_certificate'],True));self.assertGreater(o.reuses,0)
    def test_constant_shift_reuses(self):
        p=core.problem([0,2],[[1]],[[-2,2]])
        o=point.Oracle(p,F(1,10000),orbits=True);a=o.records[o.point([-2])];calls=o.solves;b=o.records[o.point([2])]
        self.assertEqual(calls,o.solves);self.assertEqual(F(b['rayleigh_quotient'])-F(a['rayleigh_quotient']),4)
        self.assertEqual(F(b['spectral_certificate']['lower'])-F(a['spectral_certificate']['lower']),4)
    def test_even_change_not_equivalent(self):
        p=core.problem([0,1],[[0,0,1]],[[-1,1]])
        o=point.Oracle(p,F(1,100),orbits=True);o.point([-1]);calls=o.solves;o.point([1]);self.assertGreater(o.solves,calls)
    def test_tampered_reflection_bound_rejected(self):
        o=point.Oracle(tilt(),F(1,100),orbits=True);r=deepcopy(o.records[o.point([3])]['spectral_certificate']);r['q_coefficients']=['0','4'];self.assertFalse(core.v3.verify(r,True))
    def test_orbit_path_same_global_interval(self):
        a=core.search(tilt(),version=13)['certificate'];b=core.search(tilt(),version=14)['certificate']
        self.assertEqual(a['global_lower'],b['global_lower']);self.assertEqual(a['candidate_upper'],b['candidate_upper'])


class V15SimplexTests(unittest.TestCase):
    def test_triangle_interior_minimum(self):
        value,x,_=simplex.minimum(F(1,2),[-1,-1],[[2,0],[0,2]],[[0,0],[2,0],[0,2]])
        self.assertEqual(value,0);self.assertEqual(x,[F(1,2)]*2)
    def test_singular_minimum_on_boundary(self):
        value,_,_=simplex.minimum(1,[-2,-2],[[2,2],[2,2]],[[0,0],[2,0],[0,2]])
        self.assertEqual(value,0)
    def test_indefinite_minimum_vertex(self):
        value,_,_=simplex.minimum(0,[0,0],[[-2,0],[0,-2]],[[0,0],[2,0],[0,2]])
        self.assertEqual(value,-4)
    def test_simplex_uses_fewer_points_on_trap(self):
        a=core.search(tilt(),version=14);b=core.search(tilt(),version=15)
        self.assertLess(b['statistics']['points'],a['statistics']['points']);self.assertEqual(b['certificate']['status'],'epsilon_global')
    def test_center_coordinate_must_match(self):
        c=core.search(tilt(),version=15,max_leaves=1)['certificate'];n=c['tree'];n['indices'][-1]=n['indices'][0];self.assertFalse(core.verify(c))
    def test_cell_bound_not_weaker(self):
        p=tilt();box=[list(map(F,r)) for r in p['box']];o=point.Oracle(p,F(1,1000));ids=[o.point(x) for x in core.positions(box)]
        self.assertGreaterEqual(core.cell_bound(p,box,o.records,ids,'center_simplex'),core.cell_bound(p,box,o.records,ids,'quadratic'))


class V16MomentTests(unittest.TestCase):
    def model(self):return core.all_real_problem([0,8],[[0,1],[0,0,1]],[0,0],[[F(1,5),0],[0,1]])
    def test_impossible_variance_cell(self):
        p=self.model();cut=core.moment_exclusion(p,[[F(4),F(5)],[F(-1,10),F(0)]])
        self.assertIsNotNone(cut);self.assertLess(F(cut['variance_upper']),0)
    def test_realizable_moments_not_excluded(self):
        p=self.model();self.assertIsNone(core.moment_exclusion(p,[[F(-1),F(1)],[F(-1,2),F(-1,4)]]))
    def test_cuts_are_actually_used(self):
        c=core.search(self.model(),F(1,10000),version=16,domain='all_real')['certificate']
        self.assertGreater(c['moment_excluded_leaves'],0);self.assertTrue(core.verify(c))
    def test_bounded_domain_cannot_use_stationarity(self):
        with self.assertRaises(ValueError):core.search(self.model(),version=16,moment_cuts=True)
    def test_mislabeled_domain_rejected(self):
        c=core.search(self.model(),F(1,10000),version=16,domain='all_real')['certificate'];c['domain']='box';self.assertFalse(core.verify(c))
    def test_wrong_moment_polynomial_gets_no_cut(self):
        p=self.model();p['directions'][1]=['0','0','2'];self.assertIsNone(core.moment_exclusion(p,[[F(4),F(5)],[F(0),F(1)]]))
    def test_nonpositive_hessian_rejected(self):
        with self.assertRaises(ValueError):core.all_real_problem([0],[[0,1]],[0],[[-1]])


class V17TensorTests(unittest.TestCase):
    def test_tensor_reproduces_polynomial(self):
        poly={(0,0,0):F(1),(1,1,1):F(-3),(2,0,0):F(2)};coeff=bernstein.coefficients(poly,3)
        for x in product((F(0),F(1,3),F(1)),repeat=3):
            value=sum(v*bernstein.product_value([F(__import__('math').comb(2,k))*z**k*(1-z)**(2-k) for z,k in zip(x,index)]) for index,v in coeff.items())
            self.assertEqual(value,bernstein.evaluate(poly,x))
    def test_lower_bound_covers_samples(self):
        poly={(2,0,0,0):F(2),(1,1,1,1):F(-4),(0,0,0,0):F(1)};low=min(bernstein.coefficients(poly,4).values())
        for x in product((F(0),F(1,2),F(1)),repeat=4):self.assertLessEqual(low,bernstein.evaluate(poly,x))
    def test_three_parameter_global_certificate(self):
        p=core.problem([0],[[0,1],[0,0,1],[0,0,0,1]],[[-1,1]]*3,[0]*3,[[2*int(i==j) for j in range(3)] for i in range(3)])
        c=core.search(p,F(1,5),version=17,max_leaves=64)['certificate'];self.assertEqual(c['status'],'epsilon_global');self.assertTrue(core.verify(c))
    def test_old_version_rejects_dimension(self):
        p=core.problem([0],[[0,1]]*3,[[-1,1]]*3)
        with self.assertRaises(ValueError):core.search(p,version=13)
    def test_forged_bernstein_lower_rejected(self):
        p=core.problem([0],[[0]]*3,[[-1,1]]*3)
        c=core.search(p,version=17,max_leaves=1)['certificate'];c['tree']['lower']='1';self.assertFalse(core.verify(c))
    def test_dimension_above_budget_rejected(self):
        with self.assertRaises(ValueError):core.problem([0],[[0]]*5,[[-1,1]]*5)


class V18SubdivisionTests(unittest.TestCase):
    def test_subdivision_removes_artificial_negative_bound(self):
        poly={(0,):F(1,4),(1,):F(-1),(2,):F(1)}
        raw=min(bernstein.coefficients(poly,1).values());low,tree=bernstein.refine(poly,1,2)
        self.assertEqual(raw,F(-1,4));self.assertEqual(low,0);self.assertEqual(bernstein.verify_refinement(poly,1,tree),0)
    def test_direct_expansion_matches_subdivision(self):
        poly={(2,1,0):F(-3,7),(1,0,2):F(5,2),(0,0,0):F(1,9)}
        low,tree=bernstein.refine(poly,3,11);self.assertEqual(low,bernstein.verify_refinement(poly,3,tree))
    def test_missing_child_rejected(self):
        _,tree=bernstein.refine({(2,):F(1)},1,2);del tree['right']
        with self.assertRaises(ValueError):bernstein.verify_refinement({(2,):F(1)},1,tree)
    def test_tampered_coefficient_rejected(self):
        poly={(2,):F(1)};_,tree=bernstein.refine(poly,1,2);leaf(tree)['lower']='20'
        with self.assertRaises(ValueError):bernstein.verify_refinement(poly,1,tree)
    def test_zero_polynomial_stops(self):
        low,tree=bernstein.refine({},4,128);self.assertEqual(low,0);self.assertEqual(tree,{'lower':'0'})
    def test_nested_refinement_monotone(self):
        poly={(2,0):F(1),(1,1):F(-2),(0,2):F(1)}
        lows=[bernstein.refine(poly,2,n)[0] for n in (1,2,4,8)];self.assertEqual(lows,sorted(lows))


class V19EqualityTests(unittest.TestCase):
    def model(self):return core.problem([0],[[0,1],[0,0,1]],[[-1,1]]*2,[0,0],[[1,0],[0,1]])
    def test_affine_map_satisfies_all_rows(self):
        matrix=[[1,2,3],[2,4,6]];rhs=[4,8];a,b=equality.affine_space(matrix,rhs,3)
        for z in ([F(0),F(0)],[F(3,7),F(-4)]):
            x=[t+sum(v*w for v,w in zip(row,z)) for t,row in zip(a,b)]
            for row,y in zip(matrix,rhs):self.assertEqual(sum(v*w for v,w in zip(row,x)),y)
    def test_inconsistent_equalities_rejected(self):
        with self.assertRaises(ValueError):equality.affine_space([[1,1],[2,2]],[1,3],2)
    def test_quadratic_and_operator_transform_exact(self):
        p=self.model();reduced,a,b=equality.transform(p,[[1,1]],[F(1,2)])
        for z in [F(-1),F(2,3)]:
            x=[v+row[0]*z for v,row in zip(a,b)]
            self.assertEqual(core.at(p,x),core.at(reduced,[z]));self.assertEqual(core.penalty(p,x),core.penalty(reduced,[z]))
    def test_global_certificate_and_lift(self):
        c=equality.search(self.model(),[[1,1]],[0])['certificate'];self.assertEqual(c['status'],'epsilon_global');self.assertTrue(equality.verify(c));self.assertEqual(sum(map(F,c['candidate_parameters'])),0)
    def test_corrupt_basis_rejected(self):
        c=equality.search(self.model(),[[1,1]],[0])['certificate'];c['basis'][0][0]='2';self.assertFalse(equality.verify(c))
    def test_full_rank_has_explicit_limit(self):
        with self.assertRaises(ValueError):equality.affine_space([[1,0],[0,1]],[0,0],2)


class V20RobustTests(unittest.TestCase):
    def model(self):return core.problem([0],[[0,1]],[[-1,1]],[0],[[1]])
    def test_asymmetric_robust_global_certificate(self):
        c=robust.search(self.model(),[[0,1]],[[-1,2]])['certificate'];self.assertEqual(c['status'],'epsilon_global');self.assertTrue(robust.verify(c));self.assertLess(F(c['candidate_parameters'][0]),0)
    def test_symmetric_optimum_at_zero(self):
        c=robust.search(self.model(),[[0,1]],[[-1,1]])['certificate'];self.assertEqual(c['status'],'epsilon_global');self.assertEqual(c['candidate_parameters'],['0'])
    def test_negative_dual_weight_rejected(self):
        c=robust.search(self.model(),[[0,1]],[[-1,2]])['certificate'];c['weights'][0]='-1';self.assertFalse(robust.verify(c))
    def test_missing_uncertainty_vertex_rejected(self):
        c=robust.search(self.model(),[[0,1]],[[-1,2]])['certificate'];c['scenario_records'].pop();self.assertFalse(robust.verify(c))
    def test_forged_rayleigh_plane_rejected(self):
        c=robust.search(self.model(),[[0,1]],[[-1,2]])['certificate'];c['cuts'][0]['constant']='100';self.assertFalse(robust.verify(c))
    def test_budget_exhaustion_stays_open(self):
        c=robust.search(self.model(),[[0,1]],[[-1,2]],F(1,1000000),max_steps=1)['certificate'];self.assertEqual(c['status'],'global_gap_open');self.assertTrue(robust.verify(c))


class V21ConstrainedTests(unittest.TestCase):
    def test_full_function_constraint_global_gap(self):
        c=constrained.search([0],[0,1],F(1,2))['certificate'];self.assertEqual(c['status'],'epsilon_global');self.assertTrue(constrained.verify(c));self.assertGreaterEqual(F(c['candidate_moment']),F(1,2));self.assertLess(F(c['candidate_upper']),F(1,2))
    def test_impossible_target_proved(self):
        c=constrained.search([0],[0,1],2)['certificate'];self.assertEqual(c['status'],'certified_infeasible');self.assertTrue(constrained.verify(c))
    def test_no_feasible_trial_not_called_infeasible(self):
        c=constrained.search([0],[0,1],F(1,2),max_steps=1)['certificate'];self.assertEqual(c['status'],'feasibility_not_found');self.assertTrue(constrained.verify(c))
    def test_negative_multiplier_rejected(self):
        c=constrained.search([0],[0,1],F(1,2))['certificate'];c['multiplier']='-1';self.assertFalse(constrained.verify(c))
    def test_infeasible_trial_rejected(self):
        c=constrained.search([0],[0,1],F(1,2))['certificate'];c['trial_legendre']=['1'];self.assertFalse(constrained.verify(c))
    def test_wrong_tilted_operator_rejected(self):
        c=constrained.search([0],[0,1],F(1,2))['certificate'];c['target']='2';self.assertFalse(constrained.verify(c))
    def test_blended_trial_has_exact_feasibility(self):
        p=[F(0),F(1)];c=constrained.blend_feasible(p,F(1,3),[F(1)],[F(1),F(1)]);self.assertGreaterEqual(constrained.expectation(p,c),F(1,3))


class V22SchedulerTests(unittest.TestCase):
    def test_distinct_objective_requires_its_entry_point(self):
        for version in (20,21):
            with self.assertRaises(ValueError):core.search(tilt(),version=version)
    def test_scheduler_selects_spectral_action(self):
        p=core.problem([0],[[0,1]],[[-4,4]],[0],[[F(1,5)]])
        r=core.search(p,version=22)
        self.assertGreater(r['actions']['point_refinements'],0);self.assertTrue(core.verify(r['certificate']))
    def test_four_parameter_refined_certificate(self):
        p=core.problem([0],[[0]*i+[1] for i in range(1,5)],[[F(-1,4),F(1,4)]]*4,
                       [0]*4,[[int(i==j) for j in range(4)] for i in range(4)])
        r=core.search(p,F(1,4),version=22,max_leaves=4)
        self.assertEqual(r['certificate']['status'],'epsilon_global');self.assertTrue(core.verify(r['certificate']))
    def test_scheduled_global_trap(self):
        r=core.search(tilt(),version=22);self.assertEqual(r['certificate']['status'],'epsilon_global');self.assertTrue(core.verify(r['certificate']))
    def test_polynomial_refinement_action_used(self):
        p=core.problem([0],[[0]]*3,[[-1,1]]*3,[0]*3,[[2*int(i==j) for j in range(3)] for i in range(3)])
        r=core.search(p,F(1,100),version=22,max_leaves=1,refinement_budget=8)
        self.assertGreater(r['actions']['polynomial_refinements'],0);self.assertEqual(r['certificate']['status'],'epsilon_global')
    def test_action_budget_not_a_proof(self):
        r=core.search(tilt(),F(1,10**8),version=22,max_actions=1)
        self.assertEqual(r['certificate']['status'],'global_gap_open');self.assertTrue(core.verify(r['certificate']))
    def test_split_coverage_tamper_rejected(self):
        c=core.search(tilt(),version=22)['certificate'];c['tree']['split']='7';self.assertFalse(core.verify(c))
    def test_actual_point_refinement_keeps_validity(self):
        o=point.Oracle(tilt(),F(1,10),max_modes=16);idx=o.point([6]);old=point.width(o.records[idx]);self.assertTrue(o.refine(idx,F(1,10**8)));self.assertLess(point.width(o.records[idx]),old);self.assertTrue(core.v3.verify(o.records[idx]['spectral_certificate'],True))
    def test_invalid_scheduler_budget_rejected(self):
        with self.assertRaises(ValueError):core.search(tilt(),version=22,max_actions=0)

if __name__=='__main__':unittest.main()
