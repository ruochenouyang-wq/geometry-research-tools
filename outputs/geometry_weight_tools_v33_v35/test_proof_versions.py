from fractions import Fraction as F
from copy import deepcopy
from itertools import product
import unittest
import compiler_v23 as compiler
import reductions as reduction
import algebra
import family_v29 as family
import matrix_v30 as matrix
import transport_v31 as transport
import decision_v28 as decision
import planner_v32 as planner
import error_bounds as e


def goal(potential='a*t',penalty='a*a/4',parameters=None,box=None,threshold='0',scope='full_sphere'):
    return {'geometry':'unit_sphere','function_space':scope,'eigenvalue_index':1,'parameters':['a'] if parameters is None else parameters,
            'potential':potential,'penalty':penalty,'domain':{'kind':'all_real'} if box is None else {'kind':'box','box':box},'threshold':threshold}


def problem(*args,**kwargs):return compiler.compile_goal(goal(*args,**kwargs))['problem']


def energy(p,x,trial):
    q0,dirs,c,l,h=reduction.data(p)
    for v,q in zip(x,dirs):q0=e.add(q0,e.scale(q,v))
    return e.residual_statistics(q0,trial)[0]+algebra.value(c,l,h,x)


class V23Compiler(unittest.TestCase):
    def test_exact_expansion_and_cancellation(self):
        c=compiler.compile_goal(goal('(a+1)*t-a*t+t*t','(a+1)**2-a*a'))
        self.assertEqual(c['problem']['q0'],['0','1','1']);self.assertEqual(c['problem']['directions'],[['0']]);self.assertEqual(c['problem']['penalty']['linear'],['2']);self.assertTrue(compiler.verify(c))
    def test_quadratic_cross_coefficient(self):
        p=problem('a*t+b*t*t','a*a+3*a*b+2*b*b',['a','b'])
        self.assertEqual(p['penalty']['hessian'],[['2','3'],['3','4']])
    def test_decimal_rejected(self):
        with self.assertRaises(ValueError):compiler.compile_goal(goal('0.1*a*t'))
    def test_code_injection_rejected(self):
        with self.assertRaises(ValueError):compiler.compile_goal(goal("__import__('os').system('true')"))
    def test_nonlinear_potential_rejected(self):
        with self.assertRaises(ValueError):compiler.compile_goal(goal('a*a*t'))
    def test_parameter_denominator_rejected(self):
        with self.assertRaises(ValueError):compiler.compile_goal(goal('t/a'))
    def test_higher_index_rejected(self):
        raw=goal();raw['eigenvalue_index']=2
        with self.assertRaises(ValueError):compiler.compile_goal(raw)
    def test_phi_dependence_rejected(self):
        with self.assertRaises(ValueError):compiler.compile_goal(goal('a*t+phi'))
    def test_wrong_geometry_rejected(self):
        raw=goal();raw['geometry']='arbitrary_surface'
        with self.assertRaises(ValueError):compiler.compile_goal(raw)
    def test_bound_tampering_rejected(self):
        c=compiler.compile_goal(goal());c['threshold']='1';self.assertFalse(compiler.verify(c))
    def test_zero_parameter_fixed_problem(self):
        p=problem('t*t','3',[],[]);self.assertEqual(p['directions'],[])


class V24SphereReduction(unittest.TestCase):
    def test_fixed_solver_cannot_skip_scope_rule(self):
        with self.assertRaises(ValueError):decision.search(problem('t*t','0',[],[]),0)
    def test_boolean_spectral_index_rejected(self):
        p=problem();p['eigenvalue_index']=True
        with self.assertRaises(ValueError):compiler.normalized(p)
    def test_scope_reduction_is_explicit(self):
        r=reduction.sphere_ground(problem());self.assertTrue(reduction.verify(r));self.assertEqual(r['reduced']['function_space'],'axisymmetric')
    def test_other_coefficients_unchanged(self):
        p=problem('a*t+2*t*t','a*a/4');r=reduction.sphere_ground(p)
        for k in ('q0','directions','penalty','domain'):self.assertEqual(p[k],r['reduced'][k])
    def test_forged_scope_reduction_rejected(self):
        r=reduction.sphere_ground(problem());r['reduced']['q0']=['2'];self.assertFalse(reduction.verify(r))
    def test_axisymmetric_precondition(self):
        with self.assertRaises(ValueError):reduction.sphere_ground(problem(scope='axisymmetric'))
    def test_lifted_trial_energy(self):
        p=problem();r=reduction.sphere_ground(p);self.assertEqual(energy(p,[F(2)],[F(1),F(-1)]),energy(r['reduced'],[F(2)],[F(1),F(-1)]))


class V25PassiveElimination(unittest.TestCase):
    def test_coupled_real_elimination_identity(self):
        p=problem('a*t+2*b+3','a*a+2*a*b+3*b*b',['a','b']);r=reduction.absorb_and_eliminate(p)
        self.assertTrue(reduction.verify(r));self.assertEqual(len(r['reduced']['directions']),1)
        for z in (F(-2),F(1,3),F(3)):
            x=reduction.lift(r,[z]);self.assertEqual(energy(p,x,[F(1),F(1,3)]),energy(r['reduced'],[z],[F(1),F(1,3)]))
    def test_conditional_minimum(self):
        p=problem('a*t+b','a*a+a*b+b*b',['a','b']);r=reduction.absorb_and_eliminate(p);x=reduction.lift(r,[F(1)])
        for delta in (F(-1),F(1)):
            y=x[:];y[1]+=delta;self.assertGreaterEqual(energy(p,y,[F(1)]),energy(p,x,[F(1)]))
    def test_box_indefinite_passive_coordinate(self):
        p=problem('a*t','a*a-b*b',['a','b'],[[-1,1],[-2,2]]);r=reduction.absorb_and_eliminate(p);self.assertEqual(r['reduced']['penalty']['constant'],'-4');self.assertTrue(reduction.verify(r))
    def test_box_coupling_not_silently_eliminated(self):
        p=problem('a*t','a*a+a*b+b*b',['a','b'],[[-1,1]]*2)
        with self.assertRaises(ValueError):reduction.absorb_and_eliminate(p)
    def test_singular_passive_real_block_rejected(self):
        p=problem('a*t','a*a',['a','b'])
        with self.assertRaises(ValueError):reduction.absorb_and_eliminate(p)
    def test_all_parameters_eliminated(self):
        p=problem('t*t+a','a*a');r=reduction.absorb_and_eliminate(p);self.assertEqual(r['reduced']['directions'],[]);self.assertEqual(r['reduced']['penalty']['constant'],'-1/4')
    def test_lift_tampering_rejected(self):
        r=reduction.absorb_and_eliminate(problem('a*t+b','a*a+b*b',['a','b']));r['lift_offset'][1]='5';self.assertFalse(reduction.verify(r))


class V26OrbitReduction(unittest.TestCase):
    def test_half_box(self):
        r=reduction.reflect(problem(box=[[-4,4]]));self.assertEqual(r['reduced']['domain']['box'],[['0','4']]);self.assertTrue(reduction.verify(r))
    def test_simultaneous_reflection_of_multiple_parameters(self):
        p=problem('a*t+b*t**3','a*a+a*b+b*b',['a','b'],[[-2,2]]*2);r=reduction.reflect(p)
        self.assertEqual(r['signs'],[-1,-1]);x=[F(1,3),F(1,2)];trial=[F(1),F(1,4)]
        self.assertEqual(energy(p,x,trial),energy(p,[-z for z in x],[trial[0],-trial[1]]))
    def test_nonsymmetric_domain_rejected(self):
        with self.assertRaises(ValueError):reduction.reflect(problem(box=[[-1,2]]))
    def test_linear_penalty_breaks_symmetry(self):
        with self.assertRaises(ValueError):reduction.reflect(problem(penalty='a*a+a',box=[[-2,2]]))
    def test_mixed_parity_rejected(self):
        with self.assertRaises(ValueError):reduction.reflect(problem(potential='a*(t+t*t)',box=[[-2,2]]))
    def test_even_odd_coupling_rejected(self):
        with self.assertRaises(ValueError):reduction.reflect(problem('a*t+b*t*t','a*a+a*b+b*b',['a','b'],[[-2,2]]*2))


class V27RankCompression(unittest.TestCase):
    def model(self):return problem('(a+2*b-c+3*d)*t','(a*a+b*b+c*c+d*d)/2',['a','b','c','d'])
    def test_four_to_one_exact_hessian(self):
        r=reduction.compress(self.model());self.assertEqual(r['reduced']['penalty']['hessian'],[['1/15']]);self.assertTrue(reduction.verify(r))
    def test_lift_preserves_energy(self):
        p=self.model();r=reduction.compress(p)
        for z in (F(-3),F(1,7),F(2)):
            self.assertEqual(energy(p,reduction.lift(r,[z]),[F(1),F(-1,5)]),energy(r['reduced'],[z],[F(1),F(-1,5)]))
    def test_fiber_perturbation_not_better(self):
        p=self.model();r=reduction.compress(p);x=reduction.lift(r,[F(2)]);y=x[:];y[0]+=2;y[1]-=1
        self.assertGreater(energy(p,y,[F(1)]),energy(p,x,[F(1)]))
    def test_two_independent_directions(self):
        p=problem('(a+b)*t+(b+c)*t*t','a*a+b*b+c*c',['a','b','c']);r=reduction.compress(p);self.assertEqual(len(r['reduced']['directions']),2);self.assertTrue(reduction.verify(r))
    def test_full_rank_rejected(self):
        with self.assertRaises(ValueError):reduction.compress(problem('a*t+b*t*t','a*a+b*b',['a','b']))
    def test_box_projection_cannot_reuse_real_fiber_formula(self):
        with self.assertRaises(ValueError):reduction.compress(problem('a*t+b*t','a*a+b*b',['a','b'],[[-1,1]]*2))
    def test_indefinite_penalty_rejected(self):
        with self.assertRaises(ValueError):reduction.compress(problem('(a+b)*t','a*a-b*b',['a','b']))
    def test_bad_map_rejected(self):
        r=reduction.compress(self.model());r['coefficient_map'][0][0]='2';self.assertFalse(reduction.verify(r))


class V28Decisions(unittest.TestCase):
    def test_prove_without_optimizing_gap(self):
        p=problem('a*t','a*a/4',box=[[-1,1]],scope='axisymmetric');r=decision.search(p,F(-2));self.assertEqual(r['certificate']['status'],'proved');self.assertEqual(r['statistics']['points'],0)
    def test_counterexample_exact(self):
        p=problem('a*t','a*a/10',box=[[-4,4]],scope='axisymmetric');r=decision.search(p,F(-1,10));self.assertEqual(r['certificate']['status'],'refuted');self.assertTrue(decision.verify(r['certificate']))
    def test_boundary_budget_unresolved(self):
        p=problem(box=[[-2,2]],scope='axisymmetric');r=decision.search(p,0,max_leaves=1);self.assertEqual(r['certificate']['status'],'unresolved');self.assertTrue(decision.verify(r['certificate']))
    def test_change_query_invalidates_decision(self):
        p=problem(box=[[-2,2]],scope='axisymmetric');r=decision.search(p,F(-10))['certificate'];r['threshold']='10';self.assertFalse(decision.verify(r))
    def test_fixed_operator(self):
        p=problem('t*t','3',[],[],scope='axisymmetric');r=decision.search(p,2);self.assertEqual(r['certificate']['status'],'proved');self.assertTrue(decision.verify(r['certificate']))
    def test_quadratic_singular_real_min(self):
        low,x=algebra.quadratic_min(1,[-2,-2],[[2,2],[2,2]]);self.assertEqual(low,0)
    def test_quadratic_unbounded_linear_direction(self):
        with self.assertRaises(ValueError):algebra.quadratic_min(0,[0,1],[[1,0],[0,0]])


class V29PoissonTemplates(unittest.TestCase):
    def test_uniform_linear_family(self):
        c=family.synthesize([0,1]);self.assertEqual(c['scale'],'1/4');self.assertEqual(c['poisson_solutions'],[['0','-1/2']]);self.assertTrue(family.verify(c))
    def test_centered_quadratic_family(self):
        c=family.synthesize([F(-1,3),0,1]);self.assertEqual(c['scale'],'1/36');self.assertEqual(c['means'],['0']);self.assertTrue(family.verify(c))
    def test_poisson_identity_degrees_one_to_six(self):
        for n in range(1,7):
            p=[F(0)]*n+[F(1)];r,mean=family.poisson(p);self.assertEqual(family.laplacian(r),e.add(e.scale(p,-1),[mean]))
    def test_local_energy_identity(self):
        p=[F(1),F(2),F(-3),F(1)];r,mean=family.poisson(p)
        for amp in (F(-10),F(1,3),F(7)):
            got=e.local_energy(e.scale(p,amp),e.scale(r,amp));expected=e.add([amp*mean],e.scale(e.mul([1,0,-1],e.mul(e.derivative(r),e.derivative(r))),-amp*amp));self.assertEqual(got,expected)
    def test_wrong_poisson_solution_rejected(self):
        c=family.synthesize([0,1]);c['poisson_solutions'][0][1]='-1';self.assertFalse(family.verify(c))
    def test_changed_direction_rejected(self):
        c=family.synthesize([0,1]);c['directions'][0][1]='2';self.assertFalse(family.verify(c))
    def test_constant_family(self):
        c=family.synthesize([3]);self.assertEqual(c['scale'],'0');self.assertEqual(c['means'],['3'])


class V30MatrixTemplates(unittest.TestCase):
    def test_combined_bound_beats_separate_cauchy_bound(self):
        c=matrix.synthesize([[0,1],[F(-1,3),0,1]]);self.assertEqual(c['scale'],'1/4');self.assertLess(F(c['scale']),F(1,2));self.assertTrue(matrix.verify(c))
    def test_nonidentity_weight(self):
        c=matrix.synthesize([[0,1],[0,0,1]],[[2,F(1,2)],[F(1,2),1]]);self.assertTrue(matrix.verify(c))
    def test_spd_required(self):
        for w in ([[1,2],[2,1]],[[1,0],[0,0]],[[1,1],[0,1]]):
            with self.assertRaises(ValueError):matrix.synthesize([[0,1],[0,0,1]],w)
    def test_degree_twelve_envelope(self):
        c=matrix.synthesize([[0,1],[0,0,0,0,0,0,1]]);self.assertLessEqual(len(c['range_proof']['polynomial']),13);self.assertTrue(matrix.verify(c))
    def test_false_scale_rejected(self):
        c=matrix.synthesize([[0,1],[0,0,1]]);c['scale']='0';self.assertFalse(matrix.verify(c))
    def test_subdivision_cover_rejected(self):
        c=matrix.synthesize([[0,1],[0,0,0,1]]);proof=c['range_proof']
        # Use a degree-six polynomial that requires a real partition.
        proof=family.maximum([0,1,0,0,0,0,1],max_leaves=4);proof['intervals'].pop();self.assertFalse(family.verify_range(proof))
    def test_open_scale_budget_still_valid(self):
        c=matrix.synthesize([[0,0,0,0,0,0,1]],tolerance=F(1,10**12),max_leaves=1);self.assertEqual(c['range_proof']['status'],'scale_gap_open');self.assertTrue(matrix.verify(c))
    def test_pointwise_rank_one_domination(self):
        c=matrix.synthesize([[0,1],[0,0,1]]);g=[list(map(F,row)) for row in c['quadratic_bound']];rs=[list(map(F,r)) for r in c['poisson_solutions']]
        for amp in product((F(-2),F(1,3)),repeat=2):
            upper=sum(amp[i]*g[i][j]*amp[j] for i in range(2) for j in range(2))
            for t in (F(-1),F(-1,3),F(0),F(1,2),F(1)):
                v=sum(z*e.evaluate(e.derivative(r),t) for z,r in zip(amp,rs));self.assertLessEqual((1-t*t)*v*v,upper)


class V31Transport(unittest.TestCase):
    def test_affine_shift_and_positive_remainder(self):
        p=problem('(2*a+1)*t+t**4','(2*a+1)**2/4');c=transport.apply(p,family.synthesize([0,1]));self.assertEqual(c['lower_bound'],'0');self.assertTrue(transport.verify(c))
    def test_constant_parameter_terms(self):
        p=problem('(2*a+1)*t+3*a+4','(2*a+1)**2/4-3*a-4');c=transport.apply(p,family.synthesize([0,1]));self.assertEqual(c['lower_bound'],'0')
    def test_multivariate_affine_reuse(self):
        p=problem('(a+2*b)*t+(b-a)*(t*t-1/3)','((a+2*b)**2+(b-a)**2)/4',['a','b']);c=transport.apply(p,matrix.synthesize([[0,1],[F(-1,3),0,1]]));self.assertEqual(c['lower_bound'],'0');self.assertTrue(transport.verify(c))
    def test_nonmatching_direction_rejected(self):
        with self.assertRaises(ValueError):transport.apply(problem('a*t**3'),family.synthesize([0,1]))
    def test_negative_remainder_not_ignored(self):
        p=problem('a*t-t**4');c=transport.apply(p,family.synthesize([0,1]));self.assertEqual(c['lower_bound'],'-1')
    def test_too_weak_all_real_lemma_not_a_proof(self):
        with self.assertRaises(ValueError):transport.apply(problem(penalty='a*a/5'),family.synthesize([0,1]))
    def test_forged_template_identifier_rejected(self):
        t=family.synthesize([0,1]);t['template_id']='fake'
        with self.assertRaises(ValueError):transport.apply(problem(),t)
    def test_wrong_parameter_map_rejected(self):
        c=transport.apply(problem(),family.synthesize([0,1]));c['amplitude_matrix'][0][0]='2';self.assertFalse(transport.verify(c))


class V32ProofSearch(unittest.TestCase):
    def test_full_sphere_uniform_goal_zero_spectral_calls(self):
        r=planner.solve(goal());self.assertEqual(r['certificate']['status'],'proved');self.assertEqual(r['statistics']['points'],0);self.assertTrue(planner.verify(r['certificate']))
    def test_counterexample_lifted_to_four_parameters(self):
        raw=goal('(a+b+c+d)*t','(a*a+b*b+c*c+d*d)/2',['a','b','c','d']);r=planner.solve(raw,max_leaves=8)
        self.assertEqual(r['certificate']['status'],'refuted');self.assertEqual(len(r['certificate']['candidate_parameters']),4);self.assertTrue(planner.verify(r['certificate']))
    def test_too_weak_template_remains_unresolved(self):
        r=planner.solve(goal(penalty='a*a/5',box=[[-2,2]]),max_leaves=4);self.assertEqual(r['certificate']['status'],'unresolved');self.assertTrue(planner.verify(r['certificate']))
    def test_missing_scope_edge_rejected(self):
        r=planner.solve(goal(penalty='a*a/10',box=[[-2,2]],threshold='-1/100'),synthesize=False)['certificate'];r['reductions']=[x for x in r['reductions'] if x['rule']!='sphere_ground_v24'];self.assertFalse(planner.verify(r))
    def test_corrupt_witness_lift_rejected(self):
        raw=goal('(a+b)*t','(a*a+b*b)/2',['a','b'],threshold='1');c=planner.solve(raw,synthesize=False,max_leaves=4)['certificate'];c['candidate_parameters'][0]='100';self.assertFalse(planner.verify(c))
    def test_library_reuse_without_synthesis(self):
        library=[family.synthesize([0,1])];r=planner.solve(goal('(3*a-2)*t+t**6','(3*a-2)**2/4'),library=library,synthesize=False)
        self.assertEqual(r['certificate']['status'],'proved');self.assertEqual(r['statistics']['template_syntheses'],0)
    def test_no_false_model_comparison(self):
        r=planner.solve(goal());self.assertEqual(r['statistics']['model_calls'],0)
    def test_bad_query_binding_rejected(self):
        c=planner.solve(goal())['certificate'];c['compiled_goal']['source']['potential']='2*a*t';self.assertFalse(planner.verify(c))
    def test_fixed_passive_pipeline(self):
        r=planner.solve(goal('t*t+a','a*a',threshold='-1'),synthesize=False)
        self.assertEqual(r['certificate']['status'],'proved');self.assertTrue(planner.verify(r['certificate']))

if __name__=='__main__':unittest.main()
