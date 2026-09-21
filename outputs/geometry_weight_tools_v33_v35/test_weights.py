import copy
import unittest
from unittest.mock import patch
from fractions import Fraction as F
import algebra as a
import dual_v33 as dual
import optimize_v34 as optimizer
import family_v29 as family
import matrix_v30
import planner_v35 as planner
import planner_v32 as baseline
import compiler_v23 as compiler
import run


DIRS=[[0,1],[F(-1,3),0,1]]


def goal(penalty='a*a/4+b*b/9'):
    return {'geometry':'unit_sphere','function_space':'full_sphere','eigenvalue_index':1,
            'parameters':['a','b'],'potential':'a*t+b*(t*t-1/3)','penalty':penalty,
            'domain':{'kind':'all_real'},'threshold':'0'}


class V33DualTests(unittest.TestCase):
    def setUp(self):
        self.lower=dual.lower_certificate([[0,1]],[[1]],[{'t':'0','matrix':[['1']]}])
        self.template=family.synthesize([0,1])

    def test_exact_scalar_global_optimum(self):
        cert=dual.certificate(self.template,self.lower,F(1,100000))
        self.assertEqual(cert['gap'],'0');self.assertEqual(cert['lower_bound'],'1/4');self.assertTrue(dual.verify(cert))

    def test_nonidentity_cost(self):
        lower=dual.lower_certificate([[0,1]],[[2]],[{'t':'0','matrix':[['2']]}])
        self.assertEqual(dual.certificate(self.template,lower,F(1,100))['upper_bound'],'1/2')

    def test_negative_atom_rejected(self):
        with self.assertRaises(ValueError):dual.lower_certificate([[0,1]],[[1]],[{'t':'0','matrix':[['-1']]}])

    def test_indefinite_atom_with_positive_diagonal_rejected(self):
        with self.assertRaises(ValueError):dual.lower_certificate(DIRS,[[3,0],[0,3]],[{'t':'0','matrix':[[1,2],[2,1]]}])

    def test_aggregate_cost_violation_rejected(self):
        with self.assertRaises(ValueError):dual.lower_certificate([[0,1]],[[1]],[{'t':'0','matrix':[['3/5']]},{'t':'1/2','matrix':[['3/5']]}])

    def test_outside_coordinate_rejected(self):
        with self.assertRaises(ValueError):dual.lower_certificate([[0,1]],[[1]],[{'t':'2','matrix':[['1']]}])

    def test_zero_slack_not_required(self):
        cert=dual.lower_certificate([[0,1]],[[1]],[{'t':'0','matrix':[['1/2']]}])
        self.assertEqual(cert['cost_slack'],[['1/2']]);self.assertEqual(cert['lower_bound'],'1/8')

    def test_singular_atom_allowed(self):
        cert=dual.lower_certificate(DIRS,[[1,0],[0,1]],[{'t':'0','matrix':[[1,0],[0,0]]}])
        self.assertEqual(cert['lower_bound'],'1/4');self.assertTrue(dual.verify_lower(cert))

    def test_repair_oversized_atom(self):
        cert=dual.repair_dual([[0,1]],[[1]],[{'t':'0','matrix':[[3]]}])
        self.assertTrue(dual.verify_lower(cert));self.assertLessEqual(F(cert['lower_bound']),F(1,4))
        self.assertGreater(F(cert['lower_bound']),F(2499999,10000000))

    def test_noncommuting_matrix_repair(self):
        cost=[[2,F(1,2)],[F(1,2),1]]
        cert=dual.repair_dual(DIRS,cost,[{'t':'1/3','matrix':[[5,2],[2,3]]},{'t':'-1/2','matrix':[[2,-1],[-1,2]]}])
        self.assertTrue(dual.verify_lower(cert));self.assertGreater(F(cert['lower_bound']),0)

    def test_lower_value_tamper(self):
        cert=copy.deepcopy(self.lower);cert['lower_bound']='100';self.assertFalse(dual.verify_lower(cert))

    def test_slack_tamper(self):
        cert=copy.deepcopy(self.lower);cert['cost_slack']=[['1']];self.assertFalse(dual.verify_lower(cert))

    def test_sample_tamper(self):
        cert=copy.deepcopy(self.lower);cert['atoms'][0]['t']='1';self.assertFalse(dual.verify_lower(cert))

    def test_dual_direction_binding(self):
        with self.assertRaises(ValueError):dual.certificate(family.synthesize([0,2]),self.lower,F(1,100))

    def test_gap_and_status_tamper(self):
        cert=dual.certificate(self.template,self.lower,F(1,100));cert['status']='global_gap_open';self.assertFalse(dual.verify(cert))

    def test_open_gap_stays_open(self):
        lower=dual.lower_certificate([[0,1]],[[1]],[{'t':'1','matrix':[['1']]}])
        cert=dual.certificate(self.template,lower,F(1,100));self.assertEqual(cert['status'],'global_gap_open');self.assertTrue(dual.verify(cert))

    def test_nonsymmetric_cost_rejected(self):
        with self.assertRaises(ValueError):dual.setup(DIRS,[[2,1],[0,1]])

    def test_semidefinite_cost_rejected(self):
        with self.assertRaises(ValueError):dual.setup(DIRS,[[1,0],[0,0]])

    def test_dependent_directions_rejected(self):
        with self.assertRaises(ValueError):dual.setup([[0,1],[2,2]])

    def test_constant_direction_rejected(self):
        with self.assertRaises(ValueError):dual.setup([[1]])

    def test_degree_limit(self):
        with self.assertRaises(ValueError):dual.setup([[0]*7+[1]])

    def test_empty_atoms_rejected(self):
        with self.assertRaises(ValueError):dual.lower_certificate([[0,1]],[[1]],[])

    def test_matrix_dimension_rejected(self):
        with self.assertRaises(ValueError):dual.lower_certificate(DIRS,[[1,0],[0,1]],[{'t':'0','matrix':[[1]]}])

    def test_weak_duality_trace_identity(self):
        dirs,v,c=dual.setup(DIRS,[[2,F(1,2)],[F(1,2),1]])
        cert=dual.repair_dual(dirs,c,[{'t':'-1/3','matrix':[[1,F(1,2)],[F(1,2),2]]},{'t':'1/2','matrix':[[2,0],[0,1]]}])
        g=[[F(10),F(1)],[F(1),F(10)]];slack=[[F(z) for z in row] for row in cert['cost_slack']]
        rhs=dual.inner(slack,g)
        for atom in cert['atoms']:
            y=[[F(z) for z in row] for row in atom['matrix']];aa=dual.at(v,atom['t'])
            rhs+=dual.inner(y,[[g[i][j]-aa[i][j] for j in range(2)] for i in range(2)])
        self.assertEqual(rhs,dual.inner(c,g)-F(cert['lower_bound']));self.assertGreaterEqual(rhs,0)


class V34OptimizerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=optimizer.solve(DIRS)

    def test_dual_certified_gap(self):
        c=self.result['certificate'];self.assertTrue(dual.verify(c));self.assertEqual(c['status'],'global_gap_closed');self.assertLessEqual(F(c['gap']),F(1,100000))

    def test_improves_identity_weight(self):
        self.assertLess(F(self.result['certificate']['upper_bound']),F(34,100))
        original=matrix_v30.synthesize(DIRS);self.assertEqual(sum(F(original['quadratic_bound'][i][i]) for i in range(2)),F(1,2))

    def test_full_interval_proof_retained(self):
        self.assertTrue(family.verify(self.result['certificate']['template']))

    def test_offdiagonal_geometry(self):
        result=optimizer.solve([[0,1,1],[0,2,-1]])
        self.assertTrue(dual.verify(result['certificate']));self.assertEqual(result['certificate']['status'],'global_gap_closed')
        self.assertNotEqual(F(result['certificate']['template']['quadratic_bound'][0][1]),0)

    def test_weighted_objective(self):
        result=optimizer.solve(DIRS,[[2,F(1,2)],[F(1,2),1]])
        self.assertTrue(dual.verify(result['certificate']));self.assertEqual(result['certificate']['status'],'global_gap_closed')

    def test_three_directions(self):
        result=optimizer.solve([[0,1],[0,0,1],[0,0,0,1]])
        self.assertTrue(dual.verify(result['certificate']));self.assertEqual(result['certificate']['status'],'global_gap_closed')

    def test_four_directions(self):
        result=optimizer.solve([[0,1],[0,0,1],[0,0,0,1],[0,0,0,0,1]])
        self.assertTrue(dual.verify(result['certificate']));self.assertEqual(result['certificate']['status'],'global_gap_closed')

    def test_exchange_budget_is_honest(self):
        result=optimizer.solve(DIRS,max_rounds=1)
        self.assertEqual(result['certificate']['status'],'global_gap_open');self.assertTrue(dual.verify(result['certificate']))

    def test_proposal_failure_returns_verified_open_gap(self):
        with patch('float_sdp.solve',side_effect=ValueError('injected numerical failure')):result=optimizer.solve(DIRS)
        self.assertEqual(result['certificate']['status'],'global_gap_open');self.assertTrue(dual.verify(result['certificate']))
        self.assertEqual(result['trace'][0]['outcome'],'proposal_failed')

    def test_no_model_calls(self):self.assertEqual(self.result['statistics']['model_calls'],0)

    def test_invalid_round_budget(self):
        for v in (0,25,True):
            with self.assertRaises(ValueError):optimizer.solve(DIRS,max_rounds=v)

    def test_invalid_accuracy_budget(self):
        for value in (0,F(1,10**9),2):
            with self.assertRaises(ValueError):optimizer.solve(DIRS,tolerance=value)

    def test_bad_range_budget(self):
        with self.assertRaises(ValueError):optimizer.solve(DIRS,max_range_leaves=0)

    def test_optimized_lemma_can_be_reused(self):
        loaded=run.load_library(self.result)
        self.assertEqual(loaded,[self.result['certificate']['template']]);self.assertTrue(family.verify(loaded[0]))

    def test_tampered_optimization_not_admitted_as_library(self):
        raw=copy.deepcopy(self.result);raw['certificate']['lower_bound']='100'
        with self.assertRaises(ValueError):run.load_library(raw)


class V35PlannerTests(unittest.TestCase):
    def test_anisotropic_goal_proved_without_spectral_search(self):
        result=planner.solve(goal());self.assertEqual(result['certificate']['status'],'proved');self.assertEqual(result['strategy'],'penalty_weighted_lemma')
        self.assertEqual(result['statistics']['spectral_searches'],0);self.assertTrue(baseline.verify(result['certificate']))

    def test_original_full_sphere_scope(self):
        result=planner.solve(goal());self.assertEqual(result['certificate']['original_function_space'],'full_sphere')

    def test_nonorthogonal_parameter_map(self):
        raw=goal();raw['potential']='(a+2*b)*t+(b-a)*(t*t-1/3)';raw['penalty']='(a+2*b)**2/4+(b-a)**2/9'
        result=planner.solve(raw);self.assertEqual(result['certificate']['status'],'proved');self.assertTrue(baseline.verify(result['certificate']))

    def test_positive_remainder(self):
        raw=goal();raw['potential']+=' + t**6';result=planner.solve(raw);self.assertEqual(result['certificate']['status'],'proved')

    def test_failed_ansatz_not_false_theorem(self):
        raw=goal();raw.update({'parameters':['a'],'potential':'a*t','penalty':'a*a/5','domain':{'kind':'box','box':[[-2,2]]}})
        result=planner.solve(raw,max_leaves=8);self.assertEqual(result['certificate']['status'],'unresolved')
        self.assertEqual(result['diagnostics'][0]['status'],'poisson_envelope_incompatible');self.assertTrue(planner.verify_budget(result['diagnostics'][0]))
        self.assertEqual(F(result['diagnostics'][0]['scale_lower']),F(5,4))

    def test_linear_term_can_prevent_proof(self):
        raw=goal('a*a/4+b*b/9+a/100');result=planner.solve(raw,max_leaves=4)
        self.assertNotEqual(result['strategy'],'penalty_weighted_lemma');self.assertTrue(baseline.verify(result['certificate']))

    def test_bad_library_rejected_before_fast_return(self):
        with self.assertRaises(ValueError):planner.solve(goal(),library=[{'fake':True}])

    def test_disabled_weights_fall_back(self):
        result=planner.solve(goal(),use_penalty_weights=False,max_leaves=4)
        self.assertEqual(result['strategy'],'v32_fallback');self.assertEqual(result['statistics']['weighted_attempts'],0)

    def test_disabled_synthesis_includes_weights(self):
        result=planner.solve(goal(),synthesize=False,max_leaves=4)
        self.assertEqual(result['statistics']['weighted_attempts'],0)

    def test_budget_claim_tamper(self):
        result=planner.solve(goal());cert=copy.deepcopy(result['diagnostics'][0]);cert['status']='poisson_envelope_incompatible'
        self.assertFalse(planner.verify_budget(cert))

    def test_budget_weight_binding(self):
        p=compiler.compile_goal(goal())['problem'];template=matrix_v30.synthesize(DIRS)
        with self.assertRaises(ValueError):planner.budget_certificate(p,template)

    def test_budget_direction_binding(self):
        result=planner.solve(goal());cert=copy.deepcopy(result['diagnostics'][0]);cert['problem']['directions'][0]=['0','2']
        self.assertFalse(planner.verify_budget(cert))

    def test_high_eigenvalue_rejected(self):
        raw=goal();raw['eigenvalue_index']=2
        with self.assertRaises(ValueError):planner.solve(raw)

    def test_counterexample_is_lifted_and_verified(self):
        raw=goal();raw.update({'parameters':['a'],'potential':'a*t','penalty':'a*a/10','domain':{'kind':'box','box':[[-4,4]]},'threshold':'-1/10'})
        result=planner.solve(raw,max_leaves=16);self.assertEqual(result['certificate']['status'],'refuted');self.assertTrue(baseline.verify(result['certificate']))

    def test_no_parameter_goal(self):
        raw=goal();raw.update({'parameters':[],'potential':'2*t*t','penalty':'0','domain':{'kind':'box','box':[]}})
        result=planner.solve(raw);self.assertEqual(result['certificate']['status'],'proved');self.assertTrue(baseline.verify(result['certificate']))


if __name__=='__main__':unittest.main()
