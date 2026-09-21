"""V90: explicit certificate registry; no verifier supplied by untrusted JSON.

The analytic log-Sobolev dependency is deliberately visible in decision scope.
Older certificate formats continue to use their unchanged verifiers.
"""
from fractions import Fraction as F
import legacy_cli
import exact_extrema as extrema
import positive_search as positive
import uniform_refine as uniform


def backend(name):
    if name in (None, positive.DEFAULT_RANGE_BACKEND.name):
        return positive.DEFAULT_RANGE_BACKEND
    if name == 'exact_sturm_minimum_v1':
        from precision_bridge import ExactRangeBackend
        return ExactRangeBackend()
    raise ValueError('Unknown trusted range backend')


def verify(c):
    if not isinstance(c, dict):
        return False
    try:
        method = c.get('method', c.get('format'))
        if method == 'exact_sturm_extrema_v1':
            return extrema.verify(c)
        if method == 'precise_poisson_matrix_v88':
            from precision_bridge import verify_family
            return verify_family(c)
        if method == 'precise_global_matrix_v88':
            from precision_bridge import verify_matrix
            return verify_matrix(c)
        if method == 'original_spectral_goal_v89':
            from original_goal import verify as verify_goal
            return verify_goal(c)
        if method in (uniform.METHOD, uniform.LOG_SOBOLEV_METHOD):
            return uniform.verify(c)
        if method == 'positive_exponential_lower_v1':
            return positive.verify_lower(c, backend(c['range_backend']))
        if method == 'positive_exponential_ansatz_upper_v1':
            return positive.verify_ansatz_upper(c)
        if method == 'positive_search_rayleigh_upper_v1':
            return positive.verify_rayleigh(c)
        if method == 'positive_exponential_even_projection_v1':
            return positive.verify_parity(c)
        if method == 'positive_exponential_global_optimization_v1':
            return positive.verify(c, backend(c['lower_certificate']['range_backend']))
        return legacy_cli.verify(c)
    except (ValueError, TypeError, KeyError, IndexError, ZeroDivisionError, ArithmeticError, AttributeError, RecursionError):
        return False


# The V35 lemma transport accepts only its original Poisson format.
load_library = legacy_cli.load_library


def summary_fields(c):
    """Return verified NEW-format semantics; protocol attaches content references."""
    if not verify(c):
        raise ValueError('Invalid certificate')
    method = c.get('method', c.get('format'))
    if method == 'exact_sturm_extrema_v1':
        return dict(kind='polynomial_maximum', status=c['status'],
                    task={'research_kind':'extrema','polynomial':c['polynomial'],'interval':c['interval']},
                    scope={'objective':'global polynomial maximum','interval':c['interval']},
                    target=c['tolerance'], bounds={'lower':c['maximum_lower'],'upper':c['maximum_upper'],'gap':c['error_bound']})
    if method in (uniform.METHOD, uniform.LOG_SOBOLEV_METHOD):
        scope={'geometry':'unit_sphere','space':'full_sphere','eigenvalue_index':1,
               'amplitudes':'all_real' if method==uniform.LOG_SOBOLEV_METHOD or c['all_real'] else c['amplitude_box'],
               'claim':c['claim']}
        bounds={}
        if method==uniform.LOG_SOBOLEV_METHOD:
            scope['analytic_dependency']='classical sharp S2 log-Sobolev theorem (not proved by this program)'
            bounds={'sharp_coefficient_lower':c['sharp_uniform_quadratic_coefficient'],
                    'sharp_coefficient_upper':c['sharp_uniform_quadratic_coefficient']}
            if c['counterexample']:
                bounds['counterexample_objective_upper']=c['counterexample']['spectral_objective_upper']
        elif c['status']=='proved':
            bounds={'coefficient_upper':c['coefficient']}
        return dict(kind='uniform_spectral_inequality',status=c['status'],
                    task={'research_kind':'uniform','direction':c['direction'],'coefficient':c['coefficient'],
                          'amplitudes':scope['amplitudes'],'rule':'log_sobolev' if method==uniform.LOG_SOBOLEV_METHOD else 'barta'},
                    scope=scope,target=c['coefficient'],bounds=bounds)
    if method=='positive_exponential_global_optimization_v1':
        lower=c['lower_certificate']
        bounds={'ansatz_lower':c['ansatz_lower'],'ansatz_upper':c['ansatz_upper'],'ansatz_gap':c['ansatz_gap'],
                'spectral_lower':lower['spectral_lower'],'spectral_upper':c.get('spectral_upper')}
        if 'spectral_gap' in c:bounds['spectral_gap']=c['spectral_gap']
        return dict(kind='positive_function_optimization',status=c['status'],
                    task={'research_kind':'positive','potential':lower['q_coefficients'],'basis':lower['basis']},
                    scope={'geometry':'unit_sphere','space':'full_sphere','eigenvalue_index':1,
                           'objective':'best Barta lower in specified log-polynomial basis',
                           'ansatz_upper_is_spectral_upper':False},
                    target=c['tolerance'],bounds=bounds)
    if method=='original_spectral_goal_v89':
        source=c['compiled_goal']['source']
        scope={'geometry':source['geometry'],'space':source['function_space'],'eigenvalue_index':source['eigenvalue_index'],
               'domain':source['domain'],'objective':'inf_parameters(lambda_1(-Delta+q)+penalty)'}
        if c['rule']=='log_sobolev':
            scope['analytic_dependency']='classical sharp S2 log-Sobolev theorem (not proved by this program)'
        return dict(kind='ground_goal',status=c['status'],task=source,scope=scope,
                    target=c['compiled_goal']['threshold'],bounds={'lower':c['lower_bound'],'upper':c['upper_bound']})
    if method=='precise_global_matrix_v88':
        # The bridge's old seed binds the same cost/directions and dual objective.
        return dict(kind='matrix_envelope',status=c['status'],
                    task={'directions':c['directions'],'cost':c['cost']},
                    scope={'geometry':'unit_sphere','amplitudes':'all_real','interval':['-1','1'],
                           'ansatz':'poisson_exponential','objective':'trace(C G)'},
                    target=c['tolerance'],bounds={'lower':c['lower_bound'],'upper':c['upper_bound'],'gap':c['gap']})
    if method=='precise_poisson_matrix_v88':
        return dict(kind='family_lemma',status=c['range_proof']['status'],
                    task={'directions':c['directions'],'weight':c['weight']},
                    scope={'geometry':'unit_sphere','space':'full_sphere','amplitudes':'all_real','ansatz':'poisson_exponential',
                           'objective':'maximum Poisson gradient quotient'},
                    target=c['range_proof']['tolerance'],bounds={'lower':c['range_proof']['maximum_lower'],'upper':c['range_proof']['maximum_upper']})
    raise ValueError('Component certificate has no decision summary; use exact inspect/verify')


def validate_task(raw):
    """Canonical new task schemas; no implicit domain or basis default in refs."""
    if not isinstance(raw,dict):raise ValueError('Task object required')
    kind=raw.get('research_kind')
    if kind=='extrema' and set(raw)=={'research_kind','polynomial','interval'}:
        # Zero refinement still verifies a finite, canonical polynomial task.
        c=extrema.maximize(raw['polynomial'],raw['interval'],max_nodes=0)
        if c['polynomial']!=raw['polynomial'] or c['interval']!=raw['interval']:raise ValueError('Noncanonical extrema task')
    elif kind=='positive' and set(raw)=={'research_kind','potential','basis'}:
        q,basis=positive.setup(raw['potential'],raw['basis'])
        if list(map(str,q))!=raw['potential'] or [list(map(str,b)) for b in basis]!=raw['basis']:raise ValueError('Noncanonical positive task')
    elif kind=='uniform' and set(raw)=={'research_kind','direction','coefficient','amplitudes','rule'}:
        direction=raw['direction']
        if not isinstance(direction,list) or not 1<=len(direction)<=7:raise ValueError('Uniform direction')
        values=[F(x) for x in direction];F(raw['coefficient'])
        if list(map(str,values))!=direction:raise ValueError('Noncanonical uniform direction')
        if raw['rule'] not in ('barta','log_sobolev'):raise ValueError('Uniform rule')
        domain=raw['amplitudes']
        if domain!='all_real' and (not isinstance(domain,list) or len(domain)!=2 or F(domain[0])>=F(domain[1])):raise ValueError('Uniform domain')
        if raw['rule']=='log_sobolev' and (len(values)>2 or domain!='all_real'):raise ValueError('Log-Sobolev task scope')
    else:raise ValueError('Unsupported research task reference')
    return raw
