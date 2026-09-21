"""V185--V194: exact geometric scale transfer and certified radius searches.

q(t) uses t=x3/R in normalized coordinates, or q(x3) in ambient coordinates.
All spectral claims concern the compression to all real mean-zero H1(S2_R).
Natural area is used: masses scale by R^2 and Dirichlet energies do not.
The min--max perturbation rule and sphere dilation are analytical dependencies;
these rational software certificates are not formal-assistant proofs.
"""
from copy import deepcopy
from common import F, exact, same, digest, save, projected, sectors
import wide_potential as wide
import precision
import ordered_spectrum as ordered
import witness

POINT = 'sphere_radius_ground_v185'
PULLBACK = 'ambient_polynomial_pullback_v186'
TRIAL = 'sphere_radius_trial_v187'
ORDERED = 'sphere_radius_ordered_v188'
CELL = 'uniform_radius_cell_v189'
GLOBAL = 'radius_global_minimum_v190'
CHECKPOINT = 'radius_search_checkpoint_v191'
VIOLATION = 'sphere_radius_counterexample_v192'
JOB = 'sphere_radius_research_job_v194'
SCOPE = 'all_real_mean_zero_H1_on_S2_radius_R'


def _q(raw):
    return wide.potential(raw)


def _radius(raw):
    r = exact(raw)
    if r <= 0 or r > 32 or r.denominator > 2**128:
        raise ValueError('Radius must be exact, positive, <=32 with bounded denominator')
    return r


def _coordinates(value):
    if value not in ('normalized_t', 'ambient_x3'):
        raise ValueError('Coordinates must be normalized_t or ambient_x3')
    return value


def _problem(q, coordinates):
    return {'geometry': 'round_two_sphere', 'scope': SCOPE, 'mean_zero': True,
            'q_coefficients': list(map(str, _q(q))), 'coordinates': _coordinates(coordinates),
            'measure': 'natural_area_on_S2_R'}


def _interval(raw):
    if not isinstance(raw, (list, tuple)) or len(raw) != 2:
        raise ValueError('A closed radius interval is required')
    a, b = map(_radius, raw)
    if a >= b:
        raise ValueError('The radius interval must have positive length')
    return a, b


def _tolerance(raw):
    t = exact(raw)
    if not F(1,10**25) <= t <= 1:
        raise ValueError('Tolerance must be exact and in [1e-25,1]')
    return t


def _unit_q(q, r, coordinates):
    q, r, coordinates = _q(q), _radius(r), _coordinates(coordinates)
    return _q([a*r**(j+2 if coordinates == 'ambient_x3' else 2) for j,a in enumerate(q)])


def physical_pullback(q, radius):
    """V186: pull q(x3) back under x=R*y and rescale the whole operator."""
    q, r = _q(q), _radius(radius)
    physical = [a*r**j for j,a in enumerate(q)]
    return {'format': PULLBACK, 'original_q': list(map(str,q)), 'radius': str(r),
            'original_coordinate': 'ambient_x3', 'unit_coordinate': 't=y3=x3/R',
            'pulled_back_potential': list(map(str,physical)),
            'unit_operator_potential': list(map(str,_unit_q(q,r,'ambient_x3'))),
            'operator_rule': 'H_R is unitarily equivalent to R^-2*(-Delta_unit+R^2*q(R*t))'}


def _unit_verifier(name, cert, q):
    if name == 'projected_v94':
        return projected.verify_full(cert, expected_q=q, expected_mean_zero=True)
    if name == 'precision_v114':
        return precision.verify_full(cert, expected_q=q, expected_mean_zero=True)
    if name == 'wide_v124':
        return wide.verify_full(cert, expected_q=q, expected_mean_zero=True)
    return False


def _backend_name(cert):
    kind = cert.get('format')
    if kind == projected.FULL_FORMAT:
        return 'projected_v94'
    if kind == wide.FULL_FORMAT:
        return 'wide_v124'
    # Precision intentionally exposes its verifier as authority, not a hardcoded tag.
    if precision.verify_full(cert):
        return 'precision_v114'
    raise ValueError('Unsupported or invalid unit full-sphere evidence')


def scale_ground(q, radius, unit_certificate, coordinates='normalized_t'):
    """V185: exact full-space ground enclosure from unit-sphere evidence."""
    problem, r = _problem(q, coordinates), _radius(radius)
    p = _unit_q(q, r, coordinates)
    name = _backend_name(unit_certificate)
    if not _unit_verifier(name, unit_certificate, p):
        raise ValueError('Unit evidence has wrong potential or mean-zero scope')
    lo, hi = F(unit_certificate['lower'])/r**2, F(unit_certificate['upper'])/r**2
    return {'format': POINT, 'problem': problem, 'radius': str(r),
            'unit_q_coefficients': list(map(str,p)), 'unit_backend': name,
            'unit_certificate': deepcopy(unit_certificate), 'lower': str(lo),
            'upper': str(hi), 'width': str(hi-lo), 'eigenvalue_index': 1,
            'scale_factor': str(1/r**2),
            'rule': 'C_R(q)=R^-2*C_unit(R^2*q(R*t)) for ambient; q(t) for normalized',
            'formal_assistant_checked': False}


def scale_trial(q, radius, coefficients, m=0, degrees=None, coordinates='normalized_t'):
    """V187: exact natural-area mass, Dirichlet and potential energies.

    Wide-degree potential multiplication is exact; the unit trial's kinetic
    proof uses q=0 so it does not inherit V94's degree-six potential limit.
    """
    problem, r = _problem(q, coordinates), _radius(radius)
    p = _unit_q(q, r, coordinates)
    unit = witness.rayleigh_certificate([0], coefficients, m=m, mean_zero=True, degrees=degrees)
    u = dict(zip(unit['degrees'], map(F,unit['coefficients'])))
    mass = F(unit['norm_squared_divided_by_azimuth_factor'])
    kinetic = F(unit['kinetic_divided_by_azimuth_factor'])
    potential = F(0)
    for l,c in u.items():
        acted = wide.q_times_basis(p,m,l)
        potential += sum((c*v*u.get(k,F(0))*wide.basis_mass(m,k) for k,v in acted.items()),F(0))
    physical_mass = r*r*mass
    energy = kinetic+potential
    return {'format': TRIAL, 'problem': problem, 'radius': str(r),
            'unit_q_coefficients': list(map(str,p)), 'unit_kinetic_trial': unit,
            'function_rule': 'u_R(x)=sum c_l P_l^m(x3/R)*cos(m*phi)',
            'azimuth_integral_factor': unit['azimuth_integral_factor'],
            'mass_divided_by_azimuth_factor': str(physical_mass),
            'dirichlet_energy_divided_by_azimuth_factor': str(kinetic),
            'potential_energy_divided_by_azimuth_factor': str(potential),
            'total_energy_divided_by_azimuth_factor': str(energy),
            'mean_zero_verified': True, 'rayleigh_quotient': str(energy/physical_mass),
            'full_space_upper_bound': str(energy/physical_mass),
            'formal_assistant_checked': False}


def scale_ordered(q, radius, unit_ordered, coordinates='normalized_t'):
    """V188: transfer ordered eigenvalues, strict negative count and trace."""
    problem, r = _problem(q, coordinates), _radius(radius)
    p = _unit_q(q,r,coordinates)
    if unit_ordered.get('format') != ordered.ORDERED or not ordered.verify(unit_ordered,expected_q=p):
        raise ValueError('Expected matching complete mean-zero ordered spectrum')
    trace = ordered.negative_trace(unit_ordered['assembly'])
    count = ordered.count_below(unit_ordered['assembly'],0)
    intervals = [{'index':row['index'],'lower':str(F(row['lower'])/r**2),
                  'upper':str(F(row['upper'])/r**2)} for row in unit_ordered['ordered_intervals']]
    return {'format': ORDERED, 'problem': problem, 'radius': str(r),
            'unit_ordered': deepcopy(unit_ordered), 'unit_negative_trace': trace,
            'unit_strict_negative_count': count, 'ordered_intervals': intervals,
            'negative_trace_lower': str(F(trace['lower'])/r**2),
            'negative_trace_upper': str(F(trace['upper'])/r**2),
            'negative_count_lower': count['count_lower'], 'negative_count_upper': count['count_upper'],
            'zero_eigenvalues_are_not_negative': True,
            'multiplicities_unchanged_by_positive_dilation': True,
            'formal_assistant_checked': False}


def radius_cell(q, interval, anchor, coordinates='normalized_t'):
    """V189: pointwise uniform radius bound, including negative eigenvalues.

    For each coefficient, R^k is monotone on positive radii.  The l1 sum of
    coefficient deviations bounds the multiplication-operator norm. Multiplying
    all four interval corners prevents the common negative-eigenvalue sign bug.
    """
    problem, (a,b) = _problem(q,coordinates), _interval(interval)
    if anchor.get('format') != POINT or not verify(anchor, expected_q=q, expected_coordinates=coordinates):
        raise ValueError('A matching full mean-zero radius anchor is required')
    r0 = F(anchor['radius'])
    if not a <= r0 <= b:
        raise ValueError('Anchor lies outside the radius cell')
    deviations = []
    for j,c in enumerate(_q(q)):
        power = j+2 if coordinates == 'ambient_x3' else 2
        deviations.append(abs(c)*max(abs(a**power-r0**power),abs(b**power-r0**power)))
    delta = sum(deviations,F(0))
    source = anchor['unit_certificate']
    unit = [F(source['lower'])-delta,F(source['upper'])+delta]
    inverse_square = [1/b**2,1/a**2]
    products = [x*y for x in unit for y in inverse_square]
    lo, hi = min(products),max(products)
    analytic = None
    if len(_q(q)) == 1:
        c = _q(q)[0]
        analytic = [2/b**2+c,2/a**2+c]
        lo,hi=max(lo,analytic[0]),min(hi,analytic[1])
    if lo > hi:
        raise ArithmeticError('Inconsistent radius bounds')
    return {'format': CELL, 'problem': problem, 'interval':list(map(str,(a,b))),
            'anchor':deepcopy(anchor), 'coefficient_deviations':list(map(str,deviations)),
            'uniform_potential_perturbation':str(delta),
            'unit_eigenvalue_interval':list(map(str,unit)),
            'inverse_radius_squared_interval':list(map(str,inverse_square)),
            'four_signed_products':list(map(str,products)),
            'constant_potential_analytic_intersection':None if analytic is None else list(map(str,analytic)),
            'lower':str(lo),'upper':str(hi),
            'quantifier':'for_every_R_in_closed_interval',
            'lower_bound_inequality_quantifier':'for_every_nonzero_real_mean_zero_H1_function',
            'bound_quantity':'ground_eigenvalue_at_each_R', 'formal_assistant_checked':False}


def _budget(modes,max_modes,max_m,bits):
    if any(type(x) is not int for x in [modes,max_modes,max_m,bits]):
        raise ValueError('Budgets must be integers')
    if not 1 <= modes <= max_modes <= 32 or not 0 <= max_m <= 16 or not 8 <= bits <= 160:
        raise ValueError('Invalid spectral budget')
    return dict(modes=modes,max_modes=max_modes,max_m=max_m,bits=bits)


def backend_ground(q, radius, coordinates='normalized_t', tolerance=F(1,10**8),
                   modes=4,max_modes=12,max_m=4,bits=36):
    """V193: geometry-aware exact backend dispatch, including high degree.

    Selects precision at unit tolerance <1e-12 when V94 input restrictions
    permit it; otherwise selects the broad exact-polynomial kernel as needed.
    Returning a bound does not imply the requested width has been achieved.
    """
    r, tol = _radius(radius), _tolerance(tolerance)
    options = _budget(modes,max_modes,max_m,bits)
    p = _unit_q(q,r,coordinates)
    unit_tol = min(F(1),max(F(1,10**30),tol*r*r))
    narrow = True
    try:
        sectors.potential(p)
    except ValueError:
        narrow = False
    if narrow and unit_tol < F(1,10**12):
        cert = precision.full_ground(p,tolerance=max(F(1,10**25),unit_tol),
                    modes=modes,max_modes=max_modes,max_m=max_m,max_steps=16)
    elif narrow:
        cert = projected.full_ground(p,tolerance=unit_tol,**options)
    else:
        cert = wide.full_ground(p,tolerance=unit_tol,**options)
    return scale_ground(q,r,cert,coordinates)


class AnchorCache:
    """Verified content cache. Keys include geometry, coordinates, q and radius."""
    def __init__(self, q, coordinates='normalized_t', entries=None, tolerance=F(1,10**8), **options):
        self.problem=_problem(q,coordinates)
        self.q=_q(q); self.coordinates=coordinates; self.tolerance=_tolerance(tolerance)
        self.options=_budget(options.get('modes',4),options.get('max_modes',12),
                            options.get('max_m',4),options.get('bits',36))
        self.entries={}; self.computed=0; self.reused=0
        for cert in entries or []:
            self.add(cert)

    def add(self,cert):
        if cert.get('format') != POINT or not verify(cert,expected_q=self.q,expected_coordinates=self.coordinates):
            raise ValueError('Cached proof has a different problem or failed replay')
        r=str(_radius(cert['radius']))
        if r in self.entries and not same(self.entries[r],cert):
            # Keep the tighter interval by exact intersection would require a
            # new proof format. Cache scheduling selects a narrower proof;
            # continuation separately retains every old feasible incumbent.
            if F(self.entries[r]['width']) <= F(cert['width']):
                return
        self.entries[r]=deepcopy(cert)

    def get(self,radius):
        r=str(_radius(radius))
        if r in self.entries:
            self.reused+=1
            return deepcopy(self.entries[r])
        cert=backend_ground(self.q,r,self.coordinates,self.tolerance,**self.options)
        self.computed+=1
        self.entries[r]=deepcopy(cert)
        return cert

    def export(self):
        return [deepcopy(self.entries[r]) for r in sorted(self.entries,key=F)]


def _global_from_evidence(q,domain,cells,points,tolerance,coordinates):
    problem,(a,b),tol=_problem(q,coordinates),_interval(domain),_tolerance(tolerance)
    if not isinstance(cells,list) or not cells or len(cells)>65:
        raise ValueError('Require one to 65 cells covering the entire domain')
    previous=a
    for cell in cells:
        if cell.get('format') not in (CELL,'radius_cell_intersection_v191') or not verify(cell,expected_q=q,expected_coordinates=coordinates):
            raise ValueError('Invalid radius cell')
        x,y=map(F,cell['interval'])
        if x!=previous or y>b:
            raise ValueError('Radius cells must cover the domain without gaps or overlaps')
        previous=y
    if previous!=b:
        raise ValueError('Radius partition does not reach the domain endpoint')
    if not isinstance(points,list) or not points or len(points)>260:
        raise ValueError('Require finite verified feasible radius points')
    point_intervals=[]
    for point in points:
        if point.get('format')!=POINT or not verify(point,expected_q=q,expected_coordinates=coordinates):
            raise ValueError('Invalid feasible radius proof')
        if not a<=F(point['radius'])<=b:
            raise ValueError('Feasible radius is outside the domain')
        point_intervals.append((F(point['lower']),F(point['upper'])))
    lo=min(F(c['lower']) for c in cells)
    winner=min(range(len(points)),key=lambda i:(F(points[i]['upper']),F(points[i]['radius'])))
    hi=F(points[winner]['upper'])
    if lo>hi:
        raise ArithmeticError('Global lower exceeds feasible upper')
    maxlo=max(pair[0] for pair in point_intervals)
    maxhi=max(F(c['upper']) for c in cells)
    candidates=[c['interval'] for c in cells if F(c['lower'])<=hi]
    return {'format':GLOBAL,'problem':problem,'domain':list(map(str,(a,b))),
            'cells':deepcopy(cells),'feasible_points':deepcopy(points),
            'lower':str(lo),'upper':str(hi),'width':str(hi-lo),'tolerance':str(tol),
            'status':'certified_target_met' if hi-lo<=tol else 'certified_global_bound_open_gap',
            'best_feasible_radius':points[winner]['radius'],
            'possible_minimizer_cells':deepcopy(candidates),
            'maximum_lower':str(maxlo),'maximum_upper':str(maxhi),
            'uniqueness_claimed':False,'full_domain_coverage':True,
            'formal_assistant_checked':False}


def adaptive_minimum(q,domain,coordinates='normalized_t',tolerance=F(1,10000),
                     max_splits=8,cache=None,saved=None,**options):
    """V190: deterministic best-lower-first interval branch-and-bound.

    Every leaf is retained, even when irrelevant to the current incumbent.
    Thus coverage replay does not trust an unrecorded pruning decision.
    """
    a,b=_interval(domain);tol=_tolerance(tolerance)
    if type(max_splits)is not int or not 0<=max_splits<=64:
        raise ValueError('max_splits must be in 0..64')
    if cache is None:
        cache=AnchorCache(q,coordinates,tolerance=min(tol/16,F(1,10**8)),**options)
    if not same(cache.problem,_problem(q,coordinates)):
        raise ValueError('Cache belongs to a different geometry problem')
    previous_points=[]
    if saved is not None:
        if not verify_checkpoint(saved,expected_q=q,expected_coordinates=coordinates,expected_domain=domain):
            raise ValueError('Checkpoint is not bound to this exact radius domain')
        old=saved['certificate']
        previous_points=deepcopy(old['feasible_points'])
        cells=deepcopy(old['cells'])
        for point in old['feasible_points']:
            cache.add(point)
    else:
        cells=[radius_cell(q,[a,b],cache.get((a+b)/2),coordinates)]
    cache.get(a);cache.get(b)
    history=[]
    for step in range(max_splits+1):
        # Old and new proofs at the same radius can be incomparable. Keep
        # both: width-only replacement must never erase a better incumbent.
        combined={digest(p):p for p in previous_points+cache.export() if a<=F(p['radius'])<=b}
        points=[combined[k] for k in sorted(combined,key=lambda k:(F(combined[k]['radius']),k))]
        cert=_global_from_evidence(q,[a,b],cells,points,tol,coordinates)
        history.append({'step':step,'cells':len(cells),'lower':cert['lower'],'upper':cert['upper'],'width':cert['width']})
        if cert['status']=='certified_target_met' or step==max_splits or len(cells)==65:
            break
        idx=min(range(len(cells)),key=lambda i:(F(cells[i]['lower']),F(cells[i]['interval'][0])))
        oldcell=cells[idx];x,y=map(F,oldcell['interval']);mid=(x+y)/2
        # A shared boundary is a useful feasible point as well as a cache hit
        # when it was the parent anchor.
        cache.get(mid)
        children=[]
        for interval in [(x,mid),(mid,y)]:
            anchor=cache.get(sum(interval,F(0))/2)
            child=radius_cell(q,interval,anchor,coordinates)
            # Intersection with the parent's uniform enclosure is itself a
            # rigorous bound, retained explicitly as a nested parent proof.
            child=_restrict_cell(child,oldcell)
            children.append(child)
        cells[idx:idx+1]=children
    return {'certificate':cert,'search_history':history,
            'work':{'new_spectral_anchor_computations':cache.computed,'verified_anchor_cache_hits':cache.reused},
            'budget_exhausted_before_target':cert['status']!='certified_target_met'}


def _restrict_cell(child,parent):
    if parent.get('format') not in (CELL,'radius_cell_intersection_v191') or not verify(parent):
        raise ValueError('Invalid parent uniform enclosure')
    if child.get('format')!=CELL or not verify(child):
        raise ValueError('Invalid child uniform enclosure')
    if not same(child['problem'],parent['problem']):
        raise ValueError('Parent and child use different problems')
    a,b=map(F,child['interval']);x,y=map(F,parent['interval'])
    if not x<=a<b<=y:
        raise ValueError('Child interval outside parent')
    return {'format':'radius_cell_intersection_v191','problem':child['problem'],
            'interval':child['interval'],'child':deepcopy(child),'parent':deepcopy(parent),
            'lower':str(max(F(child['lower']),F(parent['lower']))),
            'upper':str(min(F(child['upper']),F(parent['upper']))),
            'monotone_parent_intersection':True}


def checkpoint(certificate):
    """V191: replayable continuation retaining exact geometry and all leaves."""
    if certificate.get('format')!=GLOBAL or not verify(certificate):
        raise ValueError('Can only checkpoint a verified complete radius search')
    return {'format':CHECKPOINT,'problem':deepcopy(certificate['problem']),
            'domain':certificate['domain'],'certificate_digest':digest(certificate),
            'certificate':deepcopy(certificate)}


def verify_checkpoint(saved,expected_q=None,expected_coordinates=None,expected_domain=None):
    try:
        if not same(saved,checkpoint(saved['certificate'])):
            return False
        if not verify(saved['certificate'],expected_q=expected_q,expected_coordinates=expected_coordinates):
            return False
        return expected_domain is None or saved['domain']==list(map(str,_interval(expected_domain)))
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError,RecursionError):
        return False


def resume_minimum(saved,max_splits=4,tolerance=None,**options):
    """V191: use verified old anchors and intersect descendants with parents."""
    if not verify_checkpoint(saved):
        raise ValueError('Invalid checkpoint')
    old=saved['certificate'];problem=old['problem']
    cache=AnchorCache(problem['q_coefficients'],problem['coordinates'],old['feasible_points'],
                     tolerance=min(F(tolerance or old['tolerance'])/16,F(1,10**8)),**options)
    return adaptive_minimum(problem['q_coefficients'],old['domain'],problem['coordinates'],
                            tolerance or old['tolerance'],max_splits,cache,saved,**options)


def transfer_counterexample(trial,threshold):
    """V192: retain an actual natural-area function violating a radius bound."""
    if trial.get('format')!=TRIAL or not verify(trial):
        raise ValueError('A valid explicit radius trial is required')
    threshold=exact(threshold);mu=F(trial['rayleigh_quotient'])
    return {'format':VIOLATION,'problem':deepcopy(trial['problem']),'radius':trial['radius'],
            'trial':deepcopy(trial),'threshold':str(threshold),
            'strict_margin':str(threshold-mu),'rayleigh_quotient':str(mu),
            'status':'explicit_counterexample' if mu<threshold else 'trial_does_not_refute',
            'assertion':'for_every_nonzero_real_mean_zero_u: energy(u)>=threshold*mass(u)',
            'formal_assistant_checked':False}


def _job_from_proofs(problem,evidence,threshold,counterexample=None):
    if evidence.get('format') not in (POINT,GLOBAL) or not verify(evidence):
        raise ValueError('A verified point or global radius bound is required')
    if not same(problem,evidence['problem']):
        raise ValueError('Research request and mathematical evidence differ')
    threshold=exact(threshold)
    status=('proved' if F(evidence['lower'])>=threshold else
            'refuted_by_spectral_existence' if F(evidence['upper'])<threshold else 'undetermined')
    if counterexample is not None:
        if counterexample.get('format')!=VIOLATION or not verify(counterexample,expected_threshold=threshold):
            raise ValueError('Invalid explicit counterexample or changed threshold')
        if not same(counterexample['problem'],problem):
            raise ValueError('Explicit function belongs to a different problem')
        r=F(counterexample['radius'])
        if evidence['format']==POINT and r!=F(evidence['radius']):
            raise ValueError('Explicit function belongs to a different radius')
        if evidence['format']==GLOBAL and not F(evidence['domain'][0])<=r<=F(evidence['domain'][1]):
            raise ValueError('Explicit function lies outside radius domain')
        if counterexample['status']=='explicit_counterexample':
            if status=='proved':
                raise ArithmeticError('Explicit violation contradicts lower certificate')
            status='refuted_by_explicit_function'
    # Content-address the nested mathematical authority. Work counters remain
    # diagnostic output of research_job, outside this proof's authority.
    nodes={digest(evidence):deepcopy(evidence)}
    roots={'spectral':digest(evidence)}
    if counterexample is not None:
        nodes[digest(counterexample)]=deepcopy(counterexample)
        roots['explicit_function']=digest(counterexample)
    return {'format':JOB,'problem':deepcopy(problem),'threshold':str(threshold),
            'assertion_quantifier':'every_R_in_request_and_every_nonzero_real_mean_zero_H1_function',
            'radius_request':({'radius':evidence['radius']} if evidence['format']==POINT else {'domain':evidence['domain']}),
            'status':status,'proof_nodes':nodes,'proof_roots':roots,
            'full_replay_required':True,'formal_assistant_checked':False}


def research_job(q,threshold,radius=None,domain=None,coordinates='normalized_t',
                 tolerance=F(1,10000),max_splits=4,trial=None,saved=None,**options):
    """V194: one full-scope prove/refute/open research entry point with replay."""
    if (radius is None)==(domain is None):
        raise ValueError('Choose exactly one fixed radius or radius interval')
    if radius is not None:
        if saved is not None:
            raise ValueError('Search checkpoints apply only to radius intervals')
        evidence=backend_ground(q,radius,coordinates,tolerance,**options)
        work={'new_spectral_anchor_computations':1,'verified_anchor_cache_hits':0}
    else:
        result=adaptive_minimum(q,domain,coordinates,tolerance,max_splits,saved=saved,**options)
        evidence,work=result['certificate'],result['work']
    ce=None
    if trial is not None:
        if not isinstance(trial,dict) or not set(trial)<= {'coefficients','m','degrees','radius'}:
            raise ValueError('Unknown trial specification')
        r=radius if radius is not None else trial.get('radius',evidence['best_feasible_radius'])
        if radius is not None and 'radius' in trial and _radius(trial['radius'])!=_radius(radius):
            raise ValueError('Trial radius conflicts with request')
        ce=transfer_counterexample(scale_trial(q,r,trial['coefficients'],trial.get('m',0),
                                  trial.get('degrees'),coordinates),threshold)
    cert=_job_from_proofs(_problem(q,coordinates),evidence,threshold,ce)
    return {'certificate':cert,'work':work,'model_token_measurement':None,
            'work_metric_scope':'local exact spectral anchor computations, not model tokens'}


def verify(cert,expected_q=None,expected_radius=None,expected_coordinates=None,expected_threshold=None):
    """Replay every rational transformation and every nested spectral proof."""
    try:
        kind=cert['format']
        p=cert.get('problem')
        if p is not None:
            if not same(p,_problem(p['q_coefficients'],p['coordinates'])):
                return False
            if expected_q is not None and p['q_coefficients']!=list(map(str,_q(expected_q))):
                return False
            if expected_coordinates is not None and p['coordinates']!=_coordinates(expected_coordinates):
                return False
        elif kind==PULLBACK:
            if expected_q is not None and cert['original_q']!=list(map(str,_q(expected_q))):
                return False
            if expected_coordinates is not None and expected_coordinates!='ambient_x3':
                return False
        elif expected_q is not None or expected_coordinates is not None:
            return False
        bound_radius=(cert.get('radius_request',{}).get('radius') if kind==JOB else cert.get('radius'))
        if expected_radius is not None and bound_radius!=str(_radius(expected_radius)):
            return False
        if expected_threshold is not None and cert.get('threshold')!=str(exact(expected_threshold)):
            return False
        if kind==PULLBACK:
            expected=physical_pullback(cert['original_q'],cert['radius'])
        elif kind==POINT:
            expected=scale_ground(p['q_coefficients'],cert['radius'],cert['unit_certificate'],p['coordinates'])
        elif kind==TRIAL:
            t=cert['unit_kinetic_trial']
            expected=scale_trial(p['q_coefficients'],cert['radius'],t['coefficients'],t['azimuth_m'],t['degrees'],p['coordinates'])
        elif kind==ORDERED:
            expected=scale_ordered(p['q_coefficients'],cert['radius'],cert['unit_ordered'],p['coordinates'])
        elif kind==CELL:
            expected=radius_cell(p['q_coefficients'],cert['interval'],cert['anchor'],p['coordinates'])
        elif kind=='radius_cell_intersection_v191':
            expected=_restrict_cell(cert['child'],cert['parent'])
        elif kind==GLOBAL:
            expected=_global_from_evidence(p['q_coefficients'],cert['domain'],cert['cells'],
                                          cert['feasible_points'],cert['tolerance'],p['coordinates'])
        elif kind==CHECKPOINT:
            expected=checkpoint(cert['certificate'])
        elif kind==VIOLATION:
            expected=transfer_counterexample(cert['trial'],cert['threshold'])
        elif kind==JOB:
            nodes,roots=cert['proof_nodes'],cert['proof_roots']
            if set(roots) not in ({'spectral'},{'spectral','explicit_function'}):
                return False
            if set(nodes)!=set(roots.values()) or any(digest(v)!=k for k,v in nodes.items()):
                return False
            expected=_job_from_proofs(p,nodes[roots['spectral']],cert['threshold'],
                                      nodes[roots['explicit_function']] if 'explicit_function'in roots else None)
        else:
            return False
        return same(cert,expected)
    except (ValueError,TypeError,KeyError,ArithmeticError,IndexError,RecursionError,OverflowError,AttributeError):
        return False
