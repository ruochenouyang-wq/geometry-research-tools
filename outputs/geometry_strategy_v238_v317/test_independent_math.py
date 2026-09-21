"""Independent equation-level checks, not counted as campaign increments."""
from copy import deepcopy
from fractions import Fraction as F
import unittest

import adaptive_singular as a


def integral_power(r):
    """Direct integral of |z|^r against (1-z^2) dz."""
    r = F(r)
    if r <= -1:
        raise ValueError('divergent moment')
    return F(2)/(r+1)-F(2)/(r+3)


def independent_matrices(function, powers):
    alpha, amplitude, offset = (F(function[k]) for k in ('exponent', 'amplitude', 'offset'))
    ss = tuple(map(F, powers))
    mass, energy = [], []
    for s in ss:
        mr, hr = [], []
        for t in ss:
            m = integral_power(s+t)
            kinetic = 2*m
            if s and t:
                r = s+t-2
                kinetic += s*t*(F(2)/(r+1)-F(4)/(r+3)+F(2)/(r+5))
            mr.append(m)
            hr.append(kinetic+amplitude*integral_power(s+t+alpha)+offset*m)
        mass.append(tuple(mr)); energy.append(tuple(hr))
    return tuple(mass), tuple(energy)


def strong_statistics(function, powers, coefficients):
    """Differentiate the radial polynomial, then directly integrate its image."""
    alpha, amplitude, offset = (F(function[k]) for k in ('exponent', 'amplitude', 'offset'))
    g = dict(zip(map(F, powers), map(F, coefficients)))
    image = {}
    def add(power, coefficient):
        if coefficient:
            image[power] = image.get(power, F(0))+coefficient
    for s, c in g.items():
        # -(1-z^2)g'' + 4zg' + 2g + qg.
        add(s-2, -c*s*(s-1))
        add(s, c*s*(s-1))
        add(s, 4*c*s)
        add(s, (2+offset)*c)
        add(s+alpha, amplitude*c)
    image = {s:c for s,c in image.items() if c}
    def inner(v, w):
        return sum((x*y*integral_power(r+s) for r,x in v.items()
                    for s,y in w.items()), F(0))
    mass = inner(g,g)
    mu = inner(g,image)/mass
    residual = image.copy()
    for s,c in g.items(): residual[s] = residual.get(s,F(0))-mu*c
    residual = {s:c for s,c in residual.items() if c}
    return mass, mu, inner(residual,residual)/mass, residual


def independent_solve(matrix, rhs):
    n = len(rhs)
    rows = [list(map(F,row))+[F(rhs[i])] for i,row in enumerate(matrix)]
    for j in range(n):
        pivot = next(i for i in range(j,n) if rows[i][j])
        rows[j],rows[pivot] = rows[pivot],rows[j]
        scale = rows[j][j]
        rows[j] = [x/scale for x in rows[j]]
        for i in range(n):
            if i != j:
                factor = rows[i][j]
                rows[i] = [x-factor*y for x,y in zip(rows[i],rows[j])]
    return [row[-1] for row in rows]


def independent_projection(function, powers, coefficients, transform):
    mass,_,complete,residual = strong_statistics(function,powers,coefficients)
    ss = tuple(map(F,powers))
    columns = [{s:F(transform[i][j]) for i,s in enumerate(ss) if transform[i][j]}
               for j in range(len(transform[0]))]
    def inner(v,w):
        return sum((x*y*integral_power(r+s) for r,x in v.items()
                    for s,y in w.items()),F(0))
    gram = [[inner(v,w) for w in columns] for v in columns]
    rhs = [inner(v,residual) for v in columns]
    solved = independent_solve(gram,rhs)
    projected = {}
    for c,column in zip(solved,columns):
        for s,v in column.items(): projected[s] = projected.get(s,F(0))+c*v
    outside = dict(residual)
    for s,v in projected.items(): outside[s] = outside.get(s,F(0))-v
    return inner(projected,projected)/mass,inner(outside,outside)/mass,complete


class IndependentMathTests(unittest.TestCase):
    def test_original_weak_gradient_gram_matches_operator_matrix(self):
        ss = ['0','7/4','2','7/2','15/4']
        M,H,_ = a.trial_matrices(a.ORIGINAL,ss)
        expected_M,expected_H = independent_matrices(a.ORIGINAL,ss)
        self.assertEqual(tuple(map(tuple,M)),expected_M)
        self.assertEqual(tuple(map(tuple,H)),expected_H)

    def test_original_full_strong_residual_by_direct_differentiation(self):
        ss,cs = ['0','7/4','2'],['1','-16/21','1/5']
        mass,mu,residual,terms = strong_statistics(a.ORIGINAL,ss,cs)
        result = a.statistics(a.ORIGINAL,ss,cs)
        self.assertEqual(F(result['mass']),mass)
        self.assertEqual(F(result['rayleigh']),mu)
        self.assertEqual(F(result['residual_squared']),residual)
        self.assertEqual(a.residual_terms(a.ORIGINAL,ss,cs),terms)
        self.assertNotIn(F(-1,4),terms)

    def test_signed_rescaling_preserves_rayleigh_and_full_residual(self):
        ss,cs = ['0','7/4','2'],[F(1),F(-4,7),F(1,8)]
        original = a.statistics(a.ORIGINAL,ss,cs)
        scaled = a.statistics(a.ORIGINAL,ss,[-F(13,11)*c for c in cs])
        for key in ('rayleigh','residual_squared'):
            self.assertEqual(F(original[key]),F(scaled[key]))

    def test_distributional_and_borderline_powers_are_rejected(self):
        for ss in ([0,1],[0,F(3,2)],[0,F(7,5)],[0,F(7,4),F(7,4)]):
            with self.subTest(powers=ss), self.assertRaises(ValueError):
                a.trial_matrices(a.ORIGINAL,ss)

    def test_residual_frontier_is_smoothing_not_raw_residual(self):
        raw = a.residual_terms(a.ORIGINAL,[0],[1])
        self.assertIn(F(-1,4),raw)
        proposed = a.residual_frontier(a.ORIGINAL,[0],[1])
        self.assertIn(F(7,4),proposed)
        self.assertTrue(all(s>F(3,2) for s in proposed))

    def test_complete_source_binding_rejects_another_potential(self):
        source = a.source_certificate(a.ORIGINAL)
        bad = deepcopy(source); bad['function']['amplitude'] = '-9/10'
        with self.assertRaises(ValueError):
            a.certificate(a.ORIGINAL,[0],[1],bad)

    def test_residual_pythagoras_matches_direct_function_projection(self):
        ss,cs = ['0','7/4','2'],['1','-4/7','1/5']
        for transform in ([[1,0,0],[0,1,0],[0,0,1]],
                          [[1,2,0],[0,1,3],[0,0,1]],
                          [[1,0],[F(-16,21),0],[0,1]]):
            with self.subTest(transform=transform):
                inside,outside,complete = independent_projection(a.ORIGINAL,ss,cs,transform)
                got = a.residual_diagnosis(a.ORIGINAL,ss,cs,transform)
                self.assertEqual(F(got['inside_squared']),inside)
                self.assertEqual(F(got['outside_squared']),outside)
                self.assertEqual(F(got['complete_squared']),complete)
                self.assertEqual(inside+outside,complete)
                self.assertGreater(outside,0)

    def test_exact_leading_cancellation_really_removes_singular_action(self):
        ss = ['0','7/4','2']
        C = a.cancellation_transform(a.ORIGINAL,ss)
        v = [sum((F(x)*y for x,y in zip(row,[F(1),F(1,5)])),F(0)) for row in C]
        _,_,_,terms = strong_statistics(a.ORIGINAL,ss,v)
        self.assertNotIn(F(-1,4),terms)

    def test_mass_orthogonalization_by_integrating_transformed_functions(self):
        ss = [F(0),F(7,4),F(2),F(7,2)]
        for base in (None,a.cancellation_transform(a.ORIGINAL,ss)):
            C = a.mass_orthogonal_transform(a.ORIGINAL,ss,base)
            n = len(C[0])
            for i in range(n):
                for j in range(n):
                    product = sum((F(C[k][i])*F(C[l][j])*integral_power(s+t)
                                   for k,s in enumerate(ss) for l,t in enumerate(ss)),F(0))
                    if i == j: self.assertGreater(product,0)
                    else: self.assertEqual(product,0)
            if base is not None:
                for j in range(n):
                    self.assertEqual(F(C[1][j]),-F(16,21)*F(C[0][j]))

    def test_shift_is_strictly_below_verified_m1_lower_bound(self):
        source = a.source_certificate(a.ORIGINAL)
        ss = ['0','7/4','2']
        p = a.proposal(a.ORIGINAL,ss,iterations=2,shift_policy='certified',source=source)
        self.assertLess(F(p['inverse_shift']),F(source['sectors'][1]['lower']))
        _,mu,residual,_ = strong_statistics(a.ORIGINAL,ss,p['coefficients'])
        stats = a.statistics(a.ORIGINAL,ss,p['coefficients'])
        self.assertEqual(mu,F(stats['rayleigh']))
        self.assertEqual(residual,F(stats['residual_squared']))
        self.assertEqual(residual,F(p['diagnosis']['complete_squared']))

    def test_certified_shift_rejects_tampered_spectral_source(self):
        source = deepcopy(a.source_certificate(a.ORIGINAL))
        source['sectors'][1]['lower'] = '10'
        with self.assertRaises(ValueError):
            a.proposal(a.ORIGINAL,['0','7/4'],iterations=2,shift_policy='certified',source=source)

    def test_residual_refinement_respects_complete_temple_objective(self):
        ss = ['0','7/4','2']
        initial = a.proposal(a.ORIGINAL,ss,iterations=2)['coefficients']
        refined = a.residual_refinement(a.ORIGINAL,ss,initial,steps=2)
        _,before_mu,before_r,_ = strong_statistics(a.ORIGINAL,ss,initial)
        _,after_mu,after_r,_ = strong_statistics(a.ORIGINAL,ss,refined['coefficients'])
        beta = F(a.gap_bound(a.ORIGINAL)['beta'])
        self.assertLess(before_mu,beta)
        self.assertLess(after_mu,beta)
        self.assertLessEqual(after_r/(beta-after_mu),before_r/(beta-before_mu))

    def test_zero_step_refinement_preserves_exact_function(self):
        ss,cs = ['0','7/4','2'],['1','-16/21','1/5']
        result = a.residual_refinement(a.ORIGINAL,ss,cs,steps=0,
                                      transform=a.cancellation_transform(a.ORIGINAL,ss))
        self.assertEqual(list(map(F,result['coefficients'])),list(map(F,cs)))
        self.assertEqual(result['trace'],[])

    def test_parity_gap_uses_12_only_for_even_and_retains_odd_6(self):
        gap = a.gap_bound(a.ORIGINAL,'even_odd')
        eta = F(gap['form']['eta_upper'])
        self.assertGreaterEqual(eta*eta,F(2,9))
        self.assertEqual(F(gap['beta']),F(32,3)-13*eta)
        self.assertEqual(F(gap['odd_component_lower']),F(14,3)-7*eta)
        self.assertEqual(F(gap['free_second_eigenvalue']),12)
        self.assertEqual(F(gap['odd_component_free_ground']),6)

    def test_parity_certificate_recombines_full_m1_and_global_sectors(self):
        source = a.source_certificate(a.ORIGINAL)
        ss = ['0','7/4','2']
        v = a.proposal(a.ORIGINAL,ss,iterations=4)['coefficients']
        c = a.certificate(a.ORIGINAL,ss,v,source,gap_policy='even_odd')
        _,mu,r,_ = strong_statistics(a.ORIGINAL,ss,v)
        beta,odd = F(c['gap']['beta']),F(c['gap']['odd_component_lower'])
        complete_m1 = max(min(mu-r/(beta-mu),odd),F(source['sectors'][1]['lower']))
        global_low = min([complete_m1,F(source['angular_tail_lower'])]+
                         [F(s['lower']) for s in source['sectors'] if s['azimuth_m']!=1])
        self.assertEqual(F(c['m1_lower']),complete_m1)
        self.assertEqual(F(c['lower']),global_low)
        self.assertTrue(a.verify(c))
        for key,value in [('odd_component_lower','100'),
                          ('even_and_odd_cover_full_radial_m1',False),
                          ('free_second_eigenvalue','6')]:
            bad = deepcopy(c);bad['gap'][key] = value
            self.assertFalse(a.verify(bad))

    def test_affine_potential_matrices_and_residuals_from_original_equation(self):
        for amplitude,offset in [('-1/2','-7'),('0','-3'),('1/3','5/7')]:
            q = dict(a.ORIGINAL,amplitude=amplitude,offset=offset)
            ss,cs = ['0','7/4','2'],['1','-2/7','3/11']
            with self.subTest(amplitude=amplitude,offset=offset):
                M,H,_ = a.trial_matrices(q,ss)
                em,eh = independent_matrices(q,ss)
                self.assertEqual(M,em);self.assertEqual(H,eh)
                mass,mu,r,terms = strong_statistics(q,ss,cs)
                got = a.statistics(q,ss,cs)
                self.assertEqual(F(got['mass']),mass)
                self.assertEqual(F(got['rayleigh']),mu)
                self.assertEqual(F(got['residual_squared']),r)
                self.assertEqual(a.residual_terms(q,ss,cs),terms)

    def test_offset_transport_preserves_candidates_residuals_and_width(self):
        q = dict(a.ORIGINAL,amplitude='-1/2')
        shift = -F(17,3)
        translated = dict(q,offset=str(shift))
        ss = ['0','7/4','2']
        s1,s2 = (a.source_certificate(f,modes=4,bits=16) for f in (q,translated))
        p1 = a.proposal(q,ss,iterations=3,source=s1)
        p2 = a.proposal(translated,ss,iterations=3,source=s2)
        self.assertEqual(p1['coefficients'],p2['coefficients'])
        c1 = a.certificate(q,ss,p1['coefficients'],s1,gap_policy='even_odd')
        c2 = a.certificate(translated,ss,p2['coefficients'],s2,gap_policy='even_odd')
        for key in ('lower','upper'):
            self.assertEqual(F(c2[key])-F(c1[key]),shift)
        self.assertEqual(c1['exact_width'],c2['exact_width'])
        self.assertEqual(c1['statistics']['residual_squared'],c2['statistics']['residual_squared'])
        self.assertTrue(a.verify(c2,translated,True))

    def test_exact_free_even_eigenfunctions_distinguish_ground_and_excited(self):
        q = dict(a.ORIGINAL,amplitude='0')
        source = a.source_certificate(q,modes=3,bits=16)
        ground = a.statistics(q,[0],[1])
        self.assertEqual(F(ground['rayleigh']),2)
        self.assertEqual(F(ground['residual_squared']),0)
        excited = a.statistics(q,[0,2],[1,-5])
        self.assertEqual(F(excited['rayleigh']),12)
        self.assertEqual(F(excited['residual_squared']),0)
        with self.assertRaises(ValueError):
            a.certificate(q,[0,2],[1,-5],source,gap_policy='even_odd')
        c = a.certificate(q,[0],[1],source,gap_policy='even_odd')
        self.assertLessEqual(F(c['lower']),2)
        self.assertGreaterEqual(F(c['upper']),2)
        self.assertTrue(a.verify(c))

    def test_parameter_family_moments_residual_and_semigroup(self):
        cases = [('-1/5','-2/3','1/7'),('-1/3','-1/2','-2'),
                 ('-2/5','1/3','0'),('-49/100','-1/10','3/11')]
        for alpha,amplitude,offset in cases:
            q = dict(a.ORIGINAL,exponent=alpha,amplitude=amplitude,offset=offset)
            delta = 2+F(alpha)
            ss = [F(0),delta,F(2)]
            cs = [F(1),F(amplitude)/(delta*(delta-1)),F(1,5)]
            with self.subTest(alpha=alpha,amplitude=amplitude):
                expected = tuple(sorted({i*delta+2*j for i in range(17)
                                         for j in range(17)})[:16])
                self.assertEqual(a.generated_powers(q,16),expected)
                self.assertTrue(all(s == 0 or s>F(3,2) for s in expected))
                M,H,_ = a.trial_matrices(q,ss)
                em,eh = independent_matrices(q,ss)
                self.assertEqual(M,em);self.assertEqual(H,eh)
                mass,mu,r,terms = strong_statistics(q,ss,cs)
                got = a.statistics(q,ss,cs)
                self.assertEqual(F(got['mass']),mass)
                self.assertEqual(F(got['rayleigh']),mu)
                self.assertEqual(F(got['residual_squared']),r)
                self.assertEqual(a.residual_terms(q,ss,cs),terms)
                self.assertNotIn(F(alpha),terms)
                gap = a.gap_bound(q,'even_odd')
                eta = F(gap['form']['eta_upper'])
                variance = F(amplitude)**2*(1/(2*F(alpha)+1)-1/(F(alpha)+1)**2)
                center = F(amplitude)/(F(alpha)+1)+F(offset)
                self.assertGreaterEqual(eta*eta,variance)
                self.assertLess(eta,1)
                self.assertEqual(F(gap['beta']),12*(1-eta)+center-eta)
                self.assertEqual(F(gap['odd_component_lower']),6*(1-eta)+center-eta)

    def test_family_rejects_l2_boundary_and_unavailable_coercivity(self):
        for exponent in ('-1/2','-3/4','0','1/10'):
            with self.subTest(exponent=exponent),self.assertRaises(ValueError):
                a.normalize(dict(a.ORIGINAL,exponent=exponent))
        with self.assertRaises(ValueError):
            a.normalize(dict(a.ORIGINAL,exponent='-2/5',amplitude='-1'))

    def test_family_full_solver_returns_bound_when_temple_gap_is_unavailable(self):
        q = dict(a.ORIGINAL,exponent='-1/3')
        for policy in ('nested','residual'):
            with self.subTest(policy=policy):
                result = a.solve(q,max_terms=3,iterations=2,basis_policy=policy,
                                 source_modes=3,source_bits=16)
                c = result['certificate']
                self.assertEqual(result['status'],'certified_open')
                self.assertTrue(a.verify(c,q,True,'1/100000000'))
                self.assertEqual(c['format'],a.direct.FULL)
                self.assertTrue(any(t['status']=='temple_gap_unavailable'
                                    for t in result['attempts']))

    def test_nonoriginal_parameter_point_has_replayable_full_certificate(self):
        q = dict(a.ORIGINAL,exponent='-1/3',amplitude='-1/2',offset='-2',axis=0)
        result = a.solve(q,max_terms=3,iterations=3,source_modes=3,source_bits=16,
                         leading_cancellation=True,orthogonalize=True,
                         adaptive_iterations=True,residual_steps=1)
        self.assertTrue(a.verify(result['certificate'],q,True,'1/100000000'))
        bad = deepcopy(result['certificate']);bad['function']['axis'] = 1
        self.assertFalse(a.verify(bad,expected_function=q))


if __name__ == '__main__':
    unittest.main()
