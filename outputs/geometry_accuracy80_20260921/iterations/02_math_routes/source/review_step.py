"""Independent public-case audit of step_solver; does not inspect holdout data."""
from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import time
import unittest
from unittest.mock import patch
import step_solver as step


def add(*polynomials):
    n = max((len(p) for p in polynomials), default=0)
    return [sum((p[i] if i < len(p) else F(0) for p in polynomials), F(0)) for i in range(n)]


def scale(p, a):
    return [F(a)*v for v in p]


def multiply(p, q):
    out = [F(0)]*max(0, len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): out[i+j] += a*b
    return out


def derivative(p):
    return [F(i)*p[i] for i in range(1,len(p))]


def integral(p, side):
    # Integrate at endpoints directly, with no calls to the solver's moments.
    lo, hi = (F(-1), F(0)) if side < 0 else (F(0), F(1))
    return sum((a*(hi**(i+1)-lo**(i+1))/F(i+1) for i,a in enumerate(p)), F(0))


def independent_statistics(function, terms, coefficients, m):
    w = [F(1),F(0),F(-1)]
    mass = form = energy = norm = F(0)
    traces = []
    for side in (-1,1):
        n=max(k for _,k in terms)+1
        r=[F(0)]*n
        for (support,k),a in zip(terms,coefficients):
            if support in (0,side): r[k] += F(a)
        dr=derivative(r); ddr=derivative(dr)
        q=F(function.get('offset',0))+(F(function.get('amplitude',1)) if side>0 else 0)
        # H[(1-z²)^(m/2) r] divided by (1-z²)^(m/2).
        hr=add(scale(multiply(w,ddr),-1), scale([F(0)]+dr,2*(m+1)),
               scale(r,m*(m+1)+q))
        weight=w if m else [F(1)]
        mass += integral(multiply(weight,multiply(r,r)),side)
        form += integral(multiply(weight,multiply(r,hr)),side)
        norm += integral(multiply(weight,multiply(hr,hr)),side)
        if m == 0:
            gradient=multiply(w,multiply(dr,dr))
        else:
            # Compute the original spherical gradient before integration by parts.
            gradient=add(multiply(multiply(w,w),multiply(dr,dr)),
                         scale(multiply([F(0)]+w,multiply(r,dr)),-2),
                         multiply([F(1),F(0),F(1)],multiply(r,r)))
        energy += integral(add(gradient,scale(multiply(weight,multiply(r,r)),q)),side)
        traces.append((r[0],dr[0] if dr else F(0)))
    if traces[0] != traces[1]:
        raise ValueError('Independent audit detects a value or first derivative jump')
    if energy != form:
        raise AssertionError('Strong operator and independent spherical energy disagree')
    return {'mass':str(mass),'operator_form':str(form),
            'operator_norm_squared':str(norm),'rayleigh':str(form/mass),
            'residual_squared':str(norm/mass-(form/mass)**2)}


class StepReview(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.q={'kind':'axis_profile','profile':'step','axis':2,'amplitude':'1','offset':'0'}
        cls.runs={}
        for mean_zero in (False,True):
            started=time.perf_counter()
            result=step.solve(cls.q,degree=8,mean_zero=mean_zero,tolerance='1/100000000')
            cls.runs[mean_zero]=result
            print(json.dumps({'public_function':cls.q,'mean_zero':mean_zero,
                              'degree':8,'width':result['certificate']['exact_width'],
                              'status':result['status'],'audit_seconds':time.perf_counter()-started}))

    def test_independent_energy_action_and_full_residual(self):
        terms=step.basis(5)
        coefficients=[str(F((-1)**i*(i+1),i+3)) for i in range(len(terms))]
        for m in (0,1):
            expected=independent_statistics(self.q,terms,coefficients,m)
            self.assertEqual(expected,step.trial_statistics(self.q,terms,coefficients,m))
        for mean_zero,run in self.runs.items():
            cert=run['certificate']; trial=cert['trial']
            expected=independent_statistics(cert['function'],trial['basis'],trial['coefficients'],
                                            trial['azimuth_m'])
            self.assertEqual(expected,cert['statistics'])

    def test_constant_potentials_include_known_complete_spectrum(self):
        for mean_zero in (False,True):
            q={**self.q,'amplitude':'0','offset':'3/2'}
            cert=step.solve(q,degree=2,mean_zero=mean_zero)['certificate']
            exact=F(7,2) if mean_zero else F(3,2)
            self.assertLessEqual(F(cert['lower']),exact)
            self.assertGreaterEqual(F(cert['upper']),exact)
            self.assertEqual(F(cert['statistics']['residual_squared']),0)
            self.assertEqual(F(cert['statistics']['rayleigh']),exact)
            self.assertTrue(step.verify(cert,expected_function=q,expected_mean_zero=mean_zero))

    def test_unshifted_step_variational_sanity(self):
        # Constant trial on full sphere, and x/y/z on zero-mean sphere, all see mean potential 1/2.
        for mean_zero,run in self.runs.items():
            cert=run['certificate']; lower=F(2) if mean_zero else F(0)
            upper=F(5,2) if mean_zero else F(1,2)
            self.assertGreaterEqual(F(cert['lower']),lower)
            self.assertLessEqual(F(cert['upper']),upper)

    def test_C1_violations_and_pole_endpoint_basis(self):
        source=self.runs[True]['certificate']['source']
        for terms in (((0,0),(1,1)),((0,0),(-1,0)),((0,0),(1,True)),((0,0),(1,-2))):
            with self.assertRaises(ValueError):
                step.trial_certificate(self.q,terms,['1','1'],source)
        # Distinct piecewise second derivatives are allowed, first derivative is continuous.
        terms=((0,0),(0,1),(-1,2),(1,2))
        coeff=('1','1/7','1/10','-1/10')
        for m in (0,1):
            self.assertEqual(independent_statistics(self.q,terms,coeff,m),
                             step.trial_statistics(self.q,terms,coeff,m))

    def test_verifier_does_not_depend_on_search_or_Gram_assembly(self):
        with patch.object(step,'_matrices',side_effect=AssertionError('Gram disabled')), \
             patch.object(step,'trial_statistics',side_effect=AssertionError('trial assembly disabled')), \
             patch.object(step,'propose_trial',side_effect=AssertionError('search disabled')):
            for run in self.runs.values(): self.assertTrue(step.verify(run['certificate']))

    def test_wrong_bounds_residuals_gap_and_source_are_rejected(self):
        cert=self.runs[True]['certificate']
        changes=[(['lower'],'999'),(['upper'],'-999'),(['exact_width'],'0'),
                 (['status'],'made_up_success'),(['temple_lower'],'999'),
                 (['statistics','residual_squared'],'0'),(['statistics','mass'],'1'),
                 (['gap','lower'],'999'),(['gap','azimuth_m'],0),
                 (['source','angular_tail_lower'],'999'),(['source_digest'],'0'*64),
                 (['scope'],'all_real_H1_on_unit_S2'),(['eigenvalue_index'],2),
                 (['trial','azimuth_m'],0),(['trial','strong_operator_domain'],'no_matching_needed')]
        for path,value in changes:
            altered=deepcopy(cert); target=altered
            for key in path[:-1]: target=target[key]
            target[path[-1]]=value
            with self.subTest(path=path): self.assertFalse(step.verify(altered))
        altered=deepcopy(cert); altered['statistics']['residual_squared']='NaN'
        self.assertFalse(step.verify(altered))

    def test_original_problem_space_and_target_binding(self):
        for mean_zero,run in self.runs.items():
            cert=run['certificate']
            self.assertFalse(step.verify(cert,expected_mean_zero=not mean_zero))
            self.assertFalse(step.verify(cert,expected_mean_zero=int(mean_zero)))
            self.assertFalse(step.verify(cert,expected_function={**self.q,'amplitude':'2'}))
            self.assertFalse(step.verify(cert,expected_tolerance='1/1000000000'))
            altered=deepcopy(cert); altered['mean_zero']=not mean_zero
            self.assertFalse(step.verify(altered))
        for bad in (None, [], {}, {'format':step.FORMAT}): self.assertFalse(step.verify(bad))

    def test_valid_custom_second_eigenvalue_gap_and_tampering(self):
        import spectral_gap
        for mean_zero in (False,True):
            m=int(mean_zero)
            gap=spectral_gap.certified_gap(self.q,m=m,mean_zero=mean_zero,modes=4,near_tail=4,bits=12)
            run=self.runs[mean_zero]; cert=run['certificate']; trial=cert['trial']
            new=step.trial_certificate(self.q,trial['basis'],trial['coefficients'],cert['source'],
                                       gap=gap,mean_zero=mean_zero)
            self.assertTrue(step.verify(new))
            altered=deepcopy(new); altered['gap']['lower']='1000'; altered['gap']['beta']='1000'
            self.assertFalse(step.verify(altered))

    def test_public_reflection_and_coordinate_isometry(self):
        for mean_zero,run in self.runs.items():
            cert=run['certificate']
            for q in ({**self.q,'amplitude':'-1','offset':'1'}, {**self.q,'axis':0}):
                transformed=step.solve(q,degree=8,mean_zero=mean_zero)['certificate']
                self.assertTrue(step.verify(transformed,expected_function=q,expected_mean_zero=mean_zero))
                self.assertEqual(F(cert['lower']),F(transformed['lower']))
                self.assertEqual(F(cert['upper']),F(transformed['upper']))

    def test_documented_public_H14_development_case(self):
        # Parameters disclosed in STEP_METHOD.md, not a held-out input or seed.
        q={**self.q,'amplitude':'-3/5'}
        for degree,limit in ((8,F(1,10**10)),(10,F(1,10**12))):
            run=step.solve(q,degree=degree,mean_zero=True,tolerance='1/10000000000')
            cert=run['certificate']; trial=cert['trial']
            self.assertTrue(step.verify(cert,expected_function=q,expected_mean_zero=True,
                                        expected_tolerance='1/10000000000'))
            self.assertLess(F(cert['exact_width']),limit)
            self.assertEqual(independent_statistics(q,trial['basis'],trial['coefficients'],1),
                             cert['statistics'])
            print(json.dumps({'public_development_case':'H14','degree':degree,
                              'width_display':float(F(cert['exact_width'])),
                              'independent_full_residual_checked':True}))


if __name__ == '__main__':
    path=Path(step.__file__)
    before=hashlib.sha256(path.read_bytes()).hexdigest()
    print('reviewed_step_sha256='+before)
    program=unittest.main(verbosity=2,exit=False)
    unchanged=hashlib.sha256(path.read_bytes()).hexdigest()==before
    print('source_unchanged_during_review='+str(unchanged))
    raise SystemExit(0 if program.result.wasSuccessful() and unchanged else 1)
