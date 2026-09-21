"""Reproduce round 09; all outputs stay inside round09/."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from fractions import Fraction as F
import copy
import io
import json
import time
import unittest
import anisotropic as a
from common import save, projected


def jsonable(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        if any(isinstance(k, tuple) for k in value):
            return a.encode(value)
        return {k: jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(v) for v in value]
    return value


def record(version, capability, callable_name, value, outcome):
    evidence = save(f'round09/v{version}_evidence.json', jsonable(value))
    return {'version': version, 'capability': capability, 'callable': callable_name,
            'evidence': [evidence], 'outcome': outcome}


def main():
    q = {'1,1,0': '1'}
    stages = []
    try:
        projected.full_ground(q)
        baseline = {'rejected': False}
    except (ValueError, TypeError) as error:
        baseline = {'rejected': True, 'exception': type(error).__name__, 'reason': str(error)}
    baseline.update({'problem': 'Ground of -Delta_S2+xy compressed to all mean-zero H1 functions',
                     'weakness': 'An axisymmetric q(t) interface cannot encode xy and its cross-sector couplings',
                     'input': q, 'external_mathematical_basis': 'https://dlmf.nist.gov/14.30'})
    save('round09/baseline.json', baseline)
    p = a.polynomial({'1,0,0': '1/2', '0,1,0': '-1/3'})
    cheap = a.cheap_ground_enclosure(q)
    cheap_path = save('round09/certificates/xy_cheap.json', cheap)
    stages.append(record(175, 'Immediate full-space enclosure with exact linear trial upper bound',
                         'anisotropic.cheap_ground_enclosure / verify_cheap',
                         {'p': p, 'square': a.multiply(p, p), 'dx': a.derivative(p, 0),
                          'certificate': cheap_path, 'lower': cheap['lower'], 'upper': cheap['upper'],
                          'replayed': a.verify_cheap(cheap, expected_q=q)},
                         'Sparse xyz potential q=xy now immediately yields rigorous full-space [1,9/5]'))
    moments = {str(e): a.sphere_moment(e) for e in ((0, 0, 0), (2, 0, 0), (2, 2, 0), (2, 2, 2), (1, 2, 0))}
    stages.append(record(176, 'Exact normalized spherical moments', 'anisotropic.sphere_moment / integrate',
                         {'moments': moments, 'mean_r8': a.integrate(a.power(a.R2, 4))},
                         'All polynomial integrals use closed rational moments, not sampling'))
    stages.append(record(177, 'Rational real solid harmonic construction', 'anisotropic.solid_harmonic',
                         a.harmonic_basis(1, 3), '15 nonconstant real harmonics constructed and Laplacian-checked'))
    gram = a.gram_energy(3)
    stages.append(record(178, 'Exact Gram and tangential gradient forms', 'anisotropic.gram_energy',
                         {'mass': gram['mass'], 'gram': gram['gram'], 'energy': gram['energy']},
                         'All 15 functions orthogonal with independently integrated l(l+1) energy'))
    finite = a.coupled_form(q, 3)
    cross = [{'i': i, 'j': j, 'value': finite['potential_form'][i][j]}
             for i, p in enumerate(finite['basis']) for j, b in enumerate(finite['basis'])
             if i < j and p['m'] != b['m'] and finite['potential_form'][i][j]]
    stages.append(record(179, 'Coupled finite form includes cross-azimuth entries', 'anisotropic.coupled_form',
                         {'q': q, 'A': finite['A'], 'mass': finite['mass'], 'cross_m_entries': cross},
                         f'{len(cross)} off-diagonal entries between different azimuthal indices retained'))
    columns = a.tail_couplings(q, 3, finite)
    stages.append(record(180, 'Complete finite-to-tail coupling columns', 'anisotropic.tail_couplings',
                         {'retained_degree': 3, 'tail_degrees': [4, 5], 'columns': columns},
                         'All 20 harmonics of degrees 4 and 5 included with exact masses'))
    kernel = a.CoupledKernel(q, 2)
    stages.append(record(181, 'Degree-graded lower and upper Schur comparisons', 'anisotropic.CoupledKernel.matrix',
                         {'at': '9/5', 'tail_lower': kernel.beta,
                          'lower_matrix': kernel.matrix(F(9, 5), 'lower'),
                          'upper_tail_matrix': kernel.matrix(F(9, 5), 'upper_tail'),
                          'upper_ritz_matrix': kernel.matrix(F(9, 5), 'upper_ritz')},
                         'Positive infinite tail controlled with non-unit harmonic mass denominators'))
    cases = [
        ('zero', {}, 2), ('xy_L2', q, 2), ('xy_L3', q, 3),
        ('x_L3', {'1,0,0': '1'}, 3),
        ('degree4_mixed', {'1,1,0': '1', '1,1,1': '1/5', '2,0,2': '1/7'}, 3),
    ]
    certificates, results = {}, []
    for name, potential, L in cases:
        began = time.perf_counter()
        certificate = a.certify(potential, L=L, bits=28)
        seconds = time.perf_counter()-began
        path = save(f'round09/certificates/{name}.json', certificate)
        replay = a.verify(json.loads(json.dumps(certificate)), expected_q=potential)
        certificates[name] = certificate
        results.append({'case': name, 'q': potential, 'L': L, 'dimension': len(certificate['mass']),
                        'lower': certificate['lower'], 'upper': certificate['upper'],
                        'width': certificate['exact_width'], 'lower_decimal': float(F(certificate['lower'])),
                        'upper_decimal': float(F(certificate['upper'])), 'seconds': seconds,
                        'verified': replay, 'certificate': path})
    trial = a.rayleigh(q, {'1,0,0': '1', '0,1,0': '-1'})
    stages.append(record(182, 'Whole mean-zero sphere ground certification and exact replay',
                         'anisotropic.certify / verify', {'cases': results, 'independent_trial': trial},
                         'Five full-space enclosures, including genuinely mixed degree-4 xyz potential'))
    radius = {'2,0,0': '1', '0,2,0': '1', '0,0,2': '1'}
    radius_cert = a.certify(radius, L=2)
    radius_path = save('round09/certificates/radius_identity.json', radius_cert)
    stages.append(record(183, 'Sphere-identity potential reduction improves global range',
                         'anisotropic.sphere_reduce / potential_range',
                         {'identity': a.sphere_reduce(radius), 'before': a.potential_range(radius, False),
                          'after': a.potential_range(radius), 'certificate': radius_path,
                          'spectral_lower': radius_cert['lower'], 'spectral_upper': radius_cert['upper']},
                         'q=x²+y²+z² reduces to 1 and the certified mean-zero ground is exactly 3'))
    rotations = []
    for name, potential, coeffs, rotation in [
        ('rotated_x', {'1,0,0': '1'}, [0, 1], [[0, 0, -1], [0, 1, 0], [1, 0, 0]]),
        ('rotated_oblique_square', {'2,0,0': '9/25', '1,1,0': '24/25', '0,2,0': '16/25'},
         [0, 0, 1], [[0, 0, 1], ['4/5', '-3/5', 0], ['3/5', '4/5', 0]])]:
        began = time.perf_counter()
        cert = a.transfer_rotation(potential, coeffs, rotation)
        seconds = time.perf_counter()-began
        path = save(f'round09/certificates/{name}.json', cert)
        rotations.append({'case': name, 'certificate': path, 'lower': cert['lower'], 'upper': cert['upper'],
                          'width': cert['exact_width'], 'seconds': seconds,
                          'verified': a.verify_rotation(json.loads(json.dumps(cert)), expected_q=potential)})
    stages.append(record(184, 'Exact proper-rotation certificate transfer',
                         'anisotropic.transfer_rotation / verify_rotation', rotations,
                         'Two rational rotations reuse the fast axisymmetric full-sphere solver'))
    suite = unittest.defaultTestLoader.loadTestsFromName('test_anisotropic')
    stream = io.StringIO()
    tested = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    (Path(__file__).resolve().parent/'tests.log').write_text(stream.getvalue())
    validation = {'all_passed': tested.wasSuccessful() and all(r['verified'] for r in results+rotations)
                  and a.verify(radius_cert) and a.verify_cheap(cheap, expected_q=q) and baseline['rejected'],
                  'tests_run': tested.testsRun, 'failures': len(tested.failures), 'errors': len(tested.errors),
                  'certificates_replayed': len(results)+len(rotations)+2,
                  'baseline_rejected': baseline['rejected'], 'formal_assistant_checked': False}
    save('round09/RESULTS.json', {'cases': results, 'rotations': rotations, 'validation': validation})
    save('round09/VALIDATION.json', validation)
    save('round09/STAGES.json', stages)
    print(json.dumps({'validation': validation, 'spectral_cases': results}, indent=2))
    if not validation['all_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
