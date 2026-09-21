"""Reproduce baseline, ten executable stages, and round 6 acceptance evidence."""
from datetime import datetime, timezone
from fractions import Fraction as F
from pathlib import Path
import io
import json
import sys
import tempfile
import time
import unittest

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from common import save, digest
import constraints as c
import research_tool as old
import test_constraints


def safe(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, list):
        return [safe(x) for x in value]
    if isinstance(value, dict):
        return {k: safe(v) for k, v in value.items()}
    return value


def main():
    request = {'op': 'mean_zero_ground', 'q': [0], 'exclude_harmonic_degrees_leq': 1}
    with tempfile.TemporaryDirectory() as directory:
        try:
            response = old.Service(directory).call(request)
        except (ValueError, TypeError) as error:
            response = old.error_response(error)
    save('round06/baseline.json', {
        'request': request, 'observed_response': response,
        'external_problem': 'sphere interpolation with low spherical harmonic orthogonality',
        'exact_calibration': 'q=0; orthogonality to all degrees <=1 gives sharp value 6',
        'interpretation': 'V94 correctly rejects the unsupported additional constraint; its ordinary mean-zero value 2 is not a wrong answer to its own different problem',
        'timestamp_utc': datetime.now(timezone.utc).isoformat()})
    stages, certificates, runs = [], [], []

    def stage(version, capability, call, thunk, outcome):
        start = time.perf_counter()
        result = safe(thunk())
        seconds = time.perf_counter()-start
        path = save('round06/evidence/v%d.json' % version,
                    {'version': version, 'callable': call, 'wall_seconds': seconds,
                     'output': result, 'outcome': outcome})
        stages.append({'version': version, 'capability': capability, 'callable': call,
                       'evidence': [path], 'outcome': outcome})
        return result

    def cert_case(name, thunk):
        start = time.perf_counter()
        result = thunk()
        replay = json.loads(json.dumps(result))
        accepted = c.verify(replay)
        if not accepted:
            raise AssertionError('Saved proof replay failed: '+name)
        path = save('round06/certificates/'+name+'.json', result)
        certificates.append(path)
        spectral = result.get('evidence', result)
        runs.append({'case': name, 'certificate': path, 'verified': accepted,
                     'wall_seconds_including_replay': time.perf_counter()-start,
                     'lower': spectral['lower'], 'upper': spectral['upper'],
                     'decimal_enclosure': spectral['decimal_enclosure'],
                     'exact_width': spectral['exact_width'],
                     'scope': spectral['scope']})
        return result

    stage(145, 'Analytic full-sphere harmonic-exclusion lower bound and actual harmonic Ritz upper',
          'constraints.analytic_exclusion_bound', lambda: cert_case('analytic_q2t2_cutoff1',
              lambda: c.analytic_exclusion_bound([0, 0, 2], 1)),
          'For q=2t² and cutoff 1 the full-space analytic enclosure is [6,148/21], with an actual P2 Ritz witness; q=0 calibrates exactly to 6 and 12.')
    stage(146, 'Assemble constrained sectors with degree start max(m,cutoff+1)',
          'constraints.sector_assembly',
          lambda: [c.sector_assembly([0, 0, 2], m, 1, 4) for m in (0, 1, 3)],
          'Starts are 2,2,3; projection follows complete multiplication by q.')
    stage(147, 'Constraint-aware lower bounds for radial and omitted angular tails',
          'constraints.tail_floors', lambda: c.tail_floors([0], 2, 0, 4, 1),
          'The omitted m>=1 floor is 12 under cutoff 2, and radial tail floor is 56.')
    sector = stage(148, 'Exact Schur eigenvalue certificate for a low-harmonic constrained sector',
          'constraints.certify_exclusion_sector',
          lambda: cert_case('sector_q2t2_cutoff1_m0', lambda: c.certify_exclusion_sector(
              [0, 0, 2], m=0, cutoff=1, modes=8, bits=44)),
          'Radial tail included; dense exact inertia independently replays the endpoints.')
    full = stage(149, 'Full sphere aggregation with complete harmonic exclusion',
          'constraints.full_exclusion_ground',
          lambda: cert_case('full_q2t2_cutoff1', lambda: c.full_exclusion_ground(
              [0, 0, 2], cutoff=1, modes=8, bits=44, tolerance=F(1, 10**10))),
          'All nonaxisymmetric sectors and infinite angular tail are controlled.')
    stage(150, 'Infinite-radial-space analytic variational bracket for arbitrary finite moments',
          'constraints.analytic_moment_bound', lambda: cert_case('analytic_moment_q2t2_gP0P1P2',
              lambda: c.analytic_moment_bound([0, 0, 2], [[1, 1, 1]], 0, 6)),
          'A codimension bound leaves a nonzero trial vector in the first two harmonic modes, giving the rigorous pre-Schur bracket [0,4].')
    stage(151, 'Mass-weighted exact constraint nullspace with redundancy removal',
          'constraints.mass_nullspace', lambda: c.mass_nullspace(
              [[1, 1, 1], [2, 2, 2], [0]], [F(2), F(2, 3), F(2, 5), F(2, 7)]),
          'Duplicate/scaled rows have rank 1; free coordinates use exact L2 mass weighting.')
    stage(152, 'Generalized reduced energy, non-diagonal mass, and transformed tail couplings',
          'constraints.constrained_forms', lambda: c.constrained_forms([0, 0, 2], 0, 6, [[1, 1, 1]]),
          'Reduced mass contains off-diagonal 2/15; all entries participate in A-xM and the Schur correction.')
    stage(153, 'Infinite-dimensional arbitrary finite-moment single-sector spectrum',
          'constraints.certify_moment_sector',
          lambda: cert_case('moment_q2t2_gP0P1P2', lambda: c.certify_moment_sector(
              [0, 0, 2], [[1, 1, 1]], modes=8, bits=44)),
          'The moment-constrained quotient is certified with its full radial tail and single-sector scope.')
    stage(154, 'Three-way constrained weighted inequality with bound quantifier',
          'constraints.constrained_inequality',
          lambda: cert_case('inequality_cutoff1_q0_threshold6', lambda: c.constrained_inequality(
              [0], 6, cutoff=1)),
          'The sharp q=0 threshold 6 is proved for all functions orthogonal to degrees 0 and 1.')

    for name, q, cutoff in [('full_q0_cutoff1', [0], 1), ('full_q0_cutoff2', [0], 2),
                            ('full_shift_cutoff2', ['-7/3'], 2)]:
        cert_case(name, lambda q=q, cutoff=cutoff: c.full_exclusion_ground(q, cutoff))
    for cutoff in (1, 2):
        cert_case('analytic_q0_cutoff%d' % cutoff,
                  lambda cutoff=cutoff: c.analytic_exclusion_bound([0], cutoff))
    for name, rows, k in [('moment_q0_gP0P2_k2', [[1, 0, 1]], 2),
                         ('moment_q0_removeP0P1', [[1], [0, 1]], 1),
                         ('moment_q0_scaled_redundant', [[1, 1, 1], [-3, -3, -3]], 1)]:
        cert_case(name, lambda rows=rows, k=k: c.certify_moment_sector([0], rows, k=k, modes=6, bits=44))
    cert_case('inequality_cutoff1_q2t2_threshold25over4',
              lambda: c.inequality_certificate(full, F(25, 4)))
    cert_case('inequality_cutoff1_q2t2_threshold63over10',
              lambda: c.inequality_certificate(full, F(63, 10)))
    cert_case('inequality_cutoff1_q2t2_midpoint', lambda: c.inequality_certificate(
              full, (F(full['lower'])+F(full['upper']))/2))
    log = io.StringIO()
    suite = unittest.defaultTestLoader.loadTestsFromModule(test_constraints)
    result = unittest.TextTestRunner(stream=log, verbosity=2).run(suite)
    (ROOT/'round06/tests.log').write_text(log.getvalue())
    validation = {'all_passed': result.wasSuccessful() and all(x['verified'] for x in runs),
                  'tests_run': result.testsRun, 'failures': len(result.failures),
                  'errors': len(result.errors), 'certificates_replayed': len(certificates),
                  'baseline_rejects_unsupported_constraint': response.get('status') == 'invalid_request',
                  'no_prior_release_modified': 'Agent edits are restricted to constraints.py, test_constraints.py and round06; root performs complete frozen manifest audit',
                  'proof_engine': 'exact rational arithmetic, search banded inertia, replay dense inertia',
                  'formal_assistant_checked': False}
    save('round06/STAGES.json', stages)
    save('round06/RESULTS.json', {'cases': runs, 'certificates': certificates,
                                'validation': validation})
    save('round06/VALIDATION.json', validation)
    print(json.dumps(validation, indent=2))
    if not validation['all_passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
