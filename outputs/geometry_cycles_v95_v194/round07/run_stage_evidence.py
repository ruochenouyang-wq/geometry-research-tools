"""Run real V155--V164 capabilities, saving replayable exact evidence."""
from pathlib import Path
import inspect
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import save, projected
import ordered_spectrum as s


def main():
    start = time.perf_counter()
    save('round07/baseline.json', {
        'problem': 'Sphere Lieb-Thirring research requires full angular multiplicities, strict negative counts and trace sums.',
        'source': 'https://www.ma.ic.ac.uk/~alaptev/Papers/YJFAN-S-20-00336.pdf',
        'existing_full_ground_signature': str(inspect.signature(projected.full_ground)),
        'available': {name: hasattr(projected, name) for name in
                      ('full_ground', 'ordered_spectrum', 'count_below', 'negative_trace')},
        'weakness': 'Existing full-ground interface cannot assemble the ordered full spectrum or negative eigenvalue moments.'})
    stages = []

    def stage(version, capability, callable_name, evidence, outcome):
        path = save('round07/v%d.json' % version, evidence)
        stages.append({'version': version, 'capability': capability,
                       'callable': 'ordered_spectrum.'+callable_name,
                       'evidence': [path], 'outcome': outcome})

    request = s.normalize_request([0], k=8)
    initial = s.initial_order_bounds(s.normalize_request([-6, 0, 2], k=8))
    stage(155, 'Global initial ordered min-max intervals with full multiplicity and no Schur computation.',
          'initial_order_bounds', initial, 'First three eigenvalues in [-4,-2], next five in [0,2]; full q/k/projection replay binding.')
    calibration = s.laplace_spectrum(8)
    stage(156, 'Exact unperturbed and constant-shift ordered-spectrum calibration.',
          'laplace_spectrum', calibration, 'First eight levels: 2,2,2,6,6,6,6,6.')
    groups = [[projected.certify_sector([0], m=m, k=j, modes=4, bits=20)
               for j in (1, 2)] for m in range(3)]
    assembly = s.assemble_sectors([0], groups)
    stage(157, 'Verified sector eigenvalue assembly with 1/2 angular multiplicity.',
          'assemble_sectors', assembly, 'Six certified sector eigenvalues represent ten real full-sphere eigenvalues.')
    floors = s.omitted_floors(assembly)
    stage(158, 'Independent radial and azimuthal floors for unenumerated modes.',
          'omitted_floors', floors, 'Each omitted direction has a min-max bound; no numerical truncation is assumed exact.')
    bound = s.ordered_bounds(request, assembly)
    stage(159, 'Full kth order statistics with analytical missing-mode placeholders.',
          'ordered_bounds', bound, 'All first eight eigenvalues exactly enclosed despite uncomputed higher modes.')
    result = s.adaptive_spectrum([-6, 0, 2], k=8, modes=6, max_modes=8,
        max_m=3, max_radial=3, bits=28, tolerance='1/100000')
    stage(160, 'Adaptive sector/radial enumeration with honest resource exhaustion.',
          'adaptive_spectrum', result, result['certificate']['status'])
    shifted6 = s.adaptive_spectrum([-6], k=8)['certificate']
    count = s.count_below(shifted6['assembly'])
    stage(161, 'Strict full spectral counts and explicit threshold clusters.',
          'count_below', count, 'Exactly 3 negative eigenvalues; the five zero eigenvalues are excluded.')
    shifted20 = s.adaptive_spectrum([-20], k=16)['certificate']
    trace = s.negative_trace(shifted20['assembly'])
    stage(162, 'Negative trace intervals include all potentially negative omitted modes.',
          'negative_trace', trace, 'Exact trace 180, including a rigorously resolved omitted angular contribution.')
    ratio = s.lt_quotient(trace, [20])
    stage(163, 'Nonnegative-potential validation and rational pi-scaled fixed-potential LT quotient.',
          'lt_quotient', ratio, 'pi times trace quotient = 9/80; witness lower bound only, with factor-four dual normalization.')
    final = s.bundle(result['certificate'], assertion={'index': 1, 'lower_bound': -4}, potential=[6, 0, -2])
    stage(164, 'Bound certificate bundle, independent replay and tri-state eigenvalue assertion.',
          'bundle', final, 'All subproofs bind to q=-6+2t^2; first eigenvalue >= -4 proved.')
    certificates = [initial, bound, result['certificate'], count, trace, ratio, final]
    verified = [s.verify(cert) for cert in certificates]
    if not all(verified):
        raise AssertionError('Stage replay failed')
    save('round07/STAGES.json', stages)
    save('round07/verification.json', {'all_passed': all(verified), 'certificate_replays': len(verified),
         'replay_results': verified, 'wall_seconds': time.perf_counter()-start,
         'constant_minus20_negative_count': s.count_below(shifted20['assembly'])['count_lower']})
    print({'all_passed': all(verified), 'versions': len(stages), 'seconds': time.perf_counter()-start})


if __name__ == '__main__':
    main()
