"""JSONL access to the ten research cycles; packaging is not a version increment.

Only registered operations run. Exact proofs stay on disk and terse replies keep
their mathematical scope. This interface makes no model calls or token claims.
"""
import hashlib
import inspect
import json
import os
from pathlib import Path
import re
import sys
import tempfile
from common import F, ROOT, exact, canonical, same, projected
import witness
import precision
import wide_potential as wide
import parameter_family as family
import global_parameter as global_search
import constraints
import ordered_spectrum as spectrum
import gn_variation as gn
import anisotropic
import sphere_scale


def verify(c):
    """Dispatch proof checking by explicit format, never by a user module name."""
    if not isinstance(c, dict) or not isinstance(c.get('format'), str):
        return False
    fmt = c['format']
    if fmt == 'reused_sphere_threshold':
        try:
            return same(c, threshold(c['source'], c['q_coefficients'], c['threshold']))
        except (ValueError, TypeError, KeyError, ArithmeticError):
            return False
    if fmt == projected.FULL_FORMAT:
        return projected.verify_full(c)
    if fmt == projected.SECTOR_FORMAT:
        return projected.verify_sector(c)
    if fmt == 'precision_full_ground_v114':
        return precision.verify_full(c)
    if fmt == 'precision_checkpoint_v112':
        return precision.verify_checkpoint(c)
    if fmt == 'same_sector_temple_v107':
        return precision.verify_temple(c)
    if fmt.startswith('precision_'):
        return precision.verify_sector(c)
    wide_checks = {wide.CHEAP_FORMAT: wide.verify_cheap_ground,
                   wide.FULL_FORMAT: wide.verify_full,
                   wide.SECTOR_FORMAT: wide.verify_sector,
                   wide.RANGE_FORMAT: wide.verify_range,
                   wide.KERNEL_FORMAT: wide.verify_kernel,
                   wide.TAYLOR_FORMAT: wide.verify_approximation,
                   wide.EXP_FORMAT: wide.verify_exponential}
    if fmt in wide_checks:
        return wide_checks[fmt](c)
    if fmt == global_search.FORMAT:
        return global_search.verify(c)
    if fmt == 'cheap_global_affine_sphere_v135':
        return global_search.verify_cheap_objective(c)
    if fmt == 'global_affine_threshold_v144':
        return global_search.verify_threshold(c)
    if fmt == anisotropic.FORMAT:
        return anisotropic.verify(c)
    if fmt == anisotropic.CHEAP_FORMAT:
        return anisotropic.verify_cheap(c)
    if fmt == anisotropic.ROTATION_FORMAT:
        return anisotropic.verify_rotation(c)
    if fmt.startswith('gn_'):
        return gn.verify(c)
    if fmt in (sphere_scale.POINT, sphere_scale.PULLBACK, sphere_scale.TRIAL,
               sphere_scale.ORDERED, sphere_scale.CELL, sphere_scale.GLOBAL,
               sphere_scale.CHECKPOINT, sphere_scale.VIOLATION, sphere_scale.JOB,
               'radius_cell_intersection_v191'):
        return sphere_scale.verify(c)
    # These dispatchers reject every format outside their own finite registry.
    return any(check(c) for check in (witness.verify, family.verify,
                                      constraints.verify, spectrum.verify))


def threshold(source, q, threshold):
    """Use a prior full-space proof for a new threshold, with no spectral solve."""
    q = wide.potential(q)
    if source.get('format') not in (projected.FULL_FORMAT,
                                   'precision_full_ground_v114', wide.FULL_FORMAT):
        raise ValueError('A full-sphere ground certificate is required')
    if not verify(source) or source['mean_zero'] is not True:
        raise ValueError('Verified full mean-zero evidence is required')
    if wide.potential(source['q_coefficients']) != q:
        raise ValueError('The reused proof belongs to a different potential')
    value = exact(threshold)
    lo, hi = F(source['lower']), F(source['upper'])
    return {'format': 'reused_sphere_threshold', 'geometry': 'unit_S2',
            'scope': 'all_real_mean_zero_H1_S2',
            'statement': 'integral(|grad u|^2+q*u^2)>=threshold*integral(u^2)',
            'q_coefficients': list(map(str, q)), 'threshold': str(value),
            'lower': str(lo), 'upper': str(hi),
            'status': ('proved' if lo >= value else
                       'refuted_by_spectral_existence' if hi < value else 'undetermined'),
            'explicit_function_included': False, 'new_spectral_solves': 0,
            'source_sha256': hashlib.sha256(canonical(source).encode()).hexdigest(),
            'source': source}


class Store:
    def __init__(self, directory=None):
        self.directory = Path(directory) if directory is not None else ROOT/'service_certificates'

    def path(self, identifier):
        if not isinstance(identifier, str) or re.fullmatch('[0-9a-f]{64}', identifier) is None:
            raise ValueError('A complete lowercase SHA-256 certificate id is required')
        path = self.directory/(identifier+'.json')
        if path.is_symlink():
            raise ValueError('Certificate symlinks are not supported')
        return path

    def put(self, certificate):
        if not verify(certificate):
            raise ValueError('Proof replay failed before storage')
        data = canonical(certificate).encode('utf-8')
        identifier = hashlib.sha256(data).hexdigest()
        path = self.path(identifier)
        self.directory.mkdir(parents=True, exist_ok=True)
        if path.exists():
            if path.read_bytes() != data:
                raise ValueError('Existing certificate content has changed')
            return identifier
        name = None
        try:
            with tempfile.NamedTemporaryFile(dir=self.directory, delete=False) as handle:
                name = handle.name
                handle.write(data)
            os.replace(name, path)
        finally:
            if name and os.path.exists(name):
                os.unlink(name)
        return identifier

    def get(self, identifier):
        data = self.path(identifier).read_bytes()
        if hashlib.sha256(data).hexdigest() != identifier:
            raise ValueError('Certificate content hash mismatch')
        c = json.loads(data)
        if not verify(c):
            raise ValueError('Stored proof no longer verifies')
        return c


def brief(c, identifier):
    # No numerical answer is stripped of its source format/scope/convention.
    result = {'certificate_id': identifier, 'format': c['format'],
              'verified': True, 'proof_kind': 'exact_software_certificate_not_formal_proof'}
    for key in ('status', 'geometry', 'scope', 'function_space', 'mean_zero',
                'q_coefficients', 'threshold', 'lower', 'upper', 'exact_width', 'width', 'gap',
                'ratio_convention', 'pi_K_lower', 'rayleigh_quotient', 'claim',
                'quantifier', 'constraint_spec', 'minimizer_intervals', 'uniqueness_claimed',
                'full_function_space_global_optimum_claimed',
                'explicit_function_included', 'new_spectral_solves'):
        if key in c:
            result[key] = c[key]
    if 'family' in c:
        result['family'] = c['family']
    if 'objective' in c:
        result['objective'] = c['objective']
    if 'problem' in c:
        result['problem'] = c['problem']
    if 'radius_request' in c:
        result['radius_request'] = c['radius_request']
    if c['format'] == spectrum.BUNDLE:
        result['scope'] = 'all_real_mean_zero_H1_S2_with_multiplicity'
        result['q_coefficients'] = c['ordered']['request']['q_coefficients']
        result['ordered'] = c['ordered']['ordered_intervals']
        result['negative_trace'] = {k: c['negative_trace'][k]
                                    for k in ('lower', 'upper', 'status')}
    if c['format'].startswith('gn_'):
        result['universal_GN_sharp_constant_proved'] = False
    if c['format'] == 'gn_refine_v171':
        result.update(status=c['stop'], scope=gn.SCOPE, ratio_convention=gn.RATIO,
                      final_pi_K_lower=c['final']['pi_K_lower'],
                      accepted_steps=sum(row['accepted'] for row in c['steps']))
    if c['format'] == 'gn_spectral_diagnostic_v174':
        result['scope'] = 'fixed_potential_necessary_condition_diagnostic_only'
    return result


def family_bound(q0, directions, box, threshold, **options):
    return family.universal_inequality(
        family.validate_family(q0, directions, box), threshold, **options)


def ordered(q, **options):
    return spectrum.bundle(spectrum.adaptive_spectrum(q, **options)['certificate'])


def radius_job(q, threshold, **options):
    return sphere_scale.research_job(q, threshold, **options)['certificate']


OPERATIONS = {
    'ground': (precision.full_ground, {'q', 'mean_zero', 'tolerance', 'modes', 'max_modes',
                                      'max_m', 'max_steps', 'use_temple'}),
    'wide_ground': (wide.full_ground, {'q', 'mean_zero', 'modes', 'max_modes', 'bits',
                                     'max_m', 'tolerance'}),
    'exponential': (wide.exponential_ground, {'a', 'degree', 'modes', 'max_modes', 'bits',
                                             'max_m', 'tolerance'}),
    'counterexample': (witness.find_counterexample, {'q', 'threshold', 'mean_zero',
                                                    'modes', 'max_steps'}),
    'family': (family_bound, {'q0', 'directions', 'box', 'threshold', 'modes', 'bits'}),
    'global_minimum': (global_search.branch_and_bound,
                       {'objective', 'max_leaves', 'tolerance', 'strategy', 'bits', 'modes'}),
    'constraints': (constraints.constrained_inequality,
                     {'q', 'threshold', 'cutoff', 'rows', 'm', 'modes', 'bits', 'max_m', 'tolerance'}),
    'ordered_spectrum': (ordered, {'q', 'k', 'mean_zero', 'tolerance', 'modes', 'max_modes',
                                 'max_m', 'max_radial', 'bits'}),
    'gn_ratio': (gn.ratio, {'raw'}),
    'gn_refine': (gn.refine, {'raw', 'iterations', 'degree_budget', 'tolerance'}),
    'gn_diagnostic': (gn.spectral_diagnostic, {'raw', 'modes', 'max_modes', 'bits', 'max_m', 'tolerance'}),
    'anisotropic': (anisotropic.certify, {'q', 'L', 'bits', 'reduce'}),
    'rotation': (anisotropic.transfer_rotation, {'q', 'coefficients', 'rotation'}),
    'radius': (radius_job, {'q', 'threshold', 'radius', 'domain', 'coordinates',
                           'tolerance', 'max_splits', 'trial', 'modes', 'max_modes',
                           'max_m', 'bits'}),
}


class Service:
    def __init__(self, directory=None):
        self.store = Store(directory)

    def call(self, request):
        if not isinstance(request, dict):
            raise ValueError('Request must be a JSON object')
        op = request.get('op')
        if op in ('fetch', 'verify', 'threshold', 'resume'):
            allowed = {'op', 'certificate_id'}
            if op == 'threshold':
                allowed |= {'q', 'threshold'}
            if op == 'resume':
                allowed |= {'q', 'tolerance', 'modes', 'max_modes', 'max_m', 'max_steps', 'use_temple'}
            if set(request)-allowed:
                raise ValueError('Unknown fields for '+op)
            c = self.store.get(request['certificate_id'])
            if op == 'threshold':
                c = threshold(c, request['q'], request['threshold'])
            elif op == 'resume':
                checkpoint = precision.checkpoint(c)
                options = {k: v for k, v in request.items() if k not in ('op', 'certificate_id')}
                c = precision.full_ground(saved=checkpoint, **options)
            identifier = self.store.put(c)
            response = brief(c, identifier)
            if op == 'fetch':
                response['certificate'] = c
            return response
        if not isinstance(op, str) or op not in OPERATIONS:
            raise ValueError('Unknown operation')
        function, fields = OPERATIONS[op]
        if set(request)-fields-{'op'}:
            raise ValueError('Unknown fields for '+op)
        arguments = {k: v for k, v in request.items() if k != 'op'}
        inspect.signature(function).bind(**arguments)
        c = function(**arguments)
        return brief(c, self.store.put(c))


def main():
    service = Service()
    for line in sys.stdin:
        try:
            response = service.call(json.loads(line))
        except (ValueError, KeyError, TypeError, ArithmeticError, OSError, RecursionError) as error:
            response = {'status': 'rejected_or_computation_failed',
                        'exception': type(error).__name__, 'message': str(error)}
        print(canonical(response), flush=True)


if __name__ == '__main__':
    main()
