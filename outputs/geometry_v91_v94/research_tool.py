"""Small JSONL interface for certified projected sphere spectra.

No model calls, expression evaluation, arbitrary-path reads, or token claims.
JSON clients can access certificates only by a complete content SHA-256 ID.
"""
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import tempfile

sys.dont_write_bytecode = True
import projected_spectrum as spectrum


CERTIFICATE_DIRECTORY = Path(__file__).resolve().parent / 'certificates' / 'service'
FULL_OPTIONS = {'modes', 'max_modes', 'bits', 'max_m', 'tolerance'}
FIELDS = {
    'sector': {'q', 'm', 'mean_zero', 'k', 'modes', 'bits'},
    'mean_zero_ground': {'q'} | FULL_OPTIONS,
    'weighted_poincare': {'q', 'threshold'} | FULL_OPTIONS,
    'verify': {'certificate_id'},
    'fetch': {'certificate_id'},
}


def wire(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False)


def exact(value, name):
    if type(value) is Fraction:
        return value
    if type(value) is int:
        value = str(value)
    if (not isinstance(value, str) or len(value) > 2000
            or re.fullmatch(r'-?\d+(?:/[1-9]\d*)?', value) is None):
        raise ValueError(name + ' must be an integer or exact rational string, never a float or boolean')
    return Fraction(value)


def integer(value, name, minimum, maximum):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError('%s must be an integer in %d..%d' % (name, minimum, maximum))
    return value


def q_input(raw):
    if not isinstance(raw, list) or not 1 <= len(raw) <= 7:
        raise ValueError('q must contain 1..7 coefficients in ascending powers')
    return spectrum.potential([str(exact(value, 'q coefficient')) for value in raw])


def full_options(request):
    options = {
        'modes': integer(request.get('modes', 8), 'modes', 1, 64),
        'max_modes': integer(request.get('max_modes', 32), 'max_modes', 1, 64),
        'bits': integer(request.get('bits', 44), 'bits', 8, 160),
        'max_m': integer(request.get('max_m', 16), 'max_m', 0, 64),
        'tolerance': exact(request.get('tolerance', Fraction(1, 10**10)), 'tolerance'),
    }
    if options['max_modes'] < options['modes']:
        raise ValueError('max_modes must be at least modes')
    if not Fraction(1, 10**30) <= options['tolerance'] <= 1:
        raise ValueError('tolerance must lie in [1/10^30,1]')
    return options


def verify_certificate(certificate):
    if not isinstance(certificate, dict):
        return False
    kind = certificate.get('format')
    if kind == spectrum.SECTOR_FORMAT:
        return spectrum.verify_sector(certificate, independent=True)
    if kind == spectrum.FULL_FORMAT:
        return spectrum.verify_full(certificate, independent=True)
    if kind == spectrum.INEQUALITY_FORMAT:
        return spectrum.verify_inequality(certificate)
    return False


def full_request_matches(certificate, q, options):
    """Bind recorded scope, tolerance and truncations to this exact request."""
    return (isinstance(certificate, dict)
            and certificate.get('format') == spectrum.FULL_FORMAT
            and spectrum.verify_full(certificate, expected_q=q, expected_mean_zero=True)
            and certificate['tolerance'] == str(options['tolerance'])
            and all(sector['azimuth_m'] <= options['max_m']
                    and options['modes'] <= sector['modes'] <= options['max_modes']
                    for sector in certificate['sectors']))


class CertificateStore:
    """An internal directory argument allows isolated tests; JSON has no path option."""
    def __init__(self, directory=None):
        self.directory = Path(directory) if directory is not None else CERTIFICATE_DIRECTORY
        self.directory.mkdir(parents=True, exist_ok=True)

    def path(self, identifier):
        if not isinstance(identifier, str) or re.fullmatch(r'[0-9a-f]{64}', identifier) is None:
            raise ValueError('certificate_id must be a complete lowercase SHA-256 hex ID')
        path = self.directory / (identifier + '.json')
        if path.is_symlink():
            raise ValueError('Certificate symlinks are not supported')
        return path

    def put(self, certificate):
        if not verify_certificate(certificate):
            raise ValueError('Unverified certificate rejected')
        data = wire(certificate).encode('utf-8')
        identifier = hashlib.sha256(data).hexdigest()
        path = self.path(identifier)
        if path.exists():
            if path.read_bytes() != data:
                raise ValueError('Existing certificate content has changed')
            return identifier
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=self.directory, delete=False) as handle:
                temporary = handle.name
                handle.write(data)
            os.replace(temporary, path)
        finally:
            if temporary is not None and os.path.exists(temporary):
                os.unlink(temporary)
        return identifier

    def get(self, identifier):
        data = self.path(identifier).read_bytes()
        if hashlib.sha256(data).hexdigest() != identifier:
            raise ValueError('Certificate hash mismatch')
        certificate = json.loads(data)
        if not verify_certificate(certificate):
            raise ValueError('Certificate proof verification failed')
        return certificate


def summary(certificate, identifier):
    kind = certificate['format']
    inequality = kind == spectrum.INEQUALITY_FORMAT
    spectral = certificate['evidence'] if inequality else certificate
    if kind == spectrum.SECTOR_FORMAT:
        scope = {'geometry': 'unit_S2', 'space': 'single_azimuth_sector_only',
                 'q_coefficients': certificate['q_coefficients'],
                 'm': certificate['azimuth_m'], 'mean_zero': certificate['mean_zero'],
                 'k': certificate['eigenvalue_index'], 'projection': certificate['projection']}
        status = 'certified_sector_bound'
    else:
        scope = {'geometry': 'unit_S2', 'space': spectral['scope'],
                 'q_coefficients': spectral['q_coefficients'], 'mean_zero': spectral['mean_zero'],
                 'k': spectral['eigenvalue_index_in_constrained_space']}
        status = certificate['status']
        if inequality:
            scope.update(threshold=certificate['threshold'], quantifier=certificate['quantifier'],
                         bound_quantity='ground_energy_in_the_full_mean_zero_space')
    result = {'status': status, 'scope': scope, 'lower': spectral['lower'],
              'upper': spectral['upper'], 'width': spectral['exact_width'],
              'certificate_id': identifier}
    if inequality:
        result['refutation_does_not_include_an_explicit_function'] = certificate['refutation_does_not_include_an_explicit_function']
    return result


class Service:
    def __init__(self, directory=None):
        self.store = CertificateStore(directory)

    def call(self, request):
        if not isinstance(request, dict):
            raise ValueError('Request must be a JSON object')
        operation = request.get('op')
        if not isinstance(operation, str) or operation not in FIELDS:
            raise ValueError('Unknown operation')
        extra = set(request) - FIELDS[operation] - {'op'}
        if extra:
            raise ValueError('Unknown fields: ' + ', '.join(sorted(extra)))
        if operation in ('verify', 'fetch'):
            identifier = request['certificate_id']
            certificate = self.store.get(identifier)
            result = summary(certificate, identifier)
            result['verified'] = True
            if operation == 'fetch':
                result['certificate'] = certificate
            return result
        q = q_input(request['q'])
        if operation == 'sector':
            m = integer(request.get('m', 0), 'm', 0, 64)
            mean_zero = request.get('mean_zero', True)
            if type(mean_zero) is not bool:
                raise ValueError('mean_zero must be a boolean')
            modes = integer(request.get('modes', 12), 'modes', 1, 64)
            k = integer(request.get('k', 1), 'k', 1, modes)
            bits = integer(request.get('bits', 44), 'bits', 8, 160)
            certificate = spectrum.certify_sector(q, m=m, mean_zero=mean_zero, k=k, modes=modes, bits=bits)
            valid = (certificate.get('format') == spectrum.SECTOR_FORMAT
                     and spectrum.verify_sector(certificate, expected_q=q, expected_m=m,
                                                expected_mean_zero=mean_zero, expected_k=k)
                     and certificate['modes'] == modes)
        elif operation == 'mean_zero_ground':
            options = full_options(request)
            certificate = spectrum.full_ground(q, mean_zero=True, **options)
            valid = full_request_matches(certificate, q, options)
        else:
            threshold = exact(request['threshold'], 'threshold')
            options = full_options(request)
            certificate = spectrum.weighted_poincare(q, threshold, **options)
            valid = (certificate.get('format') == spectrum.INEQUALITY_FORMAT
                     and spectrum.verify_inequality(certificate, expected_q=q, expected_threshold=threshold)
                     and full_request_matches(certificate['evidence'], q, options))
        if not valid:
            raise ValueError('Generated proof does not bind the requested problem, scope, tolerance and truncation sizes')
        identifier = self.store.put(certificate)
        return summary(certificate, identifier)


def error_response(exception):
    status = 'computation_failed' if isinstance(exception, ArithmeticError) else 'invalid_request'
    return {'status': status, 'error': type(exception).__name__,
            'message': str(exception)[:2000], 'mathematical_verdict': None}


def serve(input_stream, output_stream, service=None):
    service = Service() if service is None else service
    for line in input_stream:
        if not line.strip():
            continue
        try:
            request = json.loads(line)
            result = service.call(request)
        except (ValueError, TypeError, KeyError, IndexError, OSError, ArithmeticError,
                AttributeError, RuntimeError, RecursionError) as exception:
            result = error_response(exception)
        output_stream.write(wire(result) + '\n')
        output_stream.flush()


if __name__ == '__main__':
    if len(sys.argv) != 1:
        sys.stdout.write(wire(error_response(ValueError('No command-line paths or options; use JSONL requests'))) + '\n')
        raise SystemExit(2)
    serve(sys.stdin, sys.stdout)
