"""Exact, assumption-checked reuse of unit-sphere mathematical certificates."""
from copy import deepcopy
from fractions import Fraction as F
import backend

FORMAT = 'geometry_proof_transport_v1'
MEASURE = 'd_sigma/(4*pi)_on_unit_S2'


def rational(value):
    return backend.pieces.rational(value)


def function(value):
    return backend.direct.normalize(value)


def _scope(mean_zero):
    return 'all_real_mean_zero_H1_on_unit_S2' if mean_zero else 'all_real_H1_on_unit_S2'


def _spectral_claim(q, mean_zero, lower, upper):
    if type(mean_zero) is not bool:
        raise ValueError('mean_zero must be a boolean')
    lower, upper = rational(lower), rational(upper)
    if lower > upper:
        raise ValueError('Contradictory spectral intervals')
    return {'kind': 'spectrum', 'function': function(q), 'geometry': 'unit_S2',
            'measure': MEASURE, 'scope': _scope(mean_zero), 'mean_zero': mean_zero,
            'eigenvalue_index': 1, 'lower': str(lower), 'upper': str(upper),
            'exact_width': str(upper-lower), 'full_infinite_space_covered': True,
            'formal_proof_assistant_checked': False}


def _source_claim(source, verifier=None):
    if type(source) is not dict:
        raise ValueError('Source must be a mathematical certificate')
    if source.get('format') == FORMAT:
        valid = verify(source, verifier=verifier)
    elif verifier is not None:
        valid = verifier(source) is True
    else:
        valid = backend.legacy.handle({'op': 'verify', 'certificate': source}).get('verified') is True
    if not valid:
        raise ValueError('Source mathematical evidence did not verify')
    if (source.get('geometry') == 'unit_S2' and
            source.get('full_infinite_space_covered') is True and
            type(source.get('eigenvalue_index')) is int and source['eigenvalue_index'] == 1 and
            source.get('measure') == MEASURE):
        claim = _spectral_claim(source['function'], source['mean_zero'], source['lower'], source['upper'])
        if source.get('scope') != claim['scope']:
            raise ValueError('Source spectral space mismatch')
        return claim
    raise ValueError('Unsupported source conclusion')


def _one_source(sources, verifier):
    if type(sources) is not list or len(sources) != 1:
        raise ValueError('Exactly one source is required')
    return _source_claim(sources[0], verifier)


def _derive(operation, sources, params, verifier=None):
    if operation == 'constant_spectrum':
        if sources or set(params) != {'function', 'mean_zero'}:
            raise ValueError('Invalid constant-spectrum premises')
        q = function(params['function'])
        if rational(q['amplitude']) != 0:
            raise ValueError('A constant potential requires zero amplitude')
        c = rational(q['offset'])
        value = c + (2 if params['mean_zero'] is True else 0)
        return _spectral_claim(q, params['mean_zero'], value, value)
    if operation == 'spectral_axis':
        if set(params) != {'axis'}:
            raise ValueError('Axis rotation requires exactly an axis')
        old = _one_source(sources, verifier)
        if old['kind'] != 'spectrum':
            raise ValueError('Spectral transport requires spectral evidence')
        q = dict(old['function'], axis=params['axis'])
        return _spectral_claim(q, old['mean_zero'], old['lower'], old['upper'])
    raise ValueError('Unsupported proof transport')


def _make(operation, sources, params, verifier=None):
    claim = _derive(operation, sources, params, verifier)
    return {'format': FORMAT, 'operation': operation, 'sources': deepcopy(sources),
            'parameters': deepcopy(params), **claim}


def constant_spectrum(value=0, mean_zero=True, axis=2):
    q = function({'kind': 'axis_profile', 'profile': 'step', 'axis': axis,
                  'amplitude': '0', 'offset': str(rational(value))})
    return _make('constant_spectrum', [], {'function': q, 'mean_zero': mean_zero})


def verify(certificate, verifier=None):
    try:
        if type(certificate) is not dict or certificate.get('format') != FORMAT:
            return False
        rebuilt = _make(certificate['operation'], certificate['sources'],
                        certificate['parameters'], verifier)
        return backend.canonical(certificate) == backend.canonical(rebuilt)
    except (ValueError, TypeError, KeyError, ArithmeticError, IndexError):
        return False


def spectral_axis(source, axis, verifier=None):
    return _make('spectral_axis', [source], {'axis': axis}, verifier)
