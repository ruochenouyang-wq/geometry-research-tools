"""Exact, assumption-checked reuse of unit-sphere mathematical certificates."""
from copy import deepcopy
from fractions import Fraction as F
import backend

FORMAT = 'geometry_proof_transport_v1'
MEASURE = 'probability_surface_measure'


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
