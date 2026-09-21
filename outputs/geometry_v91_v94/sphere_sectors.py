"""Exact rational forms for one azimuthal sector on the unit two-sphere.

The basis consists of the real associated Legendre functions P_l^m, with
degrees ``start <= l < start + modes``.  Integration here is over t in [-1,1];
the common azimuthal normalization cancels from the generalized eigenproblem.
``start=1,m=0`` is the compression to zero-mean axisymmetric functions.  For a
nonconstant potential this is a projected operator, not an invariant sector.

Only finite forms are assembled here: no claim about the omitted spectrum is
made.  All intermediate multiplication degrees are retained, including degrees
below ``start``; projection takes place after the whole polynomial is applied.
Python standard library only.  The frozen V90 potential validator is reused.
"""
from fractions import Fraction as F
import importlib.util
from math import factorial
from pathlib import Path
import sys


MAX_M = 64
MAX_START = 128
MAX_MODES = 64
MAX_DEGREE = 256
MAX_Q_BASIS_DEGREE = MAX_DEGREE - 6


def _load_legacy():
    path = Path(__file__).resolve().parents[1] / "geometry_v36_v90" / "spectral_certifier.py"
    spec = importlib.util.spec_from_file_location("_sphere_sectors_v90_validator", path)
    module = importlib.util.module_from_spec(spec)
    # Importing a frozen module must not create or refresh its bytecode files.
    previous = sys.dont_write_bytecode
    try:
        sys.dont_write_bytecode = True
        spec.loader.exec_module(module)
    finally:
        sys.dont_write_bytecode = previous
    return module


_legacy = _load_legacy()


def potential(q):
    """Validate degree <= 6 and V90 coefficient limits; reject floats/bools."""
    if not isinstance(q, (list, tuple)) or not 1 <= len(q) <= 7:
        raise ValueError("q must have 1..7 rational coefficients in ascending powers")
    if any(type(c) not in (int, str, F) for c in q):
        raise ValueError("q coefficients must be integers, Fractions or rational strings")
    return _legacy.potential([str(c) for c in q])


def _check_basis(m, l, limit=MAX_DEGREE):
    if type(m) is not int or not 0 <= m <= MAX_M:
        raise ValueError("m must be an integer in 0..64")
    if type(l) is not int or not m <= l <= limit:
        raise ValueError("basis degree must be an integer between m and %d" % limit)


def basis_mass(m, l):
    """Return integral (P_l^m)^2 dt = 2(l+m)!/((2l+1)(l-m)!)."""
    _check_basis(m, l)
    return F(2 * factorial(l + m), (2 * l + 1) * factorial(l - m))


def _times_t(series, m):
    result = {}
    for l, coefficient in series.items():
        result[l + 1] = result.get(l + 1, F(0)) + coefficient * F(l - m + 1, 2 * l + 1)
        if l > m:
            result[l - 1] = result.get(l - 1, F(0)) + coefficient * F(l + m, 2 * l + 1)
    return {l: coefficient for l, coefficient in result.items() if coefficient}


def _q_times_basis(q, m, l):
    result, term = {}, {l: F(1)}
    for power, coefficient in enumerate(q):
        if coefficient:
            for degree, value in term.items():
                result[degree] = result.get(degree, F(0)) + coefficient * value
        if power + 1 < len(q):
            term = _times_t(term, m)
    return {degree: coefficient for degree, coefficient in result.items() if coefficient}


def q_times_basis(q, m, l):
    """Return the complete associated-Legendre expansion of q(t) P_l^m(t).

    No window projection is applied.  Input l <= 250 guarantees every output
    degree is <= 256, so every output has a supported ``basis_mass``.
    """
    _check_basis(m, l, MAX_Q_BASIS_DEGREE)
    return _q_times_basis(potential(q), m, l)


def assemble(q, m, start, modes):
    """Return ``A`` (form matrix), ``M`` (mass diagonal), and ``degrees``.

    A represents -Delta + q in the specified window, with exactly rational
    entries.  Require 0 <= m <= 64, m <= start <= 128 and 1 <= modes <= 64.
    """
    _check_basis(m, start, MAX_START)
    if type(modes) is not int or not 1 <= modes <= MAX_MODES:
        raise ValueError("modes must be an integer in 1..64")
    q = potential(q)
    degrees = list(range(start, start + modes))
    mass = [basis_mass(m, l) for l in degrees]
    matrix = [[F(0) for _ in degrees] for _ in degrees]
    for j, l in enumerate(degrees):
        for degree, coefficient in _q_times_basis(q, m, l).items():
            if start <= degree < start + modes:
                matrix[degree - start][j] += mass[degree - start] * coefficient
        matrix[j][j] += l * (l + 1) * mass[j]
    if any(matrix[i][j] != matrix[j][i] for i in range(modes) for j in range(i)):
        raise ArithmeticError("associated Legendre form is not symmetric")
    return {"A": matrix, "M": mass, "degrees": degrees}
