"""Shared exact arithmetic and read-only imports for the ten research cycles."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
V94 = ROOT.parent/'geometry_v91_v94'
V90 = ROOT.parent/'geometry_v36_v90'
sys.path.insert(0, str(V94))
import projected_spectrum as projected
import sphere_sectors as sectors
import spectral_certifier as base
import error_bounds as arithmetic
import exact_extrema


def exact(value):
    if type(value) is F:
        return value
    return base.rational(value)


def rational_vector(raw, maximum=256):
    if not isinstance(raw, list) or not 1 <= len(raw) <= maximum:
        raise ValueError('Invalid rational vector')
    return [exact(x) for x in raw]


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=False, allow_nan=False)


def same(a, b):
    return canonical(a) == canonical(b)


def save(relative, value):
    path = ROOT/relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False)+'\n')
    return str(path.relative_to(ROOT))


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def qnorm(q):
    return sum((abs(exact(x)) for x in q), F(0))


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def quadratic(a, v):
    return sum((v[i]*a[i][j]*v[j] for i in range(len(v)) for j in range(len(v))), F(0))
