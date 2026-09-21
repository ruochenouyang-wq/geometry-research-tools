"""Read-only dependencies and exact serialization for the failed-case retests."""
from pathlib import Path
import hashlib
import json
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
PREVIOUS = ROOT.parent/'geometry_cycles_v95_v194'
sys.path.insert(0,str(PREVIOUS))
import common as previous_common
from common import F, exact, same, canonical, digest, projected, sectors, base, arithmetic, exact_extrema
import constraints as old_constraints
import ordered_spectrum as old_spectrum
import gn_variation as old_gn
import anisotropic as old_anisotropic
import portal as old_portal
sys.path.insert(0,str(ROOT))


def save(relative, value):
    path = ROOT/relative
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,sort_keys=True,allow_nan=False)+'\n')
    return str(path.relative_to(ROOT))


def source_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def original_case(identifier):
    book=json.loads((PREVIOUS/'challenges.json').read_text())
    return next(case for case in book['cases'] if case['id']==identifier)
