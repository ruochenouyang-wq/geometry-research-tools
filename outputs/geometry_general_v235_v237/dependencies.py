"""Read-only access to the frozen exact geometry engines."""
from pathlib import Path
import hashlib
import json
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
PREVIOUS = ROOT.parent/'geometry_retests_v195_v234'
sys.path.insert(0, str(PREVIOUS))
import support
from support import F, exact, same, canonical, digest, arithmetic, base, exact_extrema
from support import old_anisotropic as oldmath
import global_constraints as constraints
import parity_spectrum as parity
import wide_potential as wide
sys.path.insert(0, str(ROOT))


def save(relative, value):
    path = ROOT/relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False)+'\n')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()
