"""Read-only bridge to the frozen mathematical baseline."""
from pathlib import Path
import hashlib
import importlib
import json
import sys

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
FROZEN = ROOT.parent / 'geometry_precision_followup'
sys.dont_write_bytecode = True
if str(FROZEN) not in sys.path:
    sys.path.insert(0, str(FROZEN))
legacy = importlib.import_module('precision_tool')
direct = legacy.direct
pieces = legacy.pieces
enriched = legacy.enriched
from fractions import Fraction as F


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def source_hashes():
    names = ('geometry_precision_followup', 'geometry_general_v235_v237',
             'geometry_generalization_notes', 'geometry_retests_v195_v234',
             'geometry_cycles_v95_v194', 'geometry_v91_v94', 'geometry_v36_v90')
    result = {}
    for name in names:
        for path in sorted((ROOT.parent/name).rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts and path.suffix != '.pyc':
                result[str(path.relative_to(PROJECT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    for path in [PROJECT/'README.md', PROJECT/'.gitignore',
                 ROOT.parent/'geometry_research_strategy/研究指导报告.md']:
        if path.exists():
            result[str(path.relative_to(PROJECT))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result
