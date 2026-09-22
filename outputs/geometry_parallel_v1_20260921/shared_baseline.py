"""Shared read-only baseline loader for the three independent experiments."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
BASELINE_DIR = Path(__file__).resolve().parent.parent / "geometry_accuracy80_20260921"
if str(BASELINE_DIR) not in sys.path:
    sys.path.insert(0, str(BASELINE_DIR))
import solver as baseline

def solve(task):
    return baseline.solve(task)

def verify(certificate, task):
    return baseline.assess(certificate, baseline.normalize_task(task))
