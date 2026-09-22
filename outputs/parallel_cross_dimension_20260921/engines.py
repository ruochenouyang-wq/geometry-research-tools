"""Evaluation-only adapters for the frozen parallel mathematical prototypes.

No new mathematical search is introduced here. A and B use their specialized
routes; unsupported inputs do not silently delegate to the baseline. C retains
its existing, explicitly reported fallback policy.
"""
import importlib.util
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
ORIGINAL = ROOT.parent / 'geometry_parallel_v1_20260921'


def load(path, name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


class ExistingEngine:
    def __init__(self, system, family_spec=None):
        self.system = system
        path = ORIGINAL / ('shared_baseline.py' if system == 'baseline'
                           else 'c_adaptive/variant.py')
        self.module = load(path, '_cross_dimension_original_' + system)

    def prepare(self, seconds):
        return {'ok': True, 'no_preparation': True, 'system': self.system}

    def solve(self, task):
        return self.module.solve(task)

    def verify(self, certificate, task, full=True):
        return self.module.verify(certificate, task)


def make_engine(system, family_spec=None):
    if system in ('baseline', 'C'):
        return ExistingEngine(system, family_spec)
    if system not in ('A', 'B'):
        raise ValueError('Unknown system')
    module = load(ROOT / (system.lower() + '_adapter.py'),
                  '_cross_dimension_adapter_' + system)
    return module.Engine(family_spec=family_spec)
