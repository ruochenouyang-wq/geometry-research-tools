"""C3 access to frozen C1/C2 components; old source files are never modified."""
import importlib.util
from pathlib import Path
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
C2_ROOT = ROOT.parent / 'c_demand_extension_20260922'


def load(path, name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


_old = load(C2_ROOT / 'runtime.py', '_c3_frozen_c2_runtime')
c1, baseline, Deadline, Execution = _old.c1, _old.baseline, _old.Deadline, _old.Execution


def verify(certificate, task):
    """No cache, alternate checker, or caller-supplied verification receipt."""
    return baseline.assess(certificate, baseline.normalize_task(task))


def step_component():
    return load(C2_ROOT / 'step_extension.py', '_c3_frozen_c2_step')
