"""Read-only access to C1 and its original task-bound verifier."""
import importlib.util
from pathlib import Path
import sys
from time import perf_counter

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'geometry_parallel_v1_20260921'


def load(path, name):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


c1 = load(OLD / 'c_adaptive/variant.py', '_c2_frozen_c1')
baseline = c1.baseline
Deadline = c1.Deadline


def verify(certificate, task):
    return baseline.assess(certificate, baseline.normalize_task(task))


class Execution(c1.Execution):
    """C1-compatible instrumentation with a reserved fallback interval."""
    def __init__(self, started, limit, progress=None):
        super().__init__(started, limit)
        self.progress = progress
        self.stage_deadline = None

    def emit(self, row):
        if self.progress is not None:
            self.progress({**row, 'since_start_seconds': perf_counter()-self.started})

    def remaining(self):
        whole = super().remaining()
        return whole if self.stage_deadline is None else max(
            0.0, min(whole, self.stage_deadline-perf_counter()))

    def reserve_fallback(self, fraction=0.30):
        self.stage_deadline = perf_counter() + super().remaining()*(1-fraction)

    def release_reserve(self):
        self.stage_deadline = None

    def note(self, reason_code, **details):
        row = {'action': 'decision', 'status': 'decision',
               'reason_code': reason_code, **details}
        self.attempts.append(row)
        self.emit(row)

    def run(self, action, function, **details):
        self.check()
        begin = perf_counter()
        row = {'action': action, **details}
        self.emit({**row, 'status': 'started'})
        try:
            value = function()
            row['status'] = 'returned'
            return value
        except Exception as error:
            row.update(status='failed', error=type(error).__name__, reason=str(error))
            raise
        finally:
            row['elapsed_seconds'] = perf_counter()-begin
            self.attempts.append(row)
            self.emit(row)
