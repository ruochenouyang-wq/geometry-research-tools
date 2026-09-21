"""Read-only access to the frozen v317 implementation and its math backend."""
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parent
PREVIOUS=ROOT.parent/'geometry_strategy_v238_v317'
sys.dont_write_bytecode=True
if str(PREVIOUS) not in sys.path:
    sys.path.insert(0,str(PREVIOUS))
import backend
import research_service as previous_service
from verification import verify_any,assess

F=backend.F
