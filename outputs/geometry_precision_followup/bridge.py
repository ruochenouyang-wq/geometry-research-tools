"""Read-only frozen arithmetic and previous function definitions."""
from pathlib import Path
import sys, json, hashlib
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
PREVIOUS=ROOT.parent/'geometry_general_v235_v237'
sys.path.insert(0,str(PREVIOUS))
from dependencies import F, exact, same, canonical, digest, base, arithmetic, wide
import function_models as previous_functions
import spectral_transfer as previous_transfer
sys.path.insert(0,str(ROOT))
def save(name,value):
    p=ROOT/name;p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(value,ensure_ascii=False,sort_keys=True,indent=2,allow_nan=False)+'\n')
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
