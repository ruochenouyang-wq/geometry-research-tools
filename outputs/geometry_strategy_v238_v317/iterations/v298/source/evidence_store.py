"""Content-addressed research evidence. A stored object is not a proved claim."""
from pathlib import Path
import hashlib
import json


def canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False,
                      separators=(',', ':'), allow_nan=False).encode('utf-8')


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


class EvidenceStore:
    def __init__(self, root, verifier=None):
        self.root = Path(root)
        self.root.mkdir(parents=True, exist_ok=True)
        self.verifier = verifier

    def put(self, value):
        data = canonical(value)
        ref = hashlib.sha256(data).hexdigest()
        path = self.root/(ref+'.json')
        if path.exists():
            if path.read_bytes() != data:
                raise ValueError('Existing content does not match reference')
        else:
            path.write_bytes(data)
        return ref

    def get(self, ref):
        return json.loads((self.root/(ref+'.json')).read_bytes())
