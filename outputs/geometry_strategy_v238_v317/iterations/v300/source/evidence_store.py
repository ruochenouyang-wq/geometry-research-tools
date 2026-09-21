"""Content-addressed research evidence. A stored object is not a proved claim."""
from pathlib import Path
import hashlib
import json
import math
import os
import tempfile


def _json_value(value):
    if value is None or type(value) in (str, bool, int):
        return
    if type(value) is float and math.isfinite(value):
        return
    if type(value) is list:
        for item in value:
            _json_value(item)
        return
    if type(value) is dict and all(type(key) is str for key in value):
        for item in value.values():
            _json_value(item)
        return
    raise ValueError('Evidence must be strict finite JSON with string object keys')


def canonical(value):
    _json_value(value)
    return json.dumps(value, sort_keys=True, ensure_ascii=False,
                      separators=(',', ':'), allow_nan=False).encode('utf-8')


def digest(value):
    return hashlib.sha256(canonical(value)).hexdigest()


def summary(value):
    data = canonical(value)
    keys = sorted(value) if isinstance(value, dict) else []
    fmt = value.get('format') if isinstance(value, dict) else None
    return {'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data),
            'format_preview': fmt[:128] if isinstance(fmt, str) else None,
            'top_level_fields': len(keys), 'key_preview': [k[:128] for k in keys[:8]],
            'complete_content_hashed': True, 'summary_is_certificate': False}


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
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(dir=self.root, prefix='.pending-', delete=False) as handle:
                    temporary = Path(handle.name)
                    handle.write(data)
                    handle.flush()
                    os.fsync(handle.fileno())
                try:
                    os.link(temporary, path)
                except FileExistsError:
                    if path.read_bytes() != data:
                        raise ValueError('Concurrent content does not match reference')
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
        if self.get(ref) != value:
            raise ValueError('Read-after-write content mismatch')
        return ref

    def get(self, ref):
        data = (self.root/(ref+'.json')).read_bytes()
        if hashlib.sha256(data).hexdigest() != ref:
            raise ValueError('Evidence content hash mismatch')
        value = json.loads(data)
        if canonical(value) != data:
            raise ValueError('Stored evidence is not canonical JSON')
        return value

    def summary(self, ref):
        return summary(self.get(ref))
