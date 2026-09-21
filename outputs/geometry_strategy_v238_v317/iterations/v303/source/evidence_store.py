"""Content-addressed research evidence. A stored object is not a proved claim."""
from pathlib import Path
import hashlib
import json
import math
import os
import tempfile
import re
import stat
import sys
from copy import deepcopy

KNOWN_FORMATS = frozenset(('direct_axis_moment_full_v1',
                           'sphere_structured_piecewise_L2_v1',
                           'singular_fractional_trial_temple_v1'))


def _json_value(value, depth=0):
    if depth > 96:
        raise ValueError('Evidence nesting limit exceeded')
    if value is None or type(value) in (str, bool, int):
        return
    if type(value) is float and math.isfinite(value):
        return
    if type(value) is list:
        for item in value:
            _json_value(item, depth+1)
        return
    if type(value) is dict and all(type(key) is str for key in value):
        for item in value.values():
            _json_value(item, depth+1)
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


def default_verify(certificate):
    """Replay the frozen backend; never launch a solver or trust stored flags."""
    import backend
    return backend.legacy.handle({'op':'verify', 'certificate':certificate})['verified'] is True


def default_verifier_version():
    """Fingerprint frozen Python sources at store construction; no cross-run cache."""
    import backend
    names = ('geometry_precision_followup', 'geometry_general_v235_v237',
             'geometry_generalization_notes', 'geometry_retests_v195_v234',
             'geometry_cycles_v95_v194', 'geometry_v91_v94', 'geometry_v36_v90')
    files = {}
    for name in names:
        for path in sorted((backend.ROOT.parent/name).rglob('*.py')):
            files[str(path.relative_to(backend.ROOT.parent))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return 'frozen-python:'+digest({'files':files,'python':sys.version})


class EvidenceStore:
    def __init__(self, root, verifier=None, *, max_bytes=32*1024*1024, verifier_version=None):
        if type(max_bytes) is not int or not 64 <= max_bytes <= 256*1024*1024:
            raise ValueError('Evidence byte limit must be 64..268435456')
        supplied = Path(root)
        if supplied.is_symlink():
            raise ValueError('Store root must not be a symlink')
        self.root = supplied.resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.verifier = verifier
        self.max_bytes = max_bytes
        if verifier_version is not None and (type(verifier_version) is not str or not 1 <= len(verifier_version) <= 128):
            raise ValueError('Verifier version must be a nonempty string of at most 128 characters')
        self.verifier_version = (default_verifier_version() if verifier is None and verifier_version is None
                                 else verifier_version)
        self._verified_cache = {}
        self.verification_calls = 0
        self.cache_hits = 0

    def _path(self, ref):
        if type(ref) is not str or not re.fullmatch('[0-9a-f]{64}', ref):
            raise ValueError('Require a full SHA-256 content reference')
        if self.root.is_symlink() or self.root.resolve() != self.root:
            raise ValueError('Store root location changed')
        return self.root/(ref+'.json')

    def _read(self, path):
        flags = os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0) | getattr(os, 'O_NONBLOCK', 0)
        fd = os.open(path, flags)
        try:
            info = os.fstat(fd)
            if not stat.S_ISREG(info.st_mode) or info.st_size > self.max_bytes:
                raise ValueError('Evidence must be a bounded regular file')
            with os.fdopen(fd, 'rb', closefd=False) as handle:
                data = handle.read(self.max_bytes+1)
            if len(data) > self.max_bytes:
                raise ValueError('Evidence size limit exceeded')
            return data
        finally:
            os.close(fd)

    def put(self, value):
        data = canonical(value)
        if len(data) > self.max_bytes:
            raise ValueError('Evidence size limit exceeded')
        ref = hashlib.sha256(data).hexdigest()
        path = self._path(ref)
        if path.exists() or path.is_symlink():
            if self._read(path) != data:
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
                    if self._read(path) != data:
                        raise ValueError('Concurrent content does not match reference')
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)
        if self.get(ref) != value:
            raise ValueError('Read-after-write content mismatch')
        return ref

    def get(self, ref):
        data = self._read(self._path(ref))
        if hashlib.sha256(data).hexdigest() != ref:
            raise ValueError('Evidence content hash mismatch')
        value = json.loads(data)
        if canonical(value) != data:
            raise ValueError('Stored evidence is not canonical JSON')
        return value

    def summary(self, ref):
        return summary(self.get(ref))

    def format_status(self, ref):
        value = self.get(ref)
        fmt = value.get('format') if isinstance(value, dict) else None
        return {'format': fmt if isinstance(fmt, str) else None,
                'status': 'recognized_unverified' if isinstance(fmt,str) and fmt in KNOWN_FORMATS
                          else 'unknown_format_unverified',
                'verified': False}

    def verify(self, ref):
        # Always re-read complete bytes before looking up mathematical verification.
        value = self.get(ref)
        fmt = value.get('format') if isinstance(value, dict) else None
        if self.verifier is None and (type(fmt) is not str or fmt not in KNOWN_FORMATS):
            return {'reference':ref, 'verified':False, 'status':'unknown_format'}
        version = self.verifier_version
        cache_key = (ref, version, id(self.verifier) if self.verifier is not None else None)
        if version is not None and cache_key in self._verified_cache:
            self.cache_hits += 1
            return {'reference':ref, 'verified':True, 'status':'verified',
                    'verifier_version':version, 'cache_hit':True}
        verifier = self.verifier or default_verify
        self.verification_calls += 1
        try:
            accepted = verifier(deepcopy(value)) is True
        except Exception as error:
            return {'reference':ref, 'verified':False, 'status':'verifier_error',
                    'error_type':type(error).__name__}
        if accepted and version is not None:
            self._verified_cache[cache_key] = True
        return {'reference':ref, 'verified':accepted, 'verifier_version':version,
                'cache_hit':False, 'status':'verified' if accepted else 'rejected'}
