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
RELATION_TYPES = frozenset(('equivalent','implies','upper','lower','approximation','heuristic'))
NODE_TYPES = frozenset(('problem','claim','candidate','certificate'))
NODE_FORMAT = 'research_evidence_node_v1'
RELATION_FORMAT = 'research_evidence_relation_v1'


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


def problem_claim(certificate):
    """Extract identity metadata only; this function never proves a certificate."""
    if type(certificate) is not dict:
        raise ValueError('Certificate object required for problem binding')
    fmt = certificate.get('format')
    if fmt in ('direct_axis_moment_full_v1','singular_fractional_trial_temple_v1'):
        keys = ('function','geometry','scope','mean_zero','eigenvalue_index','measure','tolerance')
        return {'kind':'spectrum', **{k:deepcopy(certificate[k]) for k in keys}}
    if fmt == 'sphere_structured_piecewise_L2_v1':
        return {'kind':'approximation','geometry':'unit_S2',
                'function':deepcopy(certificate['function']),
                'scope':certificate['scope'],'norm':certificate['norm']}
    raise ValueError('No checked problem-binding schema for this format')


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
        self._cache_verifier_refs = {}
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

    def verify(self, ref, expected_problem=None):
        # Always re-read complete bytes before looking up mathematical verification.
        value = self.get(ref)
        if expected_problem is not None:
            try:
                matches = type(expected_problem) is dict and canonical(problem_claim(value)) == canonical(expected_problem)
            except (ValueError, KeyError, TypeError):
                matches = False
            if not matches:
                return {'reference':ref,'verified':False,'status':'problem_mismatch'}
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
            if self.verifier is not None:
                # Keep a strong reference so Python cannot recycle an old callback id.
                self._cache_verifier_refs[id(self.verifier)] = self.verifier
        return {'reference':ref, 'verified':accepted, 'verifier_version':version,
                'cache_hit':False, 'status':'verified' if accepted else 'rejected'}


class ResearchGraph:
    """Typed research dependencies; naming a relation never proves it."""
    def __init__(self, store, relation_rules=None):
        if not isinstance(store, EvidenceStore):
            raise ValueError('ResearchGraph requires an EvidenceStore')
        self.store = store
        self.nodes = {}
        self.relations = {}
        self.relation_rules = {} if relation_rules is None else dict(relation_rules)
        for name, rule in self.relation_rules.items():
            if (type(name) is not str or not name or len(name)>128 or name.startswith('certified_')
                    or type(rule) is not dict or set(rule) != {'version','kinds','check'}
                    or type(rule['version']) is not str or not rule['version']
                    or type(rule['kinds']) is not list or not rule['kinds']
                    or any(type(k) is not str or k not in RELATION_TYPES for k in rule['kinds'])
                    or not callable(rule['check'])):
                raise ValueError('Relation rule requires explicit name, version, kinds and trusted checker')

    def add_node(self, kind, statement, evidence=None):
        if type(kind) is not str or kind not in NODE_TYPES or type(statement) is not dict:
            raise ValueError('Require a supported node type and statement object')
        if evidence is not None:
            self.store.get(evidence)
        node = {'format':NODE_FORMAT, 'kind':kind, 'statement':deepcopy(statement),
                'evidence_reference':evidence}
        if digest(node) not in self.nodes and len(self.nodes) >= 512:
            raise ValueError('Research node budget exceeded')
        ref = self.store.put(node)
        self.nodes[ref] = node
        return ref

    def add_relation(self, kind, sources, target, *, conditions=None, scope=None, rule=None):
        if type(kind) is not str or kind not in RELATION_TYPES:
            raise ValueError('Unsupported research relation')
        if (type(sources) is not list or not 1 <= len(sources) <= 32
                or any(type(ref) is not str for ref in sources) or len(set(sources)) != len(sources)):
            raise ValueError('Require one to 32 distinct source nodes')
        if any(ref not in self.nodes for ref in sources) or type(target) is not str or target not in self.nodes:
            raise ValueError('Relation endpoints must exist in this graph')
        if conditions is not None and type(conditions) is not dict:
            raise ValueError('Conditions must be an object')
        if scope is not None and type(scope) is not dict:
            raise ValueError('Scope must be an object')
        if rule is not None and (type(rule) is not str or not rule or len(rule)>128):
            raise ValueError('Rule must be a short explicit identifier')
        relation = {'format':RELATION_FORMAT,'kind':kind,'sources':sorted(sources),'target':target,
                    'conditions':deepcopy(conditions) if conditions is not None else {},
                    'scope':deepcopy(scope) if scope is not None else {},'rule':rule}
        if digest(relation) not in self.relations and len(self.relations) >= 2048:
            raise ValueError('Research relation budget exceeded')
        self._topological(extra=relation)
        ref = self.store.put(relation)
        self.relations[ref] = relation
        return ref

    def _topological(self, extra=None):
        adjacency = {ref:set() for ref in self.nodes}
        indegree = {ref:0 for ref in self.nodes}
        relations = [self.store.get(ref) for ref in self.relations]
        if extra is not None:
            relations.append(extra)
        for relation in relations:
            target = relation['target']
            for source in relation['sources']:
                if source not in adjacency or target not in adjacency:
                    raise ValueError('Dangling research dependency')
                if target not in adjacency[source]:
                    adjacency[source].add(target)
                    indegree[target] += 1
        ready = sorted(ref for ref,n in indegree.items() if n == 0)
        ordered = []
        while ready:
            source = ready.pop()
            ordered.append(source)
            for target in sorted(adjacency[source]):
                indegree[target] -= 1
                if indegree[target] == 0:
                    ready.append(target)
        if len(ordered) != len(self.nodes):
            raise ValueError('Circular research justification is forbidden')
        return ordered

    def topological_order(self):
        return self._topological()

    def relation_status(self, ref):
        if ref not in self.relations:
            raise ValueError('Unknown relation reference')
        relation = self.store.get(ref)
        result = {'reference':ref,'kind':relation['kind'],'verified':False}
        kind, rule_name = relation['kind'], relation['rule']
        if kind == 'heuristic':
            return {**result,'status':'heuristic'}
        if rule_name is None:
            return {**result,'status':'unverified_no_rule'}
        builtins = {'certified_bound_v1':{'upper','lower'},
                    'certified_approximation_v1':{'approximation'}}
        rule = self.relation_rules.get(rule_name)
        if rule_name not in builtins and rule is None:
            return {**result,'status':'unverified_unknown_rule'}
        allowed = builtins[rule_name] if rule_name in builtins else rule['kinds']
        if kind not in allowed:
            return {**result,'status':'rule_kind_mismatch'}
        try:
            self._topological()
            sources = [self.store.get(source) for source in relation['sources']]
            target = self.store.get(relation['target'])
            if target['kind'] != 'claim' or type(target['statement'].get('problem')) is not dict:
                return {**result,'status':'target_claim_missing'}
            problem = target['statement']['problem']
            if canonical(relation['scope']) != canonical(problem):
                return {**result,'status':'scope_mismatch'}
            for key, expected in relation['conditions'].items():
                if key not in problem or canonical(expected) != canonical(problem[key]):
                    return {**result,'status':'condition_not_established'}
            certificates = []
            for source in sources:
                evidence = source['evidence_reference']
                if source['kind'] != 'certificate' or evidence is None:
                    return {**result,'status':'source_not_certified'}
                if canonical(source['statement']) != canonical(problem):
                    return {**result,'status':'source_problem_mismatch'}
                if not self.store.verify(evidence,expected_problem=problem)['verified']:
                    return {**result,'status':'source_verification_failed'}
                certificates.append(self.store.get(evidence))
            if rule_name in builtins:
                if len(certificates) != 1:
                    return {**result,'status':'rule_arity_mismatch'}
                import backend
                cert = certificates[0]
                statement = target['statement']
                if rule_name == 'certified_bound_v1':
                    if problem.get('kind') != 'spectrum' or set(statement) != {'problem','value'}:
                        return {**result,'status':'bound_claim_mismatch'}
                    value = backend.pieces.rational(statement['value'])
                    endpoint = backend.pieces.rational(cert[kind])
                    accepted = value >= endpoint if kind == 'upper' else value <= endpoint
                else:
                    if (problem.get('kind') != 'approximation'
                            or set(statement) != {'problem','error_upper','approximant_digest'}):
                        return {**result,'status':'approximation_claim_mismatch'}
                    accepted = (statement['approximant_digest'] == digest(cert['cells'])
                                and backend.pieces.rational(statement['error_upper'])
                                >= backend.pieces.rational(cert['error_upper']))
                version = rule_name
            else:
                context = {'problem':deepcopy(problem),'certificates':deepcopy(certificates)}
                accepted = rule['check'](deepcopy(relation),deepcopy(sources),deepcopy(target),context) is True
                version = rule['version']
            return {**result,'verified':accepted,'status':'verified' if accepted else 'relation_rejected',
                    'rule':rule_name,'rule_version':version}
        except Exception as error:
            return {**result,'status':'relation_check_error','error_type':type(error).__name__}

    def export(self):
        return {'format':'research_graph_v1','nodes':sorted(self.nodes),
                'relations':sorted(self.relations)}

    @classmethod
    def from_export(cls, store, data, relation_rules=None):
        """Read-only reconstruction: hashes, schemas, endpoints and dependencies."""
        if type(data) is not dict or set(data) != {'format','nodes','relations'} or data['format'] != 'research_graph_v1':
            raise ValueError('Invalid research graph manifest')
        for field, limit in (('nodes',512),('relations',2048)):
            values = data[field]
            if (type(values) is not list or len(values)>limit
                    or any(type(v) is not str for v in values) or len(set(values)) != len(values)):
                raise ValueError('Invalid graph reference inventory')
        graph = cls(store,relation_rules)
        for ref in data['nodes']:
            node = store.get(ref)
            if (type(node) is not dict or set(node) != {'format','kind','statement','evidence_reference'}
                    or node['format'] != NODE_FORMAT or type(node['kind']) is not str
                    or node['kind'] not in NODE_TYPES or type(node['statement']) is not dict):
                raise ValueError('Invalid stored research node')
            if node['evidence_reference'] is not None:
                store.get(node['evidence_reference'])
            graph.nodes[ref] = node
        for ref in data['relations']:
            relation = store.get(ref)
            if (type(relation) is not dict
                    or set(relation) != {'format','kind','sources','target','conditions','scope','rule'}
                    or relation['format'] != RELATION_FORMAT or type(relation['kind']) is not str
                    or relation['kind'] not in RELATION_TYPES):
                raise ValueError('Invalid stored relation schema')
            sources = relation['sources']; target = relation['target']; rule = relation['rule']
            if (type(sources) is not list or not 1<=len(sources)<=32
                    or any(type(v) is not str or v not in graph.nodes for v in sources)
                    or len(set(sources)) != len(sources) or type(target) is not str or target not in graph.nodes
                    or type(relation['conditions']) is not dict or type(relation['scope']) is not dict
                    or (rule is not None and (type(rule) is not str or not 1<=len(rule)<=128))):
                raise ValueError('Invalid stored relation endpoints or conditions')
            graph.relations[ref] = relation
        graph._topological()
        return graph

    def replay(self):
        checked = self.from_export(self.store,self.export(),self.relation_rules)
        results = [checked.relation_status(ref) for ref in sorted(checked.relations)]
        return {'structurally_valid':True,'nodes':len(checked.nodes),'relations':results,
                'verified_relations':sum(item['verified'] for item in results),
                'all_relations_verified':bool(results) and all(item['verified'] for item in results),
                'formal_proof_assistant_checked':False}
