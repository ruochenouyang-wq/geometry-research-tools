"""Read-only release integrity and saved mathematical certificate replay."""
import json
from dependencies import ROOT, F, sha, digest, wide
import function_models as f
import spectral_transfer as s
import research_tool as tool


def previous_integrity():
    snapshot=json.loads((ROOT/'PREVIOUS_AT_START.json').read_text())
    counts={}
    for directory, old in snapshot.items():
        folder=ROOT.parent/directory
        current={str(p.relative_to(folder)):sha(p) for p in folder.rglob('*')
                 if p.is_file() and '__pycache__' not in p.parts}
        differences=[p for p in sorted(set(old)|set(current)) if old.get(p)!=current.get(p)]
        if differences:raise AssertionError((directory,differences))
        counts[directory]=len(current)
    return counts


def nodes(value):
    if isinstance(value,dict):
        if 'format' in value:yield value
        for child in value.values():yield from nodes(child)
    elif isinstance(value,list):
        for child in value:yield from nodes(child)


def replay_saved():
    dispatch={f.ANALYTIC_FORMAT:f.verify_model,f.L2_FORMAT:f.verify_model,
              s.SOURCE:s.verify_source,s.LINF:s.verify,s.L2:s.verify,s.FIXED:s.verify,
              s.COERCIVITY:s.verify_fixed_reference,wide.FULL_FORMAT:wide.verify_full,
              tool.RECORD:tool.verify_record}
    seen=set();counts={};files=set()
    for p in sorted(ROOT.rglob('*.json')):
        if p.name in ('MANIFEST.json','VALIDATION.json','PREVIOUS_AT_START.json'):continue
        for node in nodes(json.loads(p.read_text())):
            fmt=node['format']
            if fmt not in dispatch:continue
            identifier=digest(node)
            if identifier in seen:continue
            if dispatch[fmt](node) is not True:raise AssertionError((str(p),fmt))
            seen.add(identifier);files.add(str(p.relative_to(ROOT)))
            counts[fmt]=counts.get(fmt,0)+1
    result=json.loads((ROOT/'RESULTS.json').read_text())
    for row in result['certificates']:
        if sha(ROOT/row['path'])!=row['sha256']:raise AssertionError(row['path'])
    norms=json.loads((ROOT/'FIXED_REFERENCE_NORMS.json').read_text())
    if [row['degree'] for row in norms]!=[0,4,8,12,16,20,24]:
        raise AssertionError('Fixed-reference norm degrees')
    for row in norms:
        if F(row['coercivity']['distance_squared'])+F(row['model']['error_squared'])!=F(2,9):
            raise AssertionError('Fixed-reference norm identity')
        if F(row['coercivity']['distance_upper'])>=F(1,2):
            raise AssertionError('Fixed-reference stability')
    return {'unique_recognized_mathematical_records':len(seen),'record_types':counts,
            'new_transfer_certificates':len(result['certificates']),
            'fixed_reference_norm_identities':len(norms),
            'files_with_first_occurrence_of_verified_records':len(files)}


def check():
    manifest=json.loads((ROOT/'MANIFEST.json').read_text())
    indexed=[row['path'] for row in manifest['entries']]
    actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*')
            if p.is_file() and '__pycache__' not in p.parts
            and p.name not in ('MANIFEST.json','VALIDATION.json')}
    if len(indexed)!=len(set(indexed)) or set(indexed)!=actual:
        raise AssertionError('Release file inventory')
    for row in manifest['entries']:
        p=ROOT/row['path']
        if p.is_symlink() or sha(p)!=row['sha256'] or p.stat().st_size!=row['bytes']:
            raise AssertionError(row['path'])
    stages=json.loads((ROOT/'STAGES.json').read_text())
    if [v['version'] for v in stages]!=[235,236,237]:raise AssertionError('Version ledger')
    for stage in stages:
        for item in stage['evidence']:
            if not (ROOT/item).is_file():raise AssertionError(item)
    return {'all_passed':True,'release':'V235–V237','executable_increments':3,
            'indexed_files':len(manifest['entries']),
            'previous_unchanged_files':previous_integrity(),'replayed':replay_saved(),
            'implementation_tests':manifest['tests_passed'],
            'optimization_search_rerun':False,'formal_proof_assistant_checked':False}


if __name__=='__main__':print(json.dumps(check(),ensure_ascii=False,indent=2,sort_keys=True))
