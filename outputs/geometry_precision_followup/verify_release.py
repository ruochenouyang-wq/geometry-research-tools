"""Read-only hash, certificate, original-case and measured-result replay."""
import json
from bridge import ROOT,F,sha,digest,previous_functions
import direct_moments as direct
import piecewise_models as pieces
import enriched_trial as enriched


def previous_integrity():
    original=json.loads((ROOT/'PREVIOUS_AT_START.json').read_text())
    counts={}
    for directory,snapshot in original.items():
        folder=ROOT.parent/directory
        current={str(p.relative_to(folder)):sha(p) for p in folder.rglob('*')
                 if p.is_file() and '__pycache__' not in p.parts}
        changed=[p for p in set(snapshot)|set(current) if snapshot.get(p)!=current.get(p)]
        if changed:raise AssertionError((directory,sorted(changed)))
        counts[directory]=len(current)
    return counts


def nodes(value):
    if isinstance(value,dict):
        if 'format' in value:yield value
        for child in value.values():yield from nodes(child)
    elif isinstance(value,list):
        for child in value:yield from nodes(child)


def replay():
    methods={direct.FULL:direct.verify_full,direct.SECTOR:direct.verify_sector,
             pieces.FORMAT:pieces.verify_piecewise,enriched.FORMAT:enriched.verify,
             previous_functions.L2_FORMAT:previous_functions.verify_model}
    seen=set();counts={}
    for path in sorted(ROOT.rglob('*.json')):
        if path.name in ('PREVIOUS_AT_START.json','MANIFEST.json','VALIDATION.json'):continue
        for value in nodes(json.loads(path.read_text())):
            fmt=value['format']
            if fmt not in methods:continue
            identifier=digest(value)
            if identifier in seen:continue
            if methods[fmt](value) is not True:raise AssertionError((str(path),fmt))
            seen.add(identifier);counts[fmt]=counts.get(fmt,0)+1
    return {'unique_mathematical_records':len(seen),'by_format':counts,
            'optimization_search_rerun':False}


def result_bindings():
    approximate=json.loads((ROOT/'APPROXIMATION_RESULTS.json').read_text())
    for row in approximate['results']:
        cert=json.loads((ROOT/row['path']).read_text())
        assert cert['function']==approximate['original_function']
        assert cert['error_upper']==row['error_upper']
        assert row['target_met']==(F(cert['error_upper'])<=F(row['target']))
        assert cert['expanded_polynomial_coefficient_slots']==row['expanded_coefficient_slots']
    results=json.loads((ROOT/'DIRECT_RESULTS.json').read_text())
    for row in results['results']:
        p=ROOT/row['path'];cert=json.loads(p.read_text())
        assert sha(p)==row['sha256']
        for name in ('function','lower','upper','exact_width','status','tolerance'):
            assert row[name]==cert[name]
    for row in results['matched_radial_budget']:
        assert sha(ROOT/row['old_path'])==row['old_sha256']
        old=json.loads((ROOT/row['old_path']).read_text())
        new=json.loads((ROOT/row['new_path']).read_text())
        assert old['function']==new['function']==row['original_function']
        assert old['exact_width']==row['old_width'] and new['exact_width']==row['new_width']
        assert F(row['width_ratio'])==F(row['old_width'])/F(row['new_width'])
    results=json.loads((ROOT/'ENRICHED_RESULTS.json').read_text())
    for row in results['results']:
        cert=json.loads((ROOT/row['path']).read_text())
        assert cert['function']==results['original_function']
        for name in ('lower','upper','exact_width','status'):assert cert[name]==row[name]
        assert cert['statistics']['residual_squared']==row['residual_squared']
        assert len(cert['trial']['powers'])==row['terms']
    controls=json.loads((ROOT/'BASIS_CONTROL_RESULTS.json').read_text())
    for row in controls['rows']:
        poly=json.loads((ROOT/row['polynomial_path']).read_text())
        fract=json.loads((ROOT/row['fractional_path']).read_text())
        assert poly['source']==fract['source'] and poly['function']==fract['function']
        assert len(poly['trial']['powers'])==len(fract['trial']['powers'])==row['terms']
        assert F(poly['exact_width'])/F(fract['exact_width'])==F(row['width_ratio'])
    return True


def check():
    manifest=json.loads((ROOT/'MANIFEST.json').read_text())
    indexed=[row['path'] for row in manifest['entries']]
    actual={str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file()
            and '__pycache__' not in p.parts and p.name not in ('MANIFEST.json','VALIDATION.json')}
    assert len(indexed)==len(set(indexed)) and set(indexed)==actual
    for row in manifest['entries']:
        p=ROOT/row['path']
        assert not p.is_symlink() and sha(p)==row['sha256'] and p.stat().st_size==row['bytes']
    for row in manifest['project_files']:
        assert sha(ROOT.parent.parent/row['path'])==row['sha256']
    return {'all_passed':True,'implementation_tests_passed':manifest['tests_passed'],
            'indexed_files':len(indexed),'project_files_checked':len(manifest['project_files']),
            'previous_unchanged_files':previous_integrity(),'certificate_replay':replay(),
            'result_bindings_verified':result_bindings(),
            'formal_proof_assistant_checked':False,'model_speedup_measured':False}


if __name__=='__main__':print(json.dumps(check(),ensure_ascii=False,indent=2,sort_keys=True))
