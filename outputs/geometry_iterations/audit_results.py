"""Recheck every saved mathematical document using dense endpoint counts."""
from pathlib import Path
from fractions import Fraction as F
import json
import platform
from datetime import datetime,timezone
from check_document import verify_document


def main():
    root=Path(__file__).resolve().parent;checked=[];failures=[];benchmarks=[]
    for path in sorted((root/'results').rglob('*.json')):
        doc=json.loads(path.read_text())
        if path.name=='benchmark.json':
            # Accuracy declarations in timing rows must agree with exact widths.
            for row in doc['rows']:
                if 'exact_width' in row and 'tolerance' in row:
                    met=F(row['exact_width'])<=F(row['tolerance'])
                    stated=row.get('target_met',row.get('status')=='target_met')
                    if met!=stated:failures.append(str(path.relative_to(root))+': inaccurate target claim')
            benchmarks.append(str(path.relative_to(root)));continue
        if verify_document(doc,independent=True):checked.append(str(path.relative_to(root)))
        else:failures.append(str(path.relative_to(root)))
    report={'checked_at_utc':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),
            'system':platform.system(),'machine':platform.machine(),
            'independent_dense_endpoint_checks':True,'documents_checked':len(checked),
            'benchmarks_checked':len(benchmarks),'checked_documents':checked,'failures':failures,
            'valid':not failures,
            'scope':'Recomputed rational certificates and target claims. Shared Legendre assembly; not a formal proof assistant.'}
    (root/'VALIDATION.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['documents_checked','benchmarks_checked','failures','valid']},ensure_ascii=False))
    return 0 if not failures else 1


if __name__=='__main__':raise SystemExit(main())
