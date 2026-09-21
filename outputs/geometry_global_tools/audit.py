"""Independent reconstruction of saved global and local-trap claims."""
from pathlib import Path
from fractions import Fraction as F
from datetime import datetime,timezone
import json
import platform
import global_search as g
from run import verify_document


def main():
    root=Path(__file__).resolve().parent;data=root/'results';checked=[]
    for path in sorted(data.glob('*.json')):
        document=json.loads(path.read_text())
        if 'certificate' in document:
            assert verify_document(document),path.name;checked.append(path.name)
    local=json.loads((data/'local_search.json').read_text());p=g.normalized(local['problem'])
    record=local['candidate_record'];a=g.rationals(local['candidate_parameters'])
    assert record['at']==local['candidate_parameters'] and g.inside(a,[g.rationals(pair) for pair in p['box']])
    q=g.at(p,a);sc=record['spectral_certificate']
    assert g.v3.verify(sc,independent=True) and sc['eigenvalue_index']==1 and g.v3.base.potential(sc['q_coefficients'])==q
    mu,_,_=g.v3.old.residual_statistics(q,g.v7.coefficients(record['trial_legendre']))
    assert record['rayleigh_quotient']==str(mu) and F(local['candidate_upper'])==mu+g.penalty(p,a)
    assert local['global_optimality_proved'] is False
    full=json.loads((data/'tilted_double_well_concavity.json').read_text())['certificate']
    positive=json.loads((data/'positive_basin_global.json').read_text())['certificate']
    trap=json.loads((data/'local_trap_evidence.json').read_text());endpoint_lowers=[]
    for a in [F(0),F(6)]:
        rec=next(x for x in positive['points'] if x['at']==[str(a)])
        lo=F(rec['spectral_certificate']['lower'])+g.penalty(positive['problem'],[a])
        assert lo>F(positive['candidate_upper']);endpoint_lowers.append(str(lo))
    difference=F(positive['global_lower'])-F(full['candidate_upper']);assert difference>0
    assert trap['positive_basin_boundary_lowers']==endpoint_lowers
    for field,value in [('local_candidate',local['candidate_parameters']),('local_upper',local['candidate_upper']),
                        ('positive_basin_lower',positive['global_lower']),('positive_basin_upper',positive['candidate_upper']),
                        ('full_domain_upper',full['candidate_upper']),('certified_improvement_over_every_positive_basin_point',str(difference))]:
        assert trap[field]==value
    bench=json.loads((data/'benchmark.json').read_text())
    for row in bench['rows']:
        assert F(row['global_gap'])==F(row['candidate_upper'])-F(row['global_lower'])
        assert (F(row['global_gap'])<=F(row['tolerance']))==(row['status']=='epsilon_global')
        assert len(row['samples_seconds'])==3 and row['median_seconds']==sorted(row['samples_seconds'])[1]
    report={'checked_at_utc':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),
            'system':platform.system(),'machine':platform.machine(),'valid':True,
            'global_certificates_checked':len(checked),'documents':checked,'benchmark_rows_checked':len(bench['rows']),
            'local_candidate_checked':True,'inferior_local_minimum_witness_checked':True,
            'dense_spectral_endpoint_reconstruction':True,'formal_assistant_checked':False}
    (root/'VALIDATION.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))


if __name__=='__main__':main()
