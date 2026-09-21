"""Real tokenizer comparisons with explicit baselines, retrieval costs and no model claims."""
from pathlib import Path
import copy,json,statistics,tempfile
from fractions import Fraction as F
import protocol as p
import token_meter as meter
import compiler_v23 as compiler
import legacy_cli
import run

HERE=Path(__file__).resolve().parent
BASELINE=HERE/'fixtures/v35'


def sufficient_baseline(result,path):
    # Same required semantics, direct file evidence, full task and exact rational bounds.
    out=p.summary(result,compact_numbers=False);out['evidence_file']=path;out.pop('evidence')
    out['diagnostics']=p.aggregate_diagnostics(result.get('diagnostics',[]))
    return out


def profile_request(raw):
    new=copy.deepcopy(raw)
    if all(new.get(k)==v for k,v in p.PROFILE_FIELDS.items()):
        for k in p.PROFILE_FIELDS:new.pop(k)
        new['profile']=p.PROFILE
    assert p.expand_request(new)==raw
    return new


def coefficient_requests(raw):
    cp=compiler.compile_goal(raw);problem=cp['problem'];base=profile_request(raw)
    poly=copy.deepcopy(base);poly.pop('potential');poly['potential_coefficients']={'base':problem['q0'],'directions':problem['directions']}
    quad=copy.deepcopy(base);quad.pop('penalty');q=problem['penalty'];d=len(problem['directions'])
    quad['quadratic']={'constant':q['constant'],'linear':q['linear'],'hessian_upper':[q['hessian'][i][j] for i in range(d) for j in range(i,d)]}
    both=copy.deepcopy(poly);both.pop('penalty');both['quadratic']=quad['quadratic']
    for candidate in (poly,quad,both):
        compiled=compiler.compile_goal(p.expand_request(candidate));assert compiled['problem']==cp['problem'] and compiled['threshold']==cp['threshold']
    return {'profile':base,'coefficients':poly,'triangular_quadratic':quad,'both':both}


def main():
    output=HERE/'results';output.mkdir(exist_ok=True);rows=[];events=[];input_probes=[]
    with tempfile.TemporaryDirectory(prefix='geometry-token-bench-') as tmp:
        service=p.Service(tmp)
        for file in sorted((BASELINE/'results/proofs').glob('*.json')):
            raw=json.loads(file.read_text());receipt=service.receipt(raw);proofref=receipt['evidence']
            retrieved=service.inspect(proofref);assert p.checked(raw)==retrieved
            budget=service.budget(receipt);assert budget['fits']
            source=service.store.get(receipt['task']);kind=receipt['kind']
            if kind=='ground_goal':
                req={'op':'submit','goal':profile_request(source)}
                variants=coefficient_requests(source)
                input_probes.append({'case':file.stem,'variants':{name:meter.counts(meter.wire(value)) for name,value in variants.items()},'all_compile_to_same_problem':True})
            else:req={'op':'optimize' if kind=='matrix_envelope' else 'lemma','spec':source}
            strong=sufficient_baseline(raw,'result.json')
            cli={'output':'result.json','method':raw['certificate']['method'],'status':raw['certificate'].get('status'),'statistics':raw.get('statistics')}
            compact_cert=meter.wire(retrieved)
            formats={'old_complete_pretty_file':file.read_text(),'same_complete_compact_json':meter.wire(raw),
                     'old_cli_brief':json.dumps(cli,ensure_ascii=False),
                     'sufficient_exact_baseline':meter.wire(strong),'decision_receipt':meter.wire(receipt),
                     'actual_tool_response':meter.wire(budget),
                     'request':meter.wire(req),'get_full_proof_request':meter.wire({'op':'inspect','ref':proofref}),
                     'get_full_proof_response':compact_cert}
            measures={name:meter.counts(text) for name,text in formats.items()}
            row={'case':file.stem,'kind':kind,'status':receipt['status'],'text_tokens':measures,
                 'receipt':receipt,'request':req,'same_exact_evidence_after_retrieval':True}
            rows.append(row);events.extend([{'role':'request','text':formats['request']},{'role':'response','text':formats['actual_tool_response']}])
        gate={}
        for encoding in meter.ENCODINGS:
            sizes=sorted(r['text_tokens']['decision_receipt'][encoding] for r in rows)
            totals=[r['text_tokens']['actual_tool_response'][encoding]+r['text_tokens']['request'][encoding] for r in rows]
            median=statistics.median(sizes);p90=sizes[min(len(sizes)-1,(9*len(sizes)+9)//10-1)]
            gate[encoding]={'median_receipt':median,'p90_receipt':p90,'median_request_and_receipt':statistics.median(totals),
                            'numeric_size_gate':median<=250 and p90<=500 and statistics.median(totals)<=600}
        schema=meter.wire(run.OPERATIONS)
        sessions={'first_use_plus_all_requests':meter.conversation(events,p.PROFILE_TEXT,schema),
                  'one_time_instructions':meter.counts(p.PROFILE_TEXT),'one_time_operation_schema':meter.counts(schema),
                  'interpretation':'Exact text counts only. Full proof retrieval adds the explicit request and proof response; it is not free and can erase savings.'}
        result={'encoding_implementation':'tiktoken 0.11.0','encodings':list(meter.ENCODINGS),'actual_model_calls':0,'actual_api_usage':None,
                'protocol':'Development certificate corpus. Required-information baseline and shorter but information-limited old CLI are both shown. Counts include explicit text only.',
                'cases':rows,'input_format_probes':input_probes,'size_gate':gate,'sessions':sessions}
        (output/'token_benchmark.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps({'cases':len(rows),'size_gate':gate,'sessions':sessions},indent=2))

if __name__=='__main__':main()
