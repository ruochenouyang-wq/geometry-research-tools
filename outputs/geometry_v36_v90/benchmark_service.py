"""Actual V90 tool calls, retained proofs and full visible response token costs."""
from pathlib import Path
from fractions import Fraction as F
import json,statistics,tempfile,time
import run,protocol,token_meter as meter

HERE=Path(__file__).resolve().parent
GOAL={'profile':'sphere_ground','parameters':['a'],'potential':'a*t','penalty':'a*a/5',
      'domain':{'kind':'all_real'},'threshold':'0'}
REQUESTS={
 'extrema_cubic':{'op':'extrema','spec':{'polynomial':[0,1,0,-1]}},
 'extrema_budget_open':{'op':'extrema','spec':{'polynomial':[0,1,0,-1]},'max_nodes':0},
 'uniform_sharp_sixth':{'op':'uniform','spec':{}},
 'uniform_counterexample_seventh':{'op':'uniform','spec':{'coefficient':'1/7'}},
 'uniform_barta_fifth':{'op':'uniform','rule':'barta','spec':{'coefficient':'1/5','all_real':True}},
 'uniform_barta_obstruction':{'op':'uniform','rule':'barta','spec':{'coefficient':'1/5','order':1}},
 'positive_quadratic':{'op':'positive','spec':{'potential':[0,0,6]},'options':{'degree':2}},
 'positive_sixth':{'op':'positive','spec':{'potential':[0,0,6]},'options':{'degree':6}},
 'positive_sturm':{'op':'positive','spec':{'potential':[0,0,6]},'options':{'degree':6,'range_backend':'exact_sturm_minimum_v1'}},
 'matrix_precise_scalar':{'op':'refine','spec':{'directions':[[0,1]]},'budget':0},
 'matrix_precise_odd_even':{'op':'refine','spec':{'directions':[[0,1],[0,0,1]]},'budget':8},
 'original_goal_sharp':{'op':'research','goal':GOAL},
 'original_goal_barta':{'op':'research','goal':GOAL,'rule':'barta'},
 'original_goal_open':{'op':'research','goal':dict(GOAL,penalty='a*a/7')},
}


def main():
    output=HERE/'results';proofs=output/'service_certificates';proofs.mkdir(parents=True,exist_ok=True)
    examples=HERE/'examples/v90';examples.mkdir(parents=True,exist_ok=True)
    rows=[];events=[]
    with tempfile.TemporaryDirectory(prefix='geometry-service-benchmark-') as tmp:
        s=protocol.Service(tmp)
        for name,request in REQUESTS.items():
            (examples/(name+'.json')).write_text(json.dumps(request,indent=2)+'\n')
            started=time.perf_counter();response=run.call(s,request);elapsed=time.perf_counter()-started
            assert response['fits'];receipt=json.loads(response['text'])
            retrieve={'op':'inspect','ref':receipt['evidence']};evidence=run.call(s,retrieve)
            assert protocol.checked(evidence)==evidence
            (proofs/(name+'.json')).write_text(json.dumps(evidence,indent=2)+'\n')
            exact=protocol.summary(evidence);exact.pop('evidence');exact['evidence_file']='result.json'
            formats={'request':meter.wire(request),'receipt':meter.wire(receipt),'actual_tool_response':meter.wire(response),
                     'sufficient_exact_baseline':meter.wire(exact),'full_compact_proof':meter.wire(evidence),
                     'retrieval_request':meter.wire(retrieve),'retrieval_response':meter.wire(evidence)}
            counts={key:meter.counts(value) for key,value in formats.items()}
            rows.append({'case':name,'request':request,'receipt':receipt,'status':receipt['status'],
                         'seconds_including_verification_and_budget':elapsed,'text_tokens':counts,
                         'certificate':'results/service_certificates/'+name+'.json','exact_evidence_retrieved':True})
            events.extend([{'role':'request','text':formats['request']},{'role':'response','text':formats['actual_tool_response']}])
        gate={}
        for enc in meter.ENCODINGS:
            sizes=sorted(r['text_tokens']['receipt'][enc] for r in rows)
            totals=[r['text_tokens']['request'][enc]+r['text_tokens']['actual_tool_response'][enc] for r in rows]
            p90=sizes[(9*len(sizes)+9)//10-1]
            gate[enc]={'median_receipt':statistics.median(sizes),'p90_receipt':p90,
                       'median_request_and_actual_response':statistics.median(totals),
                       'numeric_size_gate':statistics.median(sizes)<=250 and p90<=500 and statistics.median(totals)<=600}
        output_data={'solver_triggered_model_calls':0,'actual_model_comparison':False,
                     'development_subagents':5,'tokenizer_mapping_to_5_6_or_astra_asserted':False,
                     'cases':rows,'size_gate':gate,
                     'session':meter.conversation(events,protocol.PROFILE_TEXT,meter.wire(run.OPERATIONS)),
                     'limitations':['Timing is one local tool execution per case, not a model timing result',
                                    'Full proof retrieval and first-use capabilities cost tokens',
                                    'Stored references require their originating store; exported certificate files are portable',
                                    'All request and response fields are counted; unknown API framing and hidden reasoning are unavailable']}
        (output/'service_benchmark.json').write_text(json.dumps(output_data,indent=2)+'\n')
        print(json.dumps({'cases':len(rows),'size_gate':gate,'all_evidence_replayed':True},indent=2))

if __name__=='__main__':main()
