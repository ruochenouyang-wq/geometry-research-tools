"""Per-feature evidence V36–V55. Local protocol workloads, not model trials."""
from pathlib import Path
from fractions import Fraction as F
import copy,json,tempfile
import protocol as p
import token_meter as m
import compiler_v23 as compiler
from benchmark_tokens import coefficient_requests,profile_request

HERE=Path(__file__).resolve().parent
BASE=HERE/'fixtures/v35'


def measure(before,after):
    return {'before_text_tokens':m.counts(m.wire(before)),'after_text_tokens':m.counts(m.wire(after))}


def main():
    records=[]
    def record(v,feature,check,**data):records.append({'version':v,'feature':feature,'checked_property':check,**data})
    raw=json.loads((BASE/'results/proofs/v35_generic.json').read_text())
    simple=json.loads((BASE/'results/proofs/v35_anisotropic_goal.json').read_text())
    matrix=json.loads((BASE/'results/proofs/v34_optimized_four_directions.json').read_text())
    goal=simple['certificate']['compiled_goal']['source']
    with tempfile.TemporaryDirectory(prefix='geometry-version-probes-') as tmp:
        s=p.Service(tmp);r=s.receipt(raw);rs=s.receipt(simple);rm=s.receipt(matrix)
        record(36,'real_text_token_meter','both real encoders loaded; model mapping and API usage are explicitly unknown',counts=m.counts('Geometry proof: lower <= optimum <= upper.'))
        projected=p.summary(raw)
        record(37,'verified_decision_projection','same independently checked certificate; untrusted trace and generated-library repetitions omitted',**measure(raw,projected))
        compact=p.summary(raw,compact_numbers=True)
        for key,value in projected['bounds'].items():
            if value is not None:assert (F(compact['bounds'][key])<=F(value) if key=='lower' else F(compact['bounds'][key])>=F(value))
        record(38,'directed_bound_display','display bounds enclose exact rationals including sign; exact certificate retained',**measure(projected['bounds'],compact['bounds']))
        assert s.inspect(r['evidence'])==p.checked(raw)
        record(39,'immutable_evidence_reference','full exact evidence recovered and checked',**measure(compact,r))
        pointer='/certificate/compiled_goal/source';selected=s.inspect(r['evidence'],pointer)
        assert selected==raw['certificate']['compiled_goal']['source']
        record(40,'selective_evidence_query','selected exact field equals original field; this measures a field query, not full proof reading',**measure(p.checked(raw),selected))
        atoms=s.inspect(rm['evidence'],'/certificate/dual/atoms');first=s.inspect(rm['evidence'],'/certificate/dual/atoms',0,2)
        assert first['items']==atoms[:2] and first['total']==len(atoms) and first['next']==2
        record(41,'explicit_pagination','first-page identity and total/next preserved; reading all pages adds overhead',**measure(atoms,first))
        prof=profile_request(goal);assert p.expand_request(prof)==goal
        record(42,'explicit_geometry_profile','fixed sphere/ground/full-space declared by named profile; source restored exactly',**measure(goal,prof))
        variants=coefficient_requests(goal)
        record(43,'polynomial_coefficient_input','compiles to identical operator and original target; this simple case may cost more tokens',**measure(prof,variants['coefficients']))
        record(44,'triangular_quadratic_input','same full Hessian including doubled off-diagonal convention',**measure(prof,variants['triangular_quadratic']))
        revised=s.revise(rs['task'],{'threshold':'-1/100'});source=s.store.get(revised)
        assert source==dict(goal,threshold='-1/100')
        record(45,'goal_reference_patch','unchanged parameters, function space and domain remain bound',**measure(source,{'goal_ref':rs['task'],'changes':{'threshold':'-1/100'}}))
        refs=[rs['task']]*8
        batched=s.batch(refs,options={'max_leaves':4});assert len(batched['items'])==8
        record(46,'batch_tool_request','same ordered set of eight outcomes; response is not claimed compressed',**measure([{'op':'solve','goal_ref':ref,'options':{'max_leaves':4}} for ref in refs],{'op':'batch','goal_refs':refs,'options':{'max_leaves':4}}))
        dedup=s.batch(refs,options={'max_leaves':4},deduplicate=True)
        assert [dedup['unique'][i] for i in dedup['order']]==batched['items'];assert s.calls==1
        record(47,'deduplicated_outcome_table','lossless order reconstruction and one actual solve for eight identical requests',actual_solver_calls=s.calls,**measure(batched,dedup))
        weak=json.loads((BASE/'examples/weak_ansatz.json').read_text());weakref=s.register(weak)
        mix=refs[:4]+[weakref];full=s.batch(mix,options={'max_leaves':4});exceptions=s.batch(mix,options={'max_leaves':4},exceptions_only=True)
        assert [i['status'] for i in exceptions['inventory']]==[i['status'] for i in full['items']]
        assert len(exceptions['attention'])==1 and exceptions['attention'][0]['result']['status']=='unresolved'
        record(48,'exception_view_with_complete_inventory','every task status retained and unresolved detail present',**measure(full,exceptions))
        lemma=json.loads((BASE/'results/proofs/v34_optimized_odd_even.json').read_text());lrefs=s.add_lemma(lemma)
        fresh=p.Service(tmp);assert fresh.lemma_index()==s.lemma_index()
        record(49,'persistent_verified_lemma_index','another service restores same checked lemma; index cannot mutate stored proof',**measure(lemma,lrefs))
        grouped=p.aggregate_diagnostics(raw['diagnostics']);assert sum(x['count'] for x in grouped)==len(raw['diagnostics'])
        record(50,'scoped_diagnostic_aggregation','diagnostics validated as reachable representations of the original goal',**measure(raw['diagnostics'],grouped))
        okay=s.budget(rs,500);small=s.budget(rs,16)
        assert okay['fits'] and not small['fits'] and 'status' not in small
        record(51,'semantic_token_budget','undersized budget explicitly fails without a mathematical verdict',sufficient=okay['text_tokens'],too_small=small)
        poll=s.poll(rs,rs['evidence']);assert poll=={'unchanged':rs['evidence']}
        record(52,'checked_unchanged_poll','same evidence rechecked; previous receipt remains required context',**measure(rs,poll))
        cp=s.checkpoint([rs['task']],batched['items']);restored=p.Service(tmp).restore(cp)
        assert restored['results']==batched['items']
        record(53,'reference_backed_interaction_checkpoint','fresh service recovers exact verified receipts and library; not numerical solver resume',**measure(batched['items'],{'checkpoint':cp}))
        context=s.context(batched['items'],[weakref]);again=p.Service(tmp).restore(context['checkpoint'])
        assert again['goals']==[weakref] and again['results']==batched['items']
        record(54,'recoverable_session_context','all task bindings, unresolved work and results recoverable; arbitrary user chat is not summarized',**measure(batched['items'],context))
        submitted=s.submit(prof,options={'max_leaves':4},token_limit=500)
        packet=json.loads(submitted['text']);assert packet['status']=='proved' and s.inspect(packet['evidence'])['certificate']['status']=='proved'
        record(55,'one_call_checked_goal_workflow','registered input, same mathematical solver, exact evidence and bounded presentation',receipt_tokens=submitted['text_tokens'],fits=submitted['fits'])
    assert [r['version'] for r in records]==list(range(36,56))
    out={'actual_model_calls':0,'counting_scope':'explicit text only','experiments':records,
         'limitations':['per-feature different workloads, not one monotonic speedup curve','negative token gains retained','reference retrieval and first-use instructions cost extra','no observed model understanding or hidden reasoning usage']}
    (HERE/'results/version_token_experiments.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'versions':len(records),'all_checked':True}))

if __name__=='__main__':main()
