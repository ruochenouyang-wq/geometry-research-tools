"""V37–V55: verified, reference-backed, budgeted access to the SAME mathematics.

Compression is a presentation/interaction layer. It never alters a certificate.
Every proof acceptance and every reference read is independently reconstructed.
"""
from pathlib import Path
from fractions import Fraction as F
import copy,hashlib,json,os,re,tempfile
import compiler_v23 as compiler
import planner_v35 as planner
import planner_v32
import research_registry as verifier
import family_v29 as family
import dual_v33
import optimize_v34
import matrix_v30
import token_meter as meter


PROFILE='sphere_ground'
PROFILE_FIELDS={'geometry':'unit_sphere','function_space':'full_sphere','eigenvalue_index':1}
PROFILE_TEXT='sphere_ground = fixed unit S2, ground eigenvalue (including zero mode), full energy-form space; potentials axisymmetric. A goal reference binds the exact potential, penalty, parameters, domain and threshold. Display bounds are outward enclosures; exact evidence is retrievable. Matrix-envelope optimality is only within the stated Poisson ansatz. verified is rational replay, not a formal proof-assistant kernel.'


def checked(result):
    if not isinstance(result,dict) or not verifier.verify(result.get('certificate',result)):raise ValueError('Unverified result rejected')
    diagnostics=result.get('diagnostics',[])
    if not isinstance(diagnostics,list) or any(not planner.verify_budget(x) for x in diagnostics):raise ValueError('Unverified diagnostic rejected')
    if diagnostics:
        cert=result.get('certificate',result)
        if cert['method']!='proof_search_v32':raise ValueError('Diagnostics require an original ground-state goal')
        states,_=planner_v32.reduction_search(cert['compiled_goal']['problem'],max_states=64)
        reachable={meter.wire(p) for p,_ in states}
        if any(meter.wire(item['problem']) not in reachable for item in diagnostics):raise ValueError('Diagnostic does not concern a verified representation of this goal')
    return {'certificate':copy.deepcopy(result.get('certificate',result)),'diagnostics':copy.deepcopy(diagnostics)}


def power10(k):return F(10**k) if k>=0 else F(1,10**(-k))


def outward(value,upper=False,digits=7):
    """V38: exact directed significant-digit rounding (no binary floats)."""
    if value is None:return None
    if type(digits) is not int or not 2<=digits<=18:raise ValueError('Display precision budget')
    x=F(value)
    if not x:return '0'
    absolute=abs(x);exponent=len(str(absolute.numerator))-len(str(absolute.denominator))
    while absolute<power10(exponent):exponent-=1
    while absolute>=power10(exponent+1):exponent+=1
    shift=exponent-digits+1;scaled=x/power10(shift)
    n=(-((-scaled.numerator)//scaled.denominator)) if upper else scaled.numerator//scaled.denominator
    while n and n%10==0:n//=10;shift+=1
    if shift>=0 and len(str(abs(n)))+shift<=12:return str(n*10**shift)
    if -6<=shift<0:
        sign='-' if n<0 else '';s=str(abs(n)).rjust(1-shift,'0')
        return sign+s[:shift]+'.'+s[shift:]
    return str(n)+'e'+str(shift)


def summary(result,compact_numbers=False):
    """V37: required semantic fields + evidence, removing untrusted search verbosity."""
    verified=checked(result);c=verified['certificate'];method=c.get('method',c.get('format'))
    if method=='proof_search_v32':
        source=c['compiled_goal']['source'];kind='ground_goal'
        scope={'geometry':source['geometry'],'space':source['function_space'],'eigenvalue_index':source['eigenvalue_index'],'domain':source['domain']}
        task=source;bounds={'lower':c['lower_bound'],'upper':c['upper_bound']};target=c['threshold'];status=c['status']
    elif method=='global_matrix_envelope_v33':
        kind='matrix_envelope';task={'directions':c['template']['directions'],'cost':c['dual']['cost']}
        scope={'geometry':'unit_sphere','amplitudes':'all_real','interval':['-1','1'],'ansatz':'poisson_exponential','objective':'trace(C G)'}
        bounds={'lower':c['lower_bound'],'upper':c['upper_bound'],'gap':c['gap']};target=c['tolerance'];status=c['status']
    elif method in ('poisson_scalar_v29','poisson_matrix_v30'):
        kind='family_lemma';task={'directions':c['directions'],'weight':c['weight']}
        scope={'geometry':'unit_sphere','space':'full_sphere','amplitudes':'all_real','ansatz':'poisson_exponential'}
        bounds={'lower':c['range_proof']['maximum_lower'],'upper':c['range_proof']['maximum_upper']};target=c['range_proof']['tolerance'];status=c['range_proof']['status']
    else:
        new=verifier.summary_fields(c)
        kind,task,scope,target,bounds,status=(new[k] for k in ('kind','task','scope','target','bounds','status'))
    if compact_numbers:bounds={k:outward(v,upper=not(k=='lower' or k.endswith('_lower'))) for k,v in bounds.items()}
    out={'kind':kind,'status':status,'task':task,'scope':scope,'target':target,'bounds':bounds,
         'bound_display':'outward' if compact_numbers else 'exact','verified':True,'formal_kernel':False,'evidence':verified}
    return out


class Store:
    """V39: immutable canonical content references; never an arbitrary file reader."""
    def __init__(self,path):self.path=Path(path);self.path.mkdir(parents=True,exist_ok=True)
    def put(self,kind,value):
        if kind not in ('goal','proof','lemma','batch','checkpoint','context'):raise ValueError('Reference kind')
        data=meter.wire(value).encode();digest=hashlib.sha256(data).hexdigest();ref=kind+':'+digest[:24];path=self.path/(ref.replace(':','_')+'.json')
        if path.exists():
            if path.read_bytes()!=data:raise ValueError('Reference collision or store corruption')
            return ref
        with tempfile.NamedTemporaryFile(dir=self.path,delete=False) as f:f.write(data);temporary=f.name
        os.replace(temporary,path);return ref
    def get(self,ref):
        if not isinstance(ref,str) or not re.fullmatch(r'(goal|proof|lemma|batch|checkpoint|context):[0-9a-f]{24}',ref):raise ValueError('Invalid reference')
        path=self.path/(ref.replace(':','_')+'.json');data=path.read_bytes()
        if hashlib.sha256(data).hexdigest()[:24]!=ref.split(':')[1]:raise ValueError('Reference content changed')
        return json.loads(data)


def project(value,pointer):
    """V40: explicit JSON Pointer selection, not hidden-model or filesystem access."""
    if not isinstance(pointer,str) or (pointer and not pointer.startswith('/')):raise ValueError('JSON pointer required')
    if not pointer:return value
    for part in pointer[1:].split('/'):
        if re.search(r'~(?![01])',part):raise ValueError('Invalid JSON Pointer escape')
        part=part.replace('~1','/').replace('~0','~')
        if isinstance(value,list):
            if not re.fullmatch(r'0|[1-9][0-9]*',part):raise ValueError('Array index')
            value=value[int(part)]
        elif isinstance(value,dict):value=value[part]
        else:raise ValueError('Pointer traverses a scalar')
    return value


def page(value,offset=0,limit=8):
    """V41: deterministic pagination with total/next, never silent truncation."""
    if not isinstance(value,list) or type(offset) is not int or type(limit) is not int or offset<0 or not 1<=limit<=64:raise ValueError('Page requires array, nonnegative offset and limit 1..64')
    if offset>len(value):raise ValueError('Page offset beyond end')
    end=min(len(value),offset+limit)
    return {'items':value[offset:end],'total':len(value),'offset':offset,'next':end if end<len(value) else None}


def polynomial(raw):
    """V43: sparse exact polynomial input, without asking a model to expand it."""
    if not isinstance(raw,list) or not 1<=len(raw)<=7:raise ValueError('Polynomial needs 1..7 coefficients')
    parts=[]
    for degree,value in enumerate(raw):
        value=F(value)
        if value:parts.append('('+str(value)+')'+('' if not degree else '*t' if degree==1 else '*t**'+str(degree)))
    return '+'.join(parts) or '0'


def expand_request(raw):
    """V42–V44: explicit named scope + coefficient and triangular Hessian inputs."""
    if not isinstance(raw,dict):raise ValueError('Goal object required')
    out=copy.deepcopy(raw)
    if 'profile' in out:
        if out.pop('profile')!=PROFILE:raise ValueError('Unknown explicit geometry profile')
        for k,v in PROFILE_FIELDS.items():
            if k in out and out[k]!=v:raise ValueError('Profile conflicts with explicit scope')
            out[k]=v
    if 'potential_coefficients' in out:
        if 'potential' in out:raise ValueError('Ambiguous potential')
        packed=out.pop('potential_coefficients');names=out['parameters'];base=polynomial(packed['base']);directions=packed['directions']
        if len(directions)!=len(names):raise ValueError('Potential direction count')
        out['potential']=base+''.join('+('+name+')*('+polynomial(p)+')' for name,p in zip(names,directions))
    if 'quadratic' in out:
        if 'penalty' in out:raise ValueError('Ambiguous penalty')
        packed=out.pop('quadratic');names=out['parameters'];d=len(names);linear=packed.get('linear',[0]*d);upper=packed['hessian_upper']
        if len(linear)!=d or len(upper)!=d*(d+1)//2:raise ValueError('Packed Hessian dimension')
        pieces=['('+str(F(packed.get('constant',0)))+')']
        for name,value in zip(names,linear):
            if F(value):pieces.append('('+str(F(value))+')*'+name)
        index=0
        for i in range(d):
            for j in range(i,d):
                value=F(upper[index])/(2 if i==j else 1);index+=1
                if value:pieces.append('('+str(value)+')*'+names[i]+'*'+names[j])
        out['penalty']='+'.join(pieces)
    # Defaulting a theorem's domain or eigenvalue without an explicit profile is forbidden.
    compiler.compile_goal(out)
    return out


def aggregate_diagnostics(items):
    """V50: scoped repeated failures grouped without changing the original verdict."""
    groups={}
    for item in items:
        if not planner.verify_budget(item):raise ValueError('Invalid diagnostic')
        key=item['status'];g=groups.setdefault(key,{'status':key,'count':0,'scale_lower':None,'scale_upper':None})
        g['count']+=1;lo=F(item['scale_lower']);hi=F(item['scale_upper'])
        g['scale_lower']=lo if g['scale_lower'] is None else min(lo,g['scale_lower'])
        g['scale_upper']=hi if g['scale_upper'] is None else max(hi,g['scale_upper'])
    return [{'scope':'fixed_poisson_budget',**{k:outward(v,upper=(k=='scale_upper')) if isinstance(v,F) else v for k,v in g.items()}} for g in groups.values()]


def validate_task(raw):
    if isinstance(raw,dict) and 'geometry' in raw:return compiler.compile_goal(raw)
    if isinstance(raw,dict) and set(raw)=={'directions','cost'}:return dual_v33.setup(raw['directions'],raw['cost'])
    if isinstance(raw,dict) and set(raw)=={'directions','weight'}:return family.assemble(raw['directions'],raw['weight'])
    return verifier.validate_task(raw)


class Service:
    def __init__(self,path):
        self.store=Store(path);self.calls=0;self.cache={};self.lemmas={}
        for file in sorted(self.store.path.glob('lemma_*.json')):
            ref=file.stem.replace('_',':',1);template=self.store.get(ref)
            if not family.verify(template):raise ValueError('Unverified stored lemma')
            self.lemmas[ref]=template['directions']

    def register(self,raw):
        goal=expand_request(raw);return self.store.put('goal',goal)

    def revise(self,ref,changes):
        """V45: explicit validated top-level patches preserve all unstated conditions."""
        if not isinstance(changes,dict) or not changes or any(k not in ('potential','penalty','parameters','domain','threshold','function_space') for k in changes):raise ValueError('Unsupported goal patch')
        raw=self.store.get(ref)
        if not ref.startswith('goal:'):raise ValueError('Goal reference required')
        raw.update(copy.deepcopy(changes));return self.register(raw)

    def receipt(self,result,goal_ref=None,diagnostics='summary'):
        out=summary(result,compact_numbers=True);proof_ref=self.store.put('proof',checked(result))
        if goal_ref is None:goal_ref=self.store.put('goal',out['task'])
        elif not isinstance(goal_ref,str) or not goal_ref.startswith('goal:') or self.store.get(goal_ref)!=out['task']:raise ValueError('Result not bound to the supplied original goal')
        out['task']=goal_ref;out['evidence']=proof_ref
        if diagnostics=='summary':out['diagnostics']=aggregate_diagnostics(result.get('diagnostics',[]))
        elif diagnostics!='omit':raise ValueError('Unknown diagnostic view')
        return out

    def _checked_receipt(self,receipt):
        if not isinstance(receipt,dict) or not str(receipt.get('evidence','')).startswith('proof:'):raise ValueError('Proof receipt required')
        source=checked(self.store.get(receipt['evidence']))
        mode='summary' if 'diagnostics' in receipt else 'omit'
        truth=self.receipt(source,receipt['task'],diagnostics=mode)
        if receipt.get('diagnostics')=='available_in_evidence':truth['diagnostics']='available_in_evidence'
        if receipt!=truth:raise ValueError('Receipt changed since verification')
        return source,truth

    def inspect(self,ref,pointer='',offset=None,limit=8):
        value=self.store.get(ref)
        if ref.startswith('proof:'):checked(value)
        if ref.startswith('lemma:') and not family.verify(value):raise ValueError('Unverified lemma')
        value=project(value,pointer)
        return page(value,offset,limit) if offset is not None else copy.deepcopy(value)

    def solve(self,ref,options=None):
        if not isinstance(ref,str) or not ref.startswith('goal:'):raise ValueError('Registered goal required')
        raw=self.store.get(ref);compiler.compile_goal(raw);options=dict(options or {})
        if any(k not in ('max_leaves','max_modes','max_states','synthesize','use_penalty_weights') for k in options):raise ValueError('Unsupported solver option')
        cache_key=meter.wire([ref,options,sorted(self.lemmas)])
        # V47: identical requests with identical options reuse a CHECKED result.
        if cache_key in self.cache:
            receipt=copy.deepcopy(self.cache[cache_key]);truth=self.receipt(checked(self.store.get(receipt['evidence'])),ref)
            if receipt!=truth:raise ValueError('Cached receipt changed')
            return truth
        library=[self.store.get(x) for x in self.lemmas]
        result=planner.solve(raw,library=library,**options);self.calls+=1
        receipt=self.receipt(result,ref);self.cache[cache_key]=copy.deepcopy(receipt);return receipt

    def batch(self,refs,options=None,exceptions_only=False,deduplicate=False):
        """V46/V48: one invocation, explicit item identity, optional success manifest."""
        if not isinstance(refs,list) or not 1<=len(refs)<=32:raise ValueError('Batch size 1..32')
        receipts=[self.solve(ref,options) for ref in refs];batch_ref=self.store.put('batch',receipts)
        if deduplicate and not exceptions_only:
            unique=[];positions={};order=[]
            for receipt in receipts:
                key=meter.wire(receipt)
                if key not in positions:positions[key]=len(unique);unique.append(receipt)
                order.append(positions[key])
            return {'batch':batch_ref,'unique':unique,'order':order}
        if not exceptions_only:return {'batch':batch_ref,'items':receipts}
        inventory=[{'index':i,'task':r['task'],'status':r['status'],'evidence':r['evidence']} for i,r in enumerate(receipts)]
        return {'batch':batch_ref,'inventory':inventory,'attention':[{'index':i,'result':r} for i,r in enumerate(receipts) if r['status']!='proved']}

    def add_lemma(self,raw):
        """V49: persistent checked lemma references instead of repeated proof text."""
        templates=verifier.load_library(raw);refs=[]
        for t in templates:
            ref=self.store.put('lemma',t);self.lemmas[ref]=copy.deepcopy(t['directions']);refs.append(ref)
        return refs

    def lemma_index(self):return [{'ref':ref,'directions':copy.deepcopy(value)} for ref,value in sorted(self.lemmas.items())]

    def budget(self,receipt,limit=500,encoding='o200k_base'):
        """V51: never truncate a theorem/status to hit a text-token budget."""
        if type(limit) is not int or not 16<=limit<=100000:raise ValueError('Token budget 16..100000')
        source,truth=self._checked_receipt(receipt)
        def wrapped(value):
            text=meter.wire(value)
            out={'text':text,'text_tokens':meter.count(text,encoding),'encoding':encoding,'fits':True,'wire_tokens':0}
            for _ in range(8):
                actual=meter.count(meter.wire(out),encoding)
                if actual==out['wire_tokens']:return out
                out['wire_tokens']=actual
            raise ValueError('Token accounting did not stabilize')
        full=wrapped(receipt)
        if full['wire_tokens']<=limit:return full
        # Diagnostics remain available under the proof reference; primary semantic fields stay.
        concise=copy.deepcopy(receipt);concise['diagnostics']='available_in_evidence'
        short=wrapped(concise)
        if short['wire_tokens']<=limit:return short
        required=min((full,short),key=lambda x:x['wire_tokens'])
        return {'fits':False,'required_text_tokens':required['text_tokens'],'required_wire_tokens':required['wire_tokens'],'encoding':encoding,'error':'budget_too_small_for_required_semantics'}

    def poll(self,receipt,known_evidence):
        """V52: explicit not-modified response, only after rechecking stored evidence."""
        source,truth=self._checked_receipt(receipt)
        return {'unchanged':receipt['evidence']} if known_evidence==receipt['evidence'] else truth

    def checkpoint(self,goals,receipts):
        """V53: resume the interaction, not an invented numerical solver state."""
        if not isinstance(goals,list) or not isinstance(receipts,list):raise ValueError('Checkpoint lists required')
        for ref in goals:
            if not isinstance(ref,str) or not ref.startswith('goal:'):raise ValueError('Goal reference required')
            validate_task(self.store.get(ref))
        evidence=[]
        for receipt in receipts:
            source,truth=self._checked_receipt(receipt)
            evidence.append(receipt['evidence'])
        return self.store.put('checkpoint',{'goals':goals,'evidence':evidence,'lemmas':sorted(self.lemmas)})

    def restore(self,ref):
        if not isinstance(ref,str) or not ref.startswith('checkpoint:'):raise ValueError('Checkpoint reference required')
        saved=self.store.get(ref)
        for goal in saved['goals']:
            if not isinstance(goal,str) or not goal.startswith('goal:'):raise ValueError('Goal reference required')
            validate_task(self.store.get(goal))
        for lemma in saved['lemmas']:
            if not isinstance(lemma,str) or not lemma.startswith('lemma:'):raise ValueError('Lemma reference required')
            template=self.store.get(lemma)
            if not family.verify(template):raise ValueError('Changed checkpoint lemma')
            self.lemmas[lemma]=template['directions']
        for proof in saved['evidence']:
            if not isinstance(proof,str) or not proof.startswith('proof:'):raise ValueError('Proof reference required')
        return {'goals':saved['goals'],'results':[self.receipt(checked(self.store.get(proof))) for proof in saved['evidence']],
                'lemmas':self.lemma_index()}

    def context(self,receipts,pending):
        """V54: reference-backed mathematical session state, not arbitrary chat summarization."""
        checkpoint=self.checkpoint(pending,receipts)
        rows=[{'task':r['task'],'status':r['status'],'evidence':r['evidence']} for r in receipts]
        return {'available_input_profile':PROFILE,'scope':'per_task','checkpoint':checkpoint,'completed':rows,'pending':pending,
                'restore':'restore(checkpoint) returns verified bounds, exact task bindings and lemma references'}

    def submit(self,raw,options=None,token_limit=500):
        """V55: one goal-to-verified-receipt call, with a complete text budget check."""
        ref=self.register(raw);receipt=self.solve(ref,options);return self.budget(receipt,token_limit)

    def optimize(self,spec,tolerance='1/100000',max_rounds=8,max_range_leaves=128,token_limit=500):
        result=optimize_v34.solve(spec['directions'],spec.get('cost'),F(tolerance),max_rounds,max_range_leaves)
        self.calls+=1;return self.budget(self.receipt(result),token_limit)

    def lemma(self,spec,token_limit=500):
        result={'certificate':matrix_v30.synthesize(spec['directions'],spec.get('weight'))}
        self.calls+=1;return self.budget(self.receipt(result),token_limit)

    def extrema(self,spec,tolerance='1/1000000000000',max_nodes=32,token_limit=500):
        import exact_extrema
        if set(spec)-{'polynomial','interval'}:raise ValueError('Unknown extrema field')
        cert=exact_extrema.maximize(spec['polynomial'],spec.get('interval',[-1,1]),F(tolerance),max_nodes)
        self.calls+=1;return self.budget(self.receipt(cert),token_limit)

    def positive(self,spec,options=None,token_limit=500):
        import positive_search
        if set(spec)-{'potential','basis'}:raise ValueError('Unknown positive-function field')
        options=dict(options or {})
        allowed={'degree','symmetry','tolerance','range_tolerance','max_leaves','grid_size','max_exchanges','max_newton','range_backend'}
        if set(options)-allowed:raise ValueError('Unsupported positive search option')
        chosen=verifier.backend(options.pop('range_backend',None))
        for key in ('tolerance','range_tolerance'):
            if key in options:options[key]=F(options[key])
        result=positive_search.search(spec['potential'],basis=spec.get('basis'),range_backend=chosen,**options)
        self.calls+=1;return self.budget(self.receipt(result),token_limit)

    def uniform(self,spec,rule='log_sobolev',token_limit=500):
        import uniform_refine
        spec=dict(spec)
        if rule=='log_sobolev':
            if set(spec)-{'direction','coefficient'}:raise ValueError('Log-Sobolev rule supports only affine direction and coefficient')
            cert=uniform_refine.synthesize_log_sobolev(**spec)
        elif rule=='barta':
            allowed={'direction','coefficient','amplitude_box','order','max_leaves','symmetry','all_real'}
            if set(spec)-allowed:raise ValueError('Unknown Barta synthesis field')
            cert=uniform_refine.synthesize(**spec)
        else:raise ValueError('Unknown uniform inequality rule')
        self.calls+=1;return self.budget(self.receipt(cert),token_limit)

    def refine(self,spec,tolerance='1/1000000',budget=8,max_nodes=32,token_limit=500):
        import matrix_refine
        import precision_bridge
        if set(spec)-{'directions','cost'}:raise ValueError('Unknown matrix field')
        result=matrix_refine.refine(spec['directions'],spec.get('cost'),F(tolerance),budget=budget)
        cert=precision_bridge.precise_matrix(result,F(tolerance),max_nodes)
        self.calls+=1;return self.budget(self.receipt(cert),token_limit)

    def research(self,raw,rule='log_sobolev',coefficient=None,token_limit=500):
        """V89/V90: apply a verified uniform rule to the ORIGINAL compiled goal."""
        import original_goal
        ref=self.register(raw)
        cert=original_goal.solve(self.store.get(ref),coefficient=coefficient,rule=rule)
        self.calls+=1;return self.budget(self.receipt(cert,ref),token_limit)
