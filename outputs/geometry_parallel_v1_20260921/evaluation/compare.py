"""Public development comparison of three independent variants; no seed/model IO.

Cold task deadlines cover the entire child process. Stream phases have parent
watchdogs and every returned query proof also receives a fresh-process replay.
"""
from copy import deepcopy
from datetime import datetime, timezone
from fractions import Fraction as F
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import resource
import selectors
import signal
import statistics
import subprocess
import sys
import time
import zipfile

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'evaluation'
SYSTEMS=('baseline','A','B','C')
PATHS={'baseline':ROOT/'shared_baseline.py','A':ROOT/'a_multispectrum/variant.py',
       'B':ROOT/'b_parameter/variant.py','C':ROOT/'c_adaptive/variant.py'}
TASK_SECONDS=5.0
PREPARATION_SECONDS=30.0
EXTENSION_SECONDS=10.0
PROTOCOL_VERSION='parallel_public_development_v1.1'


def canonical(x):
    return json.dumps(x,sort_keys=True,ensure_ascii=False,separators=(',',':'),allow_nan=False)


def safe(path):
    p=Path(path)
    if any(s in p.name.lower() for s in ('private','seed')):
        raise ValueError('Private and seed files are outside this evaluator')
    return p


def read(path):
    return json.loads(safe(path).read_text())


def write(path,x):
    p=safe(path);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(canonical(x)+'\n');tmp.replace(p)


def sha(path):
    return hashlib.sha256(safe(path).read_bytes()).hexdigest()


def digest(x):
    return hashlib.sha256(canonical(x).encode()).hexdigest()


def q(profile='step',amplitude='-1',exponent=None,offset='0'):
    value={'kind':'axis_profile','profile':profile,'amplitude':str(F(amplitude)),
           'offset':str(F(offset)),'axis':2}
    if exponent is not None:value['exponent']=str(F(exponent))
    return value


def task(kind,function,tolerance='1/100000000',mean_zero=None,seconds=5):
    value={'kind':kind,'function':function,'tolerance':tolerance,'budget':{'wall_seconds':seconds}}
    if mean_zero is not None:value['mean_zero']=mean_zero
    return value


def public_inputs():
    common=[]
    for scope in (False,True):
        common.append(task('spectrum',q(),mean_zero=scope))
    for exponent,amplitude in (('-1/5','-1/3'),('-3/7','-3/5')):
        for scope in (False,True):
            common.append(task('spectrum',q('abs_power',amplitude,exponent),mean_zero=scope))
    common.extend([task('approximation',q('abs_power','-1/2','-1/3')),
                   task('approximation',q('abs_power','-2/3','-5/12'))])
    common=[{'id':'D'+str(i+1).zfill(2),'task':t,'classification':'public_development_not_holdout'} for i,t in enumerate(common)]
    family={'profile':'abs_power','exponent':'-1/4','axis':2,'amplitude_interval':['-1/2','-49/100'],
            'offset':'0','mean_zero':False,'tolerance':'1/1000000'}
    stream=[{'id':'Q'+str(j).zfill(2),'task':task('spectrum',q('abs_power',str(F(-1,2)+F(j,1700)),'-1/4'),
                tolerance='1/1000000',mean_zero=False),'classification':'public_parameter_stream_not_holdout'} for j in range(1,17)]
    extensions=[]
    for name,offset,scope,count,total in [('E01','-6',False,4,'-18'),('E02','-6',True,3,'-12'),('E03','-7',False,9,'-27'),('E04','-7',True,8,'-20')]:
        extensions.append({'id':name,'task':task('negative_spectrum',q(amplitude='0',offset=offset),
            '1/1000000',scope,10),'oracle':{'negative_count':count,'negative_sum':total},'classification':'analytic_multiplicity_regression'})
    extensions.append({'id':'E05','task':task('negative_spectrum',q(offset='-4'), '1/1000000',False,10),
                       'oracle':{'negative_count':4},'classification':'nonconstant_extension_capability_only'})
    extensions.append({'id':'E06','task':task('negative_spectrum',q('abs_power','-1/5','-1/4',offset='-6'), '1/1000000',True,10),
                       'oracle':{'negative_count':8},'classification':'singular_extension_capability_only'})
    return {'common_cases.json':common,'family_spec.json':family,'stream_cases.json':stream,'extension_cases.json':extensions}


def protocol():
    return {'version':PROTOCOL_VERSION,'classification':'public_development_only','systems':list(SYSTEMS),
        'common':{'cases':8,'repetitions':2,'tolerance':'1/100000000','entire_child_process_seconds':5},
        'stream':{'queries':16,'sessions':2,'preparation_seconds':30,'each_query_seconds':5,
                  'checkpoints':[1,4,16],'family_spec_file':'family_spec.json',
                  'order':'a=-1/2+j/1700, j=1..16','B_context':'prepare, complete family verify, prepare_context, query plus context.verify',
                  'independent_cold_replay':'every query proof in its own process, five seconds; separately measured and cumulative cost also added'},
        'extensions':{'A_cases':6,'entire_child_process_seconds':10,'negative_sum':'signed sum of negative eigenvalues',
                      'count_and_sum_reported_separately':True,'never_added_to_common_success_rate':True},
        'acceptance':'Only task-bound independent verification and exact target attainment within the declared deadline count. Fallback is identified, not credited as new ability.',
        'ordering':'serial; rotate systems per repeat/session; rotate common case order; no concurrent speed tests',
        'costs':'CPU/wall, complete JSON, proof replay, disk output, every failed attempt; B full preparation and initial/context verification charged at N=1/4/16',
        'speed_claim':'two repeats are preliminary; only same-success paired cost savings >=20% warrant further validation',
        'actual_model_tokens':None,'actual_model_cost':None,'network_or_model_calls':0,
        'holdout':'No holdout exists here. Any future unseen test needs a separate precommitted protocol, frozen code, prior-public-identity exclusion and authorized reveal. B001 is permanently public.'}


def initialize():
    DATA.mkdir(parents=True,exist_ok=True)
    for name,value in {**public_inputs(),'protocol.json':protocol()}.items():
        path=DATA/name
        if path.exists() and read(path)!=value:raise RuntimeError('Published inputs changed; require an explicit protocol revision: '+name)
        if not path.exists():write(path,value)
    return {'initialized':True,'formal_holdout':False}


def identities():
    sources={}
    for directory in sorted(ROOT.parent.glob('geometry_*')):
        if directory.is_dir():
            for path in sorted(directory.rglob('*.py')):
                sources[str(path.resolve())]=sha(path)
    artifacts={str((DATA/name).resolve()):sha(DATA/name) for name in list(public_inputs())+['protocol.json']}
    for path in sorted((DATA/'development_protocol_v1').glob('*.json')):
        artifacts[str(path.resolve())]=sha(path)
    artifacts[str((ROOT/'PLAN.md').resolve())]=sha(ROOT/'PLAN.md')
    if (ROOT/'BASELINE.json').exists():artifacts[str((ROOT/'BASELINE.json').resolve())]=sha(ROOT/'BASELINE.json')
    return {'sources':sources,'artifacts':artifacts}


def baseline_preserved():
    declared=read(ROOT/'BASELINE.json')['files'];bad=[]
    for name,expected in declared.items():
        path=ROOT.parents[1]/name
        if not path.is_file() or sha(path)!=expected:bad.append(name)
    return {'files_checked':len(declared),'all_unchanged':not bad,'changed_or_missing':bad}


def freeze():
    initialize()
    if not baseline_preserved()['all_unchanged']:raise RuntimeError('Read-only baseline changed')
    for path in PATHS.values():
        if not path.is_file():raise RuntimeError('Variant not ready: '+str(path))
    identity=identities();stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    archive=DATA/('source_snapshot_'+stamp+'.zip')
    with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED) as z:
        index={}
        for i,(path,expected) in enumerate(sorted({**identity['sources'],**identity['artifacts']}.items())):
            data=safe(path).read_bytes()
            if hashlib.sha256(data).hexdigest()!=expected:raise RuntimeError('Source changed during snapshot')
            member='files/'+str(i).zfill(6)+'.blob';z.writestr(member,data)
            index[path]={'member':member,'sha256':expected}
        z.writestr('index.json',canonical(index))
    result={'version':PROTOCOL_VERSION,'created_utc':stamp,'identity':identity,'source_archive':str(archive),
            'source_archive_sha256':sha(archive),'classification':'frozen_public_development_comparison_not_holdout'}
    path=DATA/('freeze_'+stamp+'.json');write(path,result)
    pointer={'manifest':str(path),'sha256':sha(path)};write(DATA/'LATEST_FREEZE.json',pointer)
    return pointer


def checked_freeze():
    pointer=read(DATA/'LATEST_FREEZE.json');frozen=read(pointer['manifest'])
    if sha(pointer['manifest'])!=pointer['sha256'] or frozen['identity']!=identities():
        raise RuntimeError('Code or published tasks changed after freeze; final comparison rejected')
    if sha(frozen['source_archive'])!=frozen['source_archive_sha256']:raise RuntimeError('Source archive changed')
    if not baseline_preserved()['all_unchanged']:raise RuntimeError('Read-only baseline changed')
    return pointer


def load_system(system):
    path=PATHS[system]
    sys.dont_write_bytecode=True
    sys.path[:0]=[str(path.parent),str(ROOT)]
    spec=importlib.util.spec_from_file_location('parallel_variant_'+system,path)
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    return module


class HardDeadline(BaseException):pass


def bounded(function,seconds):
    def alarm(*_):raise HardDeadline('worker deadline exceeded')
    previous=signal.signal(signal.SIGALRM,alarm)
    start=time.monotonic();cpu=time.process_time()
    try:
        signal.setitimer(signal.ITIMER_REAL,max(seconds,0.000001))
        return {'ok':True,'value':function(),'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu}
    except (Exception,HardDeadline) as error:
        return {'ok':False,'error_type':type(error).__name__,'reason':str(error),'timed_out':isinstance(error,HardDeadline),
                'wall_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu}
    finally:
        signal.setitimer(signal.ITIMER_REAL,0);signal.signal(signal.SIGALRM,previous)


def certificate(response):
    return response.get('certificate',response if 'format' in response else None) if type(response) is dict else None


def rational(value):
    if type(value) not in (int,str):raise ValueError('Exact rational required')
    return F(value)


def assess(cert,item,verifier):
    t=item['task'];invalid={'certificate_valid':False,'target_met':False,'reason':'missing_or_invalid_certificate'}
    if type(cert) is not dict:return invalid
    q0=deepcopy(t['function'])
    if cert.get('function')!=q0:return {**invalid,'reason':'original_function_mismatch'}
    if t['kind']!='approximation' and cert.get('mean_zero') is not t['mean_zero']:
        return {**invalid,'reason':'original_space_mismatch'}
    if 'tolerance' in cert and rational(cert['tolerance'])!=rational(t['tolerance']):
        return {**invalid,'reason':'original_tolerance_mismatch'}
    verdict=verifier(deepcopy(cert),deepcopy(t))
    if type(verdict) is not dict or verdict.get('certificate_valid') is not True:
        return {**invalid,'reason':'independent_verifier_rejected','verifier_result':verdict}
    if t['kind']=='negative_spectrum':
        lo=verdict.get('sum_lower',cert.get('sum_lower'));hi=verdict.get('sum_upper',cert.get('sum_upper'))
        count=verdict.get('negative_count',cert.get('negative_count'))
        count_valid=verdict.get('count_certified') is True and type(count) is int and count>=0
        width=rational(hi)-rational(lo) if lo is not None and hi is not None else None
        sum_met=width is not None and 0<=width<=rational(t['tolerance'])
        oracle=item.get('oracle',{})
        if 'negative_count' in oracle and count_valid and count!=oracle['negative_count']:
            return {**invalid,'reason':'independent_count_oracle_mismatch'}
        if 'negative_sum' in oracle and lo is not None and not rational(lo)<=rational(oracle['negative_sum'])<=rational(hi):
            return {**invalid,'reason':'independent_sum_oracle_mismatch'}
        return {'certificate_valid':True,'count_certified':count_valid,'negative_count':count,'sum_target_met':sum_met,
                'target_met':count_valid and sum_met,'sum_lower':lo,'sum_upper':hi,'value':str(width) if width is not None else None,
                'quantity':'signed_negative_spectrum_sum_width','verifier_result':verdict}
    value=rational(cert['upper'])-rational(cert['lower']) if t['kind']=='spectrum' else rational(cert['error_upper'])
    if value<0:return {**invalid,'reason':'negative_bound'}
    if t['kind']=='approximation' and 'error_squared' in cert:
        if not 0<=rational(cert['error_squared'])<=value*value:return {**invalid,'reason':'invalid_L2_enclosure'}
    return {'certificate_valid':True,'target_met':value<=rational(t['tolerance']),'value':str(value),
            'quantity':'spectral_width' if t['kind']=='spectrum' else 'function_L2_error_upper','verifier_result':verdict}


def fallback_marker(response):
    return {'explicit_fallback':bool(response.get('fallback') or response.get('baseline_fallback') or response.get('counters',{}).get('baseline_fallback_calls')),
            'counters_baseline_fallback_calls':response.get('counters',{}).get('baseline_fallback_calls'),
            'route':response.get('route'),'method':response.get('method')} if type(response) is dict else {}


def worker(spec):
    event_stream=sys.stdout
    def event(x):event_stream.write(canonical(x)+'\n');event_stream.flush()
    # Native print output from mathematical modules cannot spoof protocol events.
    log=open(spec['log'],'w');sys.stdout=log;sys.stderr=log
    if spec['mode']=='selftest_timeout':
        while True:
            try:time.sleep(5)
            except BaseException:pass
    module=load_system(spec['system'])
    if spec['mode'] in ('single','replay','attack','family_verify'):
        row={'system':spec['system'],'case_id':spec.get('case',{}).get('id'),'mode':spec['mode'],'success':False}
        def work():
            start=time.monotonic();cpu=time.process_time()
            if spec['mode']=='family_verify':
                row['verifier_result']=module.verify_family(spec['bank'],expected_spec=spec['family'])
                row['success']=row['verifier_result'].get('certificate_valid') is True
            else:
                item=spec['case']
                if spec['mode']=='single':
                    if item['task']['kind']=='negative_spectrum' and spec['system']!='A':
                        response={'certificate':None,'unsupported':True,'reason':'baseline task language does not include negative_spectrum'}
                    else:response=module.solve(deepcopy(item['task']))
                    row['response']=json.loads(canonical(response));row['fallback']=fallback_marker(response)
                    cert=certificate(row['response'])
                else:cert=spec['certificate'];row['certificate']=cert
                if spec['mode']=='attack':
                    # Test the supplied verifier itself as well as evaluator binding.
                    row['direct_verifier_result']=module.verify(deepcopy(cert),deepcopy(item['task']))
                    write(spec['output'],row)
                row['assessment']=assess(cert,item,module.verify)
                row['success']=row['assessment'].get('target_met') is True
            row['operation_wall_seconds']=time.monotonic()-start;row['operation_cpu_seconds']=time.process_time()-cpu
            write(spec['output'],row)
        execution=bounded(work,spec['limit'])
        if not execution['ok']:
            row.update(execution=execution,success=False);write(spec['output'],row)
        return
    if spec['mode']!='stream':raise ValueError('Unknown worker mode')
    prep={'system':spec['system'],'phase':'preparation','ok':False}
    state={}
    def prepare():
        if spec['system']=='B':
            prep['response']=json.loads(canonical(module.prepare_family(deepcopy(spec['family']),preparation_seconds=30)))
            bank=prep['response'].get('bank')
            if bank is None:raise ValueError('Preparation returned no proof bank')
            prep['verification']=module.verify_family(deepcopy(bank),expected_spec=deepcopy(spec['family']))
            if prep['verification'].get('certificate_valid') is not True or prep['verification'].get('target_met') is not True:
                raise ValueError('Uniform family proof did not independently verify and meet target')
            state['context']=module.prepare_context(deepcopy(bank));state['bank']=bank
        else:prep['response']={'preparation_required':False}
        prep['ok']=True;write(spec['preparation_output'],prep)
    stamp=time.monotonic();event({'event':'phase_started','name':'preparation','monotonic':stamp,'cpu_cumulative':time.process_time(),'limit':30})
    execution=bounded(prepare,30);prep['execution']=execution if not execution['ok'] else {k:v for k,v in execution.items() if k!='value'}
    write(spec['preparation_output'],prep)
    event({'event':'phase_finished','name':'preparation','monotonic':time.monotonic(),'cpu_cumulative':time.process_time(),'ok':execution['ok']})
    if not execution['ok']:return
    for item in spec['cases']:
        row={'system':spec['system'],'case_id':item['id'],'task':item['task'],'success':False}
        def query():
            if spec['system']=='B':
                raw=state['context'].query(deepcopy(item['task']));verify=state['context'].verify
                row['verification_mode']='prepared_context_plus_later_independent_cold_full_replay'
            else:
                raw=module.solve(deepcopy(item['task']));verify=module.verify
                row['verification_mode']='full_same_session_plus_later_independent_cold_full_replay'
            row['response']=json.loads(canonical(raw));row['fallback']=fallback_marker(raw)
            row['assessment']=assess(certificate(row['response']),item,verify)
            row['success']=row['assessment'].get('target_met') is True
            write(Path(spec['queries'])/(item['id']+'.json'),row)
        stamp=time.monotonic();event({'event':'phase_started','name':item['id'],'monotonic':stamp,'cpu_cumulative':time.process_time(),'limit':5})
        execution=bounded(query,5)
        if not execution['ok']:row.update(success=False,execution=execution);write(Path(spec['queries'])/(item['id']+'.json'),row)
        event({'event':'phase_finished','name':item['id'],'monotonic':time.monotonic(),'cpu_cumulative':time.process_time(),'ok':execution['ok']})
    if spec['system']=='B':
        outside=deepcopy(spec['cases'][0]['task']);outside['function']['amplitude']='-51/100'
        def outside_probe():
            raw=state['context'].query(outside)
            return {'returned':True,'response':json.loads(canonical(raw)),
                    'accepted':assess(certificate(raw),{'task':outside},state['context'].verify)['certificate_valid']}
        result=bounded(outside_probe,5)
        write(spec['outside_output'],{'probe':result,'rejected':
            (result['ok'] and result.get('value',{}).get('accepted') is not True) or
            (not result['ok'] and result.get('error_type') in ('ValueError','TypeError','KeyError')),
            'inconclusive':not result['ok'] and result.get('error_type') not in ('ValueError','TypeError','KeyError')})


def run_process(spec,destination):
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=False)
    spec={**spec,'output':str(destination/'output.json'),'log':str(destination/'worker.log'),
          'preparation_output':str(destination/'preparation.json'),'queries':str(destination/'queries'),
          'outside_output':str(destination/'outside_probe.json')}
    write(destination/'spec.json',spec)
    start=time.monotonic();usage=resource.getrusage(resource.RUSAGE_CHILDREN)
    events=[];errors=[];buffer=b'';watchdog=None;mode=spec['mode'];phase_index=0
    sequence=['preparation']+[item['id'] for item in spec.get('cases',[])]
    with open(destination/'stderr.log','wb') as log:
        child=subprocess.Popen([sys.executable,'-B',str(Path(__file__).resolve()),'worker','--spec',str(destination/'spec.json')],
            cwd=str(ROOT),stdout=subprocess.PIPE,stderr=log,start_new_session=True,
            env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',PYTHONHASHSEED='0'))
        selector=selectors.DefaultSelector();selector.register(child.stdout,selectors.EVENT_READ)
        deadline=start+(5 if mode=='stream' else spec['limit']);phase='initialization' if mode=='stream' else 'entire_process'
        while True:
            left=deadline-time.monotonic()
            if left<=0:
                watchdog=phase
                try:os.killpg(child.pid,signal.SIGKILL)
                except ProcessLookupError:pass
                break
            for key,_ in selector.select(min(left,.05)):
                chunk=os.read(key.fileobj.fileno(),65536)
                if not chunk:selector.unregister(key.fileobj);continue
                buffer+=chunk
                while b'\n' in buffer:
                    line,buffer=buffer.split(b'\n',1)
                    try:
                        msg=json.loads(line);events.append(msg)
                        if mode!='stream':raise ValueError('Unexpected event')
                        stamp=msg['monotonic']
                        if not start<=stamp<=time.monotonic():raise ValueError('Invalid timestamp')
                        if msg['event']=='phase_started':
                            if phase not in ('initialization','between_phases'):raise ValueError('Nested phase')
                            if phase_index>=len(sequence) or msg['name']!=sequence[phase_index]:raise ValueError('Changed phase order or denominator')
                            phase_index+=1
                            expected=30 if msg['name']=='preparation' else 5
                            if msg['limit']!=expected:raise ValueError('Changed phase budget')
                            phase=msg['name'];deadline=stamp+expected
                        elif msg['event']=='phase_finished':
                            if msg['name']!=phase:raise ValueError('Unexpected finished phase')
                            if stamp>deadline and msg['ok']:errors.append('completed_after_phase_deadline')
                            phase='between_phases';deadline=time.monotonic()+5
                        else:raise ValueError('Unknown event')
                    except (ValueError,KeyError,TypeError) as exc:errors.append(str(exc))
            if child.poll() is not None and not selector.get_map():break
        selector.close();child.wait();child.stdout.close()
    elapsed=time.monotonic()-start;after=resource.getrusage(resource.RUSAGE_CHILDREN)
    process={'wall_seconds':elapsed,'cpu_seconds':after.ru_utime+after.ru_stime-usage.ru_utime-usage.ru_stime,
             'returncode':child.returncode,'watchdog':watchdog,'protocol_errors':errors,'events':events,
             'unparsed_stdout_bytes':len(buffer),'parent_start_monotonic':start,'cold_process':True}
    if mode=='stream':
        prep=read(spec['preparation_output']) if Path(spec['preparation_output']).exists() else {'ok':False,'missing':True}
        rows=[]
        for item in spec['cases']:
            p=Path(spec['queries'])/(item['id']+'.json')
            row=read(p) if p.exists() else {'system':spec['system'],'case_id':item['id'],'task':item['task'],
                                           'success':False,'reason':'not_returned_or_preparation_failed'}
            begins=[e for e in events if e['event']=='phase_started' and e['name']==item['id']]
            ends=[e for e in events if e['event']=='phase_finished' and e['name']==item['id']]
            complete=len(begins)==len(ends)==1 and ends[0]['ok'] and ends[0]['monotonic']-begins[0]['monotonic']<=5 and not errors
            row['success']=row['success'] and complete
            row['phase_wall_seconds']=ends[0]['monotonic']-begins[0]['monotonic'] if ends and begins else None
            row['phase_cpu_seconds']=ends[0]['cpu_cumulative']-begins[0]['cpu_cumulative'] if ends and begins else None
            row['cumulative_wall_seconds']=ends[0]['monotonic']-start if ends else None
            row['cumulative_cpu_seconds']=ends[0]['cpu_cumulative'] if ends else None
            row['raw_path']=str(p);rows.append(row)
        result={'system':spec['system'],'preparation':prep,'rows':rows,'process':process,
                'session_complete':child.returncode==0 and watchdog is None and not errors and
                    phase_index==len(sequence) and len([e for e in events if e['event']=='phase_finished'])==len(sequence),
                'outside_probe':read(spec['outside_output']) if Path(spec['outside_output']).exists() else None}
    else:
        result=read(spec['output']) if Path(spec['output']).exists() else {'success':False,'reason':'worker_failed_without_output'}
        valid_exit=child.returncode==0 and watchdog is None and not errors and elapsed<=spec['limit']
        result.update(system=spec['system'],case_id=spec.get('case',{}).get('id'),process=process,
                      within_budget=valid_exit,success=bool(result.get('success') and valid_exit))
    result['raw_directory']=str(destination)
    for key in ('repetition','session'):
        if key in spec:result[key]=spec[key]
    write(destination/'result.json',result)
    return result


def dist(values):
    return {'raw':values,'min':min(values),'median':statistics.median(values),'max':max(values)} if values else None


def common_summary(rows,cases,repetitions,methods):
    result={}
    for system in methods:
        selected=[r for r in rows if r['system']==system];case_result=[]
        for item in cases:
            runs=[r for r in selected if r['case_id']==item['id']]
            case_result.append({'case_id':item['id'],'kind':item['task']['kind'],'runs':len(runs),
                'all_repeats_success':len(runs)==repetitions and len({r['repetition'] for r in runs})==repetitions and all(r['success'] for r in runs),
                'successful_runs':sum(r['success'] for r in runs),
                'wall_seconds':dist([r['process']['wall_seconds'] for r in runs]),
                'cpu_seconds':dist([r['process']['cpu_seconds'] for r in runs]),
                'fallback_markers':[r.get('fallback') for r in runs]})
        result[system]={'cases':len(cases),'successful_runs':sum(r['success'] for r in selected),
                        'all_repeat_successes':sum(c['all_repeats_success'] for c in case_result),'per_case':case_result,
                        'total_wall_seconds':sum(r['process']['wall_seconds'] for r in selected),
                        'total_cpu_seconds':sum(r['process']['cpu_seconds'] for r in selected),
                        'preliminary_only':True}
        result[system]['by_kind']={kind:{'cases':sum(c['kind']==kind for c in case_result),
            'all_repeat_successes':sum(c['kind']==kind and c['all_repeats_success'] for c in case_result)}
            for kind in sorted({c['kind'] for c in case_result})}
    if 'baseline' in result:
        base={r['case_id']:r for r in result['baseline']['per_case']}
        for system in methods:
            if system=='baseline':continue
            matched=[r for r in result[system]['per_case'] if r['all_repeats_success'] and base[r['case_id']]['all_repeats_success']]
            result[system]['same_success_paired']=[{'case_id':r['case_id'],
                'cpu_ratio':r['cpu_seconds']['median']/base[r['case_id']]['cpu_seconds']['median'],
                'wall_ratio':r['wall_seconds']['median']/base[r['case_id']]['wall_seconds']['median']} for r in matched]
    return result


def mutations(cert,kind):
    result=[('status_only',{'status':'target_met','certificate_valid':True,'target_met':True,
                          'function':deepcopy(cert['function']),'mean_zero':cert.get('mean_zero'),
                          'format':cert.get('format')})]
    wrong=deepcopy(cert);wrong['function']['offset']=str(F(wrong['function']['offset'])+1);result.append(('wrong_function',wrong))
    if kind!='approximation':
        wrong=deepcopy(cert);wrong['mean_zero']=not wrong['mean_zero'];result.append(('wrong_space',wrong))
    wrong=deepcopy(cert)
    if kind=='spectrum':wrong['lower']=str(F(wrong['upper'])+1)
    elif kind=='approximation':wrong['error_upper']='-1'
    else:
        def change(obj):
            if type(obj) is dict:
                for key,value in obj.items():
                    if ('count' in key or 'multiplicity' in key) and type(value) is int:obj[key]=value+1;return True
                for value in obj.values():
                    if change(value):return True
            elif type(obj) is list:
                for value in obj:
                    if change(value):return True
            return False
        if not change(wrong):wrong['negative_count']=999
    result.append(('bound_or_count_tamper',wrong))
    if 'family_certificate' in cert:
        wrong=deepcopy(cert);wrong['family_certificate']['uniform_width_upper']='0'
        result.append(('family_uniform_width_tamper',wrong))
        wrong=deepcopy(cert);wrong['family_certificate']['cells'][0]['amplitude_interval'][0]='-51/100'
        result.append(('family_coverage_tamper',wrong))
    if kind=='negative_spectrum':
        wrong=deepcopy(cert)
        def tail(obj):
            if type(obj) is dict:
                for key,value in obj.items():
                    if 'tail' in key and type(value) in (int,str):
                        try:obj[key]=str(F(value)+1);return True
                        except (ValueError,ZeroDivisionError):pass
                for value in obj.values():
                    if tail(value):return True
            elif type(obj) is list:
                for value in obj:
                    if tail(value):return True
            return False
        if tail(wrong):result.append(('tail_tamper',wrong))
    return result


def stream_summary(streams):
    result={}
    for system in SYSTEMS:
        selected=[s for s in streams if s['system']==system]
        result[system]={'expected_sessions':2,'returned_sessions':len(selected),'checkpoints':{},
            'all_session_wall_seconds':sum(s['process']['wall_seconds'] for s in selected),
            'all_session_cpu_seconds':sum(s['process']['cpu_seconds'] for s in selected),
            'complete_sessions':sum(s['session_complete'] for s in selected)}
        for n in (1,4,16):
            values=[s['checkpoints'][str(n)] for s in selected]
            result[system]['checkpoints'][str(n)]={'fixed_query_denominator':2*n,
                'successful_queries':sum(v['successful_queries'] for v in values),
                'all_sessions_all_queries_success':len(values)==2 and all(v['successful_queries']==n for v in values),
                'warm_plus_preparation_wall_seconds':dist([v['including_preparation_wall_seconds'] for v in values]),
                'warm_plus_preparation_cpu_seconds':dist([v['including_preparation_cpu_seconds'] for v in values]),
                'total_including_cold_replay_wall_seconds':dist([v['total_with_independent_cold_replays_wall_seconds'] for v in values]),
                'total_including_cold_replay_cpu_seconds':dist([v['total_with_independent_cold_replays_cpu_seconds'] for v in values]),
                'prefix_complete_each_session':[v['prefix_complete'] for v in values]}
        result[system]['query_phase_wall_seconds']=dist([r['phase_wall_seconds'] for s in selected for r in s['rows'] if r['phase_wall_seconds'] is not None])
        for label,indexes in [('first_query',range(1)),('subsequent_queries',range(1,16))]:
            for metric in ('wall','cpu'):
                key='phase_'+metric+'_seconds'
                result[system][label+'_'+metric+'_seconds']=dist([s['rows'][i][key] for s in selected for i in indexes if s['rows'][i][key] is not None])
        result[system]['outside_box_rejections']=[s.get('outside_probe') for s in selected] if system=='B' else []
    return result


def extension_summary(rows,cases):
    result={}
    for system in SYSTEMS:
        selected=[r for r in rows if r['system']==system]
        details=[]
        for item in cases:
            rs=[r for r in selected if r['case_id']==item['id']];r=rs[0] if len(rs)==1 else {}
            a=r.get('assessment',{});valid=a.get('certificate_valid') is True and r.get('within_budget') is True
            details.append({'case_id':item['id'],'classification':item['classification'],
                'certificate_valid_within_budget':valid,'count_certified':valid and a.get('count_certified') is True,
                'sum_target_met':valid and a.get('sum_target_met') is True,
                'both_success':r.get('success') is True,'negative_count':a.get('negative_count'),
                'sum_lower':a.get('sum_lower'),'sum_upper':a.get('sum_upper'),'sum_width':a.get('value'),
                'unsupported':r.get('response',{}).get('unsupported') is True})
        result[system]={'cases':len(cases),'count_successes':sum(r['count_certified'] for r in details),
                       'sum_successes':sum(r['sum_target_met'] for r in details),
                       'both_successes':sum(r['both_success'] for r in details),'per_case':details}
    return result


def run_comparison(final=False,smoke=False):
    initialize();pointer=checked_freeze() if final else None
    before=identities();stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    directory=DATA/(('final_' if final else 'smoke_')+stamp);directory.mkdir()
    cases=read(DATA/'common_cases.json');methods=SYSTEMS;repetitions=2
    if smoke:cases=cases[:1];methods=('baseline',);repetitions=1
    rows=[]
    for repeat in range(repetitions):
        order=methods[repeat%len(methods):]+methods[:repeat%len(methods)]
        for item in cases[repeat:]+cases[:repeat]:
            for system in order:
                row=run_process({'mode':'single','system':system,'case':item,'repetition':repeat,'limit':5},directory/('common_'+system+'_'+item['id']+'_r'+str(repeat)))
                row['repetition']=repeat;rows.append(row)
                write(directory/'checkpoint.json',{'phase':'common','rows':rows,'complete':False})
    streams=[];extensions=[];attacks=[]
    if not smoke:
        stream_cases=read(DATA/'stream_cases.json');family=read(DATA/'family_spec.json')
        for session in range(2):
            order=SYSTEMS[session:]+SYSTEMS[:session]
            for system in order:
                result=run_process({'mode':'stream','system':system,'session':session,'cases':stream_cases,'family':family},directory/('stream_'+system+'_s'+str(session)))
                result['session']=session;replays=[]
                for item,row in zip(stream_cases,result['rows']):
                    cert=certificate(row.get('response'))
                    if cert is None:replay={'success':False,'reason':'no_certificate','process':{'wall_seconds':0,'cpu_seconds':0},'not_executed':True}
                    else:replay=run_process({'mode':'replay','system':system,'case':item,'certificate':cert,'limit':5},directory/('cold_replay_'+system+'_s'+str(session)+'_'+item['id']))
                    row['independent_cold_replay_success']=replay['success'];replays.append(replay)
                result['cold_replays']=replays;result['checkpoints']={}
                for n in (1,4,16):
                    tail=result['rows'][n-1];rewall=sum(r['process']['wall_seconds'] for r in replays[:n]);recpu=sum(r['process']['cpu_seconds'] for r in replays[:n])
                    complete=tail['cumulative_wall_seconds'] is not None
                    wall=tail['cumulative_wall_seconds'] if complete else result['process']['wall_seconds']
                    cpu=tail['cumulative_cpu_seconds'] if complete else result['process']['cpu_seconds']
                    result['checkpoints'][str(n)]={'queries':n,'successful_queries':sum(r['success'] and r['independent_cold_replay_success'] for r in result['rows'][:n]),
                        'prefix_complete':complete,'incomplete_prefix_cost':'all_observed_session_cost_charged_if_prefix_missing',
                        'including_preparation_wall_seconds':wall,
                        'including_preparation_cpu_seconds':cpu,
                        'extra_full_cold_replay_wall_seconds':rewall,'extra_full_cold_replay_cpu_seconds':recpu,
                        'total_with_independent_cold_replays_wall_seconds':wall+rewall,
                        'total_with_independent_cold_replays_cpu_seconds':cpu+recpu}
                streams.append(result);write(directory/'streams_checkpoint.json',{'streams':streams,'complete':False})
        for item in read(DATA/'extension_cases.json'):
            for system in SYSTEMS:
                extensions.append(run_process({'mode':'single','system':system,'case':item,'limit':10},directory/('extension_'+system+'_'+item['id'])))
        representatives={}
        for row in rows+extensions:
            cert=certificate(row.get('response'));item=next((c for c in cases+read(DATA/'extension_cases.json') if c['id']==row['case_id']),None)
            if cert is not None and row.get('assessment',{}).get('certificate_valid') and item is not None:
                representatives.setdefault((row['system'],item['task']['kind']),(item,cert))
        for session in streams:
            if session['system']=='B':
                for item,row in zip(stream_cases,session['rows']):
                    cert=certificate(row.get('response'))
                    if cert is not None and row.get('assessment',{}).get('certificate_valid'):
                        representatives.setdefault(('B','family_query'),(item,cert))
        for (system,kind),(item,cert) in representatives.items():
            for name,mutation in mutations(cert,item['task']['kind']):
                response=run_process({'mode':'attack','system':system,'case':item,'certificate':mutation,'limit':5},directory/('attack_'+system+'_'+kind+'_'+name))
                direct=response.get('direct_verifier_result',{})
                attacks.append({'system':system,'kind':kind,'mutation':name,
                    'falsely_accepted':direct.get('certificate_valid') is True,
                    'conclusive':response.get('within_budget') is True and type(direct.get('certificate_valid')) is bool,
                    'result':response})
    after=identities();preservation=baseline_preserved();source_ok=before==after and preservation['all_unchanged']
    if final:
        try:source_ok=source_ok and checked_freeze()==pointer
        except (OSError,ValueError,RuntimeError):source_ok=False
    expected_groups={(s,k) for s in SYSTEMS for k in ('spectrum','approximation')}|{('A','negative_spectrum'),('B','family_query')}
    observed_groups={(a['system'],a['kind']) for a in attacks}
    outside=[s['outside_probe'] for s in streams if s['system']=='B']
    reliability_covered=expected_groups==observed_groups and len(outside)==2 and all(p and p['rejected'] and not p['inconclusive'] for p in outside)
    report={'version':PROTOCOL_VERSION,'classification':'public_development_comparison_not_holdout','final_requested':final,
        'valid_frozen_comparison':final and source_ok,'smoke_only':smoke,'source_unchanged':source_ok,
        'freeze':pointer,'identity_before':before,'identity_after':after,'baseline_preservation':preservation,'common_rows':rows,
        'common_summary':common_summary(rows,cases,repetitions,methods),'stream_sessions':streams,
        'stream_summary':stream_summary(streams) if streams else {},
        'extension_rows':extensions,'extension_summary':extension_summary(extensions,read(DATA/'extension_cases.json')) if extensions else {},
        'reliability_attacks':attacks,
        'false_acceptances':sum(a['falsely_accepted'] for a in attacks),'actual_model_tokens':None,
        'reliability_all_observed_attacks_pass':bool(attacks) and all(a['conclusive'] and not a['falsely_accepted'] for a in attacks),
        'reliability_bank_attacks_present':any(a['kind']=='family_query' for a in attacks),
        'reliability_expected_groups':[list(x) for x in sorted(expected_groups)],
        'reliability_observed_groups':[list(x) for x in sorted(observed_groups)],
        'reliability_coverage_complete':reliability_covered,
        'reliability_complete_pass':reliability_covered and all(a['conclusive'] and not a['falsely_accepted'] for a in attacks),
        'limitations':['Two repeats only; timings preliminary.','Extensions never enter common-task success denominator.',
                       'Shared mathematical dependencies remain trusted; this is not formal proof-assistant verification.',
                       'All inputs are public development data, including B001 historical data.']}
    write(directory/'report.json',report)
    if final and not source_ok:raise RuntimeError('Frozen identity changed; saved results are invalid for final comparison')
    return {'report_path':str(directory/'report.json'),'valid_frozen_comparison':report['valid_frozen_comparison'],
            'common_summary':report['common_summary'],'false_acceptances':report['false_acceptances']}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('command',choices=['initialize','smoke','freeze','run','worker','selftest'])
    parser.add_argument('--spec');args=parser.parse_args()
    if args.command=='worker':worker(read(args.spec));return
    if args.command=='initialize':result=initialize()
    elif args.command=='freeze':result=freeze()
    elif args.command=='smoke':result=run_comparison(smoke=True)
    elif args.command=='run':result=run_comparison(final=True)
    else:
        initialize();dest=DATA/('watchdog_selftest_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
        outcome=run_process({'mode':'selftest_timeout','system':'baseline','limit':.15},dest)
        item=public_inputs()['common_cases.json'][0]
        proof={'function':deepcopy(item['task']['function']),'mean_zero':False,'lower':'0','upper':'1/1000000000'}
        yes=lambda c,t:{'certificate_valid':True}
        tests={'parent_kills_BaseException_swallowing_worker':outcome['process']['watchdog']=='entire_process' and not outcome['success'],
               'exact_width_target':assess(proof,item,yes)['target_met']}
        for name,key,value in [('wrong_function','function',q(offset='1')),('wrong_space','mean_zero',True),
                               ('too_wide','upper','1/100'),('negative_width','upper','-1'),
                               ('wrong_tolerance','tolerance','1/10')]:
            changed=deepcopy(proof);changed[key]=value
            tests[name]=assess(changed,item,yes)['target_met'] is False
        denied=assess(proof,item,lambda c,t:{'status':'target_met','target_met':True})
        tests['status_only_verdict_rejected']=denied['certificate_valid'] is False
        result={'passed':all(tests.values()),'synthetic_only':True,'tests':tests,'watchdog_result':outcome}
        write(dest/'selftest_summary.json',result)
    print(json.dumps(result,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
