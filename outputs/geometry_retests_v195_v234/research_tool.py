"""JSONL interface to the four repaired problem families.

The original request is stored with every mathematical certificate. Reading a
certificate rechecks both proof and request binding, without rerunning search.
This packaging is not counted among the forty mathematical increments.
"""
from pathlib import Path
import inspect
import json
import os
import re
import sys
import tempfile
from support import ROOT, F, exact, canonical, digest, same
import global_constraints as constraints
import general_spectrum as spectrum
import compact_gn as gn
import parity_spectrum as parity

RECORD = 'geometry_retest_tool_record_v1'
SPECTRUM_FIELDS = {'q','spaces','k','quantities','threshold','tolerance','modes',
                   'max_modes','max_m','max_radial','bits','max_steps'}
OPERATIONS = {
    'constraints': (constraints.certify, {'q','constraints','L','bits','snap'}),
    'constrained_inequality': (constraints.inequality, {'q','constraints','threshold','L','bits'}),
    'spectrum': (spectrum.research_driver, SPECTRUM_FIELDS),
    'cap': (gn.evaluate_cap, {'n','original'}),
    'cap_original': (gn.bridge_original, {'raw_t_coefficients','n'}),
    'parity_ground': (parity.adaptive_ground, {'q','target','start_L','max_L','bits','wall_budget_seconds'}),
}


def arguments_for(op, arguments):
    if op not in OPERATIONS or not isinstance(arguments, dict):
        raise ValueError('Unknown operation or malformed arguments')
    function, allowed = OPERATIONS[op]
    if set(arguments)-allowed:
        raise ValueError('Unknown argument fields')
    inspect.signature(function).bind(**arguments)
    return function


def matching_proof(op, arguments, c):
    """Bind user inputs independently of the producer's search decisions."""
    arguments_for(op, arguments)
    if op in ('constraints','constrained_inequality'):
        expected_format = constraints.SPECTRAL if op == 'constraints' else constraints.INEQUALITY
        if c.get('format') != expected_format or not constraints.verify(c,
                expected_q=arguments['q'], expected_constraints=arguments['constraints']):
            return False
        if op == 'constrained_inequality':
            return (exact(c['threshold']) == exact(arguments['threshold']) and
                    c['spectral_certificate']['retained_degree'] == arguments.get('L',2))
        return c['retained_degree'] == arguments.get('L',2)
    if op == 'spectrum':
        if c.get('format') != spectrum.DRIVER or not spectrum.verify(c,expected_q=arguments['q']):
            return False
        spaces = arguments.get('spaces',('full_mean_zero','full_unprojected'))
        if not isinstance(spaces,(tuple,list)) or any(s not in ('full_mean_zero','full_unprojected') for s in spaces):
            return False
        options={k:v for k,v in arguments.items() if k not in ('q','spaces')}
        requests=[spectrum.request(arguments['q'],mean_zero=s=='full_mean_zero',**options) for s in spaces]
        return same(requests,c['requests'])
    if op == 'cap':
        return (c.get('format') == 'compact_cap_v224' and
                gn.verify(c,expected_n=arguments['n']) and
                c['original_scale'] is arguments.get('original',True))
    if op == 'cap_original':
        return (c.get('format') == 'compact_original_bridge_v224' and
                gn.verify(c,expected_n=arguments['n'],
                          expected_original_coefficients=arguments['raw_t_coefficients']))
    if op == 'parity_ground':
        return (parity.verify_result(c,expected_q=arguments['q'],
                                    expected_width=arguments.get('target',F(1,10**8))) and
                all(row.get('retained_degree',0) <= arguments.get('max_L',9)
                    for row in c['certificate']['characters']))
    return False


def verify(record):
    try:
        return (set(record) == {'format','operation','arguments','certificate'} and
                record['format'] == RECORD and
                matching_proof(record['operation'],record['arguments'],record['certificate']))
    except (ValueError,TypeError,KeyError,IndexError,ArithmeticError):
        return False


class Store:
    def __init__(self,directory=None):
        self.directory=Path(directory) if directory is not None else ROOT/'service_certificates'

    def path(self,identifier):
        if not isinstance(identifier,str) or re.fullmatch('[0-9a-f]{64}',identifier) is None:
            raise ValueError('A complete lowercase SHA-256 certificate identifier is required')
        path=self.directory/(identifier+'.json')
        if path.is_symlink():raise ValueError('Certificate symlinks are not supported')
        return path

    def put(self,record):
        if not verify(record):raise ValueError('Mathematical proof or original-request binding failed')
        identifier=digest(record);raw=canonical(record).encode();path=self.path(identifier)
        self.directory.mkdir(parents=True,exist_ok=True)
        if path.exists():
            if path.read_bytes()!=raw:raise ValueError('Existing certificate content changed')
            return identifier
        name=None
        try:
            with tempfile.NamedTemporaryFile(dir=self.directory,delete=False) as handle:
                name=handle.name;handle.write(raw)
            os.replace(name,path)
        finally:
            if name and os.path.exists(name):os.unlink(name)
        return identifier

    def get(self,identifier):
        raw=self.path(identifier).read_bytes();record=json.loads(raw)
        if digest(record)!=identifier or canonical(record).encode()!=raw:
            raise ValueError('Certificate content hash or encoding changed')
        if not verify(record):raise ValueError('Stored proof failed verification')
        return record


def brief(record,identifier):
    op,c=record['operation'],record['certificate']
    result={'operation':op,'certificate_id':identifier,'mathematical_proof_verified':True,
            'original_request_bound':True,'formal_proof_assistant_checked':False}
    if op in ('constraints','constrained_inequality'):
        evidence=c if op=='constraints' else c['spectral_certificate']
        result.update(scope=evidence['scope'],potential=evidence['potential'],
                      constraints=evidence['constraints'],lower=evidence['lower'],upper=evidence['upper'],
                      exact_width=evidence['exact_width'])
        result['status']=c.get('status','certified_ground_enclosure')
    elif op=='spectrum':
        result.update(all_requested_targets_met=c['all_requested_targets_met'],
                      status='targets_met' if c['all_requested_targets_met'] else 'certified_open')
        result['spaces']=[]
        for row in c['results']:
            entry={'space':row['request']['space'],'q':row['request']['q'],'decisions':row['decisions']}
            for quantity,proof in row['answers'].items():
                if quantity=='count':entry['count']={k:proof[k] for k in ('count_lower','count_upper','threshold')}
                elif quantity=='trace':entry['negative_trace']={k:proof[k] for k in ('lower','upper','width')}
                else:entry['ordered']={'k':proof['k'],'intervals':proof['ordered']}
            result['spaces'].append(entry)
    elif op in ('cap','cap_original'):
        cap=c if op=='cap' else c['evaluation'];ratio=cap['ratio']
        result.update(status='exact_trial_moments',n=cap['n'],function=cap['input_function'],
                      scope=cap['scope'],probability_ratio=ratio['probability_ratio'],
                      pi_times_area_K=ratio['pi_times_area_K'],
                      universal_GN_sharp_constant_proved=False)
    else:
        result.update(status=c['status'],scope=c['scope'],potential=c['potential'],
                      requested_width=c['requested_width'])
        for key in ('lower','upper','exact_width'):result[key]=c['certificate'][key]
    return result


class Service:
    def __init__(self,directory=None):self.store=Store(directory)

    def call(self,request):
        if not isinstance(request,dict):raise ValueError('Request must be a JSON object')
        op=request.get('op')
        if op in ('verify','fetch'):
            if set(request)!={'op','certificate_id'}:raise ValueError('Unknown fetch/verify fields')
            record=self.store.get(request['certificate_id']);response=brief(record,request['certificate_id'])
            if op=='fetch':response['record']=record
            return response
        if not isinstance(op,str):raise ValueError('Operation name required')
        args={k:v for k,v in request.items() if k!='op'};function=arguments_for(op,args)
        produced=function(**args)
        c=produced['certificate'] if op=='spectrum' else produced
        record={'format':RECORD,'operation':op,'arguments':args,'certificate':c}
        return brief(record,self.store.put(record))


if __name__=='__main__':
    service=Service()
    for line in sys.stdin:
        try:out=service.call(json.loads(line))
        except (ValueError,TypeError,KeyError,IndexError,ArithmeticError,OSError,RecursionError) as e:
            out={'status':'rejected_or_computation_failed','exception':type(e).__name__,'message':str(e)}
        print(canonical(out),flush=True)
