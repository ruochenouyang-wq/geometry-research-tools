"""JSONL access to certified original-function spectra and proof reuse."""
import inspect
import json
import os
from pathlib import Path
import re
import sys
import tempfile
from dependencies import ROOT, canonical, digest
import spectral_transfer as s

RECORD = 'general_function_request_record_v1'
FORMATS = {'linf':s.LINF,'l2':s.L2,'fixed':s.FIXED}
SOLVE_SIGNATURE = inspect.signature(s.solve)


def verify_record(record):
    try:
        if not isinstance(record,dict) or set(record)!={'format','arguments','certificate'}: return False
        if record['format'] != RECORD: return False
        args,c=record['arguments'],record['certificate']
        bound=SOLVE_SIGNATURE.bind(**args);bound.apply_defaults();a=bound.arguments
        options={} if a['source_options'] is None else a['source_options']
        if not isinstance(options,dict) or set(options)-{'modes','max_modes','max_m','bits','L','tolerance'}:return False
        if c['format'] != FORMATS[a['method']]:return False
        if not s.verify(c,expected_function=a['function'],expected_mean_zero=a['mean_zero'],expected_width=a['target_width']):return False
        model=c['model']
        if a['method']=='linf':
            if a['reference'] is not None:return False
            s.integer(a['order'],0,24,'order')
            return model['order']==a['order']
        s.integer(a['degree'],0,24,'degree');s.integer(a['sqrt_bits'],0,256,'sqrt_bits')
        if model['degree']!=a['degree'] or model['sqrt_bits']!=a['sqrt_bits']:return False
        if a['method']=='fixed':
            return c['coercivity']['reference']==s.poly(a['reference']) and c['coercivity']['sqrt_bits']==a['sqrt_bits']
        return a['reference'] is None
    except (ValueError,TypeError,KeyError,AttributeError,IndexError,ArithmeticError):return False


class Service:
    def __init__(self,directory=None):
        self.directory=Path(directory) if directory is not None else ROOT/'service_certificates'

    def path(self,identifier):
        if not isinstance(identifier,str) or re.fullmatch('[a-f0-9]{64}',identifier) is None:
            raise ValueError('A complete lowercase SHA256 identifier is required')
        p=self.directory/(identifier+'.json')
        if p.is_symlink():raise ValueError('Certificate symlinks are not supported')
        return p

    def put(self,record):
        if not verify_record(record):raise ValueError('Proof or original request binding failed')
        identifier=digest(record);path=self.path(identifier);raw=canonical(record).encode()
        self.directory.mkdir(parents=True,exist_ok=True)
        if path.exists():
            if path.read_bytes()!=raw:raise ValueError('Stored content changed')
            return identifier
        temporary=None
        try:
            with tempfile.NamedTemporaryFile(dir=self.directory,delete=False) as handle:
                temporary=handle.name;handle.write(raw)
            os.replace(temporary,path)
        finally:
            if temporary and os.path.exists(temporary):os.unlink(temporary)
        return identifier

    def get(self,identifier):
        raw=self.path(identifier).read_bytes();record=json.loads(raw)
        if canonical(record).encode()!=raw or digest(record)!=identifier or not verify_record(record):
            raise ValueError('Stored certificate failed hash or proof verification')
        return record

    @staticmethod
    def brief(record,identifier):
        c=record['certificate']
        return {'certificate_id':identifier,'mathematical_proof_verified':True,
                **{k:c[k] for k in ('function','scope','mean_zero','lower','upper','exact_width',
                                    'status','requested_width','error_budget')},
                'formal_proof_assistant_checked':False,
                'runtime_budget_is_not_a_mathematical_claim':True}

    def call(self,request):
        if not isinstance(request,dict):raise ValueError('JSON object required')
        op=request.get('op')
        if op in ('fetch','verify'):
            if set(request)!={'op','certificate_id'}:raise ValueError('Unknown proof request fields')
            record=self.get(request['certificate_id']);out=self.brief(record,request['certificate_id'])
            if op=='fetch':out['record']=record
            return out
        if op!='solve':raise ValueError('Supported operations: solve, verify, fetch')
        args={k:v for k,v in request.items() if k!='op'}
        SOLVE_SIGNATURE.bind(**args)
        c=s.solve(**args)
        record={'format':RECORD,'arguments':args,'certificate':c}
        identifier=self.put(record)
        return self.brief(record,identifier)


if __name__=='__main__':
    service=Service()
    for line in sys.stdin:
        try:out=service.call(json.loads(line))
        except (ValueError,TypeError,KeyError,IndexError,ArithmeticError,OSError,RecursionError) as e:
            out={'status':'not_certified','exception':type(e).__name__,'reason':str(e)}
        print(canonical(out),flush=True)
