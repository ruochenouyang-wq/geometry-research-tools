"""One JSON request per line; emit a self-contained certificate or verification.

Examples and requests are in TOOL_GUIDE.md. The tool uses no network/model call
and writes no cache automatically. Saving the JSON response preserves evidence.
"""
import json
import sys
import direct_moments as direct
import piecewise_models as pieces
import enriched_trial as enriched


def handle(request):
    if not isinstance(request,dict): raise ValueError('Request must be an object')
    op=request.get('op')
    if op=='precise_singular':
        allowed={'op','function','max_terms','tolerance','precision_bits','iterations'}
        if not set(request)<=allowed: raise ValueError('Unknown enriched request field')
        if 'function' not in request or direct.normalize(request['function'])!=enriched.FUNCTION:
            raise ValueError('This enriched basis is certified for q(z)=-|z|^(-1/4) only')
        source=direct.full_ground(request['function'],mean_zero=True,modes=4,bits=28)
        args={k:v for k,v in request.items() if k not in ('op','function')}
        run=enriched.enrich(source,**args)
        return {'ok':True,'kind':'original_singular_function_spectrum',
                'certificate':run['certificate'],'attempts':run['attempts'],
                'status':run['status']}
    if op=='spectrum':
        allowed={'op','function','mean_zero','modes','bits','max_m','tolerance','sqrt_bits','near_tail'}
        if not set(request)<=allowed: raise ValueError('Unknown spectrum request field')
        args={k:v for k,v in request.items() if k!='op'}
        if 'function' not in args: raise ValueError('Missing original function')
        certificate=direct.full_ground(**args)
        return {'ok':True,'kind':'original_function_spectrum','certificate':certificate}
    if op=='approximate':
        allowed={'op','function','levels','degree','root_ratio','sqrt_bits'}
        if not set(request)<=allowed: raise ValueError('Unknown approximation request field')
        args={k:v for k,v in request.items() if k!='op'}
        if 'function' not in args: raise ValueError('Missing original function')
        certificate=pieces.piecewise_model(**args)
        return {'ok':True,'kind':'function_L2_approximation_only','certificate':certificate}
    if op=='verify':
        if not set(request)<= {'op','certificate','function','mean_zero','tolerance'}:
            raise ValueError('Unknown verification request field')
        cert=request.get('certificate')
        if not isinstance(cert,dict): return {'ok':True,'verified':False}
        if cert.get('format')==direct.FULL:
            verified=direct.verify_full(cert,request.get('function'),request.get('mean_zero'),request.get('tolerance'))
        elif cert.get('format')==pieces.FORMAT:
            if 'mean_zero' in request or 'tolerance' in request:
                raise ValueError('A function approximation has no spectral space or target claim')
            verified=pieces.verify_piecewise(cert,request.get('function'))
        elif cert.get('format')==enriched.FORMAT:
            source_matches=('function' not in request or direct.normalize(request['function'])==enriched.FUNCTION)
            scope_matches=('mean_zero' not in request or request['mean_zero'] is True)
            verified=source_matches and scope_matches and enriched.verify(cert,expected_tolerance=request.get('tolerance'))
        else: verified=False
        return {'ok':True,'verified':verified}
    raise ValueError('Unknown operation; use spectrum, precise_singular, approximate, or verify')


if __name__=='__main__':
    for line in sys.stdin:
        if not line.strip(): continue
        try:
            result=handle(json.loads(line))
        except (ValueError,TypeError,KeyError,ArithmeticError) as e:
            result={'ok':False,'error':type(e).__name__,'reason':str(e)}
        print(json.dumps(result,ensure_ascii=False,sort_keys=True,allow_nan=False),flush=True)
