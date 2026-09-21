"""Public verification dispatch, including nested proofs and success claims."""
from fractions import Fraction as F
import v03_tail as v3
import v05_range as v5
import v07_precision as v7
import v08_cluster as v8
import v09_family as v9


def verify_certificate(cert,independent=True):
    if not isinstance(cert,dict):return False
    method=cert.get('method')
    if method in ('full_residual_temple_v1','positive_exponential_local_energy_v1'):
        return v3.old.verify(cert)
    if method is None and cert.get('model')==v3.base.MODEL:
        return v3.base.verify(cert)
    if method=='modewise_schur_v3':return v3.verify(cert,independent)
    if method=='rational_trial_refinement_v7':return v7.verify(cert)
    if method=='first_cluster_projector_v8':return v8.verify(cert)
    if method in ('uniform_affine_family_v9','affine_family_counterexample_v9'):return v9.verify(cert)
    if method=='polynomial_range_document_v5':
        try:
            q=v3.base.potential(cert['q_coefficients']);v5.verify_bounds(q,cert['range_proof'])
            return cert=={'method':method,'q_coefficients':[str(x) for x in q],'range_proof':cert['range_proof']}
        except (KeyError,TypeError,ValueError,ZeroDivisionError):return False
    return False


def nested_proofs_valid(value,independent):
    if isinstance(value,dict):
        if 'method' in value and not verify_certificate(value,independent):return False
        return all(nested_proofs_valid(x,independent) for x in value.values())
    if isinstance(value,list):return all(nested_proofs_valid(x,independent) for x in value)
    return True


def verify_document(document,independent=True):
    try:
        if not isinstance(document,dict):return False
        cert=document if 'method' in document or ('model' in document and 'lower' in document) else document['certificate']
        if not verify_certificate(cert,independent) or not nested_proofs_valid(cert,independent):return False
        if 'before_certificate' in document and not verify_certificate(document['before_certificate'],independent):return False
        status=document.get('status')
        if status in ('target_met','target_not_met'):
            tolerance=F(document['tolerance'])
            if tolerance<=0 or ((F(cert['exact_width'])<=tolerance)!=(status=='target_met')):return False
        if document.get('algorithm_version')==8:
            if document['cluster_size']!=len(cert['cluster_indices']):return False
            if any(len(col)!=document['modes'] for col in cert['trial_columns']):return False
        return True
    except (KeyError,TypeError,ValueError,ZeroDivisionError):return False
