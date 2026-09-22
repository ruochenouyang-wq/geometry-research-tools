"""Read-only projection of frozen A negative-spectrum proofs to ground intervals.

This adapter never delegates a spectrum task to the baseline solver.  It uses
only fully verified, nonempty negative spectra with individual eigenvalue
intervals.  Time budgets control execution, not the mathematical proposition.
"""
from copy import deepcopy
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path
import sys
from time import perf_counter

sys.dont_write_bytecode = True
FORMAT = 'a_negative_spectrum_ground_adapter_v1'
VERSION = 'cross-dimension-A-proof-projection-v1'
PROJECTION = {
    'method':'minimum_of_complete_negative_radial_first_eigenvalues',
    'requires_strictly_negative_spectrum_nonempty':True,
    'zero_count_and_omitted_angular_sectors_are_nonnegative':True,
    'real_multiplicity_does_not_change_minimum':True,
    'source_sum_precision_is_not_ground_precision':True,
}


def canonical(value):
    return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)


class ProjectionUnavailable(ValueError):
    def __init__(self,reason,source_assessment=None):
        super().__init__(reason)
        self.source_assessment=source_assessment


class Engine:
    def __init__(self,family_spec=None):
        # Private module name: never occupy the historical module name bridge.
        path=Path(__file__).resolve().parent.parent/'geometry_parallel_v1_20260921'/'a_multispectrum'/'variant.py'
        spec=importlib.util.spec_from_file_location('cross_dimension_frozen_a_multispectrum',path)
        self.source=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.source)
        self.source_path=str(path)

    def prepare(self,seconds):
        return {'ok':True,'no_preparation':True}

    def _task(self,task):
        if type(task) is not dict or task.get('kind')!='spectrum':
            raise ProjectionUnavailable('unsupported_kind: only spectrum tasks can be projected')
        return self.source.shared.baseline.normalize_task(task)

    def _source_task(self,task):
        return {'kind':'negative_spectrum','function':deepcopy(task['function']),
                'mean_zero':task['mean_zero'],'tolerance':task['tolerance'],
                'budget':deepcopy(task['budget'])}

    def _check_source_task(self,source_task,task):
        normalized=self.source.normalize_task(source_task)
        if canonical(normalized)!=canonical(source_task):
            raise ValueError('Source task must be canonical and contain only accepted fields')
        expected=self._source_task(task)
        for key in ('kind','function','mean_zero','tolerance'):
            if canonical(normalized[key])!=canonical(expected[key]):
                raise ValueError('Source task changed '+key)
        if normalized['budget']['wall_seconds']>task['budget']['wall_seconds']:
            raise ValueError('Source execution budget exceeds the external task budget')
        return normalized

    def _derive(self,source_cert,source_task,task):
        source_task=self._check_source_task(source_task,task)
        # Every invocation replays all frozen A evidence; no cached digest or
        # solver success flag can bypass this call.
        checked=self.source.verify(source_cert,source_task)
        if checked.get('certificate_valid') is not True or checked.get('fallback') is not False:
            raise ProjectionUnavailable('source_certificate_failed_full_replay',checked)
        if source_cert.get('format')!=self.source.FORMAT or source_cert.get('kind')!='negative_spectrum':
            raise ProjectionUnavailable('wrong_native_source_format',checked)
        if source_cert.get('geometry')!='unit_S2' or source_cert.get('full_infinite_space_covered') is not True:
            raise ProjectionUnavailable('source_does_not_cover_original_sphere_space',checked)
        if checked.get('count_certified') is not True:
            raise ProjectionUnavailable('negative_count_unresolved',checked)
        count=checked.get('negative_count')
        if type(count) is not int or count<=0:
            raise ProjectionUnavailable('no_certified_strictly_negative_eigenvalue',checked)
        if checked.get('sum_certified') is not True:
            raise ProjectionUnavailable('individual_negative_eigenvalue_intervals_missing',checked)
        eligible=[];lowers=[];uppers=[]
        for m,row in enumerate(source_cert['sectors']):
            if row['azimuth_m']!=m or row['count_certified'] is not True:
                raise ProjectionUnavailable('sector_coverage_or_count_unresolved',checked)
            if row['count_lower']==0:
                continue
            intervals=row['eigenvalue_intervals']
            if not intervals or intervals[0].get('index')!=1:
                raise ProjectionUnavailable('radial_first_eigenvalue_interval_missing',checked)
            first=intervals[0];eligible.append(m)
            lowers.append(F(first['lower']));uppers.append(F(first['upper']))
        if not eligible or F(source_cert['angular_tail_lower'])<0:
            raise ProjectionUnavailable('negative_angular_cover_not_complete',checked)
        lower,upper=min(lowers),min(uppers)
        if lower>upper or upper>0:
            raise ArithmeticError('Contradictory projected ground interval')
        width=upper-lower;met=width<=F(task['tolerance'])
        cert={'format':FORMAT,'kind':'spectrum','function':deepcopy(task['function']),
              'mean_zero':task['mean_zero'],'geometry':'unit_S2','measure':'d_sigma/(4*pi)',
              'scope':source_cert['scope'],'eigenvalue_index':1,
              'lower':str(lower),'upper':str(upper),'exact_width':str(width),
              'tolerance':task['tolerance'],'eligible_sectors':eligible,
              'projection_proof':deepcopy(PROJECTION),'source_task':deepcopy(source_task),
              'source_certificate':deepcopy(source_cert),
              'source_sum_target_met':checked['sum_target_met'],
              'full_infinite_space_covered':True,'baseline_fallback_used':False,
              'status':'target_met' if met else 'certified_open',
              'formal_proof_assistant_checked':False}
        return cert,checked

    @staticmethod
    def _assessment(cert):
        width=F(cert['upper'])-F(cert['lower'])
        met=width<=F(cert['tolerance'])
        return {'certificate_valid':True,'target_met':met,'metric':str(width),
                'lower':cert['lower'],'upper':cert['upper'],'exact_width':str(width),
                'tolerance':cert['tolerance'],'status':'target_met' if met else 'certified_open',
                'fallback':False,'source_sum_target_met':cert['source_sum_target_met'],
                'verification_mode':'full_native_source_replay'}

    def verify(self,cert,task,full=True):
        # full=False is accepted for interface compatibility, but never weakens
        # proof checking: the source is independently replayed on every call.
        try:
            t=self._task(task)
            if type(cert) is not dict or cert.get('format')!=FORMAT:
                raise ValueError('Wrong adapter certificate format')
            fresh,_=self._derive(cert['source_certificate'],cert['source_task'],t)
            if canonical(cert)!=canonical(fresh):
                raise ValueError('Projection certificate does not match exact replay')
            return self._assessment(fresh)
        except (ValueError,TypeError,KeyError,ArithmeticError,IndexError,OverflowError,RecursionError) as exc:
            return {'certificate_valid':False,'target_met':False,'metric':None,
                    'status':'verification_failed','fallback':False,
                    'error':type(exc).__name__+': '+str(exc)}

    def solve(self,task):
        started=perf_counter();source_response=None;source_task=None;source_run_task=None
        source_assessment=None;source_seconds=0.0;projection_seconds=0.0;limit=None
        adapter_extra_replay=False
        result={'version':VERSION,'route':'native_A_negative_spectrum_ground_projection',
                'certificate':None,'certificate_valid':False,'target_met':False,'metric':None,
                'fallback':False,'actual_model_calls':0,'actual_model_tokens':None}
        try:
            t=self._task(task);limit=t['budget']['wall_seconds']
            source_task=self._source_task(t)
            remaining=max(0.0,limit-(perf_counter()-started))
            source_run_task=deepcopy(source_task)
            source_run_task['budget']['wall_seconds']=remaining
            if remaining<=0:
                raise ProjectionUnavailable('solve_budget_exhausted_before_native_A')
            begin=perf_counter()
            try:
                source_response=self.source.solve(source_run_task)
            finally:
                source_seconds=perf_counter()-begin
            source_cert=source_response.get('certificate')
            if source_cert is None:
                raise ProjectionUnavailable('native_A_returned_no_negative_spectrum_certificate')
            begin=perf_counter()
            try:
                adapter_extra_replay=True
                cert,source_assessment=self._derive(source_cert,source_task,t)
            finally:
                projection_seconds=perf_counter()-begin
            result.update(certificate=cert,**self._assessment(cert))
        except ProjectionUnavailable as exc:
            source_assessment=exc.source_assessment
            reason=str(exc)
            if reason.startswith('unsupported_kind:') or reason=='no_certified_strictly_negative_eigenvalue':
                status='unsupported'
            elif reason=='solve_budget_exhausted_before_native_A':
                status='budget_exceeded'
            elif reason in ('native_A_returned_no_negative_spectrum_certificate',
                            'negative_count_unresolved',
                            'individual_negative_eigenvalue_intervals_missing',
                            'radial_first_eigenvalue_interval_missing',
                            'sector_coverage_or_count_unresolved'):
                status='no_verified_ground_certificate'
            else:
                status='verification_failed'
            result.update(status=status,reason=reason)
        except (ValueError,TypeError,KeyError,ArithmeticError,IndexError,OverflowError) as exc:
            result.update(status='invalid_or_unverified',reason=type(exc).__name__+': '+str(exc))
        elapsed=perf_counter()-started
        within=limit is None or elapsed<=limit
        result.update(source_task=source_task,source_run_task=source_run_task,
                      source_response=source_response,source_assessment=source_assessment,
                      source_solve_seconds=source_seconds,
                      source_replay_and_projection_seconds=projection_seconds,
                      total_elapsed_seconds=elapsed,within_budget=within,
                      budget_scope='whole_solve_cooperative; external parent scores the hard deadline',
                      metadata={
                          'mathematical_source':'frozen A complete negative spectrum and ordered radial intervals',
                          'adapter_solve_extra_full_source_replay':adapter_extra_replay,
                          'native_solver_policy':'A.solve also replays its returned source certificate',
                          'external_verify_policy':'each Engine.verify independently replays the entire A source',
                          'timing_interpretation':'actual adapter interface cost, not isolated Kernel timing'})
        if not within:
            result.update(status='budget_exceeded',target_met=False)
        return result
