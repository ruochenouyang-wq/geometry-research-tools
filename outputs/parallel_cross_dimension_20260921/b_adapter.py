"""Evaluation adapter for frozen B's native amplitude-family route.

No family selection or baseline fallback occurs here. The harness must charge
construction as well as prepare/query/replay, and enforce any hard deadline.
Preparation's cooperative budget includes the first owned-context replay.
"""
import importlib.util
import json
from math import isfinite
from pathlib import Path
import sys
from time import perf_counter, process_time


ORIGINAL = (Path(__file__).resolve().parent.parent /
            'geometry_parallel_v1_20260921/b_parameter/variant.py')
MODULE_NAME = '_cross_dimension_frozen_b_parameter'
INPUT_ERRORS = (ValueError, TypeError, KeyError, ArithmeticError, IndexError,
                RecursionError)


def _load():
    if MODULE_NAME in sys.modules:
        return sys.modules[MODULE_NAME]
    spec = importlib.util.spec_from_file_location(MODULE_NAME, ORIGINAL)
    module = importlib.util.module_from_spec(spec)
    sys.modules[MODULE_NAME] = module
    try:
        spec.loader.exec_module(module)
    except BaseException:
        sys.modules.pop(MODULE_NAME, None)
        raise
    return module


def _failure(status, reason):
    return {'ok': False, 'status': status, 'reason': reason,
            'certificate': None, 'certificate_valid': False,
            'target_met': False, 'fallback': False}


class Engine:
    """One explicit family and one verified, privately owned query context."""

    def __init__(self, family_spec=None):
        begin, cpu = perf_counter(), process_time()
        self.module = _load()
        self._family_supplied = family_spec is not None
        self._family_json = None
        self._family_error = 'No predeclared family was supplied'
        self._context = None
        if family_spec is not None:
            try:
                owned = json.loads(self.module.canonical(family_spec))
                self._family_json = self.module.canonical(
                    self.module.normalize_family(owned))
                self._family_error = None
            except INPUT_ERRORS as error:
                self._family_error = str(error)
        self.construction_wall_seconds = perf_counter()-begin
        self.construction_cpu_seconds = process_time()-cpu

    def prepare(self, seconds):
        """Return full native attempts and first context verification, charged."""
        begin, cpu = perf_counter(), process_time()
        self._context = None
        response, verification, context = None, None, None
        result = _failure('preparation_failed', 'Preparation has not completed')
        budget = None
        try:
            if self._family_json is None:
                result = _failure('unsupported', self._family_error)
            elif (isinstance(seconds, bool) or not isfinite(float(seconds))
                  or not 0 <= float(seconds) <= 120):
                result = _failure('invalid_budget', 'Require finite 0..120 seconds')
            else:
                budget = float(seconds)
                remaining = max(0.0, budget-(perf_counter()-begin))
                response = self.module.prepare_family(
                    json.loads(self._family_json), preparation_seconds=remaining)
                bank = response.get('bank')
                if bank is None:
                    result = _failure(response.get('status', 'preparation_failed'),
                                      'The native preparation issued no family bank')
                elif perf_counter()-begin >= budget:
                    result = _failure('preparation_budget_exceeded',
                                      'No time remains for the first context replay')
                else:
                    # Snapshot before checking the box, then replay the same copy.
                    owned = json.loads(self.module.canonical(bank))
                    actual = self.module.normalize_family(owned['family'])
                    if self.module.canonical(actual) != self._family_json:
                        raise ValueError('Prepared bank differs from the declared family')
                    context = self.module.prepare_context(owned)
                    verification = dict(context.verification)
                    if verification.get('certificate_valid') is not True:
                        result = _failure('verification_failed',
                                          'The first context replay rejected the family')
                    else:
                        result = {'ok': True, 'status': 'ready', 'fallback': False,
                                  'certificate_valid': True, 'target_met': True}
        except Exception as error:
            result = _failure('preparation_failed', type(error).__name__+': '+str(error))
        elapsed = perf_counter()-begin
        if budget is not None and elapsed > budget:
            result = _failure('preparation_budget_exceeded',
                              'Preparation plus first context replay exceeded the budget')
        if result['ok']:
            self._context = context
        result.update(preparation=response, verification=verification,
                      preparation_wall_seconds=elapsed,
                      preparation_cpu_seconds=process_time()-cpu,
                      preparation_budget_seconds=budget,
                      budget_scope='native_preparation_plus_first_context_replay; cooperative',
                      construction_wall_seconds=self.construction_wall_seconds,
                      construction_cpu_seconds=self.construction_cpu_seconds)
        return result

    def solve(self, task):
        begin, cpu = perf_counter(), process_time()
        if self._family_json is None:
            result = _failure('unsupported', self._family_error)
        elif self._context is None:
            result = _failure('not_prepared', 'A successful preparation is required')
        else:
            try:
                result = self._context.query(task)
                result['ok'] = (result.get('certificate_valid') is True and
                                result.get('status') != 'budget_exceeded')
            except INPUT_ERRORS as error:
                result = _failure('unsupported', str(error))
        result.update(solve_wall_seconds=perf_counter()-begin,
                      solve_cpu_seconds=process_time()-cpu, fallback=False)
        return result

    def verify(self, cert, task, full=True):
        """Cold full replay is self-contained; fast replay requires this context."""
        begin, cpu = perf_counter(), process_time()
        try:
            if type(full) is not bool:
                raise ValueError('full must be a boolean')
            if self._family_supplied and self._family_json is None:
                raise ValueError('The explicitly supplied family is unsupported: '+self._family_error)
            if not self.module.native(cert) or type(cert) is not dict:
                raise ValueError('Require an exact JSON-native B query certificate')
            owned = json.loads(self.module.canonical(cert))
            if owned.get('format') != self.module.QUERY:
                raise ValueError('Only native B family-query certificates are supported')
            family = self.module.normalize_family(owned['family_certificate']['family'])
            if self._family_json is not None and self.module.canonical(family) != self._family_json:
                raise ValueError('Query certificate differs from the declared family')
            # A missing family permits self-contained cold replay, never solve.
            self.module._bound_task(family, task)
            if full:
                result = self.module.verify(owned, task)
            elif self._context is None:
                result = _failure('not_prepared', 'Fast verification requires this verified context')
            else:
                result = self._context.verify(owned, task)
            result['ok'] = result.get('certificate_valid') is True
            result.setdefault('status', 'verified' if result['ok'] else 'verification_failed')
        except INPUT_ERRORS as error:
            result = _failure('verification_failed', str(error))
        result.update(verification_wall_seconds=perf_counter()-begin,
                      verification_cpu_seconds=process_time()-cpu,
                      verification_scope='full_independent_replay' if full is True else 'owned_verified_context',
                      fallback=False)
        return result
