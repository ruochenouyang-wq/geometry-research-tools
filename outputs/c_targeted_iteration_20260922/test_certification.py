"""Candidate-construction controls on one already published E09/6 trial.

These are semantic checks, not a performance comparison or fresh search.
"""
from copy import deepcopy
from dataclasses import FrozenInstanceError
from fractions import Fraction as F
import json
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import certification as c


class CertificationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        data=json.loads((Path(__file__).resolve().parent.parent/
            'c_external_comparison_20260922/numerical/RESULTS.json').read_text())
        cls.task=next(row['task'] for row in data['cases'] if row['id']=='E09')
        row=next(row for row in data['records'] if row['case_id']=='E09' and
                 row['powers_count']==6 and row['certificate'] is not None)
        cls.original=deepcopy(row['certificate']);cls.gap=cls.original['gap']
        cls.powers=cls.original['trial']['powers'];cls.coefficients=cls.original['trial']['coefficients']
        cls.baseline,cls.full,cls.sg=c._modules()
        cls.matrices=cls.full.trial_matrices(cls.task['function'],cls.powers)
        cls.stats=cls.full.statistics(cls.task['function'],cls.powers,cls.coefficients,cls.matrices)

    def prepared(self,matrices=None):
        matrices=self.matrices if matrices is None else matrices
        return SimpleNamespace(function_key=self.full.canonical(self.task['function']),
            powers=tuple(map(F,self.powers)),
            **{name:[list(row) for row in matrix]
               for name,matrix in zip(('M','H','R'),matrices)})

    def pending(self,ctx=None,**kwargs):
        ctx=c.prepare_fullspace(self.task,self.gap) if ctx is None else ctx
        return ctx.assemble(self.powers,self.coefficients,**kwargs)

    def cold_assess(self,cert,task=None):
        return self.baseline.assess(json.loads(json.dumps(cert)),self.task if task is None else task)

    def test_01_exact_legacy_equality_with_and_without_prepared(self):
        expected=self.full.certificate(self.task['function'],self.powers,
            self.coefficients,self.gap,self.task['tolerance'],mean_zero=False)
        self.assertEqual(expected,self.original)
        for prepared in (None,self.prepared()):
            result=self.pending(prepared=prepared,statistics=self.stats)
            self.assertEqual(result['certificate'],expected)
            self.assertIsNone(result['certificate_valid']);self.assertFalse(result['target_met'])
            self.assertTrue(result['verification_pending']);self.assertTrue(result['candidate_target_met'])
        verdict=self.cold_assess(result['certificate'])
        self.assertTrue(verdict['certificate_valid']);self.assertTrue(verdict['target_met'])

    def test_02_prepared_avoids_gram_and_all_local_gap_verification(self):
        with (patch.object(self.full,'trial_matrices',side_effect=AssertionError('rebuild')),
              patch.object(self.sg,'verify_gap',side_effect=AssertionError('local gap verify'))):
            result=self.pending(prepared=self.prepared(),statistics=self.stats)
        self.assertEqual(result['certificate'],self.original)
        self.assertEqual(result['construction_metadata']['local_gap_verifications'],0)

    def test_03_unprepared_rebuilds_gram_once_without_gap_verification(self):
        with (patch.object(self.full,'trial_matrices',wraps=self.full.trial_matrices) as build,
              patch.object(self.sg,'verify_gap',side_effect=AssertionError('local gap verify'))):
            self.pending(statistics=self.stats)
        self.assertEqual(build.call_count,1)

    def test_04_math_task_binding_and_budget_independence(self):
        ctx=c.prepare_fullspace(self.task,self.gap)
        for key,value in (('axis',1),('offset','1'),('amplitude','-1/3')):
            task=deepcopy(self.task);task['function'][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):
                ctx.assemble(self.powers,self.coefficients,task=task)
        for key,value in (('mean_zero',True),('tolerance','1/1000'),('kind','approximation')):
            task=deepcopy(self.task);task[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):
                ctx.assemble(self.powers,self.coefficients,task=task)
        task=deepcopy(self.task);task['budget']={'wall_seconds':1};task['id']='other-label'
        self.assertEqual(ctx.assemble(self.powers,self.coefficients,task=task)['certificate'],self.original)

    def test_05_gap_structural_binding(self):
        edits=[('mean_zero',True),('azimuth_m',1),('azimuth_m',False),
               ('geometry','unit_S3'),('scope','even_subspace'),('eigenvalue_index',1)]
        for key,value in edits:
            gap=deepcopy(self.gap);gap[key]=value
            with self.subTest(key=key,value=value),self.assertRaises(ValueError):
                c.prepare_fullspace(self.task,gap)
        gap=deepcopy(self.gap);gap['function']['axis']=1
        with self.assertRaises(ValueError):c.prepare_fullspace(self.task,gap)

    def test_06_snapshots_and_returned_copies_do_not_alias(self):
        gap=deepcopy(self.gap);task=deepcopy(self.task)
        ctx=c.prepare_fullspace(task,gap);prepared=self.prepared()
        evaluation=ctx.evaluate(self.powers,self.coefficients,prepared=prepared,statistics=self.stats)
        gap['beta']='999';task['function']['axis']=1;prepared.R[0][0]+=1
        ctx.task['function']['axis']=1;ctx.gap['beta']='999'
        evaluation.statistics['rayleigh']='999'
        first=ctx.assemble_evaluation(evaluation)
        self.assertEqual(first['certificate'],self.original)
        first['certificate']['gap']['beta']='999';first['certificate']['statistics']['rayleigh']='999'
        self.assertEqual(ctx.assemble_evaluation(evaluation)['certificate'],self.original)
        with self.assertRaises(FrozenInstanceError):evaluation.powers=()
        with self.assertRaises(FrozenInstanceError):ctx._gap_json='{}'

    def test_07_receipt_is_task_local_even_for_identical_task(self):
        a=c.prepare_fullspace(self.task,self.gap);b=c.prepare_fullspace(self.task,self.gap)
        evaluation=a.evaluate(self.powers,self.coefficients,prepared=self.prepared())
        with self.assertRaises(ValueError):b.assemble_evaluation(evaluation)

    def test_08_prepared_function_and_basis_binding(self):
        for field,value in (('function_key','wrong'),('powers',(F(0),F(2)))):
            prepared=self.prepared();setattr(prepared,field,value)
            with self.subTest(field=field),self.assertRaises(ValueError):self.pending(prepared=prepared)

    def test_09_bad_scalar_statistics_are_not_reused(self):
        for key,value in (('rayleigh','0'),('mass','1'),('verified',True)):
            stats=deepcopy(self.stats);stats[key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):
                self.pending(prepared=self.prepared(),statistics=stats)

    def test_10_symmetric_forged_matrices_remain_pending_and_fail_cold_verification(self):
        prepared=self.prepared();prepared.R[0][0]+=1
        prepared.verified=True;prepared.certificate_valid=True
        result=self.pending(prepared=prepared)
        self.assertIsNone(result['certificate_valid']);self.assertFalse(result['target_met'])
        self.assertFalse(self.cold_assess(result['certificate'])['certificate_valid'])

    def test_11_bad_gap_flag_cannot_create_validity(self):
        gap=deepcopy(self.gap);gap['verified']=True
        result=self.pending(ctx=c.prepare_fullspace(self.task,gap),prepared=self.prepared())
        self.assertIsNone(result['certificate_valid']);self.assertFalse(result['target_met'])
        self.assertFalse(self.cold_assess(result['certificate'])['certificate_valid'])

    def test_12_final_certificate_tampering_is_rejected(self):
        for key,value in (('lower','-999'),('matrix_digest','bad'),('tolerance','1/1000')):
            cert=self.pending(prepared=self.prepared())['certificate'];cert[key]=value
            with self.subTest(key=key):self.assertFalse(self.cold_assess(cert)['certificate_valid'])

    def test_13_final_verifier_rejects_wrong_original_task(self):
        cert=self.pending(prepared=self.prepared())['certificate']
        for key,value in (('mean_zero',True),('tolerance','1/1000')):
            task=deepcopy(self.task);task[key]=value
            with self.subTest(key=key):self.assertFalse(self.cold_assess(cert,task)['certificate_valid'])
        task=deepcopy(self.task);task['function']['axis']=1
        self.assertFalse(self.cold_assess(cert,task)['certificate_valid'])

    def test_14_bad_matrix_shapes_symmetry_and_inexact_input_rejected(self):
        for edit in ('shape','symmetry','float'):
            prepared=self.prepared()
            if edit=='shape':prepared.M.pop()
            elif edit=='symmetry':prepared.H[0][1]+=1
            else:prepared.R[0][0]=0.5
            with self.subTest(edit=edit),self.assertRaises(ValueError):self.pending(prepared=prepared)

    def test_15_strict_temple_condition_not_relaxed(self):
        gap=deepcopy(self.gap);gap['beta']=gap['lower']=self.stats['rayleigh']
        with self.assertRaises(ValueError):
            self.pending(ctx=c.prepare_fullspace(self.task,gap),prepared=self.prepared())


if __name__=='__main__':unittest.main(verbosity=2)
