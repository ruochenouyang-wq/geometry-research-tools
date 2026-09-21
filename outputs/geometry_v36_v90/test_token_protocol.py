"""Independent, offline V36–V55 protocol regression tests.

Run from this directory: python3 -m unittest -v test_token_protocol
All writable evidence stores are TemporaryDirectory instances.  No baseline
results are loaded and no numerical benchmark or model API is invoked.
"""
import copy
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import tempfile
import unittest
from fractions import Fraction as F
from unittest.mock import patch

import compiler_v23 as compiler
import family_v29 as family
import legacy_cli as verifier
import planner_v35 as planner
import protocol as p
import token_meter as meter


OPTIONS = {'max_leaves': 1, 'max_modes': 4, 'max_states': 4}


def goal(**changes):
    raw = {'geometry': 'unit_sphere', 'function_space': 'full_sphere',
           'eigenvalue_index': 1, 'parameters': ['a'], 'potential': 'a*t',
           'penalty': 'a*a/4', 'domain': {'kind': 'all_real'}, 'threshold': '0'}
    raw.update(copy.deepcopy(changes))
    return raw


def profiled(**changes):
    raw = goal(**changes)
    for key in p.PROFILE_FIELDS:
        raw.pop(key)
    return {'profile': p.PROFILE, **raw}


class EvidenceFixture(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.raws = {
            'proved': goal(),
            'unresolved': goal(penalty='a*a/5', domain={'kind': 'box', 'box': [[-2, 2]]}),
            'refuted': goal(penalty='a*a/10', threshold='-1/10',
                            domain={'kind': 'box', 'box': [[-4, 4]]}),
        }
        cls.results = {name: planner.solve(raw, **OPTIONS) for name, raw in cls.raws.items()}
        for name, result in cls.results.items():
            if result['certificate']['status'] != name or not verifier.verify(result['certificate']):
                raise AssertionError('Invalid independent fixture: ' + name)
        cls.template = family.synthesize([0, 1])

    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix='geometry-token-audit-')
        self.addCleanup(self.temporary.cleanup)
        self.service = p.Service(self.temporary.name)

    def receipt(self, status='proved', **kwargs):
        ref = self.service.register(self.raws[status])
        return self.service.receipt(self.results[status], ref, **kwargs)

    def file_for(self, ref):
        return Path(self.temporary.name) / (ref.replace(':', '_') + '.json')


class TokenMeterTests(unittest.TestCase):
    def test_v36_real_known_vocabulary_tokens_offline(self):
        with patch.dict(os.environ, {'TIKTOKEN_CACHE_DIR': str(Path(__file__).parent / 'token_cache')}), \
                patch.object(socket.socket, 'connect', side_effect=AssertionError('Network forbidden')):
            self.assertEqual(meter.encoder('cl100k_base').encode('hello world'), [15339, 1917])
            self.assertEqual(meter.encoder('o200k_base').encode('hello world'), [24912, 2375])
            self.assertEqual(meter.counts('hello world'), {'o200k_base': 2, 'cl100k_base': 2})

    def test_v36_canonical_wire_is_order_independent_and_lossless(self):
        left = {'汉字': [None, True, '1/3'], 'a': 'quotes " and \\ and\nline'}
        right = dict(reversed(list(left.items())))
        self.assertEqual(meter.wire(left), meter.wire(right))
        self.assertEqual(json.loads(meter.wire(left)), left)
        with self.assertRaises(ValueError):
            meter.wire({'x': float('nan')})

    def test_v36_conversation_counts_every_explicit_item_without_api_claim(self):
        events = [{'text': '题目 a*t'}, {'text': '结果 proved'}, {'text': 'context:abc'}]
        report = meter.conversation(events, instructions='instructions', schema='schema')
        for encoding in meter.ENCODINGS:
            expected = sum(meter.count(t, encoding) for t in ['instructions', 'schema'] + [e['text'] for e in events])
            self.assertEqual(report['encodings'][encoding], expected)
        self.assertEqual(report['items'], 5)
        self.assertIsNone(report['actual_api_usage'])
        self.assertFalse(report['model_mapping_asserted'])
        self.assertIn('hidden reasoning', report['basis'])

    def test_v36_unknown_encoding_is_rejected(self):
        for name in ('model_name', 'character_estimate', '', None):
            with self.subTest(name=name), self.assertRaises(ValueError):
                meter.count('test', name)


class DisplayTests(unittest.TestCase):
    def test_v38_positive_and_negative_bounds_are_outward(self):
        samples = [F(n, d) for n in (-99999999, -10000001, -7, -1, 1, 7, 10000001, 99999999)
                   for d in (3, 7, 10**12)]
        for value in samples:
            for digits in (2, 7, 18):
                with self.subTest(value=value, digits=digits):
                    self.assertLessEqual(F(p.outward(value, digits=digits)), value)
                    self.assertGreaterEqual(F(p.outward(value, upper=True, digits=digits)), value)

    def test_v38_exact_independent_rounding_examples(self):
        for value, lower, upper in [('1/3', '0.3333333', '0.3333334'),
                                     ('-1/3', '-0.3333334', '-0.3333333'),
                                     ('99999995', '99999990', '100000000'),
                                     ('-99999995', '-100000000', '-99999990')]:
            self.assertEqual(F(p.outward(value)), F(lower))
            self.assertEqual(F(p.outward(value, upper=True)), F(upper))

    def test_v38_small_and_large_magnitudes_have_tight_enclosures(self):
        for value in (F(1, 3 * 10**30), F(-1, 3 * 10**30), F(10**40 + 1, 3)):
            lower, upper = F(p.outward(value)), F(p.outward(value, upper=True))
            self.assertLessEqual(lower, value)
            self.assertGreaterEqual(upper, value)
            self.assertLessEqual(upper - lower, abs(value) / 10**5)

    def test_v38_null_zero_and_exact_powers_are_preserved(self):
        self.assertIsNone(p.outward(None))
        self.assertEqual(p.outward(0), '0')
        for value in (F(10**20), F(1, 10**20), F(-1000), F(1, 1000)):
            self.assertEqual(F(p.outward(value)), value)
            self.assertEqual(F(p.outward(value, upper=True)), value)

    def test_v38_invalid_precision_is_rejected(self):
        for digits in (True, 1, 19, 7.0, '7'):
            with self.subTest(digits=digits), self.assertRaises(ValueError):
                p.outward('1/3', digits=digits)


class InputTests(EvidenceFixture):
    def test_v42_explicit_profile_matches_full_goal(self):
        self.assertEqual(p.expand_request(profiled()), goal())
        self.assertEqual(self.service.register(profiled()), self.service.register(goal()))

    def test_v42_scope_domain_and_threshold_cannot_be_implicitly_defaulted(self):
        for field in ('geometry', 'function_space', 'eigenvalue_index', 'domain', 'threshold'):
            raw = goal(); raw.pop(field)
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.service.register(raw)

    def test_v42_profile_conflicts_and_unsupported_scope_are_rejected(self):
        for change in ({'profile': 'unknown'}, {'function_space': 'axisymmetric'},
                       {'geometry': 'other_surface'}, {'eigenvalue_index': 2}):
            raw = profiled(); raw.update(change)
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.service.register(raw)

    def test_v43_polynomial_packing_is_exact_through_degree_six(self):
        raw = profiled(parameters=['a', 'b'])
        raw.pop('potential')
        raw['potential_coefficients'] = {'base': ['-1/3', 0, 2, 0, 0, 0, '1/7'],
                                         'directions': [[0, 1], ['2/3', 0, -1]]}
        expected = goal(parameters=['a', 'b'], potential='-1/3+2*t*t+t**6/7+a*t+b*(2/3-t*t)')
        self.assertEqual(compiler.compile_goal(p.expand_request(raw))['problem'],
                         compiler.compile_goal(expected)['problem'])
        self.assertEqual(p.polynomial([0, 0]), '0')

    def test_v44_packed_hessian_matches_cross_term_convention(self):
        raw = profiled(parameters=['a', 'b', 'c'], potential='a*t+b*t*t+c*t**3')
        raw.pop('penalty')
        raw['quadratic'] = {'constant': '2/3', 'linear': [1, '-1/2', 0],
                            'hessian_upper': [2, 3, '-2/5', 4, 7, 6]}
        expected = goal(parameters=['a', 'b', 'c'], potential=raw['potential'],
                        penalty='2/3+a-b/2+a*a+3*a*b-2*a*c/5+2*b*b+7*b*c+3*c*c')
        actual = compiler.compile_goal(p.expand_request(raw))['problem']
        self.assertEqual(actual, compiler.compile_goal(expected)['problem'])
        self.assertEqual(actual['penalty']['hessian'], [['2', '3', '-2/5'], ['3', '4', '7'], ['-2/5', '7', '6']])

    def test_v43_v44_zero_parameter_packed_problem_is_exact(self):
        raw = profiled(parameters=[], domain={'kind': 'box', 'box': []})
        raw.pop('potential'); raw.pop('penalty')
        raw['potential_coefficients'] = {'base': [0, 0, '2/3'], 'directions': []}
        raw['quadratic'] = {'constant': '1/7', 'linear': [], 'hessian_upper': []}
        expected = goal(parameters=[], potential='2*t*t/3', penalty='1/7',
                        domain={'kind': 'box', 'box': []})
        self.assertEqual(compiler.compile_goal(p.expand_request(raw))['problem'],
                         compiler.compile_goal(expected)['problem'])

    def test_v43_v44_ambiguous_or_wrong_dimensions_are_rejected(self):
        invalid = []
        for coefficients in ({'base': [0], 'directions': []}, {'base': [0]*8, 'directions': [[1]]}):
            raw = profiled(); raw.pop('potential'); raw['potential_coefficients'] = coefficients; invalid.append(raw)
        raw = profiled(); raw['potential_coefficients'] = {'base': [0], 'directions': [[0, 1]]}; invalid.append(raw)
        for quadratic in ({'hessian_upper': []}, {'hessian_upper': [1], 'linear': []}):
            raw = profiled(); raw.pop('penalty'); raw['quadratic'] = quadratic; invalid.append(raw)
        raw = profiled(); raw['quadratic'] = {'hessian_upper': [1]}; invalid.append(raw)
        for raw in invalid:
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                self.service.register(raw)

    def test_v43_packed_parameter_names_do_not_enable_expression_injection(self):
        raw = profiled(parameters=['a).__class__'])
        raw.pop('potential'); raw['potential_coefficients'] = {'base': [0], 'directions': [[0, 1]]}
        with self.assertRaises((ValueError, SyntaxError)):
            self.service.register(raw)

    def test_v45_patch_preserves_original_and_every_unstated_field(self):
        original = goal(); old_ref = self.service.register(original)
        new_ref = self.service.revise(old_ref, {'threshold': '-1/10', 'penalty': 'a*a/3'})
        revised = self.service.inspect(new_ref)
        self.assertNotEqual(old_ref, new_ref)
        self.assertEqual(self.service.inspect(old_ref), original)
        self.assertEqual(revised, goal(threshold='-1/10', penalty='a*a/3'))

    def test_v45_patch_revalidates_domain_scope_and_parameter_dependencies(self):
        ref = self.service.register(goal())
        for changes in ({}, {'eigenvalue_index': 2}, {'potential': 'a*a*t'},
                        {'parameters': []}, {'domain': {'kind': 'box', 'box': []}},
                        {'function_space': 'unmentioned_subspace'}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                self.service.revise(ref, changes)


class ReferenceTests(EvidenceFixture):
    def test_v39_reference_is_canonical_persistent_and_detached_from_input(self):
        raw = goal(); ref = self.service.register(raw)
        self.assertEqual(ref, self.service.register(dict(reversed(list(raw.items())))))
        raw['parameters'][0] = 'changed'
        reopened = p.Store(self.temporary.name)
        self.assertEqual(reopened.get(ref), goal())
        received = reopened.get(ref); received['parameters'][0] = 'also_changed'
        self.assertEqual(reopened.get(ref), goal())

    def test_v39_invalid_and_traversal_references_are_rejected(self):
        for ref in ('../secret', '/etc/passwd', 'goal:../secret', 'goal:' + 'a'*23,
                    'goal:' + 'A'*24, 'unknown:' + 'a'*24, 7, None):
            with self.subTest(ref=ref), self.assertRaises(ValueError):
                self.service.store.get(ref)

    def test_v39_stored_content_change_is_rejected_on_read_and_put(self):
        raw = goal(); ref = self.service.register(raw)
        self.file_for(ref).write_text('{}')
        with self.assertRaises(ValueError): self.service.inspect(ref)
        with self.assertRaises(ValueError): self.service.register(raw)

    def test_v40_projection_supports_rfc6901_escapes_and_exact_nested_evidence(self):
        value = {'a/b': {'~c': [{'0': 'exact'}]}, '': 11}
        self.assertEqual(p.project(value, '/a~1b/~0c/0/0'), 'exact')
        self.assertEqual(p.project(value, '/'), 11)
        r = self.receipt()
        self.assertEqual(self.service.inspect(r['evidence'], '/certificate/compiled_goal/source'), goal())

    def test_v40_projection_rejects_invalid_paths_and_array_indices(self):
        for pointer in ('bad', '/x/01', '/x/-1', '/x/1.0', '/x/0/child'):
            with self.subTest(pointer=pointer), self.assertRaises((ValueError, KeyError, IndexError)):
                p.project({'x': [1]}, pointer)

    def test_v40_projection_rejects_invalid_escape_sequences(self):
        for pointer in ('/~2', '/~', '/a~3b'):
            with self.subTest(pointer=pointer), self.assertRaises(ValueError):
                p.project({'~2': 1, '~': 2, 'a~3b': 3}, pointer)

    def test_v41_pagination_reconstructs_every_item_with_explicit_completion(self):
        original = list(range(19)); ref = self.service.store.put('context', original)
        rebuilt = []; offset = 0
        while offset is not None:
            chunk = self.service.inspect(ref, offset=offset, limit=4)
            self.assertEqual(chunk['offset'], offset)
            self.assertEqual(chunk['total'], len(original))
            rebuilt.extend(chunk['items']); offset = chunk['next']
        self.assertEqual(rebuilt, original)
        self.assertEqual(p.page(original, len(original), 4)['items'], [])
        self.assertIsNone(p.page([], 0, 4)['next'])

    def test_v41_invalid_page_requests_are_not_silently_clipped(self):
        for value, offset, limit in (({}, 0, 2), ([1], -1, 2), ([1], 2, 2),
                                     ([1], 0, 0), ([1], 0, 65), ([1], True, 1), ([1], 0, True)):
            with self.subTest(offset=offset, limit=limit), self.assertRaises(ValueError):
                p.page(value, offset, limit)


class VerificationTests(EvidenceFixture):
    def test_v37_all_three_decisions_survive_compression_and_exact_retrieval(self):
        for status in self.results:
            with self.subTest(status=status):
                r = self.receipt(status); evidence = self.service.inspect(r['evidence'])
                self.assertEqual(r['status'], status)
                self.assertEqual(evidence, p.checked(self.results[status]))
                self.assertTrue(verifier.verify(evidence['certificate']))
                self.assertEqual(self.service.inspect(r['task']), self.raws[status])
                self.assertTrue(r['verified']); self.assertFalse(r['formal_kernel'])
                for key, bound in r['bounds'].items():
                    exact = evidence['certificate'][key + '_bound']
                    if exact is not None:
                        self.assertTrue(F(bound) <= F(exact) if key == 'lower' else F(bound) >= F(exact))

    def test_v37_untrusted_wrapper_verdict_and_search_text_are_discarded(self):
        result = copy.deepcopy(self.results['unresolved'])
        result.update({'status': 'proved', 'claim': 'all theorems proved', 'trace': ['fabricated']})
        compact = p.summary(result)
        self.assertEqual(compact['status'], 'unresolved')
        for field in ('claim', 'trace', 'strategy', 'seconds', 'statistics'):
            self.assertNotIn(field, compact)

    def test_v37_independent_certificate_status_bound_scope_and_proof_tampering_rejected(self):
        mutations = [lambda c: c.update(status='refuted'), lambda c: c.update(lower_bound='100'),
                     lambda c: c['compiled_goal']['source'].update(threshold='1'),
                     lambda c: c['compiled_goal']['source'].update(function_space='axisymmetric'),
                     lambda c: c['terminal']['proof']['template']['poisson_solutions'][0].__setitem__(1, '-9')]
        for mutate in mutations:
            result = copy.deepcopy(self.results['proved']); mutate(result['certificate'])
            self.assertFalse(verifier.verify(result['certificate']))
            with self.assertRaises(ValueError): p.summary(result)

    def test_v37_compilation_only_is_not_reported_as_a_decision(self):
        with self.assertRaises(ValueError):
            p.summary({'certificate': compiler.compile_goal(goal())})

    def test_v37_original_scope_and_domain_are_visible(self):
        for status in self.results:
            r = self.receipt(status)
            self.assertEqual(r['scope']['geometry'], 'unit_sphere')
            self.assertEqual(r['scope']['space'], 'full_sphere')
            self.assertEqual(r['scope']['eigenvalue_index'], 1)
            self.assertEqual(r['scope']['domain'], self.raws[status]['domain'])
            self.assertEqual(F(r['target']), F(self.raws[status]['threshold']))

    def test_v39_receipt_rejects_valid_evidence_bound_to_a_different_original_goal(self):
        for change in ({'threshold': '-1'}, {'penalty': 'a*a/3'}, {'potential': '2*a*t'},
                       {'function_space': 'axisymmetric'}, {'domain': {'kind': 'box', 'box': [[-2, 2]]}}):
            different = self.service.register(goal(**change))
            with self.subTest(change=change), self.assertRaises(ValueError):
                self.service.receipt(self.results['proved'], different)

    def test_v39_receipt_requires_goal_kind_even_if_content_matches(self):
        wrong_ref = self.service.store.put('context', goal())
        with self.assertRaises(ValueError):
            self.service.receipt(self.results['proved'], wrong_ref)

    def test_v50_diagnostic_tampering_and_unrelated_problem_transplant_are_rejected(self):
        result = copy.deepcopy(self.results['proved'])
        result['diagnostics'][0]['scale_upper'] = '0'
        with self.assertRaises(ValueError): p.checked(result)
        unrelated = planner.solve(goal(potential='a*t**3', penalty='a*a/100',
                                        domain={'kind': 'box', 'box': [[-1, 1]]}), **OPTIONS)
        result = copy.deepcopy(self.results['proved']); result['diagnostics'] = unrelated['diagnostics']
        self.assertTrue(all(planner.verify_budget(d) for d in result['diagnostics']))
        with self.assertRaises(ValueError): p.checked(result)

    def test_v50_diagnostic_group_counts_and_enclosures_do_not_decide_the_goal(self):
        items = copy.deepcopy(self.results['unresolved']['diagnostics'] * 3)
        groups = p.aggregate_diagnostics(items)
        self.assertEqual(sum(g['count'] for g in groups), len(items))
        for group in groups:
            self.assertEqual(group['scope'], 'fixed_poisson_budget')
            members = [d for d in items if d['status'] == group['status']]
            self.assertLessEqual(F(group['scale_lower']), min(F(d['scale_lower']) for d in members))
            self.assertGreaterEqual(F(group['scale_upper']), max(F(d['scale_upper']) for d in members))
        self.assertEqual(self.receipt('unresolved')['status'], 'unresolved')


class ServiceTests(EvidenceFixture):
    def test_v47_identical_calls_reuse_checked_result_without_solver_work(self):
        ref = self.service.register(goal())
        with patch.object(p.planner, 'solve', wraps=p.planner.solve) as solve:
            first = self.service.solve(ref, {'max_leaves': 1, 'max_modes': 4})
            second = self.service.solve(ref, {'max_modes': 4, 'max_leaves': 1})
            self.assertEqual(first, second); self.assertEqual(solve.call_count, 1)
            second['status'] = 'changed'
            self.assertEqual(self.service.solve(ref, {'max_leaves': 1, 'max_modes': 4}), first)
        self.assertEqual(self.service.calls, 1)

    def test_v47_cache_rechecks_stored_proof(self):
        ref = self.service.register(goal()); receipt = self.service.solve(ref, OPTIONS)
        self.file_for(receipt['evidence']).write_text('{}')
        with self.assertRaises(ValueError): self.service.solve(ref, OPTIONS)
        self.assertEqual(self.service.calls, 1)

    def test_v47_solver_options_and_reference_kind_are_validated(self):
        ref = self.service.register(goal())
        with self.assertRaises(ValueError): self.service.solve(ref, {'fake_option': True})
        wrong = self.service.store.put('context', goal())
        with self.assertRaises(ValueError): self.service.solve(wrong)
        self.assertEqual(self.service.calls, 0)

    def test_v46_v48_exception_filter_retains_complete_inventory_and_exact_evidence(self):
        refs = [self.service.register(self.raws[status]) for status in ('proved', 'unresolved', 'refuted')]
        out = self.service.batch(refs, OPTIONS, exceptions_only=True)
        self.assertEqual([r['status'] for r in out['inventory']], ['proved', 'unresolved', 'refuted'])
        self.assertEqual([r['index'] for r in out['attention']], [1, 2])
        full = self.service.inspect(out['batch'])
        self.assertEqual([r['task'] for r in full], refs)
        for row, r in zip(out['inventory'], full):
            self.assertEqual(row['evidence'], r['evidence'])
            self.assertTrue(verifier.verify(self.service.inspect(r['evidence'])['certificate']))
        self.assertEqual(self.service.batch(refs, OPTIONS)['items'], full)
        self.assertEqual(self.service.calls, 3)

    def test_v46_invalid_batch_size_is_rejected_before_any_solver_work(self):
        for refs in ([], ['x']*33, 'not a list'):
            with self.subTest(refs=refs), self.assertRaises(ValueError): self.service.batch(refs)
        self.assertEqual(self.service.calls, 0)

    def test_v49_checked_lemma_reuse_proves_without_synthesis_and_changes_cache_identity(self):
        ref = self.service.register(goal()); first = self.service.solve(ref, OPTIONS)
        lemma_refs = self.service.add_lemma([copy.deepcopy(self.template)])
        self.assertEqual(len(lemma_refs), 1)
        self.assertEqual(self.service.add_lemma([copy.deepcopy(self.template)]), lemma_refs)
        result = self.service.solve(ref, {**OPTIONS, 'synthesize': False})
        self.assertEqual(result['status'], 'proved')
        evidence = self.service.inspect(result['evidence'])
        self.assertEqual(evidence['certificate']['terminal']['kind'], 'transport')
        self.assertEqual(evidence['certificate']['terminal']['proof']['template'], self.template)
        self.service.solve(ref, OPTIONS)
        self.assertEqual(self.service.calls, 3)
        self.assertEqual(first['task'], result['task'])

    def test_v49_invalid_or_corrupted_lemmas_are_rejected(self):
        bad = copy.deepcopy(self.template); bad['scale'] = '100'
        with self.assertRaises(ValueError): self.service.add_lemma([bad])
        self.assertEqual(self.service.lemma_index(), [])
        ref = self.service.add_lemma([copy.deepcopy(self.template)])[0]
        self.file_for(ref).write_text('{}')
        with self.assertRaises(ValueError): self.service.inspect(ref)
        with self.assertRaises(ValueError): self.service.solve(self.service.register(goal()), OPTIONS)

    def test_v49_lemma_index_does_not_alias_caller_owned_input(self):
        raw = copy.deepcopy(self.template); ref = self.service.add_lemma([raw])[0]
        raw['directions'][0][1] = '999'
        self.assertEqual(self.service.lemma_index()[0]['directions'], self.service.inspect(ref)['directions'])

    def test_v49_returned_lemma_index_is_detached_from_internal_state(self):
        ref = self.service.add_lemma([copy.deepcopy(self.template)])[0]
        self.service.lemma_index()[0]['directions'][0][1] = '999'
        self.assertEqual(self.service.lemma_index()[0]['directions'], self.service.inspect(ref)['directions'])


class BudgetAndStateTests(EvidenceFixture):
    def test_v51_fitting_budget_reports_exact_text_tokens_and_complete_semantics(self):
        receipt = self.receipt()
        for encoding in meter.ENCODINGS:
            out = self.service.budget(receipt, 100000, encoding)
            self.assertTrue(out['fits']); self.assertEqual(json.loads(out['text']), receipt)
            self.assertEqual(out['text_tokens'], meter.count(out['text'], encoding))
            self.assertEqual(out['wire_tokens'], meter.count(meter.wire(out), encoding))

    def test_v51_too_small_budget_returns_error_without_truncated_theorem(self):
        receipt = self.receipt('unresolved'); before = copy.deepcopy(receipt)
        out = self.service.budget(receipt, 16)
        self.assertFalse(out['fits']); self.assertGreater(out['required_text_tokens'], 16)
        self.assertEqual(out['error'], 'budget_too_small_for_required_semantics')
        self.assertNotIn('text', out); self.assertEqual(receipt, before)
        self.assertEqual(self.service.inspect(receipt['evidence'])['certificate']['status'], 'unresolved')

    def test_v51_compacted_diagnostics_keep_status_and_reusable_receipt(self):
        receipt = self.receipt('unresolved')
        concise = copy.deepcopy(receipt); concise['diagnostics'] = 'available_in_evidence'
        limit = self.service.budget(concise, 100000)['wire_tokens']
        self.assertLess(limit, self.service.budget(receipt, 100000)['wire_tokens'])
        out = self.service.budget(receipt, limit)
        self.assertTrue(out['fits']); self.assertEqual(json.loads(out['text']), concise)
        self.assertEqual(self.service.poll(concise, concise['evidence']), {'unchanged': concise['evidence']})
        self.assertTrue(self.service.budget(concise, limit)['fits'])
        restored = self.service.restore(self.service.checkpoint([], [concise]))
        self.assertEqual(restored['results'][0]['status'], 'unresolved')
        self.assertEqual(restored['results'][0]['evidence'], receipt['evidence'])

    def test_v51_v52_forged_receipt_fields_and_invalid_budgets_are_rejected(self):
        for field, value in (('status', 'proved'), ('target', '-999'), ('verified', False),
                             ('scope', {}), ('bounds', {'lower': '999'})):
            r = self.receipt('unresolved'); r[field] = value
            with self.subTest(field=field), self.assertRaises(ValueError): self.service.budget(r)
            with self.subTest(field=field), self.assertRaises(ValueError): self.service.poll(r, None)
        for limit in (True, 15, 100001, 500.0):
            with self.subTest(limit=limit), self.assertRaises(ValueError): self.service.budget(self.receipt(), limit)

    def test_v51_v53_omitted_diagnostic_receipt_remains_valid_across_lifecycle(self):
        receipt = self.receipt(diagnostics='omit')
        self.assertNotIn('diagnostics', receipt)
        self.assertTrue(self.service.budget(receipt, 100000)['fits'])
        self.assertEqual(self.service.poll(receipt, receipt['evidence']), {'unchanged': receipt['evidence']})
        checkpoint = self.service.checkpoint([receipt['task']], [receipt])
        restored = self.service.restore(checkpoint)
        self.assertEqual(restored['results'][0]['status'], receipt['status'])

    def test_v52_poll_unchanged_is_explicit_and_rechecks_evidence(self):
        receipt = self.receipt()
        self.assertEqual(self.service.poll(receipt, receipt['evidence']), {'unchanged': receipt['evidence']})
        self.assertEqual(self.service.poll(receipt, 'unknown'), receipt)
        self.file_for(receipt['evidence']).unlink()
        with self.assertRaises(OSError): self.service.poll(receipt, receipt['evidence'])

    def test_v53_checkpoint_restores_exact_goals_decisions_bounds_and_library(self):
        receipts = [self.receipt(status) for status in ('proved', 'unresolved', 'refuted')]
        refs = [r['task'] for r in receipts]
        self.service.add_lemma([copy.deepcopy(self.template)])
        checkpoint = self.service.checkpoint(refs, receipts)
        fresh = p.Service(self.temporary.name); out = fresh.restore(checkpoint)
        self.assertEqual(out['goals'], refs); self.assertEqual(out['results'], receipts)
        self.assertEqual(out['lemmas'], self.service.lemma_index())
        self.assertEqual(fresh.calls, 0)
        self.assertEqual(fresh.solve(refs[0], {**OPTIONS, 'synthesize': False})['status'], 'proved')

    def test_v53_checkpoint_goal_references_require_correct_kind(self):
        wrong = self.service.store.put('context', goal())
        with self.assertRaises(ValueError): self.service.checkpoint([wrong], [])

    def test_v53_restore_requires_all_referenced_objects_and_valid_evidence(self):
        for target in ('goal', 'proof', 'lemma'):
            with self.subTest(target=target), tempfile.TemporaryDirectory(prefix='geometry-token-audit-') as directory:
                service = p.Service(directory); ref = service.register(goal())
                receipt = service.receipt(self.results['proved'], ref)
                lemma = service.add_lemma([copy.deepcopy(self.template)])[0]
                checkpoint = service.checkpoint([ref], [receipt])
                missing = {'goal': ref, 'proof': receipt['evidence'], 'lemma': lemma}[target]
                (Path(directory) / (missing.replace(':', '_') + '.json')).unlink()
                with self.assertRaises(OSError): p.Service(directory).restore(checkpoint)

    def test_v53_restore_rejects_wrong_proof_and_lemma_reference_kinds(self):
        receipt = self.receipt()
        for target in ('evidence', 'lemmas'):
            saved = {'goals': [receipt['task']], 'evidence': [], 'lemmas': []}
            value = p.checked(self.results['proved']) if target == 'evidence' else self.template
            saved[target] = [self.service.store.put('context', value)]
            ref = self.service.store.put('checkpoint', saved)
            with self.subTest(target=target), self.assertRaises(ValueError):
                p.Service(self.temporary.name).restore(ref)

    def test_v37_v53_family_and_matrix_envelope_receipts_keep_ansatz_scope(self):
        import dual_v33 as dual
        lower = dual.lower_certificate([[0, 1]], [[1]], [{'t': '0', 'matrix': [['1']]}])
        envelope = dual.certificate(self.template, lower, F(1, 1000))
        for certificate, kind in ((self.template, 'family_lemma'), (envelope, 'matrix_envelope')):
            receipt = self.service.receipt({'certificate': certificate})
            self.assertEqual(receipt['kind'], kind)
            self.assertEqual(receipt['scope']['ansatz'], 'poisson_exponential')
            self.assertEqual(receipt['scope']['amplitudes'], 'all_real')
            self.assertTrue(self.service.budget(receipt, 100000)['fits'])
            checkpoint = self.service.checkpoint([receipt['task']], [receipt])
            self.assertEqual(p.Service(self.temporary.name).restore(checkpoint)['results'], [receipt])

    def test_v54_v55_context_and_submit_are_reference_backed_and_budgeted(self):
        receipt = self.receipt('unresolved'); pending = self.service.register(goal(threshold='-1'))
        context = self.service.context([receipt], [pending])
        self.assertEqual(context['completed'], [{'task': receipt['task'], 'status': 'unresolved', 'evidence': receipt['evidence']}])
        restored = p.Service(self.temporary.name).restore(context['checkpoint'])
        self.assertEqual(restored['goals'], [pending]); self.assertEqual(restored['results'], [receipt])
        success = self.service.submit(profiled(), OPTIONS, 100000)
        self.assertTrue(success['fits']); self.assertEqual(json.loads(success['text'])['status'], 'proved')
        self.assertFalse(self.service.submit(profiled(), OPTIONS, 16)['fits'])


class CliTests(EvidenceFixture):
    def cli(self, requests):
        return subprocess.run([sys.executable, str(Path(__file__).parent / 'run.py'),
                               'serve', '--store', self.temporary.name],
                              input='\n'.join(meter.wire(r) for r in requests) + '\n',
                              text=True, capture_output=True, timeout=20, check=False)

    def test_v55_jsonl_solve_poll_checkpoint_and_context_chain(self):
        ref = self.service.register(goal())
        receipt = self.service.solve(ref, OPTIONS)
        checkpoint = self.service.checkpoint([ref], [receipt])
        requests = [{'op': 'solve', 'goal_ref': ref, 'options': OPTIONS},
                    {'op': 'solve', 'goal_ref': ref, 'options': OPTIONS},
                    {'op': 'poll', 'receipt': receipt, 'known_evidence': receipt['evidence']},
                    {'op': 'restore', 'ref': checkpoint},
                    {'op': 'context', 'receipts': [receipt], 'pending': [ref]}]
        process = self.cli(requests)
        self.assertEqual(process.returncode, 0, process.stderr)
        responses = [json.loads(line) for line in process.stdout.splitlines()]
        self.assertEqual(len(responses), len(requests))
        self.assertEqual(responses[0], receipt); self.assertEqual(responses[1], receipt)
        self.assertEqual(responses[2], {'unchanged': receipt['evidence']})
        self.assertEqual(responses[3]['results'], [receipt])
        self.assertEqual(responses[4]['completed'][0]['task'], ref)

    def test_v55_jsonl_invalid_expression_is_nonfatal_and_has_no_verdict(self):
        requests = [{'op': 'submit', 'goal': profiled(potential=text), 'options': OPTIONS}
                    for text in ('(', '1/0')]
        requests.append({'op': 'submit', 'goal': profiled(), 'options': OPTIONS, 'token_limit': 100000})
        process = self.cli(requests)
        self.assertEqual(process.returncode, 0, process.stderr)
        responses = [json.loads(line) for line in process.stdout.splitlines()]
        self.assertEqual(len(responses), 3)
        for response in responses[:2]:
            self.assertIn('error', response); self.assertIsNone(response['mathematical_verdict'])
        self.assertEqual(json.loads(responses[2]['text'])['status'], 'proved')

    def test_v55_single_call_persists_lemma_library_and_returns_exact_projection(self):
        lemma_ref = self.service.add_lemma([copy.deepcopy(self.template)])[0]
        request = {'op': 'inspect', 'ref': lemma_ref, 'pointer': '/poisson_solutions/0'}
        path = Path(self.temporary.name) / 'request.json'; path.write_text(meter.wire(request))
        process = subprocess.run([sys.executable, str(Path(__file__).parent / 'run.py'),
                                  'call', str(path), '--store', self.temporary.name],
                                 text=True, capture_output=True, timeout=20, check=False)
        self.assertEqual(process.returncode, 0, process.stderr)
        self.assertEqual(json.loads(process.stdout), ['0', '-1/2'])
        out = self.cli([{'op': 'lemma_index'}, {'op': 'submit', 'goal': profiled(),
                        'options': {**OPTIONS, 'synthesize': False}, 'token_limit': 100000}])
        self.assertEqual(out.returncode, 0, out.stderr)
        responses = [json.loads(line) for line in out.stdout.splitlines()]
        self.assertEqual(responses[0], [{'ref': lemma_ref, 'directions': [['0', '1']]}])
        self.assertEqual(json.loads(responses[1]['text'])['status'], 'proved')


if __name__ == '__main__':
    unittest.main()
