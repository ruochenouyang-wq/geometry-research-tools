import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import evidence_store as e

ROOT = Path(__file__).resolve().parent


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='.evidence-test-', dir=ROOT)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.store = e.EvidenceStore(self.root)

    def test_298_full_content_dedup(self):
        a = {'format':'candidate', 'nested':{'upper':'2', 'lower':'1'}}
        b = {'nested':{'lower':'1', 'upper':'2'}, 'format':'candidate'}
        self.assertEqual(self.store.put(a), self.store.put(b))
        self.assertEqual(len(list(self.root.glob('*.json'))), 1)
        c = {'format':'candidate', 'nested':{'upper':'2', 'lower':'0'}}
        self.assertNotEqual(self.store.put(c), self.store.put(a))
        self.assertEqual(self.store.get(self.store.put(a)), a)

    def test_299_strict_json_and_complete_summary(self):
        for bad in ({1:'numeric key'}, {'x':float('nan')}, {'x':float('inf')}, ('tuple',)):
            with self.subTest(bad=str(bad)):
                with self.assertRaises(ValueError): self.store.put(bad)
        a={'format':'proposal','large_evidence':list(range(1000))}
        b={'format':'proposal','large_evidence':list(range(999))+[-1]}
        sa=e.summary(a); sb=e.summary(b)
        self.assertNotEqual(sa['sha256'],sb['sha256'])
        self.assertFalse(sa['summary_is_certificate'])
        self.assertTrue(sa['complete_content_hashed'])
        self.assertEqual(sa,self.store.summary(self.store.put(a)))
        self.assertNotEqual(e.digest({'n':True}),e.digest({'n':1}))
        self.assertNotEqual(e.digest({'n':1.0}),e.digest({'n':1}))

    def test_300_tamper_and_failed_atomic_publish(self):
        value={'format':'proposal','bound':'7/8'}
        ref=self.store.put(value)
        (self.root/(ref+'.json')).write_text('{"bound":"1","format":"proposal"}')
        with self.assertRaises(ValueError): self.store.get(ref)
        with self.assertRaises(ValueError): self.store.put(value)
        new={'format':'proposal','bound':'8/9'}
        with patch.object(e.os,'link',side_effect=OSError('injected publish failure')):
            with self.assertRaises(OSError): self.store.put(new)
        self.assertFalse((self.root/(e.digest(new)+'.json')).exists())
        self.assertFalse(list(self.root.glob('.pending-*')))
        self.assertEqual(self.store.get(self.store.put(new)),new)

    def test_301_paths_sizes_and_unknown_format(self):
        for ref in ('../escape','A'*64,'0'*63,'/tmp/example',True):
            with self.subTest(ref=ref):
                with self.assertRaises(ValueError): self.store.get(ref)
        ref='0'*64
        target=self.root/'ordinary.json'; target.write_text('{}')
        (self.root/(ref+'.json')).symlink_to(target)
        with self.assertRaises((ValueError,OSError)): self.store.get(ref)
        small=e.EvidenceStore(self.root/'small',max_bytes=64)
        with self.assertRaises(ValueError): small.put({'large':'x'*100})
        oversized=small.root/(ref+'.json'); oversized.write_bytes(b'x'*65)
        with self.assertRaises(ValueError): small.get(ref)
        unknown=self.store.put({'format':'future_math_certificate','verified':True})
        self.assertEqual(self.store.format_status(unknown)['status'],'unknown_format_unverified')
        self.assertFalse(self.store.format_status(unknown)['verified'])
        chain=[]
        for _ in range(100): chain=[chain]
        with self.assertRaises(ValueError): self.store.put(chain)


if __name__ == '__main__':
    unittest.main()
