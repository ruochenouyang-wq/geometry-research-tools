import json
from pathlib import Path
import tempfile
import unittest
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


if __name__ == '__main__':
    unittest.main()
