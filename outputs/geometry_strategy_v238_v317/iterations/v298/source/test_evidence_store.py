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


if __name__ == '__main__':
    unittest.main()
