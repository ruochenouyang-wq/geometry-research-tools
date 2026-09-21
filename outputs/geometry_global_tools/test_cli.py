from pathlib import Path
from tempfile import TemporaryDirectory
from contextlib import redirect_stdout,redirect_stderr
from io import StringIO
import unittest
import run


class CLITests(unittest.TestCase):
    def test_all_three_domains_roundtrip(self):
        source=Path(__file__).resolve().parent/'examples'/'tilted.json'
        with TemporaryDirectory() as tmp,redirect_stdout(StringIO()),redirect_stderr(StringIO()):
            for i,args in enumerate([['--domain','box'],['--all-real-method','growth'],[]]):
                out=str(Path(tmp)/(str(i)+'.json'))
                self.assertEqual(run.main(['optimize','--problem',str(source),'--out',out]+args),0)
                self.assertEqual(run.main(['verify',out]),0)

    def test_budget_failure_does_not_mean_invalid_bounds(self):
        source=Path(__file__).resolve().parent/'examples'/'tilted.json'
        with TemporaryDirectory() as tmp,redirect_stdout(StringIO()),redirect_stderr(StringIO()):
            out=str(Path(tmp)/'failed.json')
            self.assertEqual(run.main(['optimize','--problem',str(source),'--domain','box','--max-leaves','1','--out',out]),2)
            self.assertEqual(run.main(['verify',out]),0)

    def test_empty_or_unknown_documents_rejected(self):
        for doc in [None,{},[],{'certificate':None},{'method':'unknown'}]:self.assertFalse(run.verify_document(doc))


if __name__=='__main__':unittest.main()
