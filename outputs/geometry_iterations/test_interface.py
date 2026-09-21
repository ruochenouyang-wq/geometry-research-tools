from copy import deepcopy
from fractions import Fraction as F
from pathlib import Path
from tempfile import TemporaryDirectory
from contextlib import redirect_stdout,redirect_stderr
from io import StringIO
import unittest
import run
from check_document import verify_document

class InterfaceTests(unittest.TestCase):
    def test_saved_spectrum_roundtrip_and_false_success(self):
        result=run.v10.adaptive([0,20,-100,0,1],max_modes=6)
        self.assertTrue(verify_document(result));bad=deepcopy(result);bad['status']='target_met'
        self.assertFalse(verify_document(bad))

    def test_all_public_routes_roundtrip(self):
        cases=[['spectrum','--q','0,1'],['range','--q','1/9,0,-2/3,0,1'],
               ['refine','--q','0,1'],['cluster','--q','0,1'],
               ['family','--q0','0','--direction','0,1','--threshold=-3/10'],
               ['function','--q','0,1']]
        with TemporaryDirectory() as tmp,redirect_stdout(StringIO()),redirect_stderr(StringIO()):
            for i,args in enumerate(cases):
                out=str(Path(tmp)/(str(i)+'.json'))
                self.assertEqual(run.main(args+['--out',out]),0,args)
                self.assertEqual(run.main(['verify',out]),0,args)

    def test_budget_failure_exit_code_and_certificate(self):
        with TemporaryDirectory() as tmp,redirect_stdout(StringIO()):
            out=str(Path(tmp)/'failed.json')
            self.assertEqual(run.main(['spectrum','--q','0,20,-100,0,1','--max-modes','6','--out',out]),2)
            self.assertEqual(run.main(['verify',out]),0)

    def test_malformed_document_rejected(self):
        for value in [None,[],{}, {'method':'unknown'}, {'certificate':None}]:self.assertFalse(verify_document(value))

    def test_retained_baseline_certificates(self):
        import previous_run
        self.assertTrue(verify_document(run.v3.base.certify([0,1],1,8,32)))
        for q in [[0,1],[0,20,-100]]:
            result=previous_run.search(q);self.assertTrue(verify_document(result))
            bad=deepcopy(result);bad['certificate']['lower']='999';self.assertFalse(verify_document(bad))

if __name__=='__main__':unittest.main()
