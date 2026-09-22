"""Analytic and mock controls only; no benchmark case or spectral search."""
from fractions import Fraction as F
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest

PATH=Path(__file__).with_name('benchmark.py')
spec=importlib.util.spec_from_file_location('numerical_adapter_under_test',PATH)
b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)


class AdapterControls(unittest.TestCase):
    def test_existing_frozen_source_path(self):
        self.assertTrue((b.C2/'fullspace_extension.py').is_file())

    def test_nontrivial_cholesky_coordinates_and_exact_lift(self):
        # L=[[2,0],[1,1]], K=diag(5,2): lowest exact x=(-1/2,1).
        # This catches both transposed whitening and missing L^-T recovery.
        p={'mr':((F(4),F(2)),(F(2),F(2))),
           'hr':((F(20),F(10)),(F(10),F(7))),
           'C':((F(1),F(0)),(F(0),F(1)),(F(1),F(1))),
           'ss':(F(0),F(1),F(2))}
        b._full=SimpleNamespace(matvec=lambda mat,v:[sum(a*x for a,x in zip(row,v)) for row in mat])
        result=b.numpy_candidate(p,b.numpy_prepare(p),{})
        values=list(map(F,result['coefficients']))
        ratio=values[0]/values[1]
        self.assertAlmostEqual(float(ratio),-0.5,places=14)
        self.assertEqual(values[2],values[0]+values[1])
        self.assertAlmostEqual(result['approximate_ritz_value_not_a_bound'],2.0,places=14)
        self.assertEqual(result['refinement_attempts'],0)

    def test_cholesky_failure_not_repaired(self):
        import numpy as np
        p={'mr':((F(1),F(2)),(F(2),F(1))),
           'hr':((F(1),F(0)),(F(0),F(1)))}
        state={}
        with self.assertRaises(np.linalg.LinAlgError):b.numpy_candidate(p,b.numpy_prepare(p),state)
        self.assertEqual(state['stage'],'cholesky')


if __name__=='__main__':unittest.main()
