import unittest,itertools,math
import numpy as np,pandas as pd,shap
from frozen_explainer import explain_row,kernel_explainer

class FrozenRegressionTests(unittest.TestCase):
    stack=None;train=None;external=None;saved=None
    def test_saved_predictions(self):
        self.assertLess(np.max(abs(self.stack.predict_full(self.external)-self.saved)),1e-10)
    def test_strict_columns(self):
        with self.assertRaises(ValueError):self.stack.predict_full(self.train.iloc[:2,::-1])
        with self.assertRaises(ValueError):self.stack.predict_full(self.train.iloc[:2,:-1])
    def test_null_function_and_reduced_interface(self):
        x=self.external.iloc[:12].copy();p=self.stack.predict_full(x)
        x.loc[:,~self.stack.mask]=np.arange((~self.stack.mask).sum())*1000
        self.assertLess(np.max(abs(p-self.stack.predict_full(x))),1e-12)
        self.assertLess(np.max(abs(p-self.stack.predict_active(x.loc[:,self.stack.active]))),1e-12)
    def test_dimensions_dummy_additivity(self):
        bg=self.train.sample(20,random_state=21);ex=kernel_explainer(self.stack,bg)
        x=self.external.iloc[:2];s=explain_row(ex,self.stack,x,512,20261021)
        self.assertEqual(s.shape,(2,30));self.assertEqual(list(s.columns),self.stack.columns)
        self.assertEqual(float(abs(s.loc[:,~self.stack.mask]).to_numpy().max()),0.)
        self.assertLess(np.max(abs(float(ex.expected_value)+s.sum(1).to_numpy()-self.stack.predict_full(x))),1e-8)
    def test_exact_nonlinear_coalition_reference(self):
        # A logistic stacking-style nonlinear function, not an additive toy.
        def f(a):return 1/(1+np.exp(-(a[:,0]*a[:,1]+.7*a[:,2]-.4)))
        bg=np.array([[0.,1.,-.5],[1.,-.5,.8]]);x=np.array([.4,.7,1.2]);v={}
        for bits in range(8):
            z=bg.copy()
            for j in range(3):
                if bits&(1<<j):z[:,j]=x[j]
            v[bits]=f(z).mean()
        exact=np.zeros(3)
        for j in range(3):
            for bits in range(8):
                if bits&(1<<j):continue
                k=bits.bit_count();w=math.factorial(k)*math.factorial(2-k)/math.factorial(3)
                exact[j]+=w*(v[bits|(1<<j)]-v[bits])
        ex=shap.KernelExplainer(f,bg,link='identity');got=np.asarray(ex.shap_values(x[None,:],nsamples=8,l1_reg=0,silent=True))[0]
        self.assertLess(np.max(abs(got-exact)),1e-10)

def run_tests(stack,train,external,saved,stream):
    FrozenRegressionTests.stack=stack;FrozenRegressionTests.train=train;FrozenRegressionTests.external=external;FrozenRegressionTests.saved=saved
    return unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(FrozenRegressionTests))
