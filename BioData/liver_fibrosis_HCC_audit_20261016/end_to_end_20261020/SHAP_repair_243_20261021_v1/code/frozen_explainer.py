import numpy as np
import pandas as pd
import shap

class FrozenStack:
    def __init__(self, models, meta, train):
        self.models=models;self.meta=meta;self.columns=list(train.columns)
        masks=[m.named_steps['RFE'].support_ for m in models.values()]
        assert all(np.array_equal(masks[0],m) for m in masks)
        self.mask=masks[0].copy();self.active=list(train.columns[self.mask]);self.anchor=train.median().to_numpy()
        assert len(self.columns)==30 and len(self.active)==15
        assert all(np.array_equal(m.classes_,[0,1]) for m in models.values())
    def validate(self, frame):
        if not isinstance(frame,pd.DataFrame) or list(frame.columns)!=self.columns:
            raise ValueError('Exact frozen feature names and order required')
        if not np.isfinite(frame.to_numpy()).all():raise ValueError('Nonfinite input')
    def predict_full(self, frame):
        self.validate(frame)
        b=pd.DataFrame({name:m.predict_proba(frame.to_numpy())[:,1] for name,m in self.models.items()})
        return self.meta.predict_proba(b)[:,1]
    def predict_active(self, x):
        x=np.asarray(x,float)
        if x.ndim!=2 or x.shape[1]!=15:raise ValueError('Expected 15 active features')
        full=np.tile(self.anchor,(len(x),1));full[:,self.mask]=x
        return self.predict_full(pd.DataFrame(full,columns=self.columns))
    def lift(self, active_values, index):
        v=np.asarray(active_values)
        if v.shape!=(len(index),15):raise ValueError('Unexpected attribution dimensions')
        # Exact dummy-feature reduction of the coalition game, not correction of old values.
        out=np.zeros((len(index),30));out[:,self.mask]=v
        return pd.DataFrame(out,index=index,columns=self.columns)

def kernel_explainer(stack,background):
    stack.validate(background)
    return shap.KernelExplainer(stack.predict_active,background.loc[:,stack.active].to_numpy(),link='identity')

def explain_row(explainer,stack,frame,budget,seed):
    stack.validate(frame);np.random.seed(seed)
    values=np.asarray(explainer.shap_values(frame.loc[:,stack.active].to_numpy(),nsamples=budget,l1_reg=0,silent=True))
    return stack.lift(values,frame.index)
