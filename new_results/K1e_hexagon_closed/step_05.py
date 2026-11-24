
import os, pickle
base='research/08_theory'
k1d=os.path.join(base,'K1d_engine_fixed')
with open(os.path.join(k1d,'E3_fixed.pkl'),'rb') as f:
    E3=pickle.load(f)
print(type(E3))
if isinstance(E3,dict):
    for k in E3: print(k, type(E3[k]))
elif isinstance(E3,(list,tuple)):
    print(len(E3))
    for e in E3[:5]: print(type(e), e if not hasattr(e,'args') else str(e)[:200])
