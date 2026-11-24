
import os, pickle
base='research/08_theory'
k1d=os.path.join(base,'K1d_engine_fixed')
with open(os.path.join(k1d,'E3_fixed.pkl'),'rb') as f:
    E3=pickle.load(f)
nl,loc=E3['nl'],E3['loc']
print(len(nl),len(loc))
print(nl[0])
print(loc[0])
