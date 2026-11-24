
import os, pickle
base='research/08_theory'
k1d=os.path.join(base,'K1d_engine_fixed')
with open(os.path.join(k1d,'Jtab_t.pkl'),'rb') as f:
    J=pickle.load(f)
print(type(J))
if isinstance(J,dict):
    ks=sorted(J.keys())
    print(len(ks))
    for kk in ks[:10]:
        print(kk, '->', str(J[kk])[:300])
