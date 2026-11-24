
import os, pickle, json
base='research/08_theory'
k1d=os.path.join(base,'K1d_engine_fixed')
with open(os.path.join(k1d,'Phi_syms.pkl'),'rb') as f:
    Phi=pickle.load(f)
print(type(Phi))
for k,v in (Phi.items() if isinstance(Phi,dict) else enumerate(Phi)):
    print(k, ':', str(v)[:400])
print('---- jloc ----')
jl=json.load(open(os.path.join(k1d,'jloc.json')))
print(len(jl))
items=sorted(jl.items())[:6]
for k,v in items: print(k,v)
