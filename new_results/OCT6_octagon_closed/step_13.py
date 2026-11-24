import pickle
K1d="research/08_theory/K1d_engine_fixed"
Jtab = pickle.load(open(K1d+'/Jtab_t.pkl','rb'))
print(type(Jtab), len(Jtab))
ks = sorted(Jtab.keys())
print("keys sample:", ks[:8])
import sympy as sp
for k in ks[:6]:
    print(k, '->', Jtab[k])
print()
print(open(K1d+'/esum_engine.py').read())
