import os, pickle, time
import sympy as sp
base = '<project>/research/08_theory'
work = base + '/TOW4_tower_depth8'
Ls = pickle.load(open(work + '/L_simplified.pkl','rb'))
print("L_simplified keys:", sorted(Ls.keys()))
d, m = sp.symbols('d m')
for k in sorted(Ls.keys()):
    v = Ls[k]
    e = sp.sympify(v, locals={'d':d,'m':m,'I':sp.I}) if isinstance(v,str) else v
    print(k, type(v).__name__, "len", len(str(v)), str(e)[:100])
print("="*60)
Lp = pickle.load(open(work + '/L_parsed.pkl','rb'))
print("L_parsed keys:", sorted(Lp.keys()))
for k in sorted(Lp.keys()):
    v = Lp[k]
    print(k, type(v).__name__, "len", len(str(v)))
