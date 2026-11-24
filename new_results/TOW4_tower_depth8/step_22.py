import os, pickle
base = '<project>/research/08_theory'
work = base + '/TOW4_tower_depth8'
ef = pickle.load(open(work + '/e_fits.pkl','rb'))
print(type(ef), list(ef.keys()) if hasattr(ef,'keys') else len(ef))
import sympy as sp
d, m = sp.symbols('d m')
for k in sorted(ef.keys()):
    v = ef[k]
    try:
        e = sp.sympify(v) if isinstance(v,str) else v
        print(k, type(v), str(e)[:120])
    except Exception as ex:
        print(k, type(v), "unparsed", str(v)[:120])
