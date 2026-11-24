# step_61: inspect brackets.pkl for edge 0 W3
import os, sys, pickle
base = "<project>"
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT9b_octagon_debug"))
sys.path.insert(0, os.path.join(base,"research/08_theory/OCT9c_octagon_close"))
import oct9c_common as oc
import sympy as sp
D = oc.load_all()
b = D['brackets']
print(type(b), list(b.keys()) if isinstance(b, dict) else len(b))
for k in (list(b.keys()) if isinstance(b,dict) else range(len(b))):
    v = b[k]
    print("KEY", k, type(v))
    if isinstance(v, dict):
        for kk in v: print("   ", kk, type(v[kk]), repr(v[kk])[:120])
    else:
        print(repr(v)[:300])
