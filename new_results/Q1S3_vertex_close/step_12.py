base2 = "<project>/research/08_theory/Q1S2_vertex_close/"
import pickle, os
for f in ["V_symbols.pkl", "symX_V3p.pkl"]:
    print("="*70); print(f)
    d = pickle.load(open(base2+f,"rb"))
    print(type(d), list(d.keys()) if isinstance(d,dict) else len(d))
    print(d if f=="symX_V3p.pkl" else "")
