import pickle
base2 = "<project>/research/08_theory/Q1S2_vertex_close/"
VS = pickle.load(open(base2+"V_symbols.pkl","rb"))
for k in VS:
    print(k, type(VS[k]), len(VS[k]) if hasattr(VS[k],'__len__') else '')
v3 = VS['V3p']
print("V3p structure:", [ (i, len(t)) for i,t in enumerate(v3)] if isinstance(v3,(list,tuple)) else v3)
print()
T0 = v3[0]
for k,c in sorted(T0.items()): print(k, c)
