import os, pickle
k1f = '<project>/research/08_theory/K1f_hexagon_order5'

# Load jloc_sym.pkl
with open(os.path.join(k1f, 'jloc_sym.pkl'), 'rb') as f:
    jloc = pickle.load(f)
print("=== jloc_sym type ===", type(jloc))
if isinstance(jloc, dict):
    print("num keys:", len(jloc))
    print("sample keys:", list(jloc.keys())[:8])
    for i,(k,v) in enumerate(jloc.items()):
        print(f"key={k}: type={type(v)}, value={v}")
        if i>=3: break
