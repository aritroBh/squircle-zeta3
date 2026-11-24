import os, pickle
k1f = '<project>/research/08_theory/K1f_hexagon_order5'

# Load jtab5.pkl
with open(os.path.join(k1f, 'jtab5.pkl'), 'rb') as f:
    jtab5 = pickle.load(f)
print("=== jtab5 type ===", type(jtab5))
if isinstance(jtab5, dict):
    print("num keys:", len(jtab5))
    print("sample keys:", list(jtab5.keys())[:5])
    # print a sample value
    for i,(k,v) in enumerate(jtab5.items()):
        print(f"key={k}: type={type(v)}")
        print(f"  value={v}")
        if i>=1: break
