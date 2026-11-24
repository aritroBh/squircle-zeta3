import os, pickle
k1f = '<project>/research/08_theory/K1f_hexagon_order5'

# Load E5_polys.pkl
with open(os.path.join(k1f, 'E5_polys.pkl'), 'rb') as f:
    E5_polys = pickle.load(f)

print("=== E5_polys type ===")
print(type(E5_polys))
if isinstance(E5_polys, dict):
    print("keys:", list(E5_polys.keys()))
    for k in list(E5_polys.keys()):
        v = E5_polys[k]
        print(f"  E5_polys[{k}]: type={type(v)}, ", end="")
        if isinstance(v, dict):
            print(f"len={len(v)}, sample keys={list(v.keys())[:3]}")
        else:
            print(f"value={v}")
