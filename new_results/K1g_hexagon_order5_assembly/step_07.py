import os, pickle
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f, 'E5_polys.pkl'), 'rb') as f:
    E5_polys = pickle.load(f)

E5 = E5_polys[5]
print("=== E5_polys[5] sample entries ===")
for i, (k, v) in enumerate(E5.items()):
    print(f"  key={k}: type(v)={type(v)}")
    print(f"    value = {v}")
    if i >= 4:
        break
