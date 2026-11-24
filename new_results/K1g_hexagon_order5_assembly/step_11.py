import os, pickle
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f, 'jloc_sym.pkl'), 'rb') as f:
    jloc = pickle.load(f)
print("all keys:", sorted(jloc.keys()))
print()
for k in sorted(jloc.keys()):
    print(f"Jloc{k} = {jloc[k]}")
