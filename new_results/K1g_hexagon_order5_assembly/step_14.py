import os, pickle
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f, 'atom_table.pkl'), 'rb') as f:
    atom_tab = pickle.load(f)
print("=== atom_table type ===", type(atom_tab))
print("num keys:", len(atom_tab))
for k in sorted(atom_tab.keys(), key=str):
    print(f"{k} -> {atom_tab[k]}")
