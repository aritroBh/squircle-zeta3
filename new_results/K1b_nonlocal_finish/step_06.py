import os, json
d = "<project>/research/08_theory/K1_nonlocal_law/"
print(sorted(os.listdir(d)))
with open(d+"atom_values.json") as f:
    print(f.read())
