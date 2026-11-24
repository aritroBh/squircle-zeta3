import os, pickle
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
with open(os.path.join(k1f, 'jtab5.pkl'), 'rb') as f:
    jtab5 = pickle.load(f)

# Look at several values to understand the atom notation
for k in [(0,0,0,1),(0,0,1,1),(0,0,0,2),(0,1,0,1),(0,0,0,5),(3,1,1,1),(5,1,1,3)]:
    if k in jtab5:
        print(f"J{k} = {jtab5[k]}")
        print()
