import sys, pickle
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
sys.path.insert(0, "<project>/research/08_theory/OCT5D2_octagon_depth_test")
base = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
sf2 = pickle.load(open(base+"symbol_fit2.pkl","rb"))
sf = pickle.load(open(base+"symbol_fit2_solution.pkl","rb"))
cols = sf2['cols']; x = sf['x']
print("n cols:", len(cols))
nz = [(i, cols[i], x[i]) for i in range(len(x)) if x[i] != 0]
print("nonzero:", len(nz))
for i,c,v in nz:
    print(i, c, "=", v)
# distinct col types
from collections import Counter
print(Counter(c[1] for c in cols))
