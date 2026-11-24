import os, pickle, sys
base = '<project>'
WTS = os.path.join(base, 'research/08_theory/WTS_weight_k_structure/theory')
TOW1 = os.path.join(base, 'research/08_theory/TOW1_tower_law')
sys.path.insert(0, TOW1)
sys.path.insert(0, WTS)

# Look at binteg and tower pickle structure
b3 = pickle.load(open(os.path.join(WTS, 'binteg_3.pkl'), 'rb'))
print("binteg_3 keys:", list(b3.keys()))
P = b3[(2,1)]
print("type P:", type(P), "len", len(P))
print("real:", P[0])
print("ipos:", P[1])
print("ineg:", P[2])
print("const:", P[3])
print()
t3 = pickle.load(open(os.path.join(WTS, 'tower_3.pkl'), 'rb'))
print("tower_3 keys:", list(t3.keys()))
print("tower_3[(2,1)]:", t3.get((2,1)))
print("tower_3[(3,0)]:", t3.get((3,0)))
