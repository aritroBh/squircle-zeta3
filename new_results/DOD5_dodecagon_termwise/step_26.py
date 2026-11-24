import pickle
base = "<project>/research/08_theory/"
D1 = base+"DOD1_dodecagon/"
sd = pickle.load(open(D1+"S_atom_decomp.pkl","rb"))
print("S_atom_decomp keys:", list(sd.keys()))
for k in sd:
    v = sd[k]
    print(" ", k, type(v), (len(v) if hasattr(v,'__len__') else v))
    if hasattr(v,'items'):
        items = list(v.items())
        print("    sample:", items[:2])
