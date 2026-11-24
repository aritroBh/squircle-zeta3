import pickle
base = "<project>/research/08_theory/DEC1_decagon"
mons = pickle.load(open(f"{base}/mons.pkl","rb"))
ks = list(mons.keys())
print("mons keys sample:", ks[:5])
print("mons value sample:", mons[ks[0]])
print()
Sa = pickle.load(open(f"{base}/S_assembly.pkl","rb"))
ks2 = list(Sa.keys())
print("S_assembly keys:", ks2[:10])
for k in ks2[:3]:
    print(k, "->", Sa[k])
