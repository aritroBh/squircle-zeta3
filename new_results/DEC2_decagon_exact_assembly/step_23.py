base = "<project>/research/08_theory/DEC1_decagon"
with open(f"{base}/P3_660.txt") as f:
    print(f.read())
mons = pickle.load(open(f"{base}/mons.pkl","rb"))
ks = list(mons.keys())
print("\nmons keys sample:", ks[:5])
print("mons value sample:", mons[ks[0]])
Sa = pickle.load(open(f"{base}/S_assembly.pkl","rb"))
ks2 = list(Sa.keys())
print("\nS_assembly keys:", ks2[:10])
print("sample value:", Sa[ks2[0]])
