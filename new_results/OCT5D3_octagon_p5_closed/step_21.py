import pickle, sys
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
base = "<project>/research/08_theory/OCT5D2_octagon_depth_test/"
for f in ["w5_430.pkl","w5_250.pkl"]:
    d = pickle.load(open(base+f,"rb"))
    print(f, type(d), list(d.keys()) if isinstance(d,dict) else "")
    for k,v in (d.items() if isinstance(d,dict) else []):
        print("  ",k, str(v)[:120])
# gtab cache info
g = pickle.load(open(base+"gtab_d430_N1440.pkl","rb"))
print(type(g), len(g))
raw, letters = g
print("letters:", letters)
print("words:", len(raw), "sample:", list(raw.items())[:2])
