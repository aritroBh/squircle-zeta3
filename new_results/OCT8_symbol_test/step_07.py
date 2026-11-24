import pickle, os
os.chdir("<project>/research/08_theory/OCT7_generic_octagon")
for f in ["verdicts700.pkl", "lam_lll_450.pkl", "geometry.pkl"]:
    try:
        d = pickle.load(open(f,"rb"))
        print("="*20, f, type(d))
        s = repr(d)
        print(s[:3000])
    except Exception as e:
        print(f, "ERR", e)
