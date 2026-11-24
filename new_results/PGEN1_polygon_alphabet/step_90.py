# step_62: find W3_e0 expression/value; check s123 W3 entries
import os, sys, pickle
base = "<project>"
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
s123 = pickle.load(open(os.path.join(O9B,"sol_B123.pkl"),"rb"))
for e in (1,):
    print("s123[1]['W3']:", type(s123[1]['W3']), repr(s123[1]['W3'])[:300])
    print("s123[1]['Wlo']:", type(s123[1]['Wlo']), repr(s123[1]['Wlo'])[:200])
for f in ["tmp_e0.pkl","Earc_edge0_gen.pkl","E3_edge0_mine.pkl"]:
    try:
        d = pickle.load(open(os.path.join(O9B,f),"rb"))
        print(f, type(d), repr(d)[:400])
    except Exception as ex:
        print(f, "ERR", ex)
