# step_55: inspect OCT9c solution structures
import pickle, os
from fractions import Fraction as Fr
base = "<project>"
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
s123 = pickle.load(open(os.path.join(O9B,"sol_B123.pkl"),"rb"))
lists = pickle.load(open(os.path.join(O9B,"Sunit_lists.pkl"),"rb"))
s0 = pickle.load(open(os.path.join(O9B,"sol_B0.pkl"),"rb"))
for e in (1,2,3):
    d = s123[e]
    print(f"s123[{e}] keys:", list(d.keys()))
    for k,v in d.items():
        print("   ", k, type(v), repr(v)[:200])
        break
print("Gargs[0:5]:", lists['Gargs'][:5], " Dargs[0:5]:", lists['Dargs'][:5])
print("s0 keys:", list(s0.keys()))
for k in s0:
    print("  s0[",k,"]:", type(s0[k]), repr(s0[k])[:250])
