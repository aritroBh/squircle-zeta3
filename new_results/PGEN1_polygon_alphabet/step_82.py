# step_54: task 3: L3(xi) per edge for OCT9c vs OCT9d bases
import pickle, os
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 100
base = "<project>"
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
O9D = os.path.join(base,"research/08_theory/OCT9d_octagon_short")
s123 = pickle.load(open(os.path.join(O9B,"sol_B123.pkl"),"rb"))
lists = pickle.load(open(os.path.join(O9B,"Sunit_lists.pkl"),"rb"))
s0 = pickle.load(open(os.path.join(O9B,"sol_B0.pkl"),"rb"))
print("s123[1] type:", type(s123[1]), "len:", len(s123[1]) if hasattr(s123[1],'__len__') else '')
print("first entries:", s123[1][:5] if hasattr(s123[1],'__getitem__') else s123[1])
Gargs = [Fr(a) for a in lists['Gargs']]
print("len Gargs:", len(Gargs), "Dargs:", len(lists['Dargs']))
print("s0 keys:", list(s0.keys()))
