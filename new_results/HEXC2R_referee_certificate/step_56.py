
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../HEXC2_level6_certificate")
sys.path.insert(0, "../HEXC1_hexagon_certificate")
from fractions import Fraction as Fr
import reflib as R
import hexc2_common as HC
import json, re

G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
def my_fit_vec(name, weight):
    D = json.load(open(G1 + "/" + name))
    out = {}
    for key, val in D:
        m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
        if key == "1": e, w = 0, ()
        elif key == "c": e, w = 1, ()
        else:
            e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
        if len(w) != weight: continue
        v = {(): Fr(1)} if not w else R.my_hpl_to_cmzv(w)
        out = R.vadd((out, 1), (v, Fr(val) * Fr(1, 2) ** e))
    return out

def prod(*vecs):
    out = {(): Fr(1)}
    for v in vecs: out = R.mul(out, v)
    return out

myW2 = R.vadd((prod(R.my_atom_log2(), R.my_atom_log2()), Fr(1,2)), (R.my_li_atom(2, Fr(1,3)), 1),
              (prod(R.my_atom_log32(), R.my_atom_log32()), Fr(1,2)))
hcW2 = HC.W2_closed_vec()
print("W2_closed equal:", myW2 == hcW2)
if myW2 != hcW2:
    for w in set(myW2) | set(hcW2):
        if myW2.get(w, 0) != hcW2.get(w, 0):
            print("  ", w, "mine", myW2.get(w, 0), "theirs", hcW2.get(w, 0))

myR32 = my_fit_vec("hpl_fit_k3_(0,0,0).json", 2)
hcR32 = HC.fit_vec("hpl_fit_k3_(0,0,0).json", 2)
print("R32 equal:", myR32 == hcR32)
if myR32 != hcR32:
    for w in set(myR32) | set(hcR32):
        if myR32.get(w, 0) != hcR32.get(w, 0):
            print("  ", w, "mine", myR32.get(w, 0), "theirs", hcR32.get(w, 0))
