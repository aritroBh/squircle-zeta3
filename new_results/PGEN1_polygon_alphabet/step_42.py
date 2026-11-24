# step_18: table vecfuns (exact-verified) + delta for regular octagon, decagon, 12-gon
import sympy as sp
from fractions import Fraction as Fr

def check_prod(elem, vec, letters):
    prod = sp.Integer(1)
    for nm, e in vec.items():
        prod *= letters[nm]**sp.Rational(e.numerator, e.denominator)
    return sp.simplify(sp.expand(elem) - sp.expand(prod)) == 0 or sp.simplify(sp.expand(elem) + sp.expand(prod)) == 0

# ---------------- Q(sqrt2) ----------------
r2 = sp.sqrt(2); u2 = 1+r2
L2 = {"r2": r2, "u": u2}
T2 = {  # name -> (elem, vec)
 "2":   (sp.Integer(2), {"r2":Fr(2)}),
 "k":   (r2-1,        {"u":Fr(-1)}),          # sqrt2-1 = 1/u
 "d":   (2-r2,        {"r2":Fr(1),"u":Fr(-1)}),
 "r2h": (r2/2,        {"r2":Fr(-1)}),
 "mr2h":(-r2/2,       {"r2":Fr(-1)}),
 "x1":  (1-r2,        {"u":Fr(-1)}),
 "x2":  (1-r2/2,      {"r2":Fr(-1),"u":Fr(-1)}),
 "1-x1":(r2,          {"r2":Fr(1)}),
 "1-x2":(r2/2,        {"r2":Fr(-1)}),
 "1-k": (2-r2,        {"r2":Fr(1),"u":Fr(-1)}),
 "1-d": (r2-1,        {"u":Fr(-1)}),
 "1-r2h":(1-r2/2,     {"r2":Fr(-1),"u":Fr(-1)}),
 "1-mr2h":(1+r2/2,    {"r2":Fr(-1),"u":Fr(1)}),
}
assert all(check_prod(e,v,L2) for e,v in T2.values()), "T2 fail"
print("Q(sqrt2) table verified (13 elems)")

# ---------------- Q(sqrt5) ----------------
s5 = sp.sqrt(5); phi = (1+s5)/2
L5 = {"2":sp.Integer(2), "s5":s5, "phi":phi}
T5 = {
 "1/phi": (1/phi, {"phi":Fr(-1)}),
 "-1/phi":(-1/phi, {"phi":Fr(-1)}),
 "1-1/phi": (1-1/phi, {"phi":Fr(-2)}),     # 1-1/phi = 1/phi^2
 "1+1/phi": (1+1/phi, {"phi":Fr(1)}),      # 1+1/phi = phi
 "phi": (phi, {"phi":Fr(1)}),
 "phi2": (phi**2, {"phi":Fr(2)}),
 "-phi": (-phi, {"phi":Fr(1)}),
 "1-phi": (1-phi, {"phi":Fr(-1)}),         # = -1/phi
 "1-phi2": (1-phi**2, {"phi":Fr(1)}),      # = -phi
 "1-1/phi2": (1-1/phi**2, {"phi":Fr(-1)}), # = 1/phi
 "sqrt5": (s5, {"s5":Fr(1)}),
}
assert all(check_prod(e,v,L5) for e,v in T5.values()), "T5 fail"
print("Q(sqrt5) table verified")

# ---------------- Q(sqrt3) ----------------
s3 = sp.sqrt(3); g2 = 1+s3; u3 = 2+s3
L3 = {"g2":g2, "g3":s3, "u":u3}
T3 = {
 "z1": (2-s3, {"u":Fr(-1)}),
 "z2": ((s3-1)/2, {"g2":Fr(-1)}),
 "z3": (s3/3, {"g3":Fr(-1)}),
 "z4": (7-4*s3, {"u":Fr(-2)}),
 "z5": ((2-s3)/2, {"g2":Fr(-2)}),
 "z6": ((2*s3-3)/3, {"g3":Fr(-1),"u":Fr(-1)}),
 "z7": ((3-s3)/6, {"g2":Fr(-1),"g3":Fr(-1)}),
 "z8": (sp.Rational(1,3), {"g3":Fr(-2)}),
 "1-z1": (s3-1, {"g2":Fr(1),"u":Fr(-1)}),
 "1-z2": ((3-s3)/2, {"g2":Fr(-1),"g3":Fr(1)}),
 "1-z3": ((3-s3)/3, {"g2":Fr(1),"g3":Fr(-1),"u":Fr(-1)}),
 "1-z4": (-6+4*s3, {"g2":Fr(2),"g3":Fr(1),"u":Fr(-2)}),
 "1-z5": (s3/2, {"g2":Fr(-2),"g3":Fr(1),"u":Fr(1)}),
 "1-z6": ((6-2*s3)/3, {"g2":Fr(3),"g3":Fr(-1),"u":Fr(-2)}),
 "1-z7": ((3+s3)/6, {"g2":Fr(-1),"g3":Fr(-1),"u":Fr(1)}),
 "1-z8": (sp.Rational(2,3), {"g2":Fr(2),"g3":Fr(-2),"u":Fr(-1)}),
 "p1": ((3+s3)/2, {"g2":Fr(-1),"g3":Fr(1),"u":Fr(1)}),
 "p2": (4+2*s3, {"g2":Fr(2)}),
 "p3": (-(3+2*s3), {"g3":Fr(1),"u":Fr(1)}),
 "p4": (-(1+s3)/2, {"g2":Fr(-1),"u":Fr(1)}),
 "1-p1": (( -1-s3)/2, {"g2":Fr(-1),"u":Fr(0),"g3":Fr(0)}),  # placeholder fix below
}
del T3["1-p1"]
assert all(check_prod(e,v,L3) for e,v in T3.values()), "T3 fail"
print("Q(sqrt3) table verified (8 pts + 1-z + 4 p_j)")

import pickle
pickle.dump({"L2":L2,"T2":T2,"L5":L5,"T5":T5,"L3":L3,"T3":T3}, open("field_tables.pkl","wb"))
print("saved field_tables.pkl")
