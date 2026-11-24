import pickle, sys
sys.path.insert(0, ".")
from fractions import Fraction as F
from math import gcd
import obstruction32 as OB
D = pickle.load(open("dod5d_symbols.pkl","rb"))
table, (n3, n2) = OB.obstruction(D['total'], 3)
def fact(v):
    p, q = v[0], v[1]
    if p == 0 and q == 0: return "0"
    den = p.denominator*q.denominator
    P, Q = p*den, q*den
    g = gcd(abs(P), abs(Q)); g = gcd(g, den)
    # extract integer gcd of P,Q only (keep denominator)
    gg = gcd(abs(P), abs(Q))
    if gg > 1 and den == 1:
        return f"{gg}*({P//gg}{Q//gg:+}r3)"
    if den == 1:
        return f"({P}{Q:+}r3)"
    return f"({P}{Q:+}r3)/{den}"
names3 = ['[X,[X,Y]]','[Y,[X,Y]]','[Z,[X,Y]]','[X,[X,Z]]','[Y,[X,Z]]','[Z,[X,Z]]','[Y,[Y,Z]]','[Z,[Y,Z]]']
names2 = ['[X,Y]','[X,Z]','[Y,Z]']
print("TOTAL obstruction table, X=log(2+r3), Y=log(1+r3), Z=log(r3):")
for i, row in enumerate(table):
    print(f"{names3[i]:12s} | " + " | ".join(f"{fact(v):>28s}" for v in row))
