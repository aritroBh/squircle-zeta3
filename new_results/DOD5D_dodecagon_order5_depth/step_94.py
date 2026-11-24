import pickle, sys
sys.path.insert(0, ".")
from fractions import Fraction as F
from math import gcd
import obstruction32 as OB
D = pickle.load(open("dod5d_symbols.pkl","rb"))
table, (n3, n2) = OB.obstruction(D['total'], 3)
def fact(v):
    p, q = F(v[0]), F(v[1])
    if p == 0 and q == 0: return "0"
    den = p.denominator*q.denominator
    P, Q = (p*den).numerator, (q*den).numerator
    gg = gcd(abs(P), abs(Q))
    core = f"{P//gg}{Q//gg:+}r3" if gg else f"{P}{Q:+}r3"
    pre = f"{gg}" if gg > 1 else ""
    if den != 1:
        return f"{pre+'*' if pre else ''}({core})/{den}"
    return f"{pre}*({core})" if pre else f"({core})"
names3 = ['[X,[X,Y]]','[Y,[X,Y]]','[Z,[X,Y]]','[X,[X,Z]]','[Y,[X,Z]]','[Z,[X,Z]]','[Y,[Y,Z]]','[Z,[Y,Z]]']
names2 = ['[X,Y]','[X,Z]','[Y,Z]']
for i, row in enumerate(table):
    print(f"{names3[i]:12s} | " + " | ".join(f"{fact(v):>26s}" for v in row))
