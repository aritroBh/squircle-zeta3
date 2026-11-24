import sys
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from cor7_pt import QuadField, QuadPoint
import mpmath as mp

QF = QuadField(17)
P = QuadPoint(QF, "dbg", (Fr(7,16), Fr(1,16)), primes=[2])
F = P.F
one = F.one()
letters = {"2": F.from_int(2), "c": P.c0, "1-c": F.sub(one, P.c0), "1+c": F.add(one, P.c0)}
for nm, el in letters.items():
    v = P.vec(el)
    print(nm, el, "N =", F.norm(el), "vec =", v)
print("gens:", P.gensE)
