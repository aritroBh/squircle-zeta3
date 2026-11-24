import sympy as sp
from fractions import Fraction
import pickle
ORD = 14
def ser_mul(a, b, n=ORD):
    out = {}
    for i, ai in a.items():
        for j, bj in b.items():
            if i+j <= n: out[i+j] = out.get(i+j, Fraction(0)) + ai*bj
    return {k: v for k, v in out.items() if v != 0}
def ser_pow(a, p, n=ORD):
    out = {0: Fraction(1)}
    for _ in range(p): out = ser_mul(out, a, n)
    return out
sech = {2*k: Fraction(sp.euler(2*k)) / sp.factorial(2*k) for k in range(ORD//2+1)}
gd = {2*k+1: sech[2*k]/(2*k+1) for k in range(ORD//2+1) if 2*k+1 <= ORD}
u = {k: v for k, v in sech.items() if k > 0}
sech12 = {0: Fraction(1)}
for j in range(1, ORD+1):
    c = Fraction(1)
    for i in range(j): c *= Fraction(1, 2) - i
    c /= sp.factorial(j)
    for k, v in ser_pow(u, j).items(): sech12[k] = sech12.get(k, Fraction(0)) + c*v
sech12 = {k: v for k, v in sech12.items() if v != 0}
sech32 = ser_mul(sech, sech12)
H32 = ser_mul({k: Fraction(1, 2**k*sp.factorial(k)) for k in range(ORD+1)}, sech32)
e = {}
for r in range(6):
    s = ser_mul(ser_pow(gd, 2*r), H32)
    for m in range(13): e[(r,m)] = s.get(m, Fraction(0))
pickle.dump({f"{r},{m}": str(e[(r,m)]) for (r,m) in e}, open('e_rm_mine.pkl','wb'))
print("e table extended to r<=5, m<=12; e[5,10] =", e[(5,10)], " e[5,11] =", e[(5,11)])
