import sympy as sp
from fractions import Fraction

ORD = 14
# dict-based rational power series: {degree: Fraction}
def ser_mul(a, b, n=ORD):
    out = {}
    for i, ai in a.items():
        for j, bj in b.items():
            if i+j <= n: out[i+j] = out.get(i+j, Fraction(0)) + ai*bj
    return {k: v for k, v in out.items() if v != 0}
def ser_pow(a, p, n=ORD):
    out = {0: Fraction(1)}
    for _ in range(p):
        out = ser_mul(out, a, n)
    return out

# sech x = sum E_{2n} x^{2n}/(2n)!
sech = {2*k: Fraction(sp.euler(2*k)) / sp.factorial(2*k) for k in range(ORD//2+1)}
# gd' = sech  =>  gd = integral
gd = {2*k+1: sech[2*k]/(2*k+1) for k in range(ORD//2+1) if 2*k+1 <= ORD}
# sech^{1/2} = (1+u)^{1/2}, u = sech - 1
u = {k: v for k, v in sech.items() if k > 0}
sech12 = {0: Fraction(1)}
term = {0: Fraction(1)}
for j in range(1, ORD+1):
    # binom(1/2, j) u^j
    c = Fraction(1)
    for i in range(j): c *= Fraction(1, 2) - i
    c /= sp.factorial(j)
    uj = ser_pow(u, j)
    for k, v in uj.items():
        sech12[k] = sech12.get(k, Fraction(0)) + c*v
sech12 = {k: v for k, v in sech12.items() if v != 0}
sech32 = ser_mul(sech, sech12)
expx2 = {k: Fraction(1, 2**k * sp.factorial(k)) for k in range(ORD+1)}
H32 = ser_mul(expx2, sech32)   # 2^{3/2} H(x)

MAXR, MAXM = 4, 12
e = {}
for r in range(MAXR+1):
    g2r = ser_pow(gd, 2*r)
    s = ser_mul(g2r, H32)
    for m in range(MAXM+1):
        e[(r,m)] = s.get(m, Fraction(0))

import pickle
theirs = {}
with open('../K3_fourier_all_orders/e_rm.pkl','rb') as f:
    for key, val in pickle.load(f).items():
        r, m = map(int, key.split(',')); theirs[(r,m)] = sp.Rational(val)
mism = [(k, e[k], theirs[k]) for k in theirs if sp.Rational(e[k]) != theirs[k]]
print("e_rm: my exact rebuild vs K3 table:", "MISMATCHES "+str(mism) if mism else "EXACT MATCH (all 50 entries)")
print("e_{r,2r}=1 all r:", all(e[(r,2*r)]==1 for r in range(MAXR+1)))
print("e_{r,2r+1}=1/2 all r:", all(e[(r,2*r+1)]==Fraction(1,2) for r in range(MAXR+1) if 2*r+1<=MAXM))
print("e_{r,2r+2}=-5/8-r/3:", all(sp.Rational(e[(r,2*r+2)])==-sp.Rational(5,8)-sp.Rational(r,3) for r in range(MAXR+1) if 2*r+2<=MAXM))
print("sample row r=1:", [str(e[(1,m)]) for m in range(2, 13)])
pickle.dump({f"{r},{m}": str(e[(r,m)]) for r in range(MAXR+1) for m in range(MAXM+1)}, open('e_rm_mine.pkl','wb'))
