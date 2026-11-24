"""Leading (k^{d-1}) coefficient of the degree-d kernels.  D_d = Delta^{d-1} c(d+1) = (d-1)! * [k^{d-1}] c(k).
Test: at the pole v=-a, D_d(key=(i,j)) / binom(d,j) = (-1)^{j+1} C_d with C_d independent of the key?  Exact (Q(sqrt2))."""
import json, os
from fractions import Fraction as Fr
from math import comb
HERE = os.path.dirname(os.path.abspath(__file__))
J = json.load(open(os.path.join(HERE, "genfun.json")))
for d in range(2, 9):
    Cs = {}
    for name, rec in J.items():
        dd, key = name.split("|"); 
        if int(dd) != d: continue
        i, j = map(int, key.split(","))
        p0, p1 = rec["-a"]["D"][d-1]
        val = (Fr(p0), Fr(p1))
        Cs[(i, j)] = ((-1)**(j+1) * val[0] / comb(d, j), (-1)**(j+1) * val[1] / comb(d, j))
    same = len(set(Cs.values())) == 1
    print(f"d={d}: C_d = D_d/((-1)^(j+1) binom(d,j)) at pole -a over keys: {sorted(set(Cs.values()))}  key-independent: {same}")

# ---- identification of C_d with the Taylor coefficients of sech(x)+tanh(x) = (1+sinh x)/cosh x  (exact, sympy) ----
import sympy as sp
X = sp.symbols('X')
ser = sp.series((1 + sp.sinh(X))/sp.cosh(X), X, 0, 10).removeO()
print("series of (1+sinh x)/cosh x:", sp.expand(ser))
allsame = True
for d in range(2, 9):
    Cs = set()
    for name, rec in J.items():
        dd, key = name.split("|")
        if int(dd) != d: continue
        i, j = map(int, key.split(","))
        p0, p1 = rec["-a"]["D"][d-1]
        Cs.add(((-1)**(j+1) * Fr(p0) / comb(d, j), (-1)**(j+1) * Fr(p1) / comb(d, j)))
    (c0, c1), = Cs
    target = ser.coeff(X, d)
    ok = (c1 == 0) and (sp.Rational(c0.numerator, c0.denominator) == target)
    allsame = allsame and ok
    print(f"d={d}: C_d = {c0}   [x^{d}] (sech x + tanh x) = {target}   equal: {ok}")
print("C_d = [x^d](sech x + tanh x) for d=2..8 (exact):", allsame)

# ---- all poles: D_d / binom(d,j) per key, to see which poles carry the leading coefficient ----
print("\nleading coefficient D_d(key,pole)/binom(d,j) (exact; pole '0','a','-a','-b'):")
for d in range(2, 7):
    for name, rec in J.items():
        dd, key = name.split("|")
        if int(dd) != d: continue
        i, j = map(int, key.split(","))
        row = {}
        for pole, v in rec.items():
            p0, p1 = v["D"][d-1]
            q0, q1 = Fr(p0)/comb(d, j), Fr(p1)/comb(d, j)
            row[pole] = f"{q0}" + (f" + ({q1})*sqrt2" if q1 != 0 else "")
        print(f"  d={d} key=({i},{j}): {row}")
