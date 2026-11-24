# 24_sunit_zeta8.py -- S-unit equation solutions z + w = 1 in Q(zeta8), S = {primes above 2, 3},
# generators: zeta8 (torsion), a = 1+sqrt2, pi2 = 1-zeta8, p1 = sqrt2+i, p2 = sqrt2-i.
import itertools, mpmath as mp, json
from fractions import Fraction as Fr
mp.mp.dps = 30
# elements as tuples (c0,c1,c2,c3) for c0 + c1 z + c2 z^2 + c3 z^3, z^4 = -1
def mul(x, y):
    r = [Fr(0)]*7
    for i in range(4):
        for j in range(4): r[i+j] += x[i]*y[j]
    return tuple(r[k] - (r[k+4] if k+4 < 7 else 0) for k in range(4))
def add(x, y): return tuple(a+b for a, b in zip(x, y))
def neg(x): return tuple(-a for a in x)
ONE = (Fr(1), Fr(0), Fr(0), Fr(0)); Z = (Fr(0), Fr(1), Fr(0), Fr(0))
def pw(x, k):
    r = ONE
    for _ in range(k): r = mul(r, x)
    return r
def emb(x, j):   # zeta -> exp(i pi j/4), j odd
    w = mp.expjpi(mp.mpf(j)/4)
    return sum(mp.mpf(c.numerator)/c.denominator*w**k for k, c in enumerate(x))
def norm(x):
    v = 1
    for j in (1, 3, 5, 7): v *= emb(x, j)
    return mp.re(v)
z = Z; z2 = pw(Z, 2); z3 = pw(Z, 3)
r2 = add(z, neg(z3))                 # sqrt2 = z - z^3
i_ = z2
a = add(ONE, r2)
pi2 = add(ONE, neg(z))
p1 = add(r2, i_); p2 = add(r2, neg(i_))
def inv_num(x):  # numeric inverse via embeddings not needed; we use negative exponents via a^-1 = sqrt2-1 etc.
    pass
ainv = add(r2, neg(ONE))
# inverses of pi2, p1, p2 up to units: pi2^-1 = (1 - z^-1)... use exact: 1/(1-z) = conj stuff / norm; do it numerically-exact by solving
def inverse(x):
    # x * y = 1: solve linear system over Q
    import sympy
    c = sympy.symbols('c0:4')
    y = tuple(c)
    prod = [0]*7
    for i_ in range(4):
        for j in range(4): prod[i_+j] += x[i_]*y[j]
    eqs = [prod[k] - (prod[k+4] if k+4 < 7 else 0) - (1 if k == 0 else 0) for k in range(4)]
    sol = sympy.solve(eqs, c)
    return tuple(Fr(str(sol[ci])) for ci in c)
gens = {"a": (a, ainv), "pi2": (pi2, inverse(pi2)), "p1": (p1, inverse(p1)), "p2": (p2, inverse(p2))}
def power(g, e):
    x, xi = gens[g]
    return pw(x, e) if e >= 0 else pw(xi, -e)
sols = {}
R = {"a": range(-3, 4), "pi2": range(-4, 5), "p1": range(-2, 3), "p2": range(-2, 3)}
for ea in R["a"]:
    for e2 in R["pi2"]:
        for e3 in R["p1"]:
            for e4 in R["p2"]:
                base = mul(mul(power("a", ea), power("pi2", e2)), mul(power("p1", e3), power("p2", e4)))
                for k in range(8):
                    x = mul(base, pw(Z, k))
                    w = add(ONE, neg(x))
                    if all(c == 0 for c in w): continue
                    N = abs(norm(w))
                    if N < mp.mpf(10)**-20: continue
                    lg2 = mp.log(N, 2)
                    ok = False
                    for j in range(-12, 13):
                        t = lg2 - j*mp.log(3, 2)
                        if abs(t - mp.nint(t)) < mp.mpf(10)**-15: ok = True; break
                    if not ok: continue
                    zv = emb(x, 1)
                    key = (mp.nstr(mp.re(zv), 15), mp.nstr(mp.im(zv), 15))
                    sols[key] = (zv, (k, ea, e2, e3, e4))
print(len(sols), "solutions")
out = [{"re": mp.nstr(mp.re(v[0]), 25), "im": mp.nstr(mp.im(v[0]), 25), "exps": v[1]} for v in sols.values()]
json.dump(out, open("sunit_zeta8.json", "w"), indent=0)
