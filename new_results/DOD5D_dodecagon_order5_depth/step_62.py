src = '''"""dod5d_grewrite.py -- exact Q(sqrt3) rewrite of the 12-gon's monomial integrals into G-functions.

Letters for G(a1..ak;1): 0, 1, p0=4+2r3, p1=-(3+2r3), p2=-(1+r3)/2, p3=(3+r3)/2  (p_i = -1/kappa_i).
I(n; a,b,c0,c1,c2,c3) = int_0^1 t^n (log t)^a (log(1-t))^b prod_i (log(1+kappa_i t))^c_i dt.
Exact rewrite via shuffle + integration by parts (identical recursions to OCT5D2's oct5d2_engine):
  prod of logs = a! b! c0!...c3! sum_{w in shuffle(0^a,1^b,p0^c0,..,p3^c3)} G(w; t)
  IG(n, w) = int_0^1 t^n G(w; t) dt
           = (1 - w1^(n+1))/(n+1) G(w;1) - (1/(n+1)) sum_{j=0}^n w1^(n-j) IG(j, w[1:]),
  IG(n, ()) = 1/(n+1).
All coefficients exact QS3.
"""
from fractions import Fraction as F
from functools import cache
from math import factorial
from qs3_symbol import q3, q3add, q3mul, q3div, ZERO, ONE

L0 = (F(0), F(0)); L1 = (F(1), F(0))
P0 = (F(4), F(2))            # 4+2r3 = V^2
P1 = (F(-3), F(-2))          # -(3+2r3) = -r3*U
P2 = (F(-1,2), F(-1,2))      # -(1+r3)/2
P3 = (F(3,2), F(1,2))        # (3+r3)/2
LETTERS = [L0, L1, P0, P1, P2, P3]

@cache
def shuf(u, v):
    if not u: return {v: ONE}
    if not v: return {u: ONE}
    out = {}
    for w, c in shuf(u[1:], v).items(): out[(u[0],)+w] = q3add(out.get((u[0],)+w, ZERO), c)
    for w, c in shuf(u, v[1:]).items(): out[(v[0],)+w] = q3add(out.get((v[0],)+w, ZERO), c)
    return out

def shuffle_words(parts):
    out = {(): ONE}
    for part in parts:
        nxt = {}
        for u, cu in out.items():
            for w, c in shuf(u, part).items():
                nxt[w] = q3add(nxt.get(w, ZERO), q3mul(cu, c))
        out = nxt
    return out

@cache
def IG(n, w):
    if not w:
        return {(): q3(F(1, n+1))}
    w1, wr = w[0], w[1:]
    out = {}
    c = q3div(q3add(ONE, tuple(-x for x in q3pow_(w1, n+1))), q3(F(n+1)))
    if c != ZERO:
        out[w] = c
    for j in range(n+1):
        cj = q3pow_(w1, n-j)
        if cj == ZERO: continue
        sub = IG(j, wr)
        for v, cv in sub.items():
            out[v] = q3add(out.get(v, ZERO), q3div(q3mul(tuple(-x for x in q3mul(cj, cv)), ONE), q3(F(n+1))))
    return {k: v for k, v in out.items() if v != ZERO}

@cache
def q3pow_(a, n):
    r = ONE
    for _ in range(n): r = q3mul(r, a)
    return r

def monomial_integral(n, a, b, c0, c1, c2, c3):
    """Exact rewrite of I(n; a,b,c0..c3) as dict {word: QS3 coeff}, |word| <= total logs."""
    cs = [c0, c1, c2, c3]
    parts = []
    if a: parts.append((L0,)*a)
    if b: parts.append((L1,)*b)
    for i in range(4):
        if cs[i]: parts.append((LETTERS[2+i],)*cs[i])
    sw = shuffle_words(parts)
    fac = q3(F(factorial(a)*factorial(b)*factorial(c0)*factorial(c1)*factorial(c2)*factorial(c3)))
    out = {}
    for w, mult in sw.items():
        ig = IG(n, w)
        for v, cv in ig.items():
            out[v] = q3add(out.get(v, ZERO), q3mul(q3mul(fac, mult), cv))
    return {k: v for k, v in out.items() if v != ZERO}
'''
with open("dod5d_grewrite.py","w") as f: f.write(src)
print("written dod5d_grewrite.py")
