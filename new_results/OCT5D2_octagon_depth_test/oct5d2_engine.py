"""oct5d2_engine.py — exact Q(sqrt2) machinery for OCT5D2.

Letters for G(a1..ak;1): 0, 1, p = 2+sqrt2, q = -(1+sqrt2).
I(n; a,b,c,d) = int_0^1 t^n (log t)^a (log(1-t))^b (log(1+k0 t))^c (log(1+k1 t))^d dt
with k0 = -(2-sqrt2)/2, k1 = sqrt2-1, so log(1+k0 t) = G(p; t), log(1+k1 t) = G(q; t).

Exact rewrite via shuffle + integration by parts:
  (log t)^a (log(1-t))^b (log(1+k0 t))^c (log(1+k1 t))^d
      = a! b! c! d! * sum_{w in 0^a ш 1^b ш p^c ш q^d} G(w; t)
  IG(n, w) := int_0^1 t^n G(w; t) dt
      = (1 - w1^(n+1))/(n+1) * G(w;1) - (1/(n+1)) * sum_{j=0}^n w1^(n-j) IG(j, w[1:]),
  IG(n, ()) = 1/(n+1).
All coefficients exact in Q(sqrt2) (QS2 from the paper's validated qsqrt2.py).
"""
import sys
from fractions import Fraction as F
from functools import cache

sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from qsqrt2 import QS2 as Q, SQ2 as R, A as AA, B as BB

ZERO, ONE = Q(0), Q(1)
P_LET = Q(2, 1)      # p = 2+sqrt2  = -1/k0
Q_LET = Q(-1, -1)    # q = -(1+sqrt2) = -1/k1

def shuffles_into(out, w, coef=1):
    """add coef * w into dict out (words -> QS2)."""
    out[w] = out.get(w, ZERO) + coef

@cache
def shuf(u, v):
    """shuffle product of two words -> dict {word: multiplicity}."""
    if not u: return {v: ONE}
    if not v: return {u: ONE}
    out = {}
    for w, c in shuf(u[1:], v).items(): out[(u[0],)+w] = out.get((u[0],)+w, ZERO) + c
    for w, c in shuf(u, v[1:]).items(): out[(v[0],)+w] = out.get((v[0],)+w, ZERO) + c
    return out

def shuffle_words(parts):
    """shuffle of a list of words -> dict {word: multiplicity}."""
    out = {(): ONE}
    for part in parts:
        nxt = {}
        for u, cu in out.items():
            for w, c in shuf(u, part).items():
                nxt[w] = nxt.get(w, ZERO) + cu*c
        out = nxt
    return out

@cache
def IG(n, w):
    """int_0^1 t^n G(w; t) dt as dict {v: QS2 coeff}, |v| <= |w|. Exact."""
    if not w:
        return {(): ONE/Q(n+1)}
    w1, wr = w[0], w[1:]
    out = {}
    c = (ONE - w1**(n+1)) / Q(n+1)
    if c != ZERO:
        out[w] = c
    pw = ONE  # w1^(n-j) for j=0..n; careful: w1^0 = 1 even for w1=0
    coeffs = []
    for j in range(n+1):
        coeffs.append(w1**(n-j))   # QS2 pow; 0**0 -> QS2(1)? check: QS2.__pow__(0) returns QS2(1) -- yes
    for j in range(n+1):
        cj = coeffs[j]
        if cj == ZERO: continue
        sub = IG(j, wr)
        for v, cv in sub.items():
            out[v] = out.get(v, ZERO) - cj*cv/Q(n+1)
    return {k: v for k, v in out.items() if v != ZERO}

def monomial_integral(n, a, b, c, d):
    """Exact rewrite of I(n; a,b,c,d) as dict {word: QS2 coeff}, |word| <= a+b+c+d."""
    import math
    parts = []
    if a: parts.append((ZERO,)*a)
    if b: parts.append((ONE,)*b)
    if c: parts.append((P_LET,)*c)
    if d: parts.append((Q_LET,)*d)
    sw = shuffle_words(parts)
    fac = Q(math.factorial(a)*math.factorial(b)*math.factorial(c)*math.factorial(d))
    out = {}
    for w, mult in sw.items():
        ig = IG(n, w)
        for v, cv in ig.items():
            out[v] = out.get(v, ZERO) + fac*mult*cv
    return {k: v for k, v in out.items() if v != ZERO}
