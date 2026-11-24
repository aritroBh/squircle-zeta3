"""admissible.py -- fast exact enumeration of admissible arguments f (f and 1-f both of the form
+-2^k prod_i (c - r_i)^{e_i} for a given finite set of rational roots r_i; letters are linear)."""
from fractions import Fraction as Fr
import itertools

def pmul(P, Q):
    out = [Fr(0)]*(len(P) + len(Q) - 1)
    for i, a in enumerate(P):
        if a:
            for j, b in enumerate(Q):
                out[i + j] += a*b
    return out

def psub(P, Q):
    n = max(len(P), len(Q)); out = [Fr(0)]*n
    for i, a in enumerate(P): out[i] += a
    for i, b in enumerate(Q): out[i] -= b
    while len(out) > 1 and out[-1] == 0: out.pop()
    return out

def peval(P, x):
    s = Fr(0)
    for a in reversed(P): s = s*x + a
    return s

def pdiv_linear(P, r):          # P / (c - r), exact
    n = len(P) - 1; q = [Fr(0)]*n; acc = Fr(0)
    for i in range(n, 0, -1):
        acc = acc*r + P[i]; q[i - 1] = acc
    return q

def factor_over(P, roots):
    """returns (const, {root: multiplicity}) if P = const * prod (c - r)^m, else None"""
    mult = {}
    P = list(P)
    while len(P) > 1 and P[-1] == 0: P.pop()
    if len(P) == 1 and P[0] == 0:
        return None
    changed = True
    while len(P) > 1 and changed:
        changed = False
        for r in roots:
            if len(P) > 1 and peval(P, r) == 0:
                P = pdiv_linear(P, r); mult[r] = mult.get(r, 0) + 1; changed = True
    if len(P) > 1:
        return None
    return P[0], mult

def is_pm_pow2(q):
    q = abs(Fr(q))
    if q == 0:
        return False
    n, d = q.numerator, q.denominator
    return (n & (n - 1)) == 0 and (d & (d - 1)) == 0

def letter_poly(r):             # the letter (c - r) as polynomial
    return [Fr(-r), Fr(1)]

def enumerate_admissible(roots, E=2, K=3):
    """f = s 2^k prod (c - r)^{e_r}; returns list of (s, k, exps tuple) with 1 - f admissible."""
    out = []
    for ex in itertools.product(range(-E, E + 1), repeat=len(roots)):
        if all(e == 0 for e in ex):
            continue
        N = [Fr(1)]; D = [Fr(1)]
        for r, e in zip(roots, ex):
            for _ in range(abs(e)):
                if e > 0: N = pmul(N, letter_poly(r))
                else: D = pmul(D, letter_poly(r))
        for k in range(-K, K + 1):
            for s in (1, -1):
                Nf = [s*Fr(2)**k*a for a in N]
                res = factor_over(psub(D, Nf), roots)     # 1 - f = (D - N)/D
                if res is not None and is_pm_pow2(res[0]):
                    out.append((s, k, ex))
    return out
