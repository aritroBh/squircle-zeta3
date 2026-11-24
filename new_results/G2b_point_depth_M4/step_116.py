
import sys
sys.path.insert(0, "<project>/research/08_theory/G2b_point_depth_M4")
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 80
from g2b_dbg import s05
C0 = Fr(1,2)

# enumerate needed J3 calls for R4 and the control (control passed, so skip its triples)
needed = set()
for w, coef in s05.signed_words():
    pts = (Fr(0),) + tuple(w) + (C0,)
    for k in (2, 3):
        for i in range(6 - k):
            gap = w[i:i + k]; ret = w[:i] + w[i + k:]
            ei, ci = pts[i], pts[i + k + 1]
            if k == 2:
                lam = s05.lam2(gap, ci, ei)
                if lam: needed.add((ret, Fr(0), C0))
            else:
                lam = s05.lam2(ret, C0, Fr(0))
                if lam: needed.add((gap, ei, ci))
print("unique J3 calls:", len(needed))

# test J2val for nan on each call's relevant (b,c,p) and probe J3 cheaply with a coarse level
def reli2(z):
    z = mp.mpf(z)
    if z == 0: return mp.mpf(0)
    if z == 1: return mp.pi**2 / 6
    if abs(z) <= mp.mpf('0.55'):
        tot = mp.mpf(0); zk = z; k = 1
        while True:
            t = zk / (k * k); tot += t
            if abs(t) < mp.mpf(10) ** (-mp.mp.dps + 10): break
            zk *= z; k += 1
        return tot
    if z > 1: return -reli2(1 / z) + mp.pi**2 / 3 - mp.log(z)**2 / 2
    if z > 0: return mp.pi**2 / 6 - mp.log(z) * mp.log(1 - z) - reli2(1 - z)
    w = z / (z - 1)
    return -mp.log(1 - z)**2 / 2 - reli2(w)

def Gfun(u, b, c):
    if u == c: return mp.mpf(0)
    if u == b: return mp.zeta(2) - mp.log(abs(b - c))**2
    return mp.log(abs(u - c)) * mp.log(abs((b - u) / (b - c))) + reli2((u - c) / (b - c))

def J2val(b, c, p, t):
    lam = mp.log(abs(p - c)) if p != c else mp.mpf(0)
    if b == c:
        F = lambda u: mp.mpf(1)/2 * mp.log(abs(u - b))**2 - lam * mp.log(abs(u - b))
        return F(t) - (F(p) if p != b else mp.mpf(0))
    def F(u):
        if u == b: return mp.zeta(2) - mp.log(abs(b - c))**2
        return Gfun(u, b, c) - lam * mp.log(abs(u - b))
    return F(t) - F(p)

# scan J2val for nan over a grid of t for each (gap/endpoints)
import itertools
bad = []
for (w3, ei, ci) in sorted(needed, key=str):
    a, b, c = w3
    eif = mp.mpf(ei.numerator)/ei.denominator; cif = mp.mpf(ci.numerator)/ci.denominator
    for tt in [eif + (cif - eif) * mp.mpf(j) / 7 for j in range(8)]:
        try:
            v = J2val(b, c, eif, tt)
            if mp.isnan(v) or mp.isinf(v):
                bad.append((w3, ei, ci, tt, 'J2 nan')); break
        except Exception as ex:
            bad.append((w3, ei, ci, tt, repr(ex))); break
print("bad J2 cases:", len(bad))
for b_ in bad[:10]: print("  ", b_)
