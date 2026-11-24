"""t1_direct.py -- referee's INDEPENDENT finite-but-large check of the sufficiency statement.

Builds the tower-law candidate U^{(k)}_{pq}(delta) = sum_s binom(k-d,s) e(delta,s) (off-diagonal),
diagonal via e-hat with the exceptional level k=2p = rho_pp, using ONLY:
  * closed forms of e(delta,s) from tow2lib (arho_simple for s=0; fitted/zigzag forms for s>=1),
  * e00 closed form  K_d (-1)^p (1+[p>q]) binom(d,p),
  * the T1 recursion exactly as stated in TOW1 (true integer binomials in the fold, not falling factorials),
and tests  U^{(k)}_{pq}(delta) == T1-RHS( U^{(k)} of parents )  at every (key, k, delta<=min(4,2d-4)),
for d up to DMAX, k = d+1 .. d+KEXTRA (diagonal: k >= 2p+1, plus child at k=2p).
All arithmetic exact in Q(i) (Fractions); the common factor 1/(2 sqrt2) of K_d is dropped (T1 is linear).
"""
import sys, time
from fractions import Fraction as Fr
from math import comb
import sympy as sp
sys.path.insert(0, '.')
import tow2lib as T

class G:
    __slots__ = ('a', 'b')
    def __init__(s, a=0, b=0): s.a = Fr(a); s.b = Fr(b)
    def __add__(s, o):
        o = o if isinstance(o, G) else G(o); return G(s.a+o.a, s.b+o.b)
    __radd__ = __add__
    def __neg__(s): return G(-s.a, -s.b)
    def __sub__(s, o):
        o = o if isinstance(o, G) else G(o); return G(s.a-o.a, s.b-o.b)
    def __rsub__(s, o): return G(o) - s
    def __mul__(s, o):
        o = o if isinstance(o, G) else G(o); return G(s.a*o.a - s.b*o.b, s.a*o.b + s.b*o.a)
    __rmul__ = __mul__
    def conj(s): return G(s.a, -s.b)
    def inv(s):
        n = s.a*s.a + s.b*s.b; return G(s.a/n, -s.b/n)
    def __truediv__(s, o):
        o = o if isinstance(o, G) else G(o); return s * o.inv()
    def __eq__(s, o):
        o = o if isinstance(o, G) else G(o); return s.a == o.a and s.b == o.b
    def iszero(s): return s.a == 0 and s.b == 0
    def __repr__(s): return f'({s.a}+{s.b}i)'

I_ = G(0, 1)
def ipow(n):
    n %= 4
    return [G(1), G(0, 1), G(-1), G(0, -1)][n]

def sp2G(v):
    v = sp.expand_complex(sp.together(v))
    re, im = v.as_real_imag()
    re = sp.Rational(re); im = sp.Rational(im)
    return G(Fr(int(re.p), int(re.q)), Fr(int(im.p), int(im.q)))

# ---- tau table (independent check against definition done in tau_check.py) ----
def tau1(n): return sp2G(T.tau1f(n))
def tau2(n): return sp2G(T.tau2f(n))
_t1 = {n: tau1(n) for n in range(0, 8)}
_t2 = {n: tau2(n) for n in range(0, 8)}

_arho_cache = {}
def arho(delta, dd, mm):
    key = (delta, dd, mm)
    if key not in _arho_cache:
        if delta == 0: _arho_cache[key] = G(1)
        else:
            e = T.arho_simple(delta)
            _arho_cache[key] = sp2G(e.subs({T.d: dd, T.m: mm}))
    return _arho_cache[key]

def dfact(n):  # (2d-3)!! helper: odd double factorial
    r = 1
    while n > 1:
        r *= n; n -= 2
    return r

def kappa(dd):
    # K_d = (-1)^d (2d-3)!! i^{d+1} / (2 sqrt2 d!)  ; drop 1/(2 sqrt2)
    import math
    return ipow(dd+1) * Fr((-1)**dd * dfact(2*dd-3), math.factorial(dd))

def e00(p, q):
    dd = p+q
    mult = 2 if p > q else 1
    return kappa(dd) * ((-1)**p * mult * comb(dd, p))

def F(delta, s, dd, mm):
    """e(delta,s)/e00 closed form for an off-diagonal key (d=dd, m=mm=first index)."""
    if s > delta//2: return G(0)
    if s == 0: return arho(delta, dd, mm)
    if (delta, s) == (2, 1): return G(Fr(-1, 2*dd-3))
    if (delta, s) == (3, 1): return I_ * Fr(dd-3, 2*dd-3)
    if (delta, s) == (4, 1):
        num = G(dd**3 - 9*dd**2) + G(Fr(113, 4)*dd) - I_*Fr(dd, 2) + I_*mm - 29
        return num * Fr(1, (2*dd-3)*(2*dd-5))
    if (delta, s) == (4, 2): return G(Fr(2, (2*dd-3)*(2*dd-5)))
    raise ValueError((delta, s))

def e_abs(p, q, delta, s):
    assert p > q
    return F(delta, s, p+q, p) * e00(p, q)

def ehat(p, delta, s):
    tot = G(0)
    for dp in range(delta+1):
        tot = tot + _t1[delta-dp] * e_abs(p+1, p, dp, s) * Fr(p+1, 4*p-dp-1)
    return tot

def rho_pp(p, delta):
    return arho(delta, 2*p, p) * e00(p, p)

def U(p, q, k, delta):
    dd = p+q
    if p == q:
        if k == 2*p: return rho_pp(p, delta)
        assert k > 2*p
        tot = G(0)
        for s in range(delta//2+1):
            tot = tot + ehat(p, delta, s) * comb(k-2*p-1, s)
        return tot
    assert k >= dd
    tot = G(0)
    for s in range(delta//2+1):
        tot = tot + e_abs(p, q, delta, s) * comb(k-dd, s)
    return tot

def T1rhs(p, q, k, delta):
    dd = p+q
    tot = G(0)
    for dp in range(delta+1):
        tot = tot + _t1[delta-dp] * U(p+1, q, k, dp) * Fr(p+1, 2*dd-dp-1)
        if p > q:
            tot = tot + _t2[delta-dp] * U(p, q+1, k, dp) * Fr(q+1, 2*dd-dp-1)
    if q == p-1:
        m = 4*p-4-delta
        ft = G(0)
        for mu in range(m, 4*p-3):
            g = G(0)
            for dp in range(0, 4*p-4-mu+1):
                g = g - U(p, p, k, dp) * _t1[4*p-4-mu-dp] * Fr(1, 4*p-3-dp)
            ft = ft + g.conj() * ipow(-mu) * comb(mu-2, mu-m) * ipow(m-2)
        tot = tot + ft * (-p)
    return tot

if __name__ == '__main__':
    DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    KEXTRA = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    bad = []; n = 0; t0 = time.time()
    per_depth = {}
    for dd in range(4, DMAX+1):
        for p in range(dd, (dd-1)//2, -1):
            q = dd-p
            ks = range(dd+1, dd+KEXTRA+1)
            if p == q: ks = range(2*p+1, 2*p+KEXTRA+1)   # diagonal: only k>=2p+1 (k=2p is initial data)
            for k in ks:
                for delta in range(0, min(4, 2*dd-4)+1):
                    n += 1
                    lhs = U(p, q, k, delta); rhs = T1rhs(p, q, k, delta)
                    per_depth[delta] = per_depth.get(delta, 0) + 1
                    if not (lhs == rhs):
                        bad.append((p, q, k, delta))
    print(f'T1 residual check: d=4..{DMAX}, k=d+1..d+{KEXTRA}: {n} cells, {len(bad)} bad ({round(time.time()-t0,1)}s)')
    print('cells per depth:', per_depth)
    print('bad sample:', bad[:10])
