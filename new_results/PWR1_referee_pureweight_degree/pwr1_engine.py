
# pwr1_engine.py -- PWR1 referee's independent exact engine for thm:pw.
# Formal series in tau with coefficients in Q[Z2,Z3,...,GAM] (Zk = zeta(k), GAM = Euler gamma).
# Weight: wt(Zk)=k, wt(GAM)=1. We verify a posteriori: coefficient of tau^k is wt-k homogeneous, GAM-free.
import sympy as sp
from functools import lru_cache

NMAX = 10
ZS = sp.symbols('Z2:'+str(NMAX+2))   # Z2..Z9
Z = {k: ZS[k-2] for k in range(2, NMAX+2)}
GAM = sp.Symbol('GAM')
TAU = sp.Symbol('TAU')

def zs(k):
    return sp.Integer(0) if k == 0 else Z[k]

class Ser:
    """dict {k: sympy expr}, series in tau up to NMAX."""
    def __init__(self, d=None):
        self.d = {k: sp.expand(v) for k, v in (d or {}).items() if v != 0}
    @staticmethod
    def zero(): return Ser({})
    @staticmethod
    def one(): return Ser({0: sp.Integer(1)})
    def copy(self): return Ser(dict(self.d))
    def __add__(a, b):
        d = dict(a.d)
        for k, v in b.d.items(): d[k] = sp.expand(d.get(k, 0) + v)
        return Ser(d)
    def __neg__(a): return Ser({k: -v for k, v in a.d.items()})
    def __sub__(a, b): return a + (-b)
    def scale(a, c):
        return Ser({k: sp.expand(c*v) for k, v in a.d.items()})
    def __mul__(a, b):
        d = {}
        for ka, va in a.d.items():
            for kb, vb in b.d.items():
                if ka+kb <= NMAX:
                    d[ka+kb] = sp.expand(d.get(ka+kb, 0) + va*vb)
        return Ser(d)
    def coeff(self, k): return sp.expand(self.d.get(k, sp.Integer(0)))

def ser_exp(a):
    # exp of series with zero constant
    out = Ser.one(); term = Ser.one()
    for m in range(1, NMAX+1):
        term = (term * a).scale(sp.Rational(1, m))
        out = out + term
    return out

def ser_log(a):
    # a has constant term 1
    out = Ser.zero(); term = Ser.one(); am1 = a - Ser.one()
    for m in range(1, NMAX+1):
        term = term * am1
        out = out + term.scale(sp.Rational(-(-1)**m, m))
    return out

def logGamma1(c):
    """log Gamma(1 + c*tau), c rational, as Ser: -GAM*c*tau + sum_{k>=2} (-1)^k Zk c^k tau^k/k."""
    d = {1: -GAM*sp.Rational(c)}
    for k in range(2, NMAX+1):
        d[k] = (-1)**k * Z[k] * sp.Rational(c)**k / k
    return Ser(d)

def partitions(s):
    # yield (tuple of parts, r, [multiplicities]) for partitions of s (nonincreasing)
    def rec(rem, mx, cur):
        if rem == 0:
            yield list(cur); return
        for x in range(min(mx, rem), 0, -1):
            cur.append(x); yield from rec(rem-x, x, cur); cur.pop()
    yield from rec(s, s, [])

def mult_fact(lam):
    from collections import Counter
    c = Counter(lam); out = sp.Integer(1)
    for v in c.values(): out *= sp.factorial(v)
    return out

def falling(n, r):
    out = sp.Integer(1)
    for i in range(r): out *= (n - i)
    return out

def moment_logM(lam, s, n):
    """log M(lambda,s): sum_i logG(c(lam_i)) + (n-r) logG(c(0)) - n logG(c_den),
       c(l) = 2l - 2s/n - (n-1)/n,  c_den = -(n-1)/n."""
    r = len(lam)
    tot = Ser.zero()
    for li in lam:
        tot = tot + logGamma1(sp.Rational(2*li) - sp.Rational(2*s, n) - sp.Rational(n-1, n))
    c0 = sp.Rational(-2*s, n) - sp.Rational(n-1, n)
    for _ in range(n - r):
        tot = tot + logGamma1(c0)
    cd = -sp.Rational(n-1, n)
    for _ in range(n):
        tot = tot - logGamma1(cd)
    return tot

def EX(j, n):
    """E[X^j] for B_n."""
    tot = Ser.zero()
    for s in range(0, j+1):
        inner = Ser.zero()
        if s == 0:
            inner = Ser.one()   # empty product -> M=1, count=1
        else:
            for lam in partitions(s):
                r = len(lam)
                cnt = sp.Rational(falling(n, r), 1) / mult_fact(lam)
                cnt *= sp.Rational(sp.factorial(s), 1)
                for li in lam: cnt /= sp.factorial(li)
                inner = inner + ser_exp(moment_logM(lam, s, n)).scale(cnt)
        tot = tot + inner.scale(sp.binomial(j, s) * (-1)**(j-s) * sp.Integer(n)**(j-s))
    return tot.scale(sp.Rational(1, n**j))

def logB(n):
    B = Ser.one()
    for j in range(1, NMAX+1):
        B = B + EX(j, n).scale(sp.binomial(sp.Rational(1, 2), j))
    return ser_log(B), B

def Lambda(n):
    return Ser({k: -Z[k] * sp.Rational((n-1)**k) * (1 - sp.Rational(1, n)**(k-1)) / k
                for k in range(2, NMAX+1)})

def logRat(n):
    out = Ser.zero()
    # (n-1) log(1-tau)
    for k in range(1, NMAX+1):
        out = out + Ser({k: -sp.Rational(n-1, k)})
    # - sum_{i=1}^{n-1} log(1 - (n-1) tau / i)
    for i in range(1, n):
        for k in range(1, NMAX+1):
            out = out + Ser({k: sp.Rational((n-1)**k, k * i**k)})
    return out

def total_log(n):
    lB, B = logB(n)
    return logRat(n) + Lambda(n) + lB, Lambda(n) + lB, B

def weight_check(expr, k):
    """True iff expr is homogeneous of weight k with wt(Zi)=i, wt(GAM)=1."""
    t = sp.Symbol('t')
    sub = {Z[i]: Z[i]*t**i for i in Z}; sub[GAM] = GAM*t
    e = sp.expand(expr.subs(sub))
    if e == 0: return True
    poly = sp.Poly(e, t)
    return all(monom[0] == k for monom in poly.monoms())
