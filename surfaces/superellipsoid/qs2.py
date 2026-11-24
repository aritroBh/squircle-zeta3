# qs2.py -- exact arithmetic in Q(sqrt2), S-unit factorisation (S={2,3}), and a weight-2
# "dilogarithm algebra": linear combinations of Li2(x) (x real S-units of Q(sqrt2), x<1),
# pi^2, pi*alpha, and products of the generator logs  g = (log sqrt2, log 3, L=log(1+sqrt2)).
# Coefficients are themselves elements of Q(sqrt2).  Everything is exact; mpmath only for checks.
from fractions import Fraction as Fr
import mpmath as mp, math, itertools

class Q2:
    __slots__ = ("p", "q")
    def __init__(self, p, q=0):
        self.p = Fr(p); self.q = Fr(q)
    def __add__(s, o): o = _c(o); return Q2(s.p+o.p, s.q+o.q)
    __radd__ = __add__
    def __neg__(s): return Q2(-s.p, -s.q)
    def __sub__(s, o): return s + (-_c(o))
    def __rsub__(s, o): return _c(o) - s
    def __mul__(s, o): o = _c(o); return Q2(s.p*o.p + 2*s.q*o.q, s.p*o.q + s.q*o.p)
    __rmul__ = __mul__
    def norm(s): return s.p*s.p - 2*s.q*s.q
    def conj(s): return Q2(s.p, -s.q)
    def inv(s): n = s.norm(); return Q2(s.p/n, -s.q/n)
    def __truediv__(s, o): return s*_c(o).inv()
    def __rtruediv__(s, o): return _c(o)*s.inv()
    def __pow__(s, k):
        r = Q2(1); b = s if k >= 0 else s.inv()
        for _ in range(abs(k)): r = r*b
        return r
    def __eq__(s, o): o = _c(o); return s.p == o.p and s.q == o.q
    def __hash__(s): return hash((s.p, s.q))
    def mp(s): return mp.mpf(s.p.numerator)/s.p.denominator + mp.sqrt(2)*mp.mpf(s.q.numerator)/s.q.denominator
    def __float__(s): return float(s.p) + math.sqrt(2)*float(s.q)
    def is_zero(s): return s.p == 0 and s.q == 0
    def __repr__(s):
        if s.q == 0: return str(s.p)
        if s.p == 0: return "%s*r2" % s.q
        return "(%s%+s*r2)" % (s.p, s.q) if s.q < 0 else "(%s+%s*r2)" % (s.p, s.q)

def _c(o): return o if isinstance(o, Q2) else Q2(o)

R2 = Q2(0, 1); A = Q2(1, 1); B = Q2(-1, 1)   # sqrt2, a = 1+sqrt2, b = sqrt2-1

def vval(n, p):
    n = abs(n); k = 0
    if n == 0: return None
    num, den = n.numerator, n.denominator
    while num % p == 0: num //= p; k += 1
    while den % p == 0: den //= p; k -= 1
    return k, Fr(num, den)

def factor(x):
    """x = sign * sqrt2^i * 3^j * a^k exactly; returns (sign, i, j, k) or None if not an S-unit."""
    x = _c(x)
    N = x.norm()
    if N == 0: return None
    i, rest = vval(N, 2)
    j2, rest = vval(rest, 3)
    if rest != 1 or j2 % 2: return None
    j = j2//2
    u = x/(R2**i * Q2(3)**j)
    fu = float(u)
    if fu == 0: return None
    k = round(math.log(abs(fu))/math.log(1+math.sqrt(2)))
    for sg in (1, -1):
        if u == sg*(A**k): return (sg, i, j, k)
    return None

# ---------------------------------------------------------------- weight-2 algebra
GEN = ("lr2", "l3", "L")   # log sqrt2, log 3, log(1+sqrt2)

def logvec(x):
    """exponent vector of log|x| in GEN."""
    f = factor(x)
    if f is None: raise ValueError("not an S-unit: %r" % (x,))
    return (Fr(f[1]), Fr(f[2]), Fr(f[3]))

class W2:
    """dict: ('Li2', Q2) | ('pi2',) | ('pia',) | ('pig', i) | ('gg', i, j) (i<=j) -> Q2 coefficient"""
    def __init__(self, d=None): self.d = dict(d or {})
    def add(self, key, c):
        c = _c(c)
        if c.is_zero(): return self
        v = self.d.get(key, Q2(0)) + c
        if v.is_zero(): self.d.pop(key, None)
        else: self.d[key] = v
        return self
    def __add__(s, o):
        r = W2(s.d)
        for k, v in o.d.items(): r.add(k, v)
        return r
    def __neg__(s): return W2({k: -v for k, v in s.d.items()})
    def __sub__(s, o): return s + (-o)
    def scale(s, c):
        c = _c(c); return W2({k: v*c for k, v in s.d.items() if not (v*c).is_zero()})
    def num(s):
        g = [mp.log(mp.sqrt(2)), mp.log(3), mp.log(1+mp.sqrt(2))]
        alpha = mp.atan(1/mp.sqrt(2))
        tot = mp.mpf(0)
        for k, c in s.d.items():
            if k[0] == 'Li2': val = mp.polylog(2, k[1].mp())
            elif k[0] == 'pi2': val = mp.pi**2
            elif k[0] == 'pia': val = mp.pi*alpha
            elif k[0] == 'pig': val = mp.pi*g[k[1]]
            elif k[0] == 'gg': val = g[k[1]]*g[k[2]]
            elif k[0] == 'g': val = g[k[1]]
            elif k[0] == '1': val = mp.mpf(1)
            tot += c.mp()*val
        return tot
    def li2keys(s): return [k[1] for k in s.d if k[0] == 'Li2']
    def __repr__(s):
        def kn(k):
            if k[0] == 'Li2': return "Li2(%r)" % (k[1],)
            if k[0] == 'pi2': return "pi^2"
            if k[0] == 'pia': return "pi*alpha"
            if k[0] == 'pig': return "pi*%s" % GEN[k[1]]
            if k[0] == 'g': return GEN[k[1]]
            if k[0] == '1': return "1"
            return "%s*%s" % (GEN[k[1]], GEN[k[2]])
        return " + ".join("%r*%s" % (v, kn(k)) for k, v in sorted(s.d.items(), key=lambda kv: str(kv[0])))

def Li2(x, c=1): return W2().add(('Li2', _c(x)), c)
def PI2(c=1): return W2().add(('pi2',), c)
def PIA(c=1): return W2().add(('pia',), c)

def logprod(x, y, c=1):
    """c * log|x| * log|y| as W2 (x, y S-units)."""
    u = logvec(x); v = logvec(y); r = W2()
    for i in range(3):
        for j in range(3):
            if u[i]*v[j] == 0: continue
            key = ('gg', min(i, j), max(i, j))
            r.add(key, _c(c)*Q2(u[i]*v[j]))
    return r

def logsq(x, c=1): return logprod(x, x, c)

def pilog(x, c=1):
    u = logvec(x); r = W2()
    for i in range(3):
        if u[i] != 0: r.add(('pig', i), _c(c)*Q2(u[i]))
    return r

# ---------------------------------------------------------------- functional equations (real, x<1)
def rel_reflection(x):      # 0<x<1:  Li2(x)+Li2(1-x)+log x log(1-x) - pi^2/6 = 0
    return Li2(x) + Li2(1-x) + logprod(x, 1-x) + PI2(Fr(-1, 6))
def rel_landen(x):          # x<1, x!=0:  Li2(x)+Li2(x/(x-1)) + (1/2) log^2(1-x) = 0
    return Li2(x) + Li2(x/(x-1)) + logsq(1-x, Fr(1, 2))
def rel_inversion(x):       # x<0:  Li2(x)+Li2(1/x) + pi^2/6 + (1/2) log^2(-x) = 0
    return Li2(x) + Li2(1/x) + PI2(Fr(1, 6)) + logsq(x, Fr(1, 2))
def rel_dup(x):             # |x|<=1: Li2(x^2) - 2Li2(x) - 2Li2(-x) = 0
    return Li2(x*x) + Li2(x, -2) + Li2(-x, -2)
def rel_abel(x, y):         # x<1,y<1,x+y<1
    return (Li2(x/(1-y)) + Li2(y/(1-x)) - Li2(x*y/((1-x)*(1-y))) - Li2(x) - Li2(y) - logprod(1-x, 1-y))

def canon_zero_li2(w):
    """drop Li2(0) and evaluate Li2(1)=pi^2/6, Li2(-1)=-pi^2/12, Li2(1/2) kept (it has a symbol? no: reduce)."""
    r = W2()
    for k, v in w.d.items():
        if k[0] == 'Li2' and k[1].is_zero(): continue
        if k[0] == 'Li2' and k[1] == Q2(1): r.add(('pi2',), v*Q2(Fr(1, 6))); continue
        if k[0] == 'Li2' and k[1] == Q2(-1): r.add(('pi2',), v*Q2(Fr(-1, 12))); continue
        r.add(k, v)
    return r

def ONEW(c=1): return W2().add(('1',), c)
def logw(x, c=1):
    """c * log|x| as weight-1 W2"""
    u = logvec(x); r = W2()
    for i in range(3):
        if u[i] != 0: r.add(('g', i), _c(c)*Q2(u[i]))
    return r
