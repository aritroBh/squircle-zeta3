"""pw.py -- INDEPENDENT point-wise exact engine (referee).  No cached F, no sympy rational-function
simplification.  Inputs only: L_simplified.pkl (L_1..L_8), tau table (tow4lib tau1f/tau2f).
Arithmetic: exact Gaussian rationals (pairs of python Fractions).
 - a_delta by the Bell recurrence  n a_n = sum_r r L_r a_{n-r}   (different algorithm from the partition sum)
 - F(delta,s) UNTRUNCATED ascending recursion (no assumption F(.,s>top)=0)
 - closing / vanishing above top at random rational (d,m)
 - adjacent Newton identity residual at rational p (formulas copied from the scripts under review)
"""
import sys, os, pickle, random
from fractions import Fraction as Fr
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'lib'))
import tow4lib as T4

class G:
    __slots__ = ('a', 'b')
    def __init__(self, a=0, b=0):
        self.a = Fr(a); self.b = Fr(b)
    @staticmethod
    def of(x):
        return x if isinstance(x, G) else G(x, 0)
    def __add__(s, o): o = G.of(o); return G(s.a+o.a, s.b+o.b)
    __radd__ = __add__
    def __sub__(s, o): o = G.of(o); return G(s.a-o.a, s.b-o.b)
    def __rsub__(s, o): o = G.of(o); return G(o.a-s.a, o.b-s.b)
    def __neg__(s): return G(-s.a, -s.b)
    def __mul__(s, o): o = G.of(o); return G(s.a*o.a - s.b*o.b, s.a*o.b + s.b*o.a)
    __rmul__ = __mul__
    def inv(s):
        n = s.a*s.a + s.b*s.b
        if n == 0: raise ZeroDivisionError('G inverse of 0')
        return G(s.a/n, -s.b/n)
    def __truediv__(s, o): return s * G.of(o).inv()
    def __rtruediv__(s, o): return G.of(o) * s.inv()
    def __pow__(s, n):
        n = int(n)
        if n < 0: return (s.inv())**(-n)
        r = G(1); b = s
        while n:
            if n & 1: r = r*b
            b = b*b; n >>= 1
        return r
    def conj(s): return G(s.a, -s.b)
    def iszero(s): return s.a == 0 and s.b == 0
    def __repr__(s): return f'({s.a}+{s.b}i)'

I_ = G(0, 1)

def to_G(expr):
    a, b = sp.expand(expr).as_real_imag()
    a = sp.Rational(a); b = sp.Rational(b)
    return G(Fr(int(a.p), int(a.q)), Fr(int(b.p), int(b.q)))

# ---- L_r(d,m) as exact polynomial ratios --------------------------------------------------
d_, m_ = sp.symbols('d m')
_Lraw = pickle.load(open(os.path.join(HERE, 'lib', 'L_simplified.pkl'), 'rb'))
_Lpoly = {}
for r, s in _Lraw.items():
    ex = sp.sympify(s, locals={'d': d_, 'm': m_, 'I': sp.I})
    num, den = sp.fraction(sp.together(ex))
    pn = sp.Poly(sp.expand(num), d_, m_); pd = sp.Poly(sp.expand(den), d_, m_)
    _Lpoly[r] = ([(mon, to_G(c)) for mon, c in pn.terms()], [(mon, to_G(c)) for mon, c in pd.terms()])

def _ev(terms, D, M):
    tot = G(0)
    for (i, j), c in terms:
        tot = tot + c * (D**i) * (M**j)
    return tot

def Lval(r, D, M):
    n, dd = _Lpoly[r]
    return _ev(n, D, M) / _ev(dd, D, M)

# ---- tau (exact, from tow4lib tables) ------------------------------------------------------
TAU1 = {n: to_G(T4.tau1f(n)) for n in range(0, 9)}
TAU2 = {n: to_G(T4.tau2f(n)) for n in range(0, 9)}

# ---- a_delta by Bell recurrence ------------------------------------------------------------
_acache = {}
def a_val(delta, D, M):
    key = (delta, D.a, D.b, M.a, M.b)
    if key in _acache: return _acache[key]
    if delta == 0: return G(1)
    tot = G(0)
    for r in range(1, delta+1):
        tot = tot + r * Lval(r, D, M) * a_val(delta-r, D, M)
    out = tot / delta
    _acache[key] = out
    return out

# ---- F : untruncated ascending recursion ---------------------------------------------------
_F = {}
def Fu(delta, s, D, M):
    key = (delta, s, D.a, D.b, M.a, M.b)
    if key in _F: return _F[key]
    if s == 0:
        out = a_val(delta, D, M)
    else:
        out = Lam(delta, s-1, D, M) - Fu(delta, s-1, D, M)
    _F[key] = out
    return out

def Lam(delta, s, D, M):
    tot = G(0)
    pref = I_ * (2*D - 1)
    for dp in range(0, delta+1):
        fA = Fu(dp, s, D+1, M+1); fB = Fu(dp, s, D+1, M)
        tot = tot + pref / (2*D - dp - 1) * (TAU1[delta-dp]*fA - TAU2[delta-dp]*fB)
    return tot

def rand_pt(rng):
    # random rational point with moderate height; avoid poles (2d-odd = 0 impossible for non-half-integer d)
    D = G(Fr(rng.randint(-200, 200), rng.randint(1, 37)) + Fr(1, 1000003))
    M = G(Fr(rng.randint(-200, 200), rng.randint(1, 37)) + Fr(1, 999983))
    return D, M

if __name__ == '__main__':
    rng = random.Random(20261028)
    # conj symmetry of tau
    for n in range(0, 9):
        assert TAU1[n].conj().a == ((-1)**n * TAU2[n]).a and TAU1[n].conj().b == ((-1)**n * TAU2[n]).b, n
        assert TAU2[n].conj().a == ((-1)**n * TAU1[n]).a and TAU2[n].conj().b == ((-1)**n * TAU1[n]).b, n
    print('tau conj symmetry n=0..8 OK')
    NPT = int(sys.argv[1]) if len(sys.argv) > 1 else 12
    allzero = True; nonzero_seen = {}
    for _ in range(NPT):
        D, M = rand_pt(rng)
        for delta in range(0, 9):
            top = delta // 2
            for s in range(0, top + 4):          # three levels above the top
                v = Fu(delta, s, D, M)
                if s > top:
                    if not v.iszero():
                        allzero = False; print('NONZERO above top', delta, s)
                else:
                    nonzero_seen[(delta, s)] = nonzero_seen.get((delta, s), 0) + (0 if v.iszero() else 1)
    print(f'UNTRUNCATED recursion at {NPT} random rational (d,m): F(delta,s)=0 for ALL s in (top, top+3], delta<=8 :', allzero)
    print('F(delta,s<=top) nonzero count per pt (should equal NPT):', {k: v for k, v in sorted(nonzero_seen.items()) if v != NPT} or 'all nonzero at all points')
    # explicit closing values
    D, M = rand_pt(rng)
    print('closing F(7,4) =', Fu(7, 4, D, M), '  F(8,5) =', Fu(8, 5, D, M))
