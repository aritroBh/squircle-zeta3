"""W2SPLIT: exact check that the reduction of Definition g1:def:red splits M_k.
Own code (python-flint, exact arithmetic in Q(sqrt2)(v)); independent of paper/ieee/proofs/scripts.
Usage: python -B split_check.py KMAX
For every level 2<=k<=KMAX and every canonical key (r,s), r>=s, d=r+s>=1, checks:
  (P) parity: W + iota(W) has no pole at +-i, where iota(W)(v)=W(-1/v)/v^2  [the lemma]
  (R) residue of W at +-i is 0 (so W - Hp(W) has no pole at +-i)            [the claim of the caveat]
  (L) W - Hp(W) has only simple poles, all in {0,+-a,+-b}, none at b, and is O(1/v) at infinity
  (K) W_10 = R_w2 (for k>=2)
  (G) for d>=2, G + G(-1/v) is a constant
Initial data from the two-variable generating function Psi(v;X,Y) (Stirling-number expansion), not from the
paper's scripts.
"""
import sys, time
from fractions import Fraction as Fr
import flint as _f
FQ = lambda a, b=1: _f.fmpq(a, b)
from math import comb, factorial
import flint
P = flint.fmpq_poly
V = P([0, 1]); ONE = P([1]); ZERO = P([])
DD = P([1, 0, 1])                 # v^2+1
Q1 = P([-1, 2, 1])                # v^2+2v-1 = (v-b)(v+a)
Q2 = P([-1, -2, 1])               # v^2-2v-1 = (v-a)(v+b)

def pgcd(a, b):
    return a.gcd(b)

class RF:
    """(n0 + sqrt2*n1)/den with n0,n1,den in Q[v]."""
    __slots__ = ('n0', 'n1', 'den')
    def __init__(self, n0, n1=None, den=None):
        self.n0 = n0; self.n1 = n1 if n1 is not None else P([]); self.den = den if den is not None else ONE
        self._norm()
    def _norm(self):
        if self.n0.is_zero() and self.n1.is_zero():
            self.den = ONE; return
        g = pgcd(pgcd(self.n0, self.n1), self.den)
        if g.degree() > 0:
            self.n0 = divexact(self.n0, g)
            self.n1 = divexact(self.n1, g)
            self.den = divexact(self.den, g)
        lc = self.den.leading_coefficient()
        if lc != 1:
            self.n0 = self.n0 / lc; self.n1 = self.n1 / lc; self.den = self.den / lc
    def iszero(self):
        return self.n0.is_zero() and self.n1.is_zero()
    def __add__(self, o):
        if isinstance(o, RF):
            g = pgcd(self.den, o.den)
            a = divexact(o.den, g); b = divexact(self.den, g)
            return RF(self.n0 * a + o.n0 * b, self.n1 * a + o.n1 * b, self.den * a)
        raise TypeError
    def __neg__(self):
        return RF(-self.n0, -self.n1, self.den)
    def __sub__(self, o):
        return self + (-o)
    def __mul__(self, o):
        if isinstance(o, RF):
            return RF(self.n0 * o.n0 + 2 * self.n1 * o.n1, self.n0 * o.n1 + self.n1 * o.n0, self.den * o.den)
        if isinstance(o, tuple):   # element x + sqrt2*y of Q(sqrt2)
            x, y = [flint.fmpq(z.numerator, z.denominator) if isinstance(z, Fr) else z for z in o]
            return RF(self.n0 * x + 2 * self.n1 * y, self.n0 * y + self.n1 * x, self.den)
        if isinstance(o, Fr): o = flint.fmpq(o.numerator, o.denominator)
        return RF(self.n0 * o, self.n1 * o, self.den)
    __rmul__ = __mul__
    def deriv(self):
        d = self.den; dd = d.derivative()
        return RF(self.n0.derivative() * d - self.n0 * dd, self.n1.derivative() * d - self.n1 * dd, d * d)
    def compose_inv(self, sign):
        """W(sign/v)/v^2."""
        n = max(self.n0.degree(), self.n1.degree(), self.den.degree())
        def rev(p):
            # v^n p(sign/v)
            c = p.coeffs(); out = [0] * (n + 1)
            for j, cj in enumerate(c):
                out[n - j] = cj * (sign ** j)
            return P(out)
        # W(sign/v) = rev(n)/rev(den);  /v^2
        return RF(rev(self.n0), rev(self.n1), rev(self.den) * V * V)

def divexact(a, b):
    q, r = divmod(a, b)
    assert r.is_zero()
    return q

def fold(W): return W.compose_inv(1)
def iota(W): return W.compose_inv(-1)

# ---- Q(sqrt2) numbers as pairs of Fractions
def qmul(x, y): return (x[0] * y[0] + 2 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])
def qadd(x, y): return (x[0] + y[0], x[1] + y[1])
def qinv(x):
    n = x[0] * x[0] - 2 * x[1] * x[1]
    return (x[0] / n, -x[1] / n)
def peval(p, x):
    acc = (FQ(0), FQ(0))
    for c in reversed(p.coeffs()):
        acc = qadd(qmul(acc, x), (c, FQ(0)))
    return acc
def rfeval(W, x):
    n0 = peval(W.n0, x); n1 = peval(W.n1, x); d = peval(W.den, x)
    num = qadd(n0, qmul((FQ(0), FQ(1)), n1))
    return qmul(num, qinv(d))
A_ = (FQ(1), FQ(1)); B_ = (FQ(-1), FQ(1))

def split_D(W):
    """Return (m, Anum0, Anum1, rest) with W = (A0+s2 A1)/DD^m + rest, deg A < 2m, rest without pole at +-i."""
    den = W.den; m = 0; rest = den
    while True:
        q, r = divmod(rest, DD)
        if r.is_zero(): rest = q; m += 1
        else: break
    if m == 0:
        return 0, ZERO, ZERO, W
    Dm = DD ** m
    g, al, be = Dm.xgcd(rest)          # al*Dm + be*rest = g = 1
    assert g == ONE
    # W = N/(Dm*rest) = N*be/Dm + N*al/rest
    A0 = (W.n0 * be) % Dm; A1 = (W.n1 * be) % Dm
    pp = RF(A0, A1, Dm)
    return m, A0, A1, W - pp

def hermite(m, A0, A1):
    """A/DD^m = (G)' + (c1 v + c0)/DD, A over Q(sqrt2). Returns G (RF) and residual numerator pair."""
    G = RF(ZERO)
    a0, a1 = A0, A1
    while m >= 2:
        # find B (deg<2) with A + 2(m-1) v B == 0 mod DD;  v^{-1} = -v mod DD
        # B = -A * inv(2(m-1)v) = A * v/(2(m-1)) mod DD
        B0 = (a0 * V) % DD * flint.fmpq(1, 2 * (m - 1)); B1 = (a1 * V) % DD * flint.fmpq(1, 2 * (m - 1))
        n0 = a0 + 2 * (m - 1) * V * B0; n1 = a1 + 2 * (m - 1) * V * B1
        n0 = divexact(n0, DD) - B0.derivative(); n1 = divexact(n1, DD) - B1.derivative()
        G = G + RF(B0, B1, DD ** (m - 1))
        a0, a1 = n0, n1; m -= 1
    if m == 1:
        return G, (a0, a1)
    return G, (ZERO, ZERO)

def Hp_parts(W):
    m, A0, A1, rest = split_D(W)
    G, (c0, c1) = hermite(m, A0, A1) if m > 0 else (RF(ZERO), (ZERO, ZERO))
    return G, (c0, c1), rest

# ---- basic functions
t = RF(Q1, None, 4 * V)                 # (v^2+2v-1)/(4v)
omt = RF(-Q2, None, 4 * V)              # 1-t = -(v^2-2v-1)/(4v)
dtdv = RF(DD, None, 4 * V * V)
D1 = RF(DD, None, V * Q1)
D2 = RF(DD, None, V * Q2)
c2 = RF(8 * V * V * Q2 * Q2, None, 16 * V * V * DD * DD)  # (1-t)^2/q^2 = c^2
cp2 = RF(8 * V * V * Q1 * Q1, None, 16 * V * V * DD * DD)  # t^2/q^2 = c'^2
# S*dt/dv = q^2/(sqrt2 v t(1-t)) = (sqrt2/2) * q^2/(v t (1-t))
SDT = RF(ZERO, -DD * DD, V * Q1 * Q2)   # S dt/dv = -sqrt2 (v^2+1)^2/(v Q1 Q2)
# check: t(1-t) = -(Q1*Q2)/(16 v^2) -> 1/(v t(1-t)) = -16 v^2/(v Q1 Q2)

def stirling2(n, k, memo={}):
    if (n, k) in memo: return memo[(n, k)]
    if n == k: r = 1
    elif k == 0 or k > n: r = 0
    else: r = k * stirling2(n - 1, k) + stirling2(n - 1, k - 1)
    memo[(n, k)] = r; return r

def ecoef(r, a):
    """[X^r](e^{2X}-1)^a = 2^r a! S(r,a)/r!"""
    return Fr(2 ** r * factorial(a) * stirling2(r, a), factorial(r))

def binhalf(j):
    num = Fr(1)
    for i in range(j): num *= (Fr(1, 2) - i)
    return num / factorial(j)

POW_c2 = {}; POW_cp2 = {}
def powc(base, store, n):
    if n not in store:
        store[n] = RF(ONE) if n == 0 else powc(base, store, n - 1) * base
    return store[n]

def Rtilde(r, s):
    """[X^r Y^s] Psi(v;X,Y), Psi = -1/2 (t^{-1}e^X + (1-t)^{-1}e^Y - S sqrt(1+u)) dt/dv."""
    d = r + s
    acc = RF(ZERO)
    for j in range(1, d + 1):
        inner = RF(ZERO)
        for aa in range(0, j + 1):
            bb = j - aa
            co = comb(j, aa) * ecoef(r, aa) * ecoef(s, bb)
            if co == 0: continue
            inner = inner + powc(c2, POW_c2, aa) * powc(cp2, POW_cp2, bb) * co
        if not inner.iszero():
            acc = acc + inner * binhalf(j)
    out = SDT * acc * Fr(1, 2)
    if s == 0:
        out = out - D1 * Fr(1, 2 * factorial(r))
    if r == 0:
        out = out + D2 * Fr(1, 2 * factorial(s))
    return out

# R_w2 = 1/(2v) + 1/(2(v-a)) - 1/(v+a) - 1/(2(v+b)); build over Q(sqrt2): use partial fractions via Q2,Q1
# 1/(2(v-a)) - 1/(2(v+b)) = (a+b)/(2 Q2) = sqrt2/Q2 ;  -1/(v+a) = -(v-b)/Q1 = -(v+1-sqrt2)/Q1
Rw2 = RF(ONE, None, 2 * V) + RF(ZERO, ONE, Q2) + RF(-P([1, 1]), ONE, Q1)

def pole_info_real(S):
    """For S with no pole at +-i: check den | v*Q1*Q2, deg num < deg den, and no pole at b."""
    ok_den = ((V * Q1 * Q2) % S.den).is_zero()
    ok_inf = max(S.n0.degree(), S.n1.degree()) < S.den.degree() or S.iszero()
    # pole at b: S*(v-b) at b; since den squarefree here, no pole iff numerator vanishes at b when den(b)=0
    nb = qadd(peval(S.n0, B_), qmul((FQ(0), FQ(1)), peval(S.n1, B_)))
    db = peval(S.den, B_)
    nopole_b = (db != (FQ(0), FQ(0))) or (nb == (FQ(0), FQ(0)))
    return ok_den, ok_inf, nopole_b

def run_level(k, log):
    keys = [(r, k_ - r) for k_ in range(1, k + 1) for r in range(k_, -1, -1) if r >= k_ - r]
    R = {}
    for r in range(0, k + 1):
        s = k - r
        if r >= s:
            Rt = Rtilde(r, s)
            R[(r, s)] = Rt + fold(Rtilde(s, r)) if r > s else Rt
    W = {}; G = {}
    allok = True; nres = 0
    for d in range(k, 0, -1):
        for r in range(d, -1, -1):
            s = d - r
            if r < s: continue
            w = R.get((r, s), RF(ZERO)) if d == k else RF(ZERO)
            if (r + 1, s) in G: w = w - G[(r + 1, s)] * D1 * (r + 1)
            if r > s and (r, s + 1) in G: w = w - G[(r, s + 1)] * D2 * (s + 1)
            if s == r - 1 and (r, r) in G: w = w - fold(G[(r, r)] * D1) * r
            W[(r, s)] = w
            Gp, (c0, c1), rest = Hp_parts(w)
            # parity check (P)
            par = w + iota(w)
            mP, A0P, A1P, _ = split_D(par)
            okP = A0P.is_zero() and A1P.is_zero()
            okR = c0.is_zero() and c1.is_zero()
            S = rest if okR else None
            if S is not None:
                okL = pole_info_real(S)
            else:
                okL = (False, False, False)
            if d >= 2:
                gb = rfeval(Gp, B_)
                G[(r, s)] = Gp - RF(P([gb[0]]), P([gb[1]]), ONE)
                Gi = iota(G[(r, s)]) * RF(V * V)          # G(-1/v)
                sm = G[(r, s)] + Gi
                okG = sm.iszero() or (sm.den.degree() == 0 and max(sm.n0.degree(), sm.n1.degree()) <= 0)
                if not okG:
                    allok = False; log(f"  k={k} key=({r},{s}) G + G(-1/v) not constant")
            nres += 1
            line_ok = okP and okR and all(okL)
            if not line_ok:
                allok = False
                log(f"  k={k} key=({r},{s}) parity={okP} residue0={okR} real-letters={okL}")
    okK = (W[(1, 0)] - Rw2).iszero()
    return allok, okK, nres, W

if __name__ == '__main__':
    kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    def log(s): print(s, flush=True)
    log(f"W2SPLIT split_check.py, levels 2..{kmax}")
    for k in range(2, kmax + 1):
        t0 = time.time()
        allok, okK, n, W = run_level(k, log)
        log(f"k={k:2d}: {n:3d} keys; (P), (R), (L), (G) for every key: {allok}; W10 - R_w2 = 0: {okK}; {time.time()-t0:.1f}s")

def controls(log):
    """Negative controls: the checks can fail."""
    # 1. a simple pole at +-i is detected
    _, (c0, c1), _ = Hp_parts(RF(ONE, None, DD))
    log(f"control: residual of 1/(v^2+1) detected (nonzero): {not (c0.is_zero() and c1.is_zero())}")
    # 2. 1/(v^2+1)^2 has a residue at i (Hermite must leave a simple part)
    _, (c0, c1), _ = Hp_parts(RF(ONE, None, DD * DD))
    log(f"control: 1/(v^2+1)^2 leaves a simple part: {not (c0.is_zero() and c1.is_zero())}")
    # 3. parity with Fold in place of iota fails on the level-4 kernels
    _, _, _, W = run_level(4, lambda s: None)
    bad = 0; tot = 0; orders = []
    for key, w in W.items():
        m, A0, A1, _ = split_D(w + fold(w)); tot += 1
        if not (A0.is_zero() and A1.is_zero()): bad += 1
        orders.append((key, split_D(w)[0]))
    log(f"control: W + Fold(W) has a pole at +-i for {bad} of {tot} level-4 keys (iota is the right map)")
    log(f"control: exponent m of (v^2+1)^m in the denominators of W at level 4: {orders}")
    # 4. the unreduced initial data at level 4 do have high poles
    w = Rtilde(3, 1)
    log(f"control: exponent of (v^2+1) in den of Rtilde_31: {split_D(w)[0]}; iota-odd mod regular: "
        f"{split_D(w + iota(w))[1].is_zero() and split_D(w + iota(w))[2].is_zero()}")

if __name__ == '__main__' and len(sys.argv) > 2 and sys.argv[2] == 'controls':
    controls(lambda s: print(s, flush=True))
