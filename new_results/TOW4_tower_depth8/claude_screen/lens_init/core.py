"""core.py -- referee's OWN implementation (lens: initial data + independent data).
Nothing here imports tow2lib/tow3lib/tow4lib/trunc_engine.  Exact arithmetic only.

Everything in the tower is sqrt(2) * (Gaussian rational) because rho = -(sqrt2/2)*(Q(i)) and the T1 map
has Q(i) coefficients (tau in Q(i)) and commutes with conj on sqrt2*Q(i).  So we divide the common
sqrt(2) out and work in Q(i) with flint.fmpq components:  u~ := u / sqrt2.
"""
import sys
sys.dont_write_bytecode = True
from flint import fmpq, fmpz
from math import comb, factorial

# ------------------------------------------------------------------ Gaussian rationals
class G:
    __slots__ = ('r', 'i')
    def __init__(self, r=0, i=0):
        self.r = r if isinstance(r, fmpq) else fmpq(r)
        self.i = i if isinstance(i, fmpq) else fmpq(i)
    def __add__(s, o):
        if not isinstance(o, G): o = G(o)
        return G(s.r + o.r, s.i + o.i)
    __radd__ = __add__
    def __sub__(s, o):
        if not isinstance(o, G): o = G(o)
        return G(s.r - o.r, s.i - o.i)
    def __rsub__(s, o):
        return G(o) - s
    def __neg__(s): return G(-s.r, -s.i)
    def __mul__(s, o):
        if not isinstance(o, G):
            o = fmpq(o) if not isinstance(o, fmpq) else o
            return G(s.r * o, s.i * o)
        return G(s.r * o.r - s.i * o.i, s.r * o.i + s.i * o.r)
    __rmul__ = __mul__
    def conj(s): return G(s.r, -s.i)
    def inv(s):
        n = s.r * s.r + s.i * s.i
        return G(s.r / n, -s.i / n)
    def __truediv__(s, o):
        if not isinstance(o, G):
            o = fmpq(o) if not isinstance(o, fmpq) else o
            return G(s.r / o, s.i / o)
        return s * o.inv()
    def iszero(s): return s.r == 0 and s.i == 0
    def __eq__(s, o):
        if not isinstance(o, G): o = G(o)
        return s.r == o.r and s.i == o.i
    def __hash__(s): return hash((int(s.r.p), int(s.r.q), int(s.i.p), int(s.i.q)))
    def __repr__(s): return f'({s.r}+{s.i}i)'
    def __pow__(s, n):
        if n < 0: return s.inv() ** (-n)
        res = G(1); b = s
        while n:
            if n & 1: res = res * b
            n >>= 1
            if n: b = b * b
        return res

ZERO = G(0); ONE = G(1); I = G(0, 1)
def ipow(n):
    return [G(1), G(0, 1), G(-1), G(0, -1)][n % 4]

def fr(x):
    """python Fraction -> fmpq"""
    return fmpq(x.numerator, x.denominator)

# ------------------------------------------------------------------ series in x (truncated), coefficient lists of G
def smul(A, B, N):
    out = [ZERO] * (N + 1)
    for a, ca in enumerate(A):
        if a > N or ca.iszero(): continue
        for b, cb in enumerate(B):
            if a + b > N: break
            if cb.iszero(): continue
            out[a + b] = out[a + b] + ca * cb
    return out

def sinv(A, N):
    a0 = A[0].inv()
    out = [a0] + [ZERO] * N
    for n in range(1, N + 1):
        acc = ZERO
        for k in range(1, n + 1):
            if k < len(A) and not A[k].iszero():
                acc = acc + A[k] * out[n - k]
        out[n] = -(acc * a0)
    return out

def gbinom(alpha, n):
    """generalized binomial coefficient binom(alpha, n), alpha integer (any sign), as fmpq"""
    r = fmpq(1)
    for t in range(n):
        r = r * fmpq(alpha - t, t + 1)
    return r

def spow_binomial(c, alpha, N):
    """(1 + c x)^alpha as series to x^N, alpha integer"""
    out = []
    cp = ONE
    for n in range(N + 1):
        out.append(cp * gbinom(alpha, n))
        cp = cp * c
    return out

# ------------------------------------------------------------------ tau_n from DL1, DL2 (own derivation)
# DL1(v) = 1/(v-b)+1/(v+a)-1/v = (v^2+1)/(v (v^2+2v-1)),  DL2(v) = 1/(v-a)+1/(v+b)-1/v = (v^2+1)/(v (v^2-2v-1)),
# a=1+sqrt2, b=sqrt2-1 (a*b=1, a-b=2).   At v=i+x:  v^2+1 = x(2i+x).
def tau_tables(nmax=9):
    N = nmax + 1                      # need x^{n+1}, n<=nmax -> up to x^{nmax+1}; DL = x*(2i+x)/den, den(0)=i*(-2)...
    xs = [ZERO, ONE]                  # series for x
    v = [I, ONE]                      # i + x
    def ser(c): return c
    v2 = smul(v, v, N + 2)
    # den1 = v*(v^2+2v-1), den2 = v*(v^2-2v-1)
    lin = [ZERO] * (N + 3)
    for sign, name in ((+1, 'd1'), (-1, 'd2')):
        q = [v2[t] + (v[t] * (2 * sign) if t < 2 else ZERO) for t in range(N + 3)]
        q[0] = q[0] - 1
        den = smul(v, q, N + 2)                       # v*(v^2 +- 2v - 1)
        num = [ZERO, G(0, 2), ONE]                    # x*(2i+x)
        # den(0) = i*(i^2 + 2i*sign -1) = i*(-2 + 2 i sign) != 0
        inv = sinv(den, N + 2)
        s = smul(num, inv, N + 2)
        if name == 'd1': t1 = s
        else: t2 = s
    T1 = [t1[n + 1] for n in range(nmax + 1)]
    T2 = [t2[n + 1] for n in range(nmax + 1)]
    return T1, T2

TAU1, TAU2 = tau_tables(9)

# ------------------------------------------------------------------ Stirling numbers (2nd kind)
_S2 = {}
def S2(n, k):
    if (n, k) in _S2: return _S2[(n, k)]
    if n == k: r = 1 if n >= 0 else 0
    elif k == 0 or k > n or n < 0: r = 0
    else: r = k * S2(n - 1, k) + S2(n - 1, k - 1)
    _S2[(n, k)] = r
    return r

def binom_half(j):
    r = fmpq(1)
    for t in range(j):
        r = r * (fmpq(1, 2) - t)
    return r / factorial(j)

# ------------------------------------------------------------------ rho engine: S_m Taylor at t*=(1+i)/2 in h, then h = x^2/(4(i+x))
TSTAR = G(fmpq(1, 2), fmpq(1, 2))
def hseries_mul(A, B, H):   # series in h truncated at h^H
    return smul(A, B, H)

def S_taylor(d, m, H=4):
    """S_m(t*+h) = sum_l c_l h^l, l<=H, from the DEFINITION
       S_m(t)=sum_{j=1}^d binom(1/2,j) q2^{d-j} sum_s binom(j,s) w(s,m) w(j-s,d-m) t^{2(j-s)} (1-t)^{2s},
       w(s,m)= s! S2(m,s) 2^m / m!.   q2(t*+h) = 2 i h + 2 h^2 (checked in checks)."""
    q = d - m
    tp = [TSTAR, ONE]                       # t*+h
    um = [ONE - TSTAR, G(-1)]               # 1-t*-h
    tp2 = smul(tp, tp, H); um2 = smul(um, um, H)
    q2 = [tp2[t] + um2[t] for t in range(H + 1)]
    # powers
    def spow(base, e):
        res = [ONE] + [ZERO] * H
        for _ in range(e): res = smul(res, base, H)
        return res
    tp2p = [[ONE] + [ZERO] * H]; um2p = [[ONE] + [ZERO] * H]
    for e in range(1, d + 1):
        tp2p.append(smul(tp2p[-1], tp2, H)); um2p.append(smul(um2p[-1], um2, H))
    out = [ZERO] * (H + 1)
    for j in range(max(1, d - H), d + 1):
        # B_j(t) series
        Bj = [ZERO] * (H + 1)
        for s in range(0, j + 1):
            # w(s,m) w(j-s, d-m)
            a = S2(m, s); b = S2(q, j - s)
            if a == 0 or b == 0: continue
            coef = fmpq(comb(j, s)) * fmpq(factorial(s) * a * 2 ** m, factorial(m)) * fmpq(factorial(j - s) * b * 2 ** q, factorial(q))
            term = smul(tp2p[j - s], um2p[s], H)
            for t in range(H + 1):
                if not term[t].iszero(): Bj[t] = Bj[t] + term[t] * coef
        qp = [ONE] + [ZERO] * H
        for _ in range(d - j): qp = smul(qp, q2, H)
        prod = smul(qp, Bj, H)
        bh = binom_half(j)
        for t in range(H + 1): out[t] = out[t] + prod[t] * bh
    return out

def a_series(d, m, N=8):
    """a_delta(d,m) = rho(delta)/rho(0), delta=0..N, and S(t*) value, from first principles:
       A(i+x)/A(i) = (1+x/i)^{2d-1} (1+x/(2i))^{2-2d} (8/G(i+x)) S(t*+h)/S(t*),  h = x^2/(4(i+x))."""
    H = N // 2
    c = S_taylor(d, m, H)
    c0 = c[0]
    # series for h(x) = x^2/(4 (i+x))
    ix = [I, ONE]
    inv_ix = sinv(ix + [ZERO] * N, N)
    hser = [ZERO, ZERO] + [inv_ix[t] * fmpq(1, 4) for t in range(N - 1)]
    hser = hser[:N + 1]
    # S ratio
    Sr = [ZERO] * (N + 1); hp = [ONE] + [ZERO] * N
    for l in range(H + 1):
        for t in range(N + 1): Sr[t] = Sr[t] + hp[t] * c[l]
        hp = smul(hp, hser, N)
    Sr = [z / c0 for z in Sr]
    # G(v) = v^4 - 6 v^2 + 1 at v=i+x, normalised by G(i)=8
    v = [I, ONE]
    v2 = smul(v, v, N); v4 = smul(v2, v2, N)
    Gs = [v4[t] - v2[t] * 6 + (ONE if t == 0 else ZERO) for t in range(N + 1)]
    assert Gs[0] == G(8)
    Ginv = sinv([z / 8 for z in Gs], N)
    f1 = spow_binomial(I.inv(), 2 * d - 1, N)                # (1 + x/i)^{2d-1}
    f2 = spow_binomial((G(0, 2)).inv(), 2 - 2 * d, N)        # (1 + x/(2i))^{2-2d}
    tot = smul(smul(smul(f1, f2, N), Ginv, N), Sr, N)
    return tot, c0

_r0cache = {}
def rho0_tilde(p, q):
    if (p, q) in _r0cache: return _r0cache[(p, q)]
    v = _rho0_tilde(p, q)
    _r0cache[(p, q)] = v
    return v

def _rho0_tilde(p, q):
    """rho_{pq}(0)/sqrt2 = (1+[p>q]) * A(i)/sqrt2,
       A(i) = -(sqrt2/2) 8^d i^{2d-1} (2i)^{2-2d} S_p(t*) / 8."""
    d = p + q
    _, c0 = a_series(d, p, 0)
    A = ipow(2 * d - 1) * ((G(0, 2)) ** (2 - 2 * d)) * c0 * fmpq(8 ** d, 8) * fmpq(-1, 2)
    return A * (2 if p > q else 1)

_aser = {}
def a_of(d, m):
    if (d, m) not in _aser:
        _aser[(d, m)] = a_series(d, m, 8)[0]
    return _aser[(d, m)]

def rho_tilde(p, q, N=8):
    return [rho0_tilde(p, q) * z for z in a_of(p + q, p)[:N + 1]]

# ------------------------------------------------------------------ T1 recursion (own implementation, spec = TOW1 proof_tower.md Sec.2)
def T1_key(i, j, delta, U):
    """u_{ij}^{(k)}(delta) from parents; U(key, dp) returns the degree-(d+1) value (G).  d=i+j."""
    d = i + j
    tot = ZERO
    for dp in range(delta + 1):
        den = 2 * d - dp - 1
        tot = tot + TAU1[delta - dp] * U((i + 1, j), dp) * fmpq(i + 1, den)
        if i > j:
            tot = tot + TAU2[delta - dp] * U((i, j + 1), dp) * fmpq(j + 1, den)
    if j == i - 1:
        tot = tot + fold(i, delta, lambda dp: U((i, i), dp))
    return tot

def fold(p, delta, Ud, binom_fn=None, sign=1):
    """Fold_p(delta;u) = -p sum_{mu=m}^{4p-4} (-i)^mu binom(mu-2,mu-m) i^{m-2} conj(g_mu),  m=4p-4-delta,
       g_mu = - sum_{dp'=0}^{4p-4-mu} u(dp') tau1_{4p-4-mu-dp'} / (4p-3-dp')."""
    bf = binom_fn or (lambda n, k: comb(n, k))
    m = 4 * p - 4 - delta
    tot = ZERO
    for mu in range(m, 4 * p - 3):
        t = 4 * p - 4 - mu
        g = ZERO
        for dp in range(t + 1):
            g = g - Ud(dp) * TAU1[t - dp] * fmpq(1, 4 * p - 3 - dp)
        gc = g.conj() if sign == 1 else g
        tot = tot + gc * ipow(-mu) * ipow(m - 2) * bf(mu - 2, mu - m)
    return tot * (-p)

# note (-i)^mu = ipow(-mu)?  (-i)^mu = (i^{-1})^mu = i^{-mu}. yes.

def keys_of_degree(n):
    return [(a, n - a) for a in range(n, (n - 1) // 2, -1)]

def tower_levelk(k, nlow, rho_fn=rho_tilde, T1=T1_key, depth=8):
    """own T1 descent at level k from the initial data of the degree-k keys down to degree nlow.
       returns dict key -> [u(0..min(depth,2n-4))]  (u~ = u/sqrt2)"""
    tower = {}
    for key in keys_of_degree(k):
        tower[key] = rho_fn(key[0], key[1], depth)
    for n in range(k - 1, nlow - 1, -1):
        dmax = min(depth, 2 * n - 4)
        for key in keys_of_degree(n):
            i, j = key
            def U(kk, dp, _t=tower):
                v = _t.get(kk)
                if v is None or dp >= len(v): return ZERO
                return v[dp]
            tower[key] = [T1(i, j, dl, U) for dl in range(dmax + 1)]
    return tower
