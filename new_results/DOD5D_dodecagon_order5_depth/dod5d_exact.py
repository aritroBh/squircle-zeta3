"""dod5d_exact.py -- exact Q(sqrt3) machinery for DOD5D.

QS3 = (Fraction p, Fraction q) meaning p + q*sqrt(3).
The overall factor u = sqrt(2-sqrt3) = (sqrt6-sqrt2)/2 is pulled out globally:
every E5 term has L+h+s0 degree exactly 1, and L/u, h/u, s0/u, alpha, beta, c0_i, kappa_i
are all in Q(sqrt3). Monomial coefficients are stored as QS3 = coefficient/u.

Edge (wedge (0,pi/6)): j0=3, j1=4, nonloc=[0,1,2,5]; alpha=beta=(2-r3)/2;
L=u(2-r3), h=u/2, s0=-u(2-r3)/2; kappas and c0s below.
Formal log slots in monomial keys: (eLT,eLU,eW0,eW1,eW2,eW3,eCA,eCC0,eCC1) where
  LT=log t, LU=log(1-t), W_i=log(1+kappa_i t),
  CA=log(alpha)=log(beta)=log(c0_2)=-2Y, CC0=log(c0_0)=-X, CC1=log(c0_1)=log(c0_5)=Z-2Y,
with X=log(2+sqrt3), Y=log(1+sqrt3), Z=log(sqrt3) (V basis).
"""
from fractions import Fraction as F
from itertools import product as iproduct, combinations_with_replacement as cwr
from math import comb, factorial
import pickle

ZERO = (F(0), F(0)); ONE = (F(1), F(0))
def q3(r, s=0): return (F(r), F(s))
def q3add(a, b): return (a[0]+b[0], a[1]+b[1])
def q3neg(a): return (-a[0], -a[1])
def q3mul(a, b): return (a[0]*b[0]+3*a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def q3div(a, b):
    n = b[0]*b[0]-3*b[1]*b[1]
    return ((a[0]*b[0]-3*a[1]*b[1])/n, (a[1]*b[0]-a[0]*b[1])/n)
def q3pow(a, n):
    r = ONE
    if n < 0: a = q3div(ONE, a); n = -n
    for _ in range(n): r = q3mul(r, a)
    return r

R3 = (F(0), F(1))            # sqrt(3)
TWOmR3 = (F(2), F(-1))       # 2-sqrt(3)
ALPHA = (F(1), F(-1,2))      # (2-sqrt3)/2
Lq = TWOmR3                  # L/u
Hq = (F(1,2), F(0))          # h/u
S0q = (F(-1), F(1,2))        # s0/u = -(2-sqrt3)/2
KAP = [(F(-1), F(1,2)), (F(-1), F(2,3)), (F(-1), F(1)), (F(-1), F(1,3))]   # kappa_0..3
C0S = [TWOmR3, (F(-3,2), F(1)), ALPHA, (F(-3,2), F(1))]                    # c0 of nonloc [0,1,2,5]
CCslot = [7, 8, 6, 8]        # const-log slot: c0_0->CC0(slot7), c0_1->CC1(8), c0_2->CA(6), c0_5->CC1(8)
# letters -1/kappa_i (QS3):
LETTERS = [q3div((F(-1),F(0)), k) for k in KAP]  # p0..p3

def binom_log_exp(a, slot_const, slot_var):
    """(CC+W)^a -> {(logexp_delta, binom)}: only slots slot_const, slot_var."""
    return [((slot_const, a-i, slot_var, i), comb(a, i)) for i in range(a+1)]

def mult_dicts(A, B):
    """multiply two dicts (n, logexp9)->QS3 (polynomial in t and formal logs)."""
    out = {}
    for (n1, le1), c1 in A.items():
        for (n2, le2), c2 in B.items():
            le = tuple(a+b for a, b in zip(le1, le2))
            key = (n1+n2, le)
            out[key] = q3add(out.get(key, ZERO), q3mul(c1, c2))
    return {k: v for k, v in out.items() if v != ZERO}

def scale_dict(A, c):
    return {k: q3mul(v, c) for k, v in A.items()}

_TAUGAM_CACHE = {}
def taugam(kind, k, e):
    """expansion of tau_k^e (kind='tau') or gam_k^e (kind='gam') as dict (n, logexp9)->QS3."""
    if (kind, k, e) in _TAUGAM_CACHE: return _TAUGAM_CACHE[(kind, k, e)]
    out = {}
    for ms in cwr(range(4), e):              # multiset of nonloc indices
        mult = factorial(e)
        for i in range(4): mult //= factorial(ms.count(i))
        # coefficient product; gam_k = sum ap_i (log a_i)^k with ap_i = c0_i*kappa_i CONSTANT (no t-factor)
        coef = ONE
        for i in ms:
            coef = q3mul(coef, C0S[i] if kind == 'tau' else q3mul(C0S[i], KAP[i]))
        # t-part: prod (1+kappa_i t) for tau (a_i(t) factor); trivial for gam
        tdict = {0: ONE}
        if kind == 'tau':
            for i in ms:
                nd = dict(tdict)
                for n, c in tdict.items(): nd[n+1] = q3add(nd.get(n+1, ZERO), q3mul(c, KAP[i]))
                tdict = nd
        # log-part: prod (CC_i+W_i)^k
        logd = {((0,)*9): ONE}
        for i in ms:
            nd = {}
            for le, c in logd.items():
                for (sc, ec, sv, ev), b in binom_log_exp(k, CCslot[i], 2+i):
                    le2 = list(le); le2[sc] += ec; le2[sv] += ev
                    le2 = tuple(le2)
                    nd[le2] = q3add(nd.get(le2, ZERO), q3mul(c, q3(b)))
            logd = nd
        for n, ct in tdict.items():
            for le, cl in logd.items():
                key = (n, le)
                out[key] = q3add(out.get(key, ZERO), q3mul(q3(mult), q3mul(coef, q3mul(ct, cl))))
    out = {k2: v for k2, v in out.items() if v != ZERO}
    _TAUGAM_CACHE[(kind, k, e)] = out
    return out

_L01_CACHE = {}
def l01(which, a):
    """(CA+LT)^a (which=0) or (CA+LU)^a (which=1) as dict (0, logexp9)->QS3."""
    if (which, a) in _L01_CACHE: return _L01_CACHE[(which, a)]
    slot_var = 0 if which == 0 else 1
    out = {}
    for i in range(a+1):
        le = [0]*9; le[slot_var] = i; le[6] = a-i
        out[(0, tuple(le))] = q3(comb(a, i))
    _L01_CACHE[(which, a)] = out
    return out

def param_factor(eL, eh, es, eab):
    """L^eL h^eh s0^es alpha^i beta^j (i+j=eab) / u  as QS3. eL+eh+es == 1 required."""
    assert eL+eh+es == 1
    r = q3pow(Lq, eL)      # handles negative powers
    r = q3mul(r, q3pow(Hq, eh))
    r = q3mul(r, q3pow(S0q, es))
    r = q3mul(r, q3pow(ALPHA, eab))
    return r

def expand_E5(terms):
    """terms: list of (exp18, rat) from E5_terms.pkl. Returns dict (n, logexp9) -> QS3 (coefficient/u)."""
    out = {}
    for exp, rat in terms:
        base = q3mul(param_factor(exp[3], exp[4], exp[5], exp[1]+exp[2]), q3(rat))
        cur = {(exp[0], (0,)*9): base}
        if exp[6]: cur = mult_dicts(cur, l01(0, exp[6]))
        if exp[7]: cur = mult_dicts(cur, l01(1, exp[7]))
        for i in range(5):
            if exp[8+i]:  cur = mult_dicts(cur, taugam('tau', i+1, exp[8+i]))
            if exp[13+i]: cur = mult_dicts(cur, taugam('gam', i+1, exp[13+i]))
        for key, c in cur.items():
            out[key] = q3add(out.get(key, ZERO), c)
    return {k: v for k, v in out.items() if v != ZERO}
