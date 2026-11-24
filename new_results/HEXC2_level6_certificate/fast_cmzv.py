"""Fast regularized level-6 CMZV evaluation via FLINT acb_poly power series (validation only).

Same conventions as numeric_cmzv.py: word first-letter-outermost,
  I(w;0->z) = int_0^z omega_{w0}(t) I(w[1:];0->t) dt,  omega_{-1}=dt/t, omega_j=zeta^j dt/(1-zeta^j t).
  I(w;0->1) = sum_k B(w[:k]) A(w[k:]), split at 1/2, geometric finite parts (log t -> log(1/2)).
Log-series: list P where P[j] is an acb_poly for the (log t)^j part (plain powers, no t^{-1});
shift_down produces transient t^{-1} coefficients handled inside integrate.
"""
from flint import acb, acb_poly, arb, ctx
from functools import lru_cache
from math import factorial

BITS = 230
N_ = 0
ZET = None
FORMS = {}   # letter -> acb_poly of Taylor coeffs (length N_+1)
AFORMS = {}  # B-side a_j forms

def setup(dps=62):
    global N_, ZET, BITS
    bits = int(dps * 3.322) + 40
    ctx.prec = bits
    BITS = bits
    N_ = int((dps + 15) / 0.30103) + 20
    ZET = (acb(0, 1) * acb.pi() / 3).exp()
    FORMS.clear(); AFORMS.clear()
    _A_series.cache_clear(); _B_series.cache_clear(); A.cache_clear(); B.cache_clear()

def _zero():
    return acb_poly([acb(0)])

def _one():
    return [acb_poly([acb(1)])]

def _mul(P, F):
    return [(pj * F).truncate(N_ + 1) for pj in P]

def _integrate_terms(c_neg, poly, j):
    """integrate (log t)^j * (c_neg/t + poly(t)) dt from 0 to s. Returns list Q."""
    Q = []
    # m >= 0 terms: Q_i += (-1)^{j-i} j!/i! T_{j-i+1}(poly)
    for i in range(j + 1):
        k = j - i + 1
        # T_k(poly): sum c_m/(m+1)^k x^{m+1}
        cs = poly.coeffs()
        scaled = [cs[m] / acb(m + 1) ** (k - 1) for m in range(len(cs))]
        Tk = acb_poly(scaled).integral().truncate(N_ + 1)
        coef = acb(((-1) ** (j - i)) * factorial(j) // factorial(i))
        while len(Q) <= i: Q.append(_zero())
        Q[i] = (Q[i] + coef * Tk).truncate(N_ + 1)
    # c_neg: Q_{j+1} += c_neg/(j+1) at x^0
    if c_neg != 0:
        while len(Q) <= j + 1: Q.append(_zero())
        Q[j + 1] = (Q[j + 1] + acb_poly([c_neg / (j + 1)])).truncate(N_ + 1)
    return Q

def _letter_step(P, letter, side):
    """apply outermost letter: I_new = int omega_letter * P dt."""
    if letter == -1:
        # dt/t: P/t then integrate: split each P_j into c_neg (t^0 coeff) + rest/t
        Q = []
        for j, pj in enumerate(P):
            cs = pj.coeffs()
            c0 = cs[0] if cs else acb(0)
            rest = acb_poly(cs[1:]) if len(cs) > 1 else _zero()
            for i, q in enumerate(_integrate_terms(c0, rest, j)):
                while len(Q) <= i: Q.append(_zero())
                Q[i] = (Q[i] + q).truncate(N_ + 1)
        return Q
    if side == 'A':
        F = FORMS[letter]
    else:
        F = AFORMS[letter]
    P2 = _mul(P, F)
    Q = []
    for j, pj in enumerate(P2):
        for i, q in enumerate(_integrate_terms(acb(0), pj, j)):
            while len(Q) <= i: Q.append(_zero())
            Q[i] = (Q[i] + q).truncate(N_ + 1)
    return Q

def _form(letter):
    z = ZET ** letter
    return acb_poly([z ** (m + 1) for m in range(N_ + 1)])

def _aform(j):
    # omega_j(1-s) = -a_j ds/(1 - a_j s), a_j = -zeta^j/(1-zeta^j): series -a_j^{m+1} s^m
    a = -ZET ** j / (1 - ZET ** j)
    return acb_poly([-a ** (m + 1) for m in range(N_ + 1)])

S_HALF = acb(1) / 2
LOG_HALF = None

def _eval(P, s, ls):
    tot = acb(0)
    for j, pj in enumerate(P):
        tot += ls ** j * pj(s)
    return tot

@lru_cache(maxsize=None)
def _A_series(w):
    if not w: return tuple([acb_poly([acb(1)])])
    return tuple(_letter_step(list(_A_series(w[1:])), w[0], 'A'))

@lru_cache(maxsize=None)
def A(w):
    if not w: return acb(1)
    return _eval(list(_A_series(w)), S_HALF, LOG_HALF)

@lru_cache(maxsize=None)
def _B_series(sw):
    if not sw: return tuple([acb_poly([acb(1)])])
    return tuple(_letter_step(list(_B_series(sw[1:])), sw[0], 'B'))

@lru_cache(maxsize=None)
def B(v):
    """I(v; 1/2->1). letter map (rev): -1 -> ds/(1-s) (0), 0 -> ds/s (-1), j -> a_j form."""
    if not v: return acb(1)
    sw = []
    for x in reversed(v):
        if x == -1: sw.append(0)
        elif x == 0: sw.append(-1)
        else: sw.append(x)
    return _eval(list(_B_series(tuple(sw))), S_HALF, LOG_HALF)

def cmzv_value(w):
    w = tuple(w)
    tot = acb(0)
    for k in range(len(w) + 1):
        tot += B(w[:k]) * A(w[k:])
    return tot

def finish_setup():
    global LOG_HALF, FORMS, AFORMS
    LOG_HALF = (S_HALF).log()
    FORMS = {j: _form(j) for j in range(6)}
    FORMS[-1] = None
    AFORMS = {j: _aform(j) for j in range(1, 6)}
    AFORMS[0] = acb_poly([acb(1)] * (N_ + 1))   # ds/(1-s)
