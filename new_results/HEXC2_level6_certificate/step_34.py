code = '''"""High-precision evaluation of regularized level-6 CMZV word values, for validation only.

Convention: word w = (w0,...,w_{n-1}), first letter outermost:
  I(w; 0->z) = int_0^z omega_{w0}(t) I(w[1:]; 0->t) dt,
  omega_{-1} = dt/t,  omega_j = zeta^j dt/(1 - zeta^j t)  (j = 0..5),  zeta = e^{i pi/3}.
Regularization: shuffle-reg with H_{dt/t}(1) = H_{dt/(1-t)}(1) = 0 (geometric finite parts).

Method: I(w;0->1) = sum_k A(w[:k]) B(w[k:]) with
  A(u) = I(u; 0->1/2)          (log-series in t at 0, log t -> log(1/2))
  B(v) = I(v; 1/2->1) = (-1)^|v| I(revmap(v); 0->1/2), s = 1-t,
       letter map: -1 -> -(ds/(1-s)),  0 -> -(ds/s),  j -> a_j ds/(1-a_j s), a_j = 1/(1-zeta^{-j}),
       log s -> log(1/2).
Convergence: worst ratio 1/2, so N ~ (dps+10)/log10(2) terms.
"""
import mpmath as mp
from functools import lru_cache

ZET = None  # set by set_dps

def set_dps(dps):
    global ZET, N
    mp.mp.dps = dps
    ZET = mp.e ** (mp.pi * mp.j / 3)
    N[0] = int((dps + 15) / 0.30103) + 20
    _clear()

N = [0]
_caches = []
def _clear():
    for c in _caches: c.cache_clear()

def _int_taylor(P, N_):
    """P: dict j -> array c[0..N_] for (log t)^j * sum c_m t^m.
    Return Q = int_0^s P(t) dt/t  as log-series (m starts at -1 internally: caller passes shifted)."""
    raise NotImplementedError

def _integrate(P, N_):
    """Q(s) = int_0^s P(t) dt  where P has powers t^m, m >= -1 encoded as array index m+1."""
    out = {}
    for j, c in P.items():
        # term t^m (log t)^j, m = idx-1
        for idx in range(len(c)):
            if c[idx] == 0: continue
            m = idx - 1
            if m == -1:
                out.setdefault(j + 1, [mp.mpc(0)] * (N_ + 2))[0] += c[idx] / (j + 1)
            else:
                # s^{m+1} sum_{i=0}^j (-1)^{j-i} j!/i! (log s)^i / (m+1)^{j-i+1}
                for i in range(j + 1):
                    coef = c[idx] * (-1) ** (j - i) * mp.factorial(j) / mp.factorial(i) / mp.mpf(m + 1) ** (j - i + 1)
                    out.setdefault(i, [mp.mpc(0)] * (N_ + 2))[idx] += coef  # idx = m+1 -> power s^{m+1}: index (m+1)+1? no: see below
    return out

def _shift_down(P):
    """divide by t: index m+1 -> m"""
    out = {}
    for j, c in P.items():
        out[j] = c[1:] + [mp.mpc(0)]
    return out

def _mul_series(P, coefs):
    """multiply each P_j by the Taylor series sum coefs[m] t^m (m>=0)."""
    N_ = len(coefs) - 1
    out = {}
    for j, c in P.items():
        new = [mp.mpc(0)] * (N_ + 2)
        for a in range(N_ + 1):
            if c[a] == 0: continue
            ca = c[a]
            for b in range(N_ + 1 - a):
                new[a + b] += ca * coefs[b]
        out[j] = new
    return out

def _form_coefs(letter, N_):
    """Taylor coefficients of omega_letter(t) = sum c_m t^m dt, m >= 0; None for dt/t."""
    if letter == -1: return None
    z = ZET ** letter
    return [z ** (m + 1) for m in range(N_ + 1)]

def _eval_logs(P, s, N_):
    tot = mp.mpc(0)
    ls = mp.log(s)
    for j, c in P.items():
        pw = sum(c[m] * s ** m for m in range(len(c)))
        tot += ls ** j * pw
    return tot

from functools import cache

def make_eval(N_):
    @lru_cache(maxsize=None)
    def A(w):
        """I(w; 0->1/2), letters -1,0..5."""
        if not w: return mp.mpc(1)
        P = _A_series(w[1:])
        if w[0] == -1:
            P = _shift_down(P)
        else:
            P = _mul_series(P, _form_coefs(w[0], N_))
        Q = _integrate(P, N_)
        return _eval_logs(Q, mp.mpf(1) / 2, N_)
    @lru_cache(maxsize=None)
    def _A_series(w):
        """log-series of I(w;0->s) around s=0, arrays length N_+2 indexed by power+1."""
        if not w: return {0: [mp.mpc(1)] + [mp.mpc(0)] * (N_ + 1)}
        P = _A_series(w[1:])
        if w[0] == -1:
            P = _shift_down(P)
        else:
            P = _mul_series(P, _form_coefs(w[0], N_))
        return _integrate(P, N_)
    @lru_cache(maxsize=None)
    def _B_series(w):
        """log-series for s-side word with letters in {-1(dt/t->ds/s), 0(ds/(1-s)), ('a',j)}."""
        if not w: return {0: [mp.mpc(1)] + [mp.mpc(0)] * (N_ + 1)}
        P = _B_series(w[1:])
        x = w[0]
        if x == -1:
            P = _shift_down(P)   # ds/s
        elif x == 0:
            P = _mul_series(P, [mp.mpc(1)] * (N_ + 1))   # ds/(1-s)
        else:
            a = 1 / (1 - ZET ** (-x))
            P = _mul_series(P, [a ** (m + 1) for m in range(N_ + 1)])
        return _integrate(P, N_)
    @lru_cache(maxsize=None)
    def B(v):
        """I(v; 1/2->1), letters -1,0..5."""
        if not v: return mp.mpc(1)
        sw = []
        for x in reversed(v):
            if x == -1: sw.append(0)      # omega_{-1}(1-s) = -ds/(1-s): coefficient folded below
            elif x == 0: sw.append(-1)    # omega_0(1-s) = -ds/s
            else: sw.append(x)
        # coefficients: -1 -> -ds/(1-s): sign -1; 0 -> -ds/s: sign -1; j -> a_j/(1-a_j s): sign +1
        sign = (-1) ** len(v)
        for x in v:
            if x in (-1, 0): sign *= -1
        P = _B_series(tuple(sw))
        return sign * _eval_logs(P, mp.mpf(1) / 2, N_)
    return A, B, _A_series, _B_series

_EVALS = {}
def cmzv_value(w, dps=None):
    """regularized I(w; 0->1)."""
    if dps is None: dps = mp.mp.dps
    if mp.mp.dps != dps: set_dps(dps)
    N_ = N[0]
    if N_ not in _EVALS: _EVALS[N_] = make_eval(N_)
    A, B, _, _ = _EVALS[N_]
    w = tuple(w)
    tot = mp.mpc(0)
    for k in range(len(w) + 1):
        tot += A(w[:k]) * B(w[k:])
    return tot
'''
open('numeric_cmzv.py','w').write(code)
print('written')
# quick test at 50 dps
import sys; sys.dont_write_bytecode=True; sys.path.insert(0,'.')
import numeric_cmzv as NC
import mpmath as mp
NC.set_dps(50)
print('(3) =', NC.cmzv_value((3,)), ' expect -log2 =', -mp.log(2))
print('(1) =', NC.cmzv_value((1,)), ' expect i*pi/3 =', mp.pi*mp.j/3)
print('(-1,0) =', NC.cmzv_value((-1,0)), ' expect zeta(2) =', mp.zeta(2))
print('(-1,-1,0) =', NC.cmzv_value((-1,-1,0)), ' expect zeta(3) =', mp.zeta(3))
