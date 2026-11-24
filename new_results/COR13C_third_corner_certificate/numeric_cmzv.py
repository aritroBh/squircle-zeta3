"""High-precision evaluation of regularized level-6 CMZV word values, for validation only.

Convention: word w = (w0,...,w_{n-1}), first letter outermost:
  I(w; 0->z) = int_0^z omega_{w0}(t) I(w[1:]; 0->t) dt,
  omega_{-1} = dt/t,  omega_j = zeta^j dt/(1 - zeta^j t)  (j = 0..5),  zeta = e^{i pi/3}.
Regularization: geometric finite part (log t at the singular end).

Method: I(w;0->1) = sum_k A(w[:k]) B(w[k:]) with
  A(u) = I(u; 0->1/2)          (log-series in t at 0, log t -> log(1/2))
  B(v) = I(v; 1/2->1) = (-1)^|v| I(revmap(v); 0->1/2), s = 1-t,
       letter map: -1 -> -(ds/(1-s)),  0 -> -(ds/s),  j -> a_j ds/(1-a_j s), a_j = -zeta^j/(1-zeta^j),
       log s -> log(1/2).
Log-series representation: {j: c[0..M]} means sum_j (log s)^j * sum_m c[m+1] s^m  (m >= -1,
index = power+1).  Convergence: worst ratio 1/2, so N ~ (dps+15)/log10(2) terms.
"""
import mpmath as mp
from functools import lru_cache

ZET = None
N = [0]
_STATE = {}

def set_dps(dps):
    mp.mp.dps = dps
    _STATE.clear()
    global ZET
    ZET = mp.e ** (mp.pi * mp.j / 3)
    N[0] = int((dps + 15) / 0.30103) + 20

def _integrate(P, M):
    """Q(s) = int_0^s P(t) dt.  P: {j: array}, array index = power+1 (power -1 allowed)."""
    out = {}
    for j, c in P.items():
        for idx in range(len(c)):
            if c[idx] == 0: continue
            m = idx - 1
            if m == -1:
                out.setdefault(j + 1, [mp.mpc(0)] * (M + 2))[1] += c[idx] / (j + 1)
            else:
                f = mp.factorial
                for i in range(j + 1):
                    coef = c[idx] * (-1) ** (j - i) * f(j) / f(i) / mp.mpf(m + 1) ** (j - i + 1)
                    out.setdefault(i, [mp.mpc(0)] * (M + 2))[idx + 1] += coef
    return out

def _shift_down(P):
    return {j: c[1:] + [mp.mpc(0)] for j, c in P.items()}

def _mul_series(P, coefs, M):
    out = {}
    for j, c in P.items():
        new = [mp.mpc(0)] * (M + 2)
        for a in range(M + 1):
            if c[a] == 0: continue
            ca = c[a]
            lim = min(M + 1 - a, len(coefs))
            for b in range(lim):
                new[a + b] += ca * coefs[b]
        out[j] = new
    return out

def _eval_logs(P, s):
    tot = mp.mpc(0)
    ls = mp.log(s)
    for j, c in P.items():
        pw = mp.mpc(0)
        sp = s ** (-1)
        for idx in range(len(c)):
            pw += c[idx] * sp
            sp *= s
        tot += ls ** j * pw
    return tot

def _one(M):
    a = [mp.mpc(0)] * (M + 2); a[1] = mp.mpc(1); return {0: a}

def make_eval(M):
    @lru_cache(maxsize=None)
    def A_series(w):
        """log-series of I(w; 0->s), letters -1,0..5."""
        if not w: return _one(M)
        P = A_series(w[1:])
        if w[0] == -1:
            P = _shift_down(P)
        else:
            z = ZET ** w[0]
            P = _mul_series(P, [z ** (m + 1) for m in range(M + 1)], M)
        return _integrate(P, M)
    @lru_cache(maxsize=None)
    def A(w):
        if not w: return mp.mpc(1)
        return _eval_logs(A_series(w), mp.mpf(1) / 2)
    @lru_cache(maxsize=None)
    def B_series(sw):
        """s-side word: -1 = ds/s (from letter 0), 0 = ds/(1-s) (from letter -1), j = a_j ds/(1-a_j s)."""
        if not sw: return _one(M)
        P = B_series(sw[1:])
        x = sw[0]
        if x == -1:
            P = _shift_down(P)
        elif x == 0:
            P = _mul_series(P, [mp.mpc(1)] * (M + 1), M)
        else:
            a = -ZET ** x / (1 - ZET ** x)
            P = _mul_series(P, [a ** (m + 1) for m in range(M + 1)], M)
        return _integrate(P, M)
    @lru_cache(maxsize=None)
    def B(v):
        """I(v; 1/2->1), v letters -1,0..5."""
        if not v: return mp.mpc(1)
        sw = []
        sign = (-1) ** len(v)
        for x in reversed(v):
            if x == 0:
                sw.append(-1); sign *= -1   # omega_0(1-s) = -ds/s
            elif x == -1:
                sw.append(0); sign *= -1    # omega_{-1}(1-s) = -ds/(1-s)
            else:
                sw.append(x)                # omega_j(1-s) = a_j ds/(1-a_j s)
        return sign * _eval_logs(B_series(tuple(sw)), mp.mpf(1) / 2)
    return A, B

def cmzv_value(w, dps=None):
    """regularized I(w; 0->1), w letters in {-1, 0..5}."""
    if dps is None: dps = mp.mp.dps
    if 'ZET' not in globals() or ZET is None or mp.mp.dps != dps:
        set_dps(dps)
    M = N[0]
    if M not in _STATE: _STATE[M] = make_eval(M)
    A, B = _STATE[M]
    w = tuple(w)
    return sum(B(w[:k]) * A(w[k:]) for k in range(len(w) + 1))
