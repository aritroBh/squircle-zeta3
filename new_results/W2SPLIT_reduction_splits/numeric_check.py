"""W2SPLIT numeric_check.py (NUMERICAL safeguard, not part of any proof).
(a) the kernels Rtilde of split_check.py against a direct Taylor expansion of -1/2 Phi_eps(t(v)) dt/dv in eps;
(b) the split: sum over keys of int_b^a S_rs l1^r l2^s dv, with S_rs = W_rs - Hp(W_rs) (simple poles at real letters only),
    against M_k = -1/2 int_0^1 [eps^k] Phi_eps dt computed directly (mpmath quadrature, Taylor in eps by mpmath.taylor).
"""
import sys
import mpmath as mp
from split_check import *
mp.mp.dps = 30
s2 = mp.sqrt(2); a = 1 + s2; b = s2 - 1
def ev(W, x):
    f = lambda p: mp.polyval([mp.mpf(int(c.p)) / int(c.q) for c in reversed(p.coeffs())], x) if not p.is_zero() else mp.mpf(0)
    return (f(W.n0) + s2 * f(W.n1)) / f(W.den)
def tof(v): return (v * v + 2 * v - 1) / (4 * v)
def Phi(e, t):
    al = t ** (e - 1); be = (1 - t) ** (e - 1)
    return al + be - mp.sqrt(al * al + be * be)
def phik(k, t):
    return mp.taylor(lambda e: Phi(e, t), 0, k)[k]
kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 5
for k in range(2, kmax + 1):
    v0 = mp.mpf('1.3'); t0 = tof(v0); l1 = mp.log(t0); l2 = mp.log(1 - t0)
    direct = -mp.mpf(1) / 2 * phik(k, t0) * (v0 * v0 + 1) / (4 * v0 * v0)
    ours = sum(ev(Rtilde(r, k - r), v0) * l1 ** r * l2 ** (k - r) for r in range(k + 1))
    print(f"k={k} (a) integrand at v=1.3: |direct - sum Rtilde l1^r l2^s| = {mp.nstr(abs(direct - ours), 3)}", flush=True)
    _, _, _, W = run_level(k, lambda s: None)
    pieces = {}
    for key, w in W.items():
        _, _, S = Hp_parts(w)
        pieces[key] = S
    def integrand(v):
        t = tof(v); L1 = mp.log(t); L2 = mp.log(1 - t)
        return sum(ev(S, v) * L1 ** r * L2 ** s for (r, s), S in pieces.items())
    split = mp.quad(integrand, [b, 1, a])
    Mk = -mp.mpf(1) / 2 * mp.quad(lambda t: phik(k, t), [0, mp.mpf(1) / 2, 1])
    print(f"k={k} (b) M_k direct = {mp.nstr(Mk, 20)}; sum of split pieces = {mp.nstr(split, 20)}; diff = {mp.nstr(abs(Mk - split), 3)}", flush=True)
