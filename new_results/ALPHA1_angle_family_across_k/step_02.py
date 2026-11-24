"""step_02.py -- exact corner reducer at c = 1/2 (r = 1/3), k = 1..4;
weight decomposition checked numerically against corner_series.py (>= 30 digits)."""
import sys, time, pickle
sys.path.insert(0, '.')
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 50
from cwfield import make_quad
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf
from corner_series import corner_coeffs

c = Fr(1,2); r = Fr(1,3)
F = make_quad(-r)
red = Reducer(F, F(r))
ci = CornerIntegrand(F, F(c), F(r))

KMAX = 4
kernels = {}
for k in range(1, KMAX+1):
    t0 = time.time()
    ipf = build_integrand_pf(ci, red, k)
    fin = red.reduce_k(ipf)
    kernels[k] = fin
    simple_th = {kk: (P[1][0] if P[1] else None) for kk, P in fin.items()}
    print(f"k={k}: keys={sorted(fin)} ({time.time()-t0:.1f}s); simple theta-pole coeffs: "
          + str({kk: (str(v)) for kk, v in simple_th.items() if v and not v.is_zero()}))
from cwbuild import pf_plain
pickle.dump({k: {kk: pf_plain(P) for kk, P in fin.items()} for k, fin in kernels.items()}, open('kern_c12.pkl','wb'))

# ---- numeric weight parts ----
rmp = mp.mpf(1)/3
th = mp.mpc(0, 1)/mp.sqrt(3)
def lt(w): return mp.log((w-rmp)*(w+1)/(2*(1-rmp)*w))
def lmt(w): return mp.log((1-w)*(w+rmp)/(2*(1-rmp)*w))
def q2mp(q):
    return mp.mpf(int(q.a.numerator))/int(q.a.denominator) + (mp.mpf(int(q.b.numerator))/int(q.b.denominator))*th

def kfun(P):
    def f(w):
        tot = mp.mpc(0)
        for rho, cf in P[0].items():
            rr = mp.mpf(int(rho.a.numerator))/int(rho.a.denominator)
            tot += q2mp(cf)/(w - rr)
        for j, e in enumerate(P[1], start=1):
            tot += q2mp(e)/(w - th)**j
        for j, e in enumerate(P[2], start=1):
            tot += q2mp(e)/(w + th)**j
        tot += q2mp(P[3])
        return tot
    return f

def weight_part(k, d):
    fin = kernels[k]
    tot = mp.mpc(0)
    for (i, j), P in fin.items():
        if i + j != d: continue
        f = kfun(P)
        g = lambda w: f(w) * lt(w)**i * lmt(w)**j
        tot += mp.quad(g, [rmp, mp.sqrt(rmp), 1])
    return tot

Mk = corner_coeffs(mp.acos(mp.mpf(1)/2)*2, KMAX)   # alpha with cos(alpha/2) = 1/2 -> alpha = 2pi/3
print("corner M_k:", [mp.nstr(x, 12) for x in Mk])
for k in range(1, KMAX+1):
    parts = {}
    tot = mp.mpc(0)
    for d in range(1, k+1):
        W = weight_part(k, d)
        parts[d] = W
        tot += W
    print(f"k={k}: M_k from kernels = {mp.nstr(tot, 30)}")
    print(f"      corner_series   = {mp.nstr(Mk[k], 30)}   diff = {mp.nstr(abs(tot - Mk[k]), 3)}")
    for d, W in parts.items():
        print(f"      W_{d+1}(k) = {mp.nstr(W, 22)}")
print("STEP02 DONE")
