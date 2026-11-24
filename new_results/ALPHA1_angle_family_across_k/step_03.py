"""step_03.py -- isolate bug: does the UNREDUCED w-integrand integrate to M_k?"""
import sys
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

Mk = corner_coeffs(mp.pi*2/3, 4)
for k in range(1, 5):
    ipf = build_integrand_pf(ci, red, k)
    tot = mp.mpc(0)
    for (i, j), P in ipf.items():
        f = kfun(P)
        tot += mp.quad(lambda w: f(w)*lt(w)**i*lmt(w)**j, [rmp, mp.sqrt(rmp), 1])
    print(f"k={k}: raw integrand = {mp.nstr(tot,25)}  corner = {mp.nstr(Mk[k],25)}  diff = {mp.nstr(abs(tot-Mk[k]),3)}")
print("STEP03 DONE")
