"""step_05.py -- single IBP step check at k=2, key (2,0) and (1,1)."""
import sys
sys.path.insert(0, '.')
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 40
from cwfield import make_quad
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf

c = Fr(1,2); r = Fr(1,3)
F = make_quad(-r)
red = Reducer(F, F(r))
ci = CornerIntegrand(F, F(c), F(r))
rmp = mp.mpf(1)/3
th = mp.mpc(0, 1)/mp.sqrt(3)
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
def lt(w): return mp.log((w-rmp)*(w+1)/(2*(1-rmp)*w))
def lmt(w): return mp.log((1-w)*(w+rmp)/(2*(1-rmp)*w))
om1 = lambda w: 1/(w-rmp) + 1/(w+1) - 1/w
om2 = lambda w: 1/(w-1) + 1/(w+rmp) - 1/w
pts = [rmp, mp.sqrt(rmp), 1]

ipf = build_integrand_pf(ci, red, 2)
for key, W in ipf.items():
    i, j = key
    simple, higher, const = red.pf_split(W)
    G = red.pf_integrate_higher(higher)
    C = red.pf_eval_at(G, red.r)
    Gm = [G[0], G[1], G[2], -C]
    fW = kfun(W); fS = kfun(simple); fG = kfun(Gm)
    LHS = mp.quad(lambda w: fW(w)*lt(w)**i*lmt(w)**j, pts)
    RHS = mp.quad(lambda w: fS(w)*lt(w)**i*lmt(w)**j, pts)
    feed1 = mp.quad(lambda w: fG(w)*om1(w)*lt(w)**(i-1)*lmt(w)**j, pts) if i >= 1 else 0
    feed2 = mp.quad(lambda w: fG(w)*om2(w)*lt(w)**i*lmt(w)**(j-1), pts) if j >= 1 else 0
    RHS = RHS - i*feed1 - j*feed2
    print(f"key {key}: LHS={mp.nstr(LHS,20)}")
    print(f"   RHS={mp.nstr(RHS,20)}  diff={mp.nstr(abs(LHS-RHS),3)}  (C={mp.nstr(q2mp(C),6)})")
print("STEP05 DONE")
