"""V1. Mellin transform of k_c = (1-c) - h_c: quadrature vs four equivalent closed forms.

  eps*M_theta(eps) = 2 sin^2(th/2) (pi eps/sin pi eps) 2F1(-eps,1+eps;2;sin^2(th/2))       [R_hyp]
                   = pi eps sin(th) P^{-1}_eps(cos th)/sin(pi eps)                       [R_ferrers]
                   = pi eps (P_{eps-1}(c) - P_{eps+1}(c))/((2eps+1) sin pi eps)           [R_legdiff]
                   = pi sin^2(th) P'_eps(c)/((1+eps) sin pi eps)                          [R_deriv]
Quadrature: M = int_0^1 [k_c(s) s^(eps-1) + k_c(1/s) s^(-eps-1)] ds, cancellation-free k_c.
Also: theta = pi/2 vs Ahmed's R(eps), and the ODE (1-c^2) M'' + eps(eps+1) M = 0 in c.
"""
import mpmath as mp
from common import *

setdps(60)
angles = [('pi/6', mp.pi/6), ('pi/3', mp.pi/3), ('pi/2', mp.pi/2), ('2pi/3', 2*mp.pi/3),
          ('3pi/4', 3*mp.pi/4), ('5pi/6', 5*mp.pi/6), ('1', mp.mpf(1)), ('2.9', mp.mpf('2.9'))]
epss = [mp.mpf('0.05'), mp.mpf('0.37'), mp.mpf('0.5'), mp.mpf('0.8'), mp.mpf('0.97')]


def M_quad_c(eps, c):
    """int_0^1 k s^(eps-1) ds = (1/eps) int_0^1 k(t^(1/eps)) dt ;
       int_1^oo k s^(eps-1) ds = int_0^1 k(1/t) t^(-eps-1) dt, t = tau^(1/(1-eps)) -> (1/(1-eps)) int_0^1 k(1/t)/t dtau."""
    f1 = lambda t: kc(t**(1/eps), c)
    def f2(tau):
        if tau == 0:
            return (1 - c*c)/2
        t = tau**(1/(1 - eps))
        return kc(1/t, c)/t
    pts = [0, mp.mpf(1)/16, mp.mpf(1)/4, mp.mpf(1)/2, 1]
    return mp.quad(f1, pts)/eps + mp.quad(f2, pts)/(1 - eps)


def M_quad(eps, th):
    return M_quad_c(eps, mp.cos(th))


worst = 0
print("V1: eps*M_theta(eps) by quadrature vs closed forms (dps=60); residuals are absolute")
print(f"{'theta':>6} {'eps':>5} {'eps*M (quadrature)':>44} {'|hyp|':>9} {'|ferr|':>9} {'|ldiff|':>9} {'|deriv|':>9}")
for name, th in angles:
    for e in epss:
        q = e*M_quad(e, th)
        r1, r2, r3, r4 = R_hyp(e, th), R_ferrers(e, th), R_legdiff(e, th), R_deriv(e, th)
        d = [abs(q - r) for r in (r1, r2, r3, r4)]
        worst = max(worst, max(d))
        print(f"{name:>6} {mp.nstr(e,3):>5} {mp.nstr(q,40):>44} " + " ".join(f"{mp.nstr(x,2):>9}" for x in d))
print("worst residual:", mp.nstr(worst, 3))

print("\nV1b: theta = pi/2 reproduces Ahmed's R(eps) (Bailey's 2F1(a,1-a;c;1/2) summation)")
w2 = 0
for e in [mp.mpf('0.01'), mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('0.75'), mp.mpf('0.99'), mp.mpf('-0.6'), mp.mpf('1.5')]:
    a, b = R_hyp(e, mp.pi/2), R_ahmed(e)
    w2 = max(w2, abs(a - b))
    print(f"  eps={mp.nstr(e,3):>5}  R_hyp={mp.nstr(a,45)}  |R_hyp - R_Ahmed|={mp.nstr(abs(a-b),3)}")
print("  worst:", mp.nstr(w2, 3))

print("\nV1c: tail coefficients gamma_m(c=0) = binom(1/2,k) for m=2k, 0 for m odd")
for m in range(2, 12):
    g = gamma_m(m, mp.mpf(0))
    ref = mp.binomial(mp.mpf(1)/2, m//2) if m % 2 == 0 else 0
    print(f"  m={m:2d}  gamma_m={mp.nstr(g,20):>24}  ref={mp.nstr(ref,20):>24}  diff={mp.nstr(abs(g-ref),2)}")

print("\nV1d: Legendre-type ODE in c:  (1-c^2) d^2M/dc^2 + eps(eps+1) M = 0  (M by quadrature at dps=80, d/dc by mp.diff)")
for c0 in [mp.mpf('0.3'), mp.mpf('-0.4')]:
    for e in [mp.mpf('0.37'), mp.mpf('0.8')]:
        Mc = lambda c: M_quad_c(e, c)
        mp.mp.dps = 80
        val = Mc(c0); d2 = mp.diff(Mc, c0, 2)
        res = (1 - c0**2)*d2 + e*(e + 1)*val
        print(f"  c={mp.nstr(c0,2):>5} eps={mp.nstr(e,3)}  ODE residual = {mp.nstr(res,3)}   (M = {mp.nstr(val,25)})")
        mp.mp.dps = 60
