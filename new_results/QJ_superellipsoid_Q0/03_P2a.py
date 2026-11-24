# 03_P2a.py -- P2a = int_0^1 arsinh(q)/sqrt(1+q^2) ds, q = sqrt(s^2+(1-s)^2).
# Real route only: q = sinh(phi), integration by parts, Fubini in an auxiliary parameter k,
# rationalising substitution tau, then two real-dilogarithm primitives (exact derivative checks).
# Output: P2a as an exact combination of real Li2 values at S-units (reduced in 05_reduce.py).
import sympy as sp, mpmath as mp, json
from common import lam, l3, L, al, PI, Li2, li2, logS, num, show, Q2, A, B, SQ2
mp.mp.dps = 90
def ok(e):
    if sp.simplify(e) == 0: return True
    return sp.radsimp(sp.simplify(sp.expand(sp.sympify(e).rewrite(sp.exp)))) == 0
ref = json.load(open("reference_100.json"))
q, ph, k, tau, w = sp.symbols('q phi k tau w', positive=True)
r2, r3 = sp.sqrt(2), sp.sqrt(3)
c = r2/r3
a_, b_ = 1 + r2, r2 - 1
h = 1/r2

print("== Step 1: q = sinh(phi)")
# P2a = 2 int_{1/sqrt2}^1 arsinh(q)/sqrt(1+q^2) * q/sqrt(2q^2-1) dq   (s -> q as in 02_P2c.py)
integ_q = sp.asinh(q)/sp.sqrt(1 + q**2)*q/sp.sqrt(2*q**2 - 1)
# with q = sinh(phi), phi > 0: arsinh q = phi, sqrt(1+q^2) = cosh(phi) (> 0), dq = cosh(phi) dphi
assert ok(1 + sp.sinh(ph)**2 - sp.cosh(ph)**2)
integ_ph = (ph/sp.cosh(ph))*sp.sinh(ph)/sp.sqrt(2*sp.sinh(ph)**2 - 1)*sp.cosh(ph)
assert ok(integ_ph - ph*sp.sinh(ph)/sp.sqrt(2*sp.sinh(ph)**2 - 1))
phi0 = sp.asinh(1/r2)
assert ok(sp.sinh(sp.log(a_)) - 1)
print("  2 q dq arsinh q/(sqrt(1+q^2) sqrt(2q^2-1)) = 2 phi sinh(phi) dphi/sqrt(2 sinh^2 phi - 1),  phi in [phi0, L],")
print("  phi0 = arsinh(1/sqrt2), L = arsinh(1) = log(1+sqrt2)   (exact)")

print("== Step 2: integration by parts with G(phi) = arcosh(c cosh phi), c = sqrt(2/3)")
G = sp.acosh(c*sp.cosh(ph))
Gp = sp.diff(G, ph)
assert ok(Gp**2 - (r2*sp.sinh(ph)/sp.sqrt(2*sp.sinh(ph)**2 - 1))**2)   # both sides positive for phi > phi0
assert ok(c**2*sp.cosh(ph)**2 - 1 - (2*sp.sinh(ph)**2 - 1)/3)
assert ok(c*sp.cosh(phi0) - 1)                       # G(phi0) = arcosh 1 = 0
assert ok(c*sp.sqrt(2) - 2/r3)                       # cosh L = sqrt2, so c cosh L = 2/sqrt3
assert ok(sp.acosh(2/r3).rewrite(sp.log) - sp.log(r3))   # arcosh(2/sqrt3) = log sqrt3
print("  G'(phi) = sqrt2 sinh phi/sqrt(2sinh^2phi-1) (>0), G(phi0)=0, G(L)=log(3)/2   (exact)")
print("  P2a = sqrt2 int phi G'(phi) dphi = sqrt2 [ L log(3)/2 - W ],  W = int_{phi0}^{L} G(phi) dphi")
Wnum = mp.quad(lambda p: mp.acosh(mp.sqrt(mp.mpf(2)/3)*mp.cosh(p)), [mp.asinh(1/mp.sqrt(2)), mp.log(1 + mp.sqrt(2))])
show("  P2a = sqrt2(L log3/2 - W) (W by quadrature) vs reference",
     mp.sqrt(2)*(mp.log(1 + mp.sqrt(2))*mp.log(3)/2 - Wnum), mp.mpf(ref["P2a"]), 80)

print("== Step 3: W by Fubini")
# arcosh(c cosh phi) = int_{1/cosh phi}^{c} cosh phi dk / sqrt(k^2 cosh^2 phi - 1)
assert ok(sp.diff(sp.acosh(k*sp.cosh(ph)), k)**2 - sp.cosh(ph)**2/(k**2*sp.cosh(ph)**2 - 1))
# inner integral over phi in [arcosh(1/k), L] for k in [1/sqrt2, c]:  (1/k) arcosh(sinh(phi)/mu), mu = sqrt(1-k^2)/k
mu = sp.sqrt(1 - k**2)/k
Hk = sp.acosh(sp.sinh(ph)/mu)/k
assert ok(sp.diff(Hk, ph)**2 - sp.cosh(ph)**2/(k**2*sp.cosh(ph)**2 - 1))
print("  d/dphi[(1/k) arcosh(sinh phi/mu)] = cosh phi/sqrt(k^2cosh^2phi-1), mu = sqrt(1-k^2)/k; = 0 at sinh phi = mu  (exact)")
print("  => W = int_{1/sqrt2}^{c} (1/k) arcosh(k/sqrt(1-k^2)) dk   (sinh L = 1)")
Wk = mp.quad(lambda kk: mp.acosh(kk/mp.sqrt(1 - kk*kk))/kk, [1/mp.sqrt(2), mp.sqrt(mp.mpf(2)/3)])
show("  W (k-form) vs W (phi-form)", Wk, Wnum, 80)

print("== Step 4: tau = sqrt(2k^2-1)/k, k = 1/sqrt(2-tau^2)")
k_t = 1/sp.sqrt(2 - tau**2)
assert ok((2*k_t**2 - 1)/k_t**2 - tau**2)               # both sides >= 0 on [0, 1/sqrt2]
assert ok(k_t**2/(1 - k_t**2) - 1/(1 - tau**2))          # k/sqrt(1-k^2) = 1/sqrt(1-tau^2) = cosh(artanh tau)
assert ok(sp.diff(k_t**2, tau)/(2*k_t**2) - tau/(2 - tau**2))
assert ok(k_t.subs(tau, 0) - 1/r2) and ok(k_t.subs(tau, h) - c)
print("  k/sqrt(1-k^2) = 1/sqrt(1-tau^2) = cosh(artanh tau);  dk/k = tau dtau/(2-tau^2);  tau: 0 -> 1/sqrt2   (exact)")
print("  => W = int_0^{1/sqrt2} artanh(tau) tau/(2-tau^2) dtau")
Wt = mp.quad(lambda t: mp.atanh(t)*t/(2 - t*t), [0, 1/mp.sqrt(2)])
show("  W (tau-form) vs W (phi-form)", Wt, Wnum, 80)
K = 1/(r2 - tau) - 1/(r2 + tau)
assert ok(K/2 - tau/(2 - tau**2))
print("  artanh tau = (log(1+tau)-log(1-tau))/2, tau/(2-tau^2) = K(tau)/2 with K odd ==>")
print("  W = (1/4) int_{-1/sqrt2}^{1/sqrt2} log(1+tau) [1/(sqrt2-tau) - 1/(sqrt2+tau)] dtau")
Ws = mp.quad(lambda t: mp.log(1 + t)*(1/(mp.sqrt(2) - t) - 1/(mp.sqrt(2) + t)), [-1/mp.sqrt(2), 0, 1/mp.sqrt(2)])/4
show("  W (symmetric form) vs W (phi-form)", Ws, Wnum, 80)

print("== Step 5: primitives in w = 1 + tau in [b/sqrt2, a/sqrt2];  sqrt2 - tau = a - w,  sqrt2 + tau = w + b")
assert ok((1 - h) - b_/r2) and ok((1 + h) - a_/r2) and ok(r2 - (a_ - 1)) and ok(r2 - (1 + b_))
Fm = -sp.log(w)*sp.log(1 - w/a_) - sp.polylog(2, w/a_)
Fp = sp.log(w)*sp.log(1 + w/b_) + sp.polylog(2, -w/b_)
dFm = sp.expand_func(sp.diff(Fm, w)); dFp = sp.expand_func(sp.diff(Fp, w))
assert ok(dFm - sp.log(w)/(a_ - w)) and ok(dFp - sp.log(w)/(w + b_))
print("  F-(w) = -log w log(1-w/a) - Li2(w/a):   F-'(w) - log w/(a-w) = 0   (exact)")
print("  F+(w) =  log w log(1+w/b) + Li2(-w/b):  F+'(w) - log w/(w+b) = 0   (exact)")
print("  on [b/sqrt2, a/sqrt2]: w>0, 0 < w/a < 1, -w/b < 0: all logs of positive numbers, Li2 real-analytic, F+- continuous")

# exact boundary values; Li2 arguments and log arguments are S-units of Q(sqrt2)
def Fm_val(wq):
    return -logS(wq)*logS(1 - wq/A) - li2(wq/A)
def Fp_val(wq):
    return logS(wq)*logS(1 + wq/B) + li2(-(wq/B))
wlo, whi = B/SQ2, A/SQ2
W = sp.Rational(1, 4)*((Fm_val(whi) - Fm_val(wlo)) - (Fp_val(whi) - Fp_val(wlo)))
W = sp.expand(W)
print("  Li2 arguments:", sorted(set(str(f.args[0]) for f in W.atoms(sp.Function))))
print("  W =", W)
show("  W exact (unreduced) vs W (phi-form)", num(W), Wnum, 80)
P2a = sp.expand(r2*(L*l3/2 - W))
print("== Result (unreduced): P2a =", P2a)
show("P2a exact (unreduced) vs 100-digit reference", num(P2a), mp.mpf(ref["P2a"]), 90)
json.dump({"P2a": sp.srepr(P2a)}, open("P2a_unreduced.json", "w"))
