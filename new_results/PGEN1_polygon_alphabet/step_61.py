import mpmath as mp
mp.mp.dps = 60
def G_int(word, z=1, dps=60):
    mp.mp.dps = dps
    if len(word) == 0: return mp.mpf(1)
    a = word[-1]; sub = word[:-1]
    pts = [mp.mpf(0)]
    if 0 < a < z: pts.append(mp.mpf(a))
    pts.append(mp.mpf(z))
    f = lambda t: G_int(sub, t, dps)/(t-a)
    return mp.quad(f, pts)
# G(2,1,3;1): predicted symbol -3(v2^3) + 2(v3 v2 v2)  [recursion B]
# candidate: -1/2 log^3 2 - 2 Li3(-2) + c2 log2 pi^2 + c3 zeta(3)
g = G_int((2,1,3))
base = -mp.mpf('0.5')*mp.log(2)**3 - 2*mp.re(mp.polylog(3,-2))
print("G(2,1,3;1) =", g)
print("base       =", base)
rem = g - base
print("remainder =", rem)
# fit rem = c2*log2*pi^2 + c3*zeta(3)
import sympy as sp
M = [[mp.log(2)*mp.pi**2, mp.zeta(3)]]
# one equation two unknowns; try small rationals via pslq over [rem, log2*pi^2, zeta(3), pi^3? no]
v = [rem, mp.log(2)*mp.pi**2, mp.zeta(3)]
rel = mp.pslq(v, tol=mp.mpf('1e-45'), maxcoeff=10000)
print("pslq [rem, log2*pi^2, zeta3]:", rel)
