"""PUR2 (C),(E): admissible kappa range, interval certificate for kappa^2 = 9/7, thresholds."""
import mpmath as mp
import sys
sys.path.insert(0, '.')
mp.mp.dps = 30
from pur2_numerics_lib import Curve, PI4

def out(*a):
    print(*a); sys.stdout.flush()

# admissible range: kappa e^{pi/(2 kappa)} = 2 + kappa^2
kmax = mp.findroot(lambda k: k * mp.exp(mp.pi / (2 * k)) - 2 - k * k, 2.5)
out("kappa_max (phi>0 on [0,V] iff kappa < kappa_max): %s ; kappa_max^2 = %s" % (mp.nstr(kmax, 20), mp.nstr(kmax**2, 20)))
out("  at kappa^2=9/7: kappa e^{pi/2kappa} - 2 - kappa^2 = %s" % mp.nstr(mp.sqrt(mp.mpf(9)/7) * mp.exp(mp.pi / (2 * mp.sqrt(mp.mpf(9)/7))) - 2 - mp.mpf(9)/7, 10))

# interval certificate for kappa^2 = 9/7:  kappa^2 <= 3 e^V cos^3(kappa V)  (then D >= 0)
iv = mp.iv
iv.dps = 30
k = iv.sqrt(iv.mpf(9) / 7)
def Phi_iv(v):
    return (iv.exp(v) * (iv.cos(k * v) + k * iv.sin(k * v)) - 1) / (1 + k * k)
C = Curve('k', kappa=mp.sqrt(mp.mpf(9) / 7))
Vlo, Vhi = C.V - mp.mpf(10)**-20, C.V + mp.mpf(10)**-20
plo, phi_ = Phi_iv(iv.mpf(Vlo)), Phi_iv(iv.mpf(Vhi))
out("\nkappa^2=9/7: V = %s ; Phi(V-1e-20) = %s (<1), Phi(V+1e-20) = %s (>1)" % (mp.nstr(C.V, 25), plo, phi_))
assert plo.b < 1 and phi_.a > 1
Vint = iv.mpf([Vlo, Vhi])
rhs = 3 * iv.exp(Vint) * iv.cos(k * Vint)**3
out("  interval: 3 e^V cos^3(kappa V) in %s ; kappa^2 = 9/7 = %s ; certified kappa^2 < rhs: %s ; kappa^2 < 3: True" % (rhs, mp.nstr(mp.mpf(9)/7, 10), (rhs.a > mp.mpf(9)/7)))
out("  kappa V = %s < pi/2 ; phi(V) = %s" % (mp.nstr(C.k * C.V, 15), mp.nstr(C.phi(C.V), 15)))

def sufficient(kap):
    Cc = Curve('k', kappa=kap)
    return min(3, 3 * mp.exp(Cc.V) * mp.cos(kap * Cc.V)**3) - kap**2

def minD(Cc, Ng=300):
    return min(Cc.D(PI4 * i / Ng)[1] for i in range(1, Ng + 1))

# kappa_s: largest kappa where the sufficient condition holds
ks = mp.findroot(sufficient, 1.2)
out("\nkappa_s (sufficient condition kappa^2 <= min(3, 3e^V cos^3 kappa V) holds for kappa <= kappa_s): %s, kappa_s^2 = %s" % (mp.nstr(ks, 15), mp.nstr(ks**2, 15)))
# kappa_D: numerically, largest kappa with min D >= 0 (bisection on grid min of D phi/M^2)
lo, hi = mp.mpf('0.5'), kmax - mp.mpf('1e-6')
for _ in range(40):
    mid = (lo + hi) / 2
    if minD(Curve('k', kappa=mid)) >= -mp.mpf(10)**-25:
        lo = mid
    else:
        hi = mid
out("kappa_D (numerical, grid 300): min D >= 0 iff kappa <= %s ; kappa_D^2 = %s" % (mp.nstr(lo, 12), mp.nstr(lo**2, 12)))

# beyond kappa_D: sharp tau* (numerical) and the proved tau_0 = 1/(1 + kappa tan(kappa V))
out("\nkappa beyond kappa_D: proved convex for tau <= tau0 = 1/(1+kappa tan(kappa V)); numerical sharp tau*")
for kap in [mp.mpf('1.4'), mp.mpf('1.6'), mp.mpf('1.8'), mp.mpf('1.95'), mp.mpf('2.05')]:
    if kap >= kmax:
        continue
    Cc = Curve('k', kappa=kap)
    tau0 = 1 / (1 + kap * mp.tan(kap * Cc.V))
    Ng = 300
    grid = [PI4 * i / Ng for i in range(1, Ng + 1)]
    def ok(tau):
        return min(Cc.hpp(u, tau) for u in grid) >= 0
    if ok(mp.mpf(1)):
        tstar = mp.mpf(1)
    else:
        a_, b_ = mp.mpf('1e-3'), mp.mpf(1)
        for _ in range(30):
            m_ = (a_ + b_) / 2
            if ok(m_):
                a_ = m_
            else:
                b_ = m_
        tstar = a_
    out("  kappa=%s kappa^2=%s: V=%s minD*phi/M^2=%s tau0(proved)=%s tau*(numerical)=%s" % (
        mp.nstr(kap, 4), mp.nstr(kap**2, 5), mp.nstr(Cc.V, 8), mp.nstr(minD(Cc), 6), mp.nstr(tau0, 6), mp.nstr(tstar, 6)))

out("\nreal c: grid min of D*phi/M^2 (D >= 0 proved for |c| < 6)")
for cc in ['0', '0.5', '1', '2', '3', '6', '10', '30', '100']:
    Cc = Curve('c', c=mp.mpf(cc))
    out("  c=%-4s V=%s min D*phi/M^2 = %s" % (cc, mp.nstr(Cc.V, 8), mp.nstr(minD(Cc), 8)))

out("\ncheck of D*phi = 6p + (2-w^2)chi - 3w int_0^p chi dp' (c=1/2, kappa^2=9/7) at u=0.5:")
for Cc in [Curve('c', c=mp.mpf(1)/2), Curve('k', kappa=mp.sqrt(mp.mpf(9)/7))]:
    u = mp.mpf('0.5'); p = 1 + mp.sin(u) - mp.cos(u); w = 1 - p
    L = Cc.Lam(u)
    chi = (lambda v: Cc.c * mp.tanh(Cc.c * v)) if Cc.c is not None else (lambda v: -Cc.k * mp.tan(Cc.k * v))
    integ = mp.quad(lambda v: chi(v) * Cc.phi(v), [0, L])   # int chi dp = int chi phi dv
    lhs = Cc.D(u)[0] * Cc.phi(L)
    out("  %s: D phi = %s ; formula = %s" % (Cc.name, mp.nstr(lhs, 20), mp.nstr(6 * p + (2 - w * w) * chi(L) - 3 * w * integ, 20)))
