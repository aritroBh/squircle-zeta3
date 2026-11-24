"""Uniform constants of Lemma H on grids (values reported; the proof only uses their finiteness, by compactness),
Lemma U inequality on a t-grid, two-consecutive-k lemma, and the measure / profile integrals."""
import numpy as np
from mpmath import mp, mpf, quad, exp, pi, inf, log, sqrt
z0 = (-np.log(2)+1j*np.pi)/2; R = abs(z0); th = np.angle(z0); rho1 = 2.0; w0 = 0.19
P = lambda w, z: 1+np.exp(-2*w*z)+np.exp(-2*z)
dP = lambda w, z: -2*w*np.exp(-2*w*z)-2*np.exp(-2*z)
def zstar(w):
    z = z0*(1+w/2)
    for _ in range(60): z = z-P(w, z)/dP(w, z)
    return z
MPsi = 0; Lip = 0; Dmin = 1e9; ReMax = -1e9; Fcirc = 0
for w in np.linspace(0, w0, 191):
    zs = zstar(w)
    # branch: track sqrt D along [0,z*]
    s = np.linspace(0, 1, 4001)[:-1]; Ds = P(w, s*zs)/(1-s); Ds = np.append(Ds, -zs*dP(w, zs))
    r = np.sqrt(Ds.astype(complex))
    for i in range(1, len(r)):
        if abs(r[i]-r[i-1]) > abs(-r[i]-r[i-1]): r[i] = -r[i]
    psi_star = np.exp((1+w)*zs)*r[-1]
    v1 = rho1/abs(zs)-1; v = np.linspace(0, v1, 4001)[1:]; zz = zs*(1+v)
    Dv = P(w, zz)/(1-zz/zs); rv = np.sqrt(Dv.astype(complex)); prev = r[-1]
    for i in range(len(rv)):
        if abs(rv[i]-prev) > abs(-rv[i]-prev): rv[i] = -rv[i]
        prev = rv[i]
    psi = np.exp((1+w)*zz)*rv
    MPsi = max(MPsi, np.abs(psi).max(), abs(psi_star)); Dmin = min(Dmin, np.abs(Dv).min(), abs(Ds[-1]))
    Lip = max(Lip, (np.abs(psi-psi_star)/(abs(zs)*v)).max()); ReMax = max(ReMax, zz.real.max(), zs.real)
    ph = np.linspace(0, 2*np.pi, 20001); zc = rho1*np.exp(1j*ph)
    Fcirc = max(Fcirc, np.abs(np.sqrt(np.abs(np.exp(2*(1+w)*zc)+np.exp(2*zc)+np.exp(2*w*zc)))).max())
print("Lemma H constants over w in [0,0.19] (cut sets {z*(w)(1+v): 0<=v<=2/|z*|-1}):")
print("  sup|Psi^h| = %.5f ; sup |Psi^h(z*(1+v))-Psi^h(z*)|/(|z*| v) = %.5f ; min|D_w| on cuts = %.5f ; max Re on cuts = %.5f (<0)" % (MPsi, Lip, Dmin, ReMax))
print("  sup_{|z|=2}|h_w| = %.4f  (bound used: sqrt3 e^{2 rho1(1+w0)}... crude = %.1f)" % (Fcirc, np.sqrt(3)*np.exp(2*rho1*(1+w0)/1)))
# Lemma U: |[tau^k] f(.;t)| <= M (max(l,l0)/(R-eps))^k with M=9^{R/l0} sqrt3 e^{2R}, eps=0.05, l0=1
from mpmath import mpf as F
mp.dps = 30
def coeffs_f(t, K):
    # f = (t1t2t3)^-tau (sum t^{2tau})^{1/2}
    Lg = [log(F(x)) for x in t]
    Q = [sum((2*L)**n for L in Lg)/mp.factorial(n) for n in range(K+1)]
    q = [sqrt(Q[0])]+[F(0)]*K
    for n in range(1, K+1): q[n] = (Q[n]-sum(q[j]*q[n-j] for j in range(1, n)))/(2*q[0])
    S = -sum(Lg); e = [S**n/mp.factorial(n) for n in range(K+1)]
    return [sum(e[j]*q[n-j] for j in range(n+1)) for n in range(K+1)]
eps, l0 = 0.05, 1.0; M = 9**(R/l0)*np.sqrt(3)*np.exp(2*R); worst = 0
rng = np.random.default_rng(1)
pts = [(1/3, 1/3, 1/3), (0.5, 0.5-1e-6, 1e-6), (1-2e-9, 1e-9, 1e-9), (0.6, 0.3, 0.1)]
for _ in range(60):
    x = rng.exponential(size=3)*np.exp(rng.uniform(-12, 0, 3)); pts.append(tuple(x/x.sum()))
for t in pts:
    c = coeffs_f(t, 40); l = np.log(max(t)/min(t))
    for k in range(41):
        worst = max(worst, float(abs(c[k]))/(M*(max(l, l0)/(R-eps))**k))
print("Lemma U: max over 64 t-points, k<=40 of |[tau^k]f| / (M (max(l,l0)/(R-eps))^k) = %.3e (must be <=1)" % worst)
# two-consecutive lemma
x = np.linspace(0, np.pi, 2000001); c = np.maximum(np.abs(np.cos(x)), np.abs(np.cos(x+th)))
print("two-consecutive lemma: min_x max(|cos x|,|cos(x+theta)|) = %.6f >= sin((pi-theta)/2) = %.6f" % (c.min(), np.sin((np.pi-th)/2)))
k = np.arange(1, 200001); cc = np.abs(np.cos((k-0.5)*th)); mm = np.maximum(cc[:-1], cc[1:])
print("  over k<=2e5: min max(|cos((k-1/2)th)|,|cos((k+1/2)th)|) = %.6f" % mm.min())
mp.dps = 30
m = quad(lambda lam: quad(lambda l: exp(-lam-l)/(1+exp(-lam)+exp(-l))**3, [lam, inf]), [0, inf])
print("measure check 6*int t1^3 e^{-lam-l} over 0<=lam<=l = %s (area of simplex 1/2)" % (6*m))
print("profile integral int_0^inf e^{-3lam/2}(1+e^{-lam})^-3 = %s ; pi/16 = %s" % (quad(lambda u: exp(-1.5*u)/(1+exp(-u))**3, [0, inf]), pi/16))
