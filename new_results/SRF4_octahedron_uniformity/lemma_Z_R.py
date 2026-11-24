"""Numerical tests of Lemma Z (zero-free disc, all w in [0,1]) and Lemma R (Rouche count, w in [0,w0]).
E_w(zeta)=1+e^{w zeta}+e^{zeta};  P_w(z)=E_w(-2z).  rho*=sqrt(pi^2+log^2 2)=2R."""
import numpy as np
rho = np.hypot(np.pi, np.log(2)); R = rho/2
def zeros_E(w, box=7.0, n=40):
    a, b = np.meshgrid(np.linspace(-box, box, n), np.linspace(-box, box, n))
    z = (a+1j*b).ravel()
    for _ in range(80):
        f = 1+np.exp(w*z)+np.exp(z); fp = w*np.exp(w*z)+np.exp(z)
        z = z - f/fp
    f = np.abs(1+np.exp(w*z)+np.exp(z)); z = z[(f < 1e-12) & (np.abs(z) < box)]
    out = []
    for x in z:
        if all(abs(x-y) > 1e-7 for y in out): out.append(x)
    return np.array(out)
def count_arg(w, r, N=20000):
    phi = np.linspace(0, 2*np.pi, N, endpoint=False); z = r*np.exp(1j*phi)
    f = 1+np.exp(w*z)+np.exp(z)
    d = np.angle(np.roll(f, -1)/f); return int(round(d.sum()/(2*np.pi))), np.abs(f).min()
print("== Lemma Z: zeros of 1+e^{w zeta}+e^{zeta}; claim |zeta|>=rho*=%.10f, |Im zeta|>=pi" % rho)
ws = np.linspace(0, 1, 1001); mins = []; minim = []; bad = 0
for w in ws:
    Z = zeros_E(w); N, fmin = count_arg(w, 6.5)
    inside = Z[np.abs(Z) < 6.5]
    if len(inside) != N: bad += 1
    mins.append(np.abs(Z).min()/rho); minim.append(np.abs(Z.imag).min()/np.pi)
mins = np.array(mins)
print(" argument-principle count == Newton count (|zeta|<6.5) on all 1001 w:", bad == 0)
print(" min over w of min|zeta|/rho* = %.12f at w=%.3f" % (mins.min(), ws[mins.argmin()]))
print(" min over w of min|Im zeta|/pi = %.12f" % min(minim))
for w in [0, 0.001, 0.01, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.999, 1]:
    i = int(round(w*1000)); print("  w=%.3f  min|zeta|/rho*=%.10f" % (w, mins[i]))
print(" interior w in [0.001,0.999]: min ratio = %.8f (>1 strictly)" % mins[1:-1].min())
# Step-3 inequality: g(beta)=log^2(2cos b)+(pi+b)^2 increasing on [0, rho*-pi]
bet = np.linspace(0, rho-np.pi, 10001); g = np.log(2*np.cos(bet))**2+(np.pi+bet)**2
print(" Step 3: g increasing on [0,%.5f]: %s ; min g'=%.6f" % (rho-np.pi, bool(np.all(np.diff(g) > 0)),
      (-2*np.log(2*np.cos(bet))*np.tan(bet)+2*(np.pi+bet)).min()))
# Step 3 consequence check: for zeros with a>=0, |b| in [pi,rho*): a >= log(2cos beta)
viol = 0; tested = 0
for w in ws:
    for x in zeros_E(w):
        if x.real >= 0 and np.pi <= abs(x.imag) < rho:
            tested += 1; viol += x.real < np.log(2*np.cos(abs(x.imag)-np.pi)) - 1e-12
print(" zeros in the Step-3 band tested: %d, violations of a>=log(2cos beta): %d" % (tested, viol))
print("== Lemma R (Rouche): P_0(z)=2+e^{-2z} on |z|=rho1")
for rho1 in [1.8, 2.0, 2.2, 2.5]:
    N = 400000; phi = np.linspace(0, 2*np.pi, N, endpoint=False); z = rho1*np.exp(1j*phi)
    m_grid = np.abs(2+np.exp(-2*z)).min(); Lip = 2*rho1*np.exp(2*rho1); m_cert = m_grid - Lip*np.pi/N
    w0 = np.log(1+m_cert)/(2*rho1)
    print("  rho1=%.2f: min|P_0| (grid)=%.6f, certified lower bound=%.6f  => w0=log(1+m)/(2rho1)=%.5f" % (rho1, m_grid, m_cert, w0))
rho1 = 2.0; w0 = 0.05
print(" chosen rho1=%.2f, w0=%.3f; check e^{2 w0 rho1}-1 = %.5f" % (rho1, w0, np.exp(2*w0*rho1)-1))
cnt = set(); argmin = []; modmin = []
for w in np.linspace(0, w0, 501):
    # zeros of P_w in |z|<rho1  <=> zeros of E_w in |zeta|<2 rho1
    N, fmin = count_arg(w, 2*rho1); cnt.add(N)
    Z = zeros_E(w); Z = Z[np.abs(Z) < 2*rho1]; zs = -Z/2; zu = zs[zs.imag > 0][0]
    argmin.append(np.angle(zu)); modmin.append(abs(zu)/R)
print(" zero counts in |z|<rho1 over 501 w in [0,w0]:", cnt)
print(" z*(w): |z*|/R from %.8f to %.8f (monotone: %s); arg z* in [%.6f, %.6f] (theta=%.6f)" % (
    modmin[0], modmin[-1], bool(np.all(np.diff(modmin) > 0)), min(argmin), max(argmin), np.angle((-np.log(2)+1j*np.pi)/2)))
z0 = (-np.log(2)+1j*np.pi)/2; h = 1e-6
Z = zeros_E(h); zs = -Z/2; zu = zs[(zs.imag > 0) & (np.abs(zs) < rho1)][0]
print(" dz*/dw at 0 (numerical) = %s ; z0/2 = %s" % ((zu-z0)/h, z0/2))
