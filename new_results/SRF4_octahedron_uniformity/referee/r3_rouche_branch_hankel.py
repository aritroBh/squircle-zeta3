# Referee item 3: Lemma R (Rouche, rho1=2, w0=0.19), Lemma B (branch), Lemma H (Hankel transfer).
import numpy as np
import mpmath as mp
R = abs(complex(-np.log(2), np.pi))/2
z0c = complex(-np.log(2), np.pi)/2
P = lambda z, w: 1 + np.exp(-2*w*z) + np.exp(-2*z)
dP = lambda z, w: -2*w*np.exp(-2*w*z) - 2*np.exp(-2*z)

# --- Lemma R: min |P0| on |z|=2, hand-bound chain, Rouche margin
phi = np.linspace(0, 2*np.pi, 2_000_001)
zc = 2*np.exp(1j*phi)
m0 = np.abs(P(zc, 0)).min()
# Lipschitz certificate: |d/dphi P0(2e^{i phi})| = |2z e^{-2z}| <= 4 e^4
lip = 4*np.exp(4)*(phi[1]-phi[0])/2
print(f"min_|z|=2 |P0| = {m0:.6f} (grid), certified >= {m0-lip:.6f}")
c1 = np.arccos(np.log(2)/4); c2 = np.arccos(-np.log(3.5)/4)
s_min = min(np.sin(c1), np.sin(c2)); psi_lo = 4*s_min
print(f"phi-window ({c1:.5f},{c2:.5f}), sin>{s_min:.6f}, |psi|>{psi_lo:.6f} in [pi,2pi]: {np.pi<=psi_lo<=4<=2*np.pi}, cos(psi)>{np.cos(psi_lo):.6f}")
hb = np.sqrt(4 - (2*abs(np.cos(psi_lo)))**2)
print(f"hand bound sqrt(4-(2*0.79173)^2) = {hb:.6f};  e^(4*0.19)-1 = {np.exp(0.76)-1:.6f};  margin {hb-(np.exp(0.76)-1):.4f}")
# direct: max over the circle and w in [0,0.19] of |P_w-P0|/|P0|
mx = 0
for w in np.linspace(0, 0.19, 39):
    mx = max(mx, (np.abs(P(zc[::20], w)-P(zc[::20], 0))/np.abs(P(zc[::20], 0))).max())
print(f"max_(w<=0.19,|z|=2) |Pw-P0|/|P0| = {mx:.4f} (<1 needed)")
# zero count in |z|<2 for w in [0,0.19] (argument principle) and |z*|/R
def wind(w, n=40000):
    th = np.linspace(0, 2*np.pi, n+1)
    return int(round(np.sum(np.diff(np.unwrap(np.angle(P(2*np.exp(1j*th), w)))))/(2*np.pi)))
cnts = {wind(w) for w in np.linspace(0, 0.19, 381)}
print("zero counts in |z|<2 over 381 w in [0,0.19]:", cnts)
def zstar(w, guess):
    z = guess
    for _ in range(60): z = z - P(z, w)/dP(z, w)
    return z
ws = np.linspace(0, 0.19, 1901); zs = []; g = z0c
for w in ws:
    g = zstar(w, g); zs.append(g)
zs = np.array(zs)
print(f"|z*|/R: monotone={np.all(np.diff(np.abs(zs))>0)}, at 0.19: {abs(zs[-1])/R:.6f}; z*'(0) fd = {(zs[1]-zs[0])/(ws[1]):.6f} vs z0/2 = {z0c/2:.6f}")
v1 = 2/np.abs(zs)-1
print(f"v1 range [{v1.min():.4f},{v1.max():.4f}] (proof: [0.121,0.244], 1+v1<=1.25)")

# --- Lemma B: continuous root of D_w along [0,z*] from +sqrt3; compare with principal root at z*
def track_root(w, zst, n=20000):
    s = np.linspace(0, 1, n+1)
    zz = s*zst
    D = np.empty(n+1, complex)
    D[:-1] = P(zz[:-1], w)/(1-s[:-1])
    D[-1] = -zst*dP(zst, w)
    ang = np.unwrap(np.angle(D))
    return np.sqrt(abs(D[-1]))*np.exp(0.5j*ang[-1]), np.sqrt(D[-1]+0j), ang
for w, zst in [(0, zs[0]), (0.05, zs[500]), (0.1, zs[1000]), (0.19, zs[-1])]:
    tr, pr, ang = track_root(w, zst)
    print(f"w={w}: tracked sqrt D(z*) = {tr:.10f}, principal = {pr:.10f}, arg range [{ang.min():.4f},{ang.max():.4f}]")
    print(f"        Psi^h_w(z*) = {np.exp((1+w)*zst)*tr:.12f}")
print(f"sqrt2*z0^(1/2) = {np.sqrt(2)*np.sqrt(z0c):.12f}")

# --- Lemma H: exact check of the leading term against series coefficients
mp.mp.dps = 60
def series_F(w, ell, t1, K):
    # F_t(zeta) = t1^{-2 zeta/ell} * sqrt(e^{2(1+w)zeta}+e^{2zeta}+e^{2w zeta})
    w = mp.mpf(w); ell = mp.mpf(ell); t1 = mp.mpf(t1)
    p = [((2*(1+w))**n + 2**n + (2*w)**n)/mp.factorial(n) for n in range(K+1)]
    s = [mp.sqrt(p[0])] + [0]*K
    for n in range(1, K+1):
        s[n] = (p[n] - sum(s[i]*s[n-i] for i in range(1, n)))/(2*s[0])
    a = [(-2*mp.log(t1)/ell)**n/mp.factorial(n) for n in range(K+1)]
    return [sum(a[j]*s[n-j] for j in range(n+1)) for n in range(K+1)]
def mpzstar(w):
    w = mp.mpf(w)
    return mp.findroot(lambda z: 1+mp.exp(-2*w*z)+mp.exp(-2*z), mp.mpc(z0c) + mp.mpc(z0c)*w/2)
def lead(w, ell, t1, k, zst, rootD):
    Psi = mp.mpf(t1)**(-2*zst/ell)*mp.exp((1+w)*zst)*rootD
    return -2/mp.pi*mp.beta(1.5, k-0.5)*mp.re(zst**(-k)*Psi)
K = 80
print("Lemma H: (coef - lead)*|z*|^k*k^{5/2}  (bounded => error O(k^-5/2 |z*|^-k)), and lead/coef")
for (w, ell) in [(0, 40), (0.05, 3), (0.19, 1)]:
    t1 = 1/(1+mp.exp(-w*ell)+mp.exp(-ell))
    c = series_F(w, ell, t1, K)
    zst = mpzstar(w)
    i = int(round(w/0.19*1900))
    tr, _, _ = track_root(w, zs[i])
    # refine the root sign with mp: principal root, after checking tracked == principal above
    rootD = mp.sqrt(-zst*(-2*w*mp.exp(-2*w*zst)-2*mp.exp(-2*zst)))
    assert abs(complex(rootD)-tr) < 1e-6
    out = []
    for k in [10, 20, 40, 80]:
        L = lead(w, ell, t1, k, zst, rootD)
        out.append(f"k={k}: {mp.nstr((c[k]-L)*abs(zst)**k*k**2.5,5)} ratio {mp.nstr(L/c[k],8)}")
    print(f"  (w,ell)=({w},{ell}): " + "; ".join(out))
# sign check [x^k](1-x)^{1/2} = -(1/pi) B(3/2,k-1/2)
for k in [2, 5, 9]:
    print(f"  k={k}: [x^k](1-x)^(1/2)={mp.nstr(mp.taylor(lambda x: mp.sqrt(1-x), 0, k)[k],12)}  -(1/pi)B={mp.nstr(-mp.beta(1.5,k-0.5)/mp.pi,12)}")
