
import mpmath as mp
mp.mp.dps = 50

def S_p(phi, m, p):
    return mp.fsum(abs(mp.cos(phi - i*mp.pi/m))**p for i in range(1, m+1))
def rho_p(phi, m, p):
    return S_p(phi, m, p)**(-mp.mpf(1)/p)
def h_p_raw(theta, m, p):
    f = lambda phi: rho_p(phi, m, p)*mp.cos(theta - phi)
    N = 2000
    phis = [2*mp.pi*k/N for k in range(N)]
    vals = [f(ph) for ph in phis]
    k0 = max(range(N), key=lambda k: vals[k])
    try:
        root = mp.findroot(lambda ph: mp.diff(f, ph), phis[k0], tol=mp.mpf('1e-42'))
        return f(root), root
    except Exception:
        return vals[k0], phis[k0]

m = 3; p = mp.mpf('1.1')
th0 = mp.pi/2 - mp.pi/(2*m)
# check several fan thetas, print argmax phi* and compare with vertex direction phi_v = pi/2
for u in ['0.05','0.15','0.28']:
    u = mp.mpf(u)
    hv, phistar = h_p_raw(th0+u, m, p)
    print(f"u={mp.nstr(u,5)}: h={mp.nstr(hv,12)} argmax phi*={mp.nstr(phistar,10)} (phi_v=pi/2={mp.nstr(mp.pi/2,8)}, edge spans roughly)")
# Lagrange check: maximize <n,x> s.t. S_p(phi)=1 param by phi is exact already (polar). 
# Cross-check h via fine dense sampling + local quadratic interpolation
def h_p_dense(theta, m, p, N=200000):
    best = mp.mpf('-1'); bestphi=None
    # only need phi near theta
    import math
    f = lambda phi: rho_p(phi, m, p)*mp.cos(theta - phi)
    lo = theta - mp.pi/2; hi = theta + mp.pi/2
    for k in range(N):
        ph = lo + (hi-lo)*k/(N-1)
        v = f(ph)
        if v > best: best=v; bestphi=ph
    return best, bestphi
for u in ['0.15']:
    u = mp.mpf(u)
    hv2, ph2 = h_p_dense(th0+u, m, p)
    print("dense check u=0.15:", mp.nstr(hv2,12), "phi*=", mp.nstr(ph2,10))
# support via h+h'' consistency: check Cauchy perimeter vs polar perimeter
def peri_polar(m,p):
    f = lambda phi: rho_p(phi,m,p)
    return None
# vertices of K_p for m=3: on kink line <n_j,x>=0
print()
print("rho at vertex direction:", mp.nstr(rho_p(mp.pi/2, m, p), 15))
print("rho at edge-mid direction (th0=pi/3):", mp.nstr(rho_p(mp.pi/3, m, p), 15), " p=1:", mp.nstr(rho_p(mp.pi/3, m, 1), 15))
