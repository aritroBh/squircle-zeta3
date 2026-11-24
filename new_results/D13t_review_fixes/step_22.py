import math, random, time
import mpmath as mp
mp.mp.dps = 60
mpr2 = mp.sqrt(2); mpp = mp.pi; ma = mpp/(2*mpr2); mbb = mp.e**(1j*mpp/4)

# --- (a) saddle sigma* in [0,1] across the whole range ---
def saddle(N):
    y = mp.sqrt(mpp/(2*N))
    f = lambda r: r**2*(1+mpr2*r*y)/(1+mpr2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-55'))
    return y, rho, (rho-1)/y**2
sigmin, sigmax = mp.mpf('9'), mp.mpf('-9')
for N in [100, 100.0001, 101, 150, 400, 1891, 1e4, 1e6, 1e8, 1e10, 1e14, 1e20]:
    y,rho,sg = saddle(mp.mpf(N))
    sigmin=min(sigmin,sg); sigmax=max(sigmax,sg)
print(f"(a) saddle sigma* over N in [100,1e20]: min {mp.nstr(sigmin,10)}, max {mp.nstr(sigmax,10)}  (need in [0,1])")

# --- (b) independent check: phi = (pi/rt2) Bh[A1'/y] exactly ---
def Phi(x,N): return (-1+1j)*ma*x + N*mp.log(mbb*x/(1+mbb*x))
random.seed(20240513)
worst = mp.mpf('0'); worstN=None
trials = []
# targeted N near 100, 1891 and random log-uniform
Ns = [100, 100.0001, 100.5, 1891.0, 1890.8, 400.0] + [10**random.uniform(2,10) for _ in range(2000)]
for Nv in Ns:
    N = mp.mpf(Nv)
    y = mp.sqrt(mpp/(2*N)); th = mpp/(mpr2*y)
    y0,rho,sg = saddle(N)
    x0 = 1/(rho*y); t = rho*y
    phi_true = mp.im(Phi(x0,N)) - th
    zz = mpr2*t/(2+mpr2*t)
    Ap = 1 - mp.atan(zz)/zz
    A1 = Ap/y
    T2n = 2*y*y*sg - mpr2*rho*y - 2*rho*A1*y
    Bh = -y*sg/(2*rho) + T2n/(2*y*(2+mpr2*t))
    r = abs(phi_true - (mpp/mpr2)*Bh)
    if r > worst: worst = r; worstN = Nv
print(f"(b) phi=(pi/rt2)Bh[A1'/y]: worst abs residual over {len(Ns)} pts (incl N=100,1891) = {mp.nstr(worst,6)} at N={worstN}")
