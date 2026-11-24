import sys
sys.path.insert(0, '../DCY2_decay_first_correction')
from mpmath import mp, mpf, pi, sqrt, exp, log, gamma, gammainc
mp.dps = 45
import dcy2_common as C
from dcy2_common import Hc, omega, wb, df2, j as jj
bR = pi/(2*sqrt(2))

def true_pieces(N):
    y = sqrt(pi/(2*N)); rho = 1/y
    s0 = rho*omega
    phi_s0 = (N-1)*log(s0) - N*log(1+s0) + jj*pi*s0/2
    def ell(m):
        if m == 1: return (N-1)*y - N*y/(1+y*wb) + jj*pi*omega/2
        return (-1)**(m+1) * y**m/m * ((N-1) - N/(1+y*wb)**m)
    l1,l2 = ell(1), ell(2)
    kap2 = mpf('2.48')**2 + 4*log(N/400)
    kap = sqrt(kap2)
    c2m = mpf('0.9')*bR*y
    sig = 1/sqrt(2*c2m)
    Delta = kap*sig
    # kept polynomial P and its coefficients (complex)
    H0 = Hc(s0)
    def Hder(m):
        return (-1)**m * mp.quad(lambda u: u**m*sqrt(1-exp(-2*u))*exp(-s0*u), [0, mp.inf])
    h = [1] + [omega**m*Hder(m)/H0/mp.factorial(m) for m in range(1,4)]
    l3,l4,l5 = ell(3), ell(4), ell(5)
    P = {}
    def add(p,c): P[p] = P.get(p,0)+c
    add(0,1); add(1,h[1]); add(2,h[2]); add(3,h[3])
    add(3,l3); add(4,l4); add(5,l5)
    add(6,l3**2/2); add(7,l3*l4); add(9,l3**3/6)
    add(4,h[1]*l3); add(5,h[1]*l4); add(5,h[2]*l3); add(7,h[1]*l3**2/2)
    # true mu_p^inf and mu_p(Delta)
    t0 = -l1/(2*l2)
    base = exp(-l2*t0**2)*sqrt(-pi/l2)
    def muinf(p):
        return base*sum(mp.binomial(p,2*m)*t0**(p-2*m)*df2(2*m-1)/(-2*l2)**m for m in range(0,p//2+1))
    def muD(p):
        f = lambda t: t**p*exp(l1*t+l2*t**2)
        return mp.quad(f, [-Delta, 0, Delta])
    drop = sum(cp*(muinf(p)-muD(p)) for p,cp in P.items())
    # true tail integral
    def Lam(t):
        u = 1+y*t
        Q = u**2+sqrt(2)*y*u+y**2
        return (N-1)*log(u) + (N/2)*log(1+sqrt(2)*y+y**2) - (N/2)*log(Q) - bR*t
    F = lambda t: Hc(s0+omega*t)/H0
    truetailR = mp.quad(lambda t: abs(F(t))*exp(Lam(t)), [Delta, 2*rho, mp.inf])
    truetailL = mp.quad(lambda t: abs(F(t))*exp(Lam(t)), [-rho+mpf('1e-6'), -Delta])
    return dict(y=y, Delta=Delta, drop=abs(drop), truetailR=truetailR, truetailL=truetailL,
                eRe=exp(phi_s0.real), N=N, th=sqrt(pi*N))

for N in [30000, 100000, 1000000]:
    tp = true_pieces(N)
    # scale: (2/pi)|e^{phi}| omega * pieces -> times N^2 e^{theta}
    fac = (2/pi)*tp['eRe']*tp['N']**2*exp(tp['th'])
    print(f"N={N}: drop-piece scaled = {float(fac*tp['drop']):.4f}, tailR scaled = {float(fac*tp['truetailR']):.4f}, tailL scaled = {float(fac*tp['truetailL']):.4f}")
