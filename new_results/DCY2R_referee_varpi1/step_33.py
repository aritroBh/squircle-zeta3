import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from mpmath import mp, mpf, pi, sqrt, exp, log
mp.dps = 60
import dcy2_audit as A

def truth_window(N):
    tq = A.true_quant(N)
    y=tq['y']; Delta=tq['Delta']; rho=1/y
    l1,l2 = tq['l'][1], tq['l'][2]
    s0 = tq['s0']
    phi_s0 = (N-1)*log(s0)-N*log(1+s0)+1j*pi*s0/2
    ellt = tq['ellt']; F = tq['F']
    # window integral
    W = mp.quad(lambda t: exp(ellt(t))*F(t), [-Delta, 0, Delta])
    # model moment sum
    P = tq['P']
    t0 = -l1/(2*l2)
    base = exp(-l2*t0**2)*sqrt(-pi/l2)
    def muinf(p):
        return base*sum(mp.binomial(p,2*m)*t0**(p-2*m)*A.df2(2*m-1)/(-2*l2)**m for m in range(0,p//2+1))
    MS = sum(cp*muinf(p) for p,cp in P.items())
    # tails
    def Lam(t):
        u = 1+y*t
        Q = u**2+sqrt(2)*y*u+y**2
        return (N-1)*log(u) + (N/2)*log(1+sqrt(2)*y+y**2) - (N/2)*log(Q) - A.bR*t
    TR = mp.quad(lambda t: exp(ellt(t))*F(t), [Delta, rho, 10*rho, 4*rho**2])
    TL = mp.quad(lambda t: exp(ellt(t))*F(t), [-rho, -Delta])
    J = (2/pi)*exp(phi_s0)*A.omega*(W + TR + TL)
    model = (2/pi)*exp(phi_s0)*A.omega*MS
    return J, model

for N in [10000, 100000, 1000000]:
    J, model = truth_window(N)
    th = sqrt(pi*N)
    d = abs(J-model)*N**2*exp(th)
    print(f"N={N:9d}: |J| = {mp.nstr(abs(J),5)}, |J-model|N^2e^th = {mp.nstr(d,8)}")
