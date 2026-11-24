import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from mpmath import mp, mpf, pi, sqrt, exp, log
import dcy2_audit as A

def truth_window2(N, dps=60):
    mp.dps = dps
    Nm = mpf(N)
    tq = A.true_quant(Nm)
    y=tq['y']; Delta=tq['Delta']; rho=1/y
    l1,l2 = tq['l'][1], tq['l'][2]
    s0 = tq['s0']
    phi_s0 = (Nm-1)*log(s0)-Nm*log(1+s0)+1j*pi*s0/2
    ellt = tq['ellt']; F = tq['F']
    H0 = A.Hc(s0)
    W = mp.quad(lambda t: exp(ellt(t))*F(t), [-Delta, 0, Delta])
    P = tq['P']
    t0 = -l1/(2*l2)
    base = exp(-l2*t0**2)*sqrt(-pi/l2)
    def muinf(p):
        return base*sum(mp.binomial(p,2*m)*t0**(p-2*m)*A.df2(2*m-1)/(-2*l2)**m for m in range(0,p//2+1))
    MS = sum(cp*muinf(p) for p,cp in P.items())
    TR = mp.quad(lambda t: exp(ellt(t))*F(t), [Delta, rho, 10*rho])
    TL = mp.quad(lambda t: exp(ellt(t))*F(t), [-rho+mpf('1e-9'), -Delta])
    diff = (2/pi)*exp(phi_s0)*A.omega*H0*(W + TR + TL - MS)
    return abs(diff)*Nm**2*exp(sqrt(pi*Nm))

for N in [400, 2000, 10000, 100000, 1000000]:
    d = truth_window2(N)
    print(f"N={N:9d}: |J-model|N^2 e^th (window, H0-fixed) = {mp.nstr(d,10)}")
