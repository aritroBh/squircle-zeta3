import sys, time
sys.path.insert(0, '../DCY2_decay_first_correction')
from mpmath import mp, mpf, pi, sqrt, exp, log
import dcy2_audit as A

def truth_window(N, dps=60):
    mp.dps = dps
    tq = A.true_quant(N)
    y=tq['y']; Delta=tq['Delta']; rho=1/y
    l1,l2 = tq['l'][1], tq['l'][2]
    s0 = tq['s0']
    phi_s0 = (N-1)*log(s0)-N*log(1+s0)+1j*pi*s0/2
    ellt = tq['ellt']; F = tq['F']
    W = mp.quad(lambda t: exp(ellt(t))*F(t), [-Delta, 0, Delta])
    P = tq['P']
    t0 = -l1/(2*l2)
    base = exp(-l2*t0**2)*sqrt(-pi/l2)
    def muinf(p):
        return base*sum(mp.binomial(p,2*m)*t0**(p-2*m)*A.df2(2*m-1)/(-2*l2)**m for m in range(0,p//2+1))
    MS = sum(cp*muinf(p) for p,cp in P.items())
    TR = mp.quad(lambda t: exp(ellt(t))*F(t), [Delta, rho, 10*rho])
    TL = mp.quad(lambda t: exp(ellt(t))*F(t), [-rho+mpf('1e-9'), -Delta])
    J = (2/pi)*exp(phi_s0)*A.omega*(W + TR + TL)
    model = (2/pi)*exp(phi_s0)*A.omega*MS
    return J, model, abs(J-model)*N**2*exp(sqrt(pi*N))

for N in [400, 2000, 10000, 100000]:
    J, model, d = truth_window(N)
    print(f"N={N:9d}: window-truth |J-model|N^2e^th = {mp.nstr(d,8)}")
