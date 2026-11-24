import os
os.chdir("<project>/research/08_theory/DCY2_decay_first_correction")
from mpmath import mp, mpf, pi, sqrt, exp, log
mp.dps = 50
bR = pi/(2*sqrt(2))
def LogM(N, t):
    N = mpf(N); y = sqrt(pi/(2*N)); kap = sqrt(2.48**2 + 4*log(N/400))
    sig = 1/sqrt(2*0.9*bR*y)
    Delta = kap*sig
    t = mpf(t)*Delta
    u = 1+y*t
    Q = u**2 + sqrt(2)*y*u + y**2
    Q0 = 1+sqrt(2)*y+y**2
    return (N-1)*log(u) + (N/2)*log(Q0) - (N/2)*log(Q) - bR*t
for N0 in [400, 3000, 10000, 30000, 100000]:
    N0 = mpf(N0)
    h = max(N0/1000, 1)
    d1 = (LogM(N0+h,-1) - LogM(N0-h,-1))/(2*h)
    d2 = (LogM(N0+h,-1) - 2*LogM(N0,-1) + LogM(N0-h,-1))/h**2
    print(f"N={int(N0):7d}  LogM={float(LogM(N0,-1)):.6f}  |d1|={float(abs(d1)):.3e}  |d2|={float(abs(d2)):.3e}")
