# Compare trapezoid M-dependence (discretization) and quadrature accuracy at individual points
from mpmath import mp, mpf, mpc, quad, sqrt, exp, log, pi, cos, sin
mp.dps = 100
def P_complex(p, dpstest=None):
    e = 1/p
    g = lambda u: exp(-u)*sqrt((e*exp(-(e-1)*u))**2 + (e*(1-exp(-u))**(e-1))**2)
    return 8*quad(g, [log(2), 4, 10, 30, mp.inf], method='tanh-sinh')
# 1) quadrature accuracy at several angles vs dps=150 reference
r = mpf('0.12')
for k in [0, 24, 48, 72]:
    p = 1 + r*mpc(cos(2*pi*k/96), sin(2*pi*k/96))
    v100 = P_complex(p)
    mp.dps = 150; v150 = P_complex(p); mp.dps = 100
    print(f"theta={k}/96*2pi  P err =", mp.nstr(abs(v100-mpc(v150)), 3))
