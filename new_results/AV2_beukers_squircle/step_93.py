import os
os.chdir('<project>/research/08_theory/AV2_beukers_squircle')
from mpmath import mp, mpf, log, log1p, exp, nstr
mp.dps = 60
p = mpf(1024); ep = 1/p
for r in [mpf('0.3'), mpf('0.7'), mpf('0.9'), mpf('0.99'), mpf('0.999'), mpf('0.9999'), mpf(1)-mpf('1e-6')]:
    lQ = ep*log1p(-r**p)
    Q = exp(lQ)
    one_Q = -mp.expm1(lQ)
    print("r=%8s  -log(one_Q)=%s   r*(-log)/1=%s" % (r, nstr(-log(one_Q),12), nstr(r*(-log(one_Q)),12)))
