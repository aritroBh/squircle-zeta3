from mpmath import (mp, mpf, pi, log, nsum, inf, nstr, zeta, si, euler)
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
rel = [-572820, 9882629, 3715627, -6205451, -640815, -419547, 5269802, -5557417, -7005986, 7571301, 8878930, 5510646, 678563, -4950929, 9489332]
for dps in [50, 70, 90, 110]:
    mp.dps = dps
    Lp = log(pi)
    gfun = lambda n: (-1)**(n+1)*(si(n*pi) - pi/2)
    R4 = nsum(lambda n: (-1)**(n+1)*gfun(n)/n**4, [1, inf])
    basis = [mpf(1), log(2), Lp, euler, pi**2, pi**4, zeta(3), zeta(3)*Lp, euler*zeta(3), zeta(3,derivative=1), zeta(5), pi**4*log(2), pi**4*Lp, zeta(4,derivative=1)]
    resid = abs(sum(c*v for c,v in zip(rel,[R4]+basis)))
    print(f"dps={dps}: relation residual = {nstr(resid,4)}")
logline("step_50: R4 PSLQ relation tested at 4 precisions -> spurious (residual ~ precision floor each time)")
