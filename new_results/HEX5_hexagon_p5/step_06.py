import re
from mpmath import mp, mpf, log, pi, zeta, sqrt, nstr, pslq

txt = open('p5_run150.txt').read()
val_str = re.search(r'P5/sqrt3 = ([0-9.]+)', txt).group(1)
print("stored P5/sqrt3 digits:", len(val_str.replace('.','')))

names = ['1','log2','pi^2','log2^2','z3','pi^2 log2','log2^3','pi^4','z3 log2',
         'pi^2 log2^2','log2^4','z5','z2 z3','pi^4 log2','z3 log2^2','pi^2 log2^3','log2^5']

for DPS in [100, 125, 145]:
    mp.dps = DPS
    L2 = log(2)
    basis = [mpf(1), L2, pi**2, L2**2, zeta(3), pi**2*L2, L2**3, pi**4, zeta(3)*L2,
             pi**2*L2**2, L2**4, zeta(5), (pi**2/6)*zeta(3), pi**4*L2, zeta(3)*L2**2,
             pi**2*L2**3, L2**5]
    x = [mpf(val_str)] + basis
    rel = pslq(x, tol=mpf(10)**(-(DPS-25)), maxcoeff=10**9, maxsteps=10000)
    if rel is None:
        print("DPS=%d: PSLQ found NO relation" % DPS); continue
    # relation: sum rel[i]*x[i] = 0 ; check
    resid = abs(sum(mpf(rel[i])*x[i] for i in range(len(x))))
    print("DPS=%d: relation found, residual = %s" % (DPS, nstr(resid, 6)))
    print("  rel =", rel)
    q0 = rel[0]
    if q0 != 0:
        print("  coefficients c_i = -rel[i+1]/rel[0]  (target: P5/sqrt3 = sum c_i b_i):")
        for i, nm in enumerate(names):
            c = mpf(-rel[i+1])/mpf(rel[0])
            if rel[i+1] != 0:
                print("    %-12s : %s   (~%s)" % (nm, str(-rel[i+1])+'/'+str(rel[0]), nstr(c, 12)))
    print("  z5 coeff exact? -rel[12]/rel[0] vs 2835/32: ",
          mpf(-rel[12])/rel[0] == mpf(2835)/32 if len(rel)>12 else '?')
