# Slope-2/e law of the notes (T3), mentioned but not stated in the paper: L_n = 1 + x* + c2/n^2 + O(log n/n^3).
from mpmath import mp, mpf, log, exp, e, pi, nstr, findroot, expm1
import importlib.util, sys
mp.dps = 50
src = open('s2_fn.py').read().split('# ---- test functions')[1].split('pr("== Thm fn')[0]
exec('from mpmath import *\nmp.dps=50\nE1=e-1\n# ---- test functions' + src)
out = []
for name in ['(e^x-1)/(e-1)', '2x/(1+x)', 'x']:
    F, Fp, a, A, B = funcs[name]
    c2 = ((a**2 - A)*pi**2/12 + (a**2 - a - A)/2)/a**3
    for n in ['1e4', '1e6', '1e8']:
        nn = mpf(n)
        g = lambda t: log(nn) + (nn-1)*log(F(1-t)) + log(Fp(1-t)) - log(2/e)
        t0 = log(nn*a*e/2)/(nn*a)
        ts = findroot(g, t0)
        D = defect(F, Fp, a, n)
        r = (1 - D - (1 - ts))  # L-1-x*
        out.append("%s n=%s (L-1-x*)*n^2=%s c2=%s  (resid)*n^3/log n=%s" % (name, n, nstr(r*nn**2, 12), nstr(c2, 12), nstr((r - c2/nn**2)*nn**3/log(nn), 6)))
print("\n".join(out)); open('s2b_slope.out','w').write("\n".join(out)+"\n")
