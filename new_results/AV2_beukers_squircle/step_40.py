
from mpmath import mp, mpf, zeta, fabs, nstr, nsum, inf
mp.dps = 100
exec(open('av2core.py').read())
from math import comb, gcd
from fractions import Fraction
from functools import reduce
def lcm(a, b): return a*b//gcd(a, b)
def L(n): return reduce(lcm, range(1, n+1), 1)
def apery(n): return sum(comb(n,k)**2 * comb(n+k,k)**2 for k in range(n+1))

def Minf(j, k):
    return nsum(lambda m: mpf(j+k+2*m+2)/((j+m+1)**2*(k+m+1)**2), [0, inf])
Mcache = {}
def Ln_inf(n):
    cs = shifted_legendre_coeffs(n)
    tot = mpf(0)
    for i, ci in enumerate(cs):
        for j, cj in enumerate(cs):
            if ci == 0 or cj == 0: continue
            key = (min(i,j), max(i,j))
            if key not in Mcache: Mcache[key] = Minf(*key)
            tot += ci*cj*Mcache[key]
    return tot
print("n | L_n(inf) | a_n(Beukers?) | a_n*LCM^3 L_n - a_n zeta(3) = b_n (rational?) | |L_n|^{1/n}")
z3 = zeta(3)
rows = []
for n in range(0, 13):
    Ln = Ln_inf(n)
    dn3 = L(n)**3
    an = apery(n)
    resid = dn3*Ln - an*z3
    # reconstruct b_n as fraction
    fr = Fraction(float(resid)).limit_denominator(10**12)
    ok = fabs(resid - mpf(fr.numerator)/mpf(fr.denominator)) < mpf('1e-30')
    rate = fabs(Ln)**(1.0/n) if n > 0 else mpf(0)
    print("n=%2d L=%s" % (n, nstr(Ln, 20)))
    print("   dn3*L - a_n*zeta3 = %s ; fraction=%s exact=%s ; |L|^(1/n)=%s ; a_n=%d" % (
        nstr(resid, 12), fr if ok else "??", ok, nstr(rate, 10), an))
    rows.append((n, mp.nstr(Ln, 50), an, str(fr) if ok else None, mp.nstr(rate, 15)))
import json
json.dump(rows, open('legendre_inf.json','w'))
print("(sqrt2-1)^4 =", nstr((2**0.5-1)**4, 15))
