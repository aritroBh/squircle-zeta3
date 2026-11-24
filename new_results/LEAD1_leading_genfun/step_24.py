import json, mpmath
from fractions import Fraction as Fr
mpmath.mp.dps = 410
from mpmath import mpf, log, pi, zeta
Dtop = json.load(open('Dtop.json'))
N = json.load(open('Nints_dps420.json'))['N']
def Nval(i,j,pole): return mpf(N[f"{i},{j}|{pole}"])
cw = {}
for w in range(2,11):
    d = w-1
    tot = mpf(0)
    for name, rec in Dtop.items():
        dd, key = name.split('|'); i, j = map(int, key.split(','))
        if int(dd) != d: continue
        for pole, pr in rec.items():
            cf = Fr(pr['Dtop'][0])  # sqrt2 part always 0
            if cf: tot += mpf(cf.numerator)/cf.denominator * Nval(i,j,pole)
    cw[w] = tot / mpmath.factorial(d-1)
    print(f"c({w}) =", mpmath.nstr(cw[w], 50))
# controls
lam = log(2)
c2t = (pi**2 + 6*lam**2)/24
c3t = -(zeta(3)/8 + pi**2*lam/24 + lam**3/12)
c4t = -(pi**4/320 + pi**2*lam**2/48 + lam**4/48 + zeta(3)*lam/8)
c6t = (168*lam**6 + 420*pi**2*lam**4 + 378*pi**4*lam**2 + 61*pi**6 + 5040*zeta(3)*lam**3
       + 2520*pi**2*zeta(3)*lam + 3780*zeta(3)**2 + 22680*zeta(5)*lam)/mpmath.factorial(9)
for w, t in [(2,c2t),(3,c3t),(4,c4t),(6,c6t)]:
    print(f"control w={w}: diff =", mpmath.nstr(cw[w]-t, 10))
json.dump({str(w): mpmath.nstr(cw[w], 405) for w in cw}, open('cw_dps420.json','w'))
