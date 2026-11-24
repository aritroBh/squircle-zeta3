import os, pickle, math
from fractions import Fraction as Fr
from math import gcd

def P(m, n, d=1):
    m, n, d = int(m), int(n), int(d)
    if d < 0: m, n, d = -m, -n, -d
    g = gcd(gcd(abs(m), abs(n)), d)
    return (m//g, n//g, d//g)
def padd(a, b):
    m1,n1,d1 = a; m2,n2,d2 = b
    return P(m1*d2 + m2*d1, n1*d2 + n2*d1, d1*d2)
def pneg(a): m,n,d = a; return P(-m,-n,d)
def psub(a,b): return padd(a, pneg(b))
def pmul(a, b):
    m1,n1,d1 = a; m2,n2,d2 = b
    return P(m1*m2 + 3*n1*n2, m1*n2 + n1*m2, d1*d2)
def pdiv(a, b):
    m2,n2,d2 = b
    den = m2*m2 - 3*n2*n2
    return pmul(a, P(m2*d2, -n2*d2, den))
def ppow(a, k):
    if k < 0: return ppow(pdiv(P(1,0), a), -k)
    r = P(1,0)
    for _ in range(k): r = pmul(r, a)
    return r
def pscale(a, q):
    m,n,d = a
    return P(m*q.numerator, n*q.numerator, d*q.denominator)
def pfloat(a):
    m,n,d = a
    return (m + n*(3.0**0.5))/d
def peq(a,b): return a==b
def psplit(a):
    m,n,d = a
    return Fr(m,d), Fr(n,d)
ONE = P(1,0); R3 = P(0,1)
P2 = P(1,1); P3 = P(0,1); EPS = P(2,1)
LN_EPS = math.log(2+3**0.5)
def pnorm(a):
    m,n,d = a
    return Fr(m*m-3*n*n, d*d)

def vp23(x):
    e2=e3=0
    while x%2==0: x//=2; e2+=1
    while x%3==0: x//=3; e3+=1
    return e2,e3,x

def factor_sunit(a):
    if peq(a, P(0,0)): return None
    N = pnorm(a)
    num, den = abs(N.numerator), N.denominator
    e2n,e3n,rn = vp23(num); e2d,e3d,rd = vp23(den)
    if rn != 1 or rd != 1: return ('NOT_SUNIT', ('residual', rn, rd))
    v2 = e2n - e2d; v3 = e3n - e3d
    u = pmul(a, pmul(ppow(P2,-v2), ppow(P3,-v3)))
    fu = abs(pfloat(u))
    k = int(round(math.log(fu)/LN_EPS)) if fu>0 else 0
    for kk in (k-1,k,k+1):
        if peq(u, ppow(EPS,kk)): return (1, v2, v3, kk)
        if peq(u, pneg(ppow(EPS,kk))): return (-1, v2, v3, kk)
    return ('NOT_SUNIT', ('unit', u))

Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
letters = {}
for i in Xs:
    letters[('x',i)] = Xs[i]
    letters[('1-x',i)] = psub(ONE, Xs[i])
for i in Xs:
    for j in Xs:
        if i!=j: letters[('x%d-x%d'%(i,j))] = psub(Xs[i], Xs[j])
print("--- letter factorizations (a = sign * p2^v2 p3^v3 eps^k) ---")
tab = {}
for nm, a in sorted(letters.items(), key=str):
    f = factor_sunit(a)
    assert f[0] != 'NOT_SUNIT', (nm, f)
    sgn,v2,v3,k = f
    tab[nm] = (sgn,v2,v3,k)
    print("   %-10s = (%+d) p2^%2d p3^%2d eps^%2d"%(str(nm), sgn,v2,v3,k))
print("ALL LETTERS ARE S-UNITS over S = {p2 = 1+sqrt3 (above 2), p3 = sqrt3 (above 3)}, units {+-1, eps = 2+sqrt3}")

# also gamma/delta and derived swap letters factorization spot check via formulas:
# gam(i,j) = x_j/(x_j-x_i);  gam-1 = x_i/(x_j-x_i);  dl = gam*(1-x_i);  dl-1 = x_i*(1-x_j)/(x_j-x_i)
gam = {}; dl = {}
for i in Xs:
    for j in Xs:
        if i==j: continue
        gam[(i,j)] = pdiv(Xs[j], psub(Xs[j], Xs[i]))
        dl[(i,j)]  = pmul(gam[(i,j)], psub(ONE, Xs[i]))
ok = True
for (i,j) in gam:
    assert peq(psub(gam[(i,j)], ONE), pdiv(Xs[i], psub(Xs[j], Xs[i])))
    assert peq(psub(dl[(i,j)], ONE), pdiv(pmul(Xs[i], psub(ONE,Xs[j])), psub(Xs[j], Xs[i])))
print("\ngamma/delta identities verified exactly: gam-1 = x_i/(x_j-x_i); dl-1 = x_i(1-x_j)/(x_j-x_i)")
for key in [('gam(1,2)',gam[(1,2)]), ('gam(1,5)',gam[(1,5)]), ('dl(1,5)',dl[(1,5)]), ('dl(1,6)',dl[(1,6)])]:
    f = factor_sunit(key[1]); print("   %-10s ="%key[0], key[1], "~%.10f"%pfloat(key[1]), "fac:", f)

# swap Li21 letters: x = gam_ik/gam_ij, y = gam_ij:
# letters: xy-1 = gam_ik-1 = x_i/(x_k-x_i); y-1 = x_i/(x_j-x_i); x-1 = x_i(x_j-x_k)/(x_j(x_k-x_i));
# x(y-1)/(x-1) = x_k/(x_j-x_k);  x; y
for (i,j,k) in [(1,2,5),(2,5,6)]:
    x = pdiv(gam[(i,k)], gam[(i,j)]); y = gam[(i,j)]
    L = {'xy-1': psub(pmul(x,y),ONE), 'y-1': psub(y,ONE), 'x-1': psub(x,ONE),
         'x(y-1)/(x-1)': pdiv(pmul(x, psub(y,ONE)), psub(x,ONE)), 'x': x, 'y': y, 'xy': pmul(x,y)}
    print("  swap letters (%d,%d,%d):"%(i,j,k))
    for nm,a in L.items():
        f = factor_sunit(a); print("      %-14s fac:"%nm, f, "~%.6f"%pfloat(a))
pickle.dump({'letters':letters,'gam':gam,'dl':dl,'tab':tab}, open('letters_exact.pkl','wb'))
print("\nsaved letters_exact.pkl")
