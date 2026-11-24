import os, pickle
from fractions import Fraction as Fr
from math import gcd
cwd = os.getcwd()
print("workdir:", cwd)

# ---------------- exact Q(sqrt3) arithmetic: element = (m,n,d) = (m + n*sqrt3)/d, d>0, gcd(|m|,|n|,d)=1
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
    # 1/(m2+n2 r3) = (m2 - n2 r3)/(m2^2-3 n2^2)
    den = m2*m2 - 3*n2*n2
    assert den != 0
    return pmul(a, P(m2*d2, -n2*d2, den))
def ppow(a, k):
    if k < 0: return ppow(pdiv(P(1,0), a), -k)
    r = P(1,0)
    for _ in range(k): r = pmul(r, a)
    return r
def pscale(a, q):  # multiply by rational q=Fr
    m,n,d = a
    return P(m*q.numerator, n*q.numerator, d*q.denominator)
def pfloat(a):
    m,n,d = a
    return (m + n*(3.0**0.5))/d
def peq(a,b): return a==b
def psplit(a):
    """return (rational part, sqrt3 part) as Fractions: a = p + q*sqrt3"""
    m,n,d = a
    return Fr(m,d), Fr(n,d)

ONE = P(1,0); R3 = P(0,1)
P2 = P(1,1)     # 1+sqrt3, norm -2, prime above 2
P3 = P(0,1)     # sqrt3,     norm -3, prime above 3
EPS = P(2,1)    # 2+sqrt3, fundamental unit, norm +1
assert pfloat(P2)>0 and pfloat(EPS)>0

# quick sanity: (1+r3)^2 = 2*(2+r3);  3 = r3^2
assert peq(ppow(P2,2), pmul(P(2,0), EPS))
assert peq(ppow(P3,2), P(3,0))
assert peq(pmul(EPS, P(2,-1)), ONE)  # (2+r3)(2-r3)=1
print("sanity OK: (1+sqrt3)^2 = 2(2+sqrt3); 3=(sqrt3)^2; (2+sqrt3)(2-sqrt3)=1")

import math
LN_EPS = math.log(2+3**0.5)
def pnorm(a):
    m,n,d = a
    return Fr(m*m-3*n*n, d*d)

def factor_sunit(a):
    """a in Q(sqrt3)^x. If a is an S-unit for S={p2,p3}: return (sign, v2, v3, k) with
       a = sign * P2^v2 * P3^v3 * EPS^k. Else return ('NOT_SUNIT', detail)."""
    if peq(a, P(0,0)): return None
    N = pnorm(a)  # rational
    num, den = abs(N.numerator), N.denominator
    def vp(x, p):
        e = 0
        while x % p == 0: x//=p; e+=1
        return e, x
    e2n, r2n = vp(num, 2); e2d, r2d = vp(den, 2)
    if r2n != 1 or r2d != 1: return ('NOT_SUNIT', ('odd part at 2', r2n, r2d))
    e3n, r3n = vp(num, 3); e3d, r3d = vp(den, 3)
    if r3n != 1 or r3d != 1: return ('NOT_SUNIT', ('odd part at 3', r3n, r3d))
    v2 = e2n - e2d; v3 = e3n - e3d
    u = pmul(a, pmul(ppow(P2,-v2), ppow(P3,-v3)))
    # u must be +- EPS^k
    fu = abs(pfloat(u))
    k = int(round(math.log(fu)/LN_EPS)) if fu>0 else 0
    for kk in (k-1,k,k+1):
        if peq(u, ppow(EPS,kk)): return (1, v2, v3, kk)
        if peq(u, pneg(ppow(EPS,kk))): return (-1, v2, v3, kk)
    return ('NOT_SUNIT', ('unit part not EPS-power', u, fu))

# define x_i
Xs = {1: P(3,-2,3), 2: P(1,-1), 5: P(3,-1,3), 6: P(2,-1,2)}
# exact check vs definitions 1-2r3/3, 1-r3, 1-r3/3, 1-r3/2
assert peq(Xs[1], psub(ONE, pmul(Fr(2,3) and P(2,0,3), R3)))
for i in Xs: print("x_%d ="%i, Xs[i], "~", pfloat(Xs[i]))

names = {1:'p2=1+sqrt3', 2:'p3=sqrt3', 3:'eps=2+sqrt3'}
print("\n--- letter factorizations ---")
letters = {}
for i in Xs:
    letters[('x',i)] = Xs[i]
    letters[('1-x',i)] = psub(ONE, Xs[i])
for i in Xs:
    for j in Xs:
        if i!=j: letters[('x%d-x%d'%(i,j),0)] = psub(Xs[i], Xs[j])
bad = []
for nm, a in sorted(letters.items(), key=str):
    f = factor_sunit(a)
    if f[0]=='NOT_SUNIT' or f is None:
        bad.append((nm,f)); print("  !!", nm, f)
    else:
        sgn,v2,v3,k = f
        print("  ", nm, "=", "(%+d) * p2^%d p3^%d eps^%d"%(sgn,v2,v3,k), " check:", abs(pfloat(a)-sgn*pfloat(ppow(P2,v2))*pfloat(ppow(P3,v3))*pfloat(ppow(EPS,k))))
print("non-S-unit letters:", bad)
pickle.dump({'letters':letters}, open('letters_exact.pkl','wb'))
