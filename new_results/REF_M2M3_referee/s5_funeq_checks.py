"""Referee step 2d: the primitives used in the weight-3 evaluation (by differentiation) and the
functional equations with their stated real domains (random points inside; a few points outside to
show that the domain matters)."""
import random
import mpmath as mp

mp.mp.dps = 50
random.seed(1)
Li2 = lambda u: mp.polylog(2, u)
Li3 = lambda u: mp.polylog(3, u)
lg = mp.log
z3 = mp.zeta(3)
pi2 = mp.pi**2


def P1(s):
    return -Li2(1 - s) if s > 0 else lg(-s)*lg(1 - s) + Li2(s)


def P2(s):
    if 0 < s < 1:
        return lg(s)**2*lg(1 - s) + 2*lg(s)*Li2(s) - 2*Li3(s)
    if s < 0:
        return lg(-s)**2*lg(1 - s) + 2*lg(-s)*Li2(s) - 2*Li3(s)
    r = 1/s
    return -lg(r)**3/3 + lg(r)**2*lg(1 - r) + 2*lg(r)*Li2(r) - 2*Li3(r)


worst = 0
for s in [mp.mpf(x) for x in ('-7.3', '-0.4', '0.2', '0.83', '1.4', '9.1')]:
    d1 = abs(mp.diff(P1, s) - lg(abs(s))/(s - 1))
    d2 = abs(mp.diff(P2, s) - lg(abs(s))**2/(s - 1))
    worst = max(worst, d1, d2)
print('primitives P1, P2: max |P\' - integrand| at 6 points in the 3 regions:', mp.nstr(worst, 3))
print('limits: P2(1-) =', mp.nstr(P2(1 - mp.mpf(10)**-30), 12), ' P2(1+) =', mp.nstr(P2(1 + mp.mpf(10)**-30), 12),
      ' -2 zeta3 =', mp.nstr(-2*z3, 12))


def K(x, y):
    lhs = (Li3(x*(1 - y)**2/(y*(1 - x)**2)) + Li3(x*y) + Li3(x/y) - 2*Li3(x*(1 - y)/(y*(1 - x)))
           - 2*Li3(x*(1 - y)/(x - 1)) - 2*Li3((1 - y)/(1 - x)) - 2*Li3((1 - y)/(y*(x - 1)))
           - 2*Li3(x) - 2*Li3(y) + 2*z3)
    rhs = lg(y)**2*lg((1 - y)/(1 - x)) - pi2/3*lg(y) - lg(y)**3/3
    return abs(mp.re(lhs) - rhs)


def Abel(x, y):
    lhs = Li2(x/(1 - y)) + Li2(y/(1 - x)) - Li2(x*y/((1 - x)*(1 - y)))
    rhs = Li2(x) + Li2(y) + lg(1 - x)*lg(1 - y)
    return abs(mp.re(lhs) - rhs)


wk = 0
for _ in range(300):
    y = mp.mpf(random.random())
    x = y - mp.mpf(random.expovariate(0.3))      # x < y, x often very negative
    wk = max(wk, K(x, y))
print('Kummer (K) real form, 300 random points of D={0<y<1, x<y}: max residual', mp.nstr(wk, 3))
print('  four instances used:', [mp.nstr(K(x, y), 3) for x, y in
      [(mp.sqrt(2) - 1, 1/mp.sqrt(2)), (2 - mp.sqrt(2), 1/mp.sqrt(2)), (1/mp.sqrt(2), 2*(mp.sqrt(2) - 1)),
       (-1, (mp.sqrt(2) - 1)**2)]])
for pt in [('0.8', '0.3'), ('1.5', '0.5'), ('0.3', '1.5'), ('-0.5', '-0.2')]:
    try:
        print('  outside D', pt, ': residual of the real-part form', mp.nstr(K(mp.mpf(pt[0]), mp.mpf(pt[1])), 5))
    except Exception as ex:
        print('  outside D', pt, ': undefined (', type(ex).__name__, ')')
wa = 0
for _ in range(300):
    x = 1 - mp.mpf(random.expovariate(0.3))
    y = min(1 - mp.mpf(random.expovariate(0.3)), 1 - x - mp.mpf(random.random())*3)
    wa = max(wa, Abel(x, y))
print('Abel (A), 300 random points of {x<1,y<1,x+y<1}: max residual', mp.nstr(wa, 3))

# two-variable sanity of the one-variable relations
wr = 0
for _ in range(100):
    z = mp.mpf(random.random()); zn = -mp.mpf(random.expovariate(0.2)); w = 2*mp.mpf(random.random()) - 1
    wr = max(wr,
             abs(Li3(z) + Li3(1 - z) + Li3(1 - 1/z) - (z3 + lg(z)**3/6 + pi2/6*lg(z) - lg(z)**2*lg(1 - z)/2)),
             abs(Li3(zn) - Li3(1/zn) + pi2/6*lg(-zn) + lg(-zn)**3/6),
             abs(Li3(w*w) - 4*Li3(w) - 4*Li3(-w)),
             abs(Li2(z) + Li2(1 - z) - pi2/6 + lg(z)*lg(1 - z)),
             abs(Li2(zn) + Li2(1/zn) + pi2/6 + lg(-zn)**2/2),
             abs(Li2(zn) + Li2(zn/(zn - 1)) + lg(1 - zn)**2/2))
print('(T3),(I3),(D3),(R),(I),(La) at 100 random points each: max residual', mp.nstr(wr, 3))
