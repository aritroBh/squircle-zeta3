import mpmath as mp
mp.mp.dps = 420
phi = (1+mp.sqrt(5))/2
s5 = mp.sqrt(5)
z3 = mp.zeta(3)
def chi5(k):
    k = k % 5
    return 0 if k==0 else (1 if k in (1,4) else -1)
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
L5 = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5))/125
# careful direct sum: block terms by period of 5, then euler-maclaurin over blocks
def block(m):
    return sum(chi5(a)/mp.mpf(5*m+a)**3 for a in range(1,5))
L5_direct = mp.nsum(block, [0, mp.inf], method='euler-maclaurin')
print('L5 Hurwitz vs block EM sum:', mp.nstr(L5-L5_direct, 5))
# N3 from DEC2 proved table coefficients
from fractions import Fraction as Fr
c1 = (Fr(-75,32), Fr(81,64)); c2 = (Fr(-35,32), Fr(53,64))
c3 = (Fr(33,64), Fr(-19,80)); c3loc = (Fr(5,32), Fr(-1,8))
def evQ(pq, s5v):
    return mp.mpf(pq[0].numerator)/pq[0].denominator + s5v*mp.mpf(pq[1].numerator)/pq[1].denominator
def qsub(x,y): return (x[0]-y[0], x[1]-y[1])
N3 = evQ(c1,s5)*L3(1/phi) + evQ(c2,s5)*L3(-1/phi) + evQ(qsub(c3,c3loc),s5)*z3
print('N3 =', mp.nstr(N3, 60))
# claimed pair
a = (Fr(1,64), Fr(31,320)); b = (Fr(875,768), Fr(-125,192))
valab = evQ(a,s5)*z3 + evQ(b,s5)*L5
print('a*z3+b*L5 =', mp.nstr(valab, 60))
print('claim C residual (420 dps):', mp.nstr(N3-valab, 5))
# PSLQ at this precision
vec = [N3, z3, s5*z3, L5, s5*L5]
rel = mp.pslq(vec, tol=mp.mpf('1e-370'), maxcoeff=10**8)
print('pslq(420dps) [N3, z3, s5*z3, L5, s5*L5]:', rel)
if rel:
    q0,q1,q2,q3,q4 = rel
    print('=> N3 = (%s/%s) z3 + (%s/%s) s5 z3 + (%s/%s) L5 + (%s/%s) s5 L5' % (-q1,q0,-q2,q0,-q3,q0,-q4,q0))
    from fractions import Fraction as Fr2
    print('a recovered:', Fr2(-q1,q0), Fr2(-q2,q0), ' b recovered:', Fr2(-q3,q0), Fr2(-q4,q0))
    print('match PGEN1 a:', (Fr2(-q1,q0),Fr2(-q2,q0))==a, ' b:', (Fr2(-q3,q0),Fr2(-q4,q0))==b)
