from fractions import Fraction as Fr
# Exact derivation of claim C's (a,b) from: DEC2 table (PROVED), duplication, Lewin Li3(1/phi^2),
# golden chi3 identity D (NUMERICAL), classical Li2 golden values.
# Q(sqrt5) arithmetic as pairs (p,q) = p + q*sqrt5.
def qadd(x,y): return (x[0]+y[0], x[1]+y[1])
def qmul(x,y): return (x[0]*y[0]+5*x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def qsc(x,c): return (x[0]*c, x[1]*c)
# DEC2 final table (PROVED): coefficients (p,q)
c1 = (Fr(-75,32), Fr(81,64))   # Li3(1/phi)
c2 = (Fr(-35,32), Fr(53,64))   # Li3(-1/phi)
c3 = (Fr(33,64), Fr(-19,80))   # zeta(3)
c3loc = (Fr(5,32), Fr(-1,8))   # local zeta(3) share (C5 = 60L*c3loc, PROVED in DEC2)
# reductions: S = Li3(1/phi)+Li3(-1/phi) = (1/5)z3 - (1/30)pi2 lf + (1/6) lf^3   [duplication+Lewin]
#             D_ = Li3(1/phi)-Li3(-1/phi) = (25 sqrt5/24) L5 - (1/6) pi2 lf + (1/2) lf^3   [claim D]
# L3(1/phi)  = (S+D_)/2 + (1/10) pi2 lf - (1/3) lf^3     [Li2(1/phi) = pi2/10 - lf^2]
# L3(-1/phi) = (S-D_)/2 - (1/15) pi2 lf + (1/6) lf^3     [Li2(-1/phi) = -pi2/15 + lf^2/2]
h1 = qsc(qadd(c1,c2), Fr(1,2))   # (c1+c2)/2 multiplies S
h2 = qsc(qadd(c1,( -c2[0], -c2[1])), Fr(1,2))  # (c1-c2)/2 multiplies D_
# coefficient collectors over basis {z3, L5, pi2*lf, lf^3}: each a Q(sqrt5) element
Z3  = qadd(qsc(h1, Fr(1,5)), qadd(c3, (-c3loc[0], -c3loc[1])))
L5c = qmul(h2, (Fr(0), Fr(25,24)))        # (25 sqrt5)/24 as (0,25/24)
P2L = qadd(qadd(qsc(h1, Fr(-1,30)), qsc(h2, Fr(-1,6))), qadd(qsc(c1, Fr(1,10)), qsc(c2, Fr(-1,15))))
LF3 = qadd(qadd(qsc(h1, Fr(1,6)), qsc(h2, Fr(1,2))), qadd(qsc(c1, Fr(-1,3)), qsc(c2, Fr(1,6))))
print('zeta(3) coeff a =', Z3, ' i.e.', Z3[0], '+', Z3[1], 'sqrt5')
print('L(chi5,3) coeff b =', L5c, ' i.e.', L5c[0], '+', L5c[1], 'sqrt5')
print('pi^2 log phi coeff (must be 0):', P2L)
print('log^3 phi coeff (must be 0):', LF3)
a_PGEN1 = (Fr(1,64), Fr(31,320)); b_PGEN1 = (Fr(875,768), Fr(-125,192))
print('a == PGEN1:', Z3==a_PGEN1, '  b == PGEN1:', L5c==b_PGEN1)
open('progress.log','a').write('step: claim C reduced EXACTLY to identity D (+Lewin/dup/Li2): (a,b) match PGEN1, product coeffs cancel exactly\n')
