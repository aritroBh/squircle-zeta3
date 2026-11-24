import cor8r_sym as S
import sympy as sp
from fractions import Fraction
cc = S.cc
TWO = sp.Integer(2); B = 2/(1+cc); Z = sp.Integer(0)
s2   = S.symG((TWO,))
s0B  = S.symG((Z,B))
s20B = S.symG((TWO,Z,B))
s02B = S.symG((Z,TWO,B))
s0B2 = S.symG((Z,B,TWO))
prod = S.S_shuffle(s2, s0B)
sum3 = S.S_add(s20B, S.S_add(s02B, s0B2))
res = S.S_add(prod, S.S_scaleLP_frac(sum3, -1))
print("residual for G(2)*G(0,B) = G(2,0,B)+G(0,2,B)+G(0,B,2): #terms:", len(res))
# numeric value check of this identity
import mpmath as mp
mp.mp.dps = 25
c0 = mp.mpf('0.31'); B0 = 2/(1+c0)
G2 = mp.log(mp.mpf('0.5')); G0B = -mp.polylog(2,1/B0)
lhs = G2*G0B
def G3(w):
    w = [mp.mpf(2) if x==2 else (mp.mpf(0) if x==0 else B0) for x in w]
    return mp.quad(lambda t1: mp.quad(lambda t2: mp.quad(lambda t3: 1/(t3-w[2]), [mp.mpf(0),t2])/(t2-w[1]), [mp.mpf(0),t1])/(t1-w[0]), [mp.mpf(0),mp.mpf(1)])
rhs = G3([2,0,B0])+G3([0,2,B0])+G3([0,B0,2])
print("numeric: lhs=", mp.nstr(lhs,15), " rhs=", mp.nstr(rhs,15), " diff:", mp.nstr(abs(lhs-rhs),3))
print("\ns0B =", {w:dict(cf) for w,cf in s0B.items()})
print("prod #terms:", len(prod), " sum3 #terms:", len(sum3))
