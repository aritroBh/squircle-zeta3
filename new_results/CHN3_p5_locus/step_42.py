exec(open('common.py').read())
from fractions import Fraction as F
import sympy.polys.rootisolation as RI
R = sp.sympify(open('R_a1a2a4.txt').read())
G = sp.expand(R / (16*a1*(a1+1)))
A_ = sp.sympify(open('subres1_A.txt').read()); B_ = sp.sympify(open('subres1_B.txt').read())
g11l = [int(c) for c in sp.Poly(G.subs({a1:1,a2:1}), a4).all_coeffs()]
res = RI.dup_refine_real_root(g11l, F(1), F(2), sp.ZZ, steps=50)
print("raw type:", type(res), res if len(str(res))<200 else str(res)[:200])
def toF(x):
    x = sp.Rational(x); return F(int(x.p), int(x.q))
lo, hi = toF(res[0]), toF(res[1])
def peval_interval(poly, iv4):
    L, U = F(0), F(0)
    lo_, hi_ = iv4
    for e,c in sp.Poly(poly, a1,a2,a4).terms():
        c = F(int(c.p), int(c.q)); e4 = e[2]
        if c >= 0: L += c*lo_**e4; U += c*hi_**e4
        else:      L += c*hi_**e4; U += c*lo_**e4
    return L, U
LA, UA = peval_interval(A_, (lo,hi)); LB, UB = peval_interval(B_, (lo,hi))
print("A in [", float(LA), ",", float(UA), "]")
print("B in [", float(LB), ",", float(UB), "]")
if LA>0: iv5 = (-UB/LA, -LB/UA)
elif UA<0: iv5 = (-LB/UA, -UB/LA)
else: iv5=None
print("a5* interval:", [float(x) for x in iv5], " positive:", iv5[0]>0)
print("a4* interval: [", float(lo), ",", float(hi), "] positive:", lo>0)
