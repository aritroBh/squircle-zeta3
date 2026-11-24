exec(open('common.py').read())
from fractions import Fraction as F
import sympy.polys.rootisolation as RI
R = sp.sympify(open('R_a1a2a4.txt').read())
G = sp.expand(R / (16*a1*(a1+1)))
A_ = sp.sympify(open('subres1_A.txt').read()); B_ = sp.sympify(open('subres1_B.txt').read())
g11l = [int(c) for c in sp.Poly(G.subs({a1:1,a2:1}), a4).all_coeffs()]
lo, hi = RI.dup_refine_real_root(g11l, F(1), F(2), sp.ZZ, steps=50)
lo, hi = F(lo), F(hi)

def peval_interval(poly, iv4):
    """interval evaluation at a1=1,a2=1,a4 in iv4=[lo,hi]; returns [L,U] Fractions."""
    L, U = F(0), F(0)
    lo_, hi_ = iv4
    for e,c in sp.Poly(poly, a1,a2,a4).terms():
        c = F(int(c.p), int(c.q))
        e4 = e[2]
        pw_lo, pw_hi = lo_**e4, hi_**e4
        if c >= 0: L += c*pw_lo; U += c*pw_hi
        else:      L += c*pw_hi; U += c*pw_lo
    return L, U

LA, UA = peval_interval(A_, (lo,hi))
LB, UB = peval_interval(B_, (lo,hi))
print("A(1,1,a4*) in [", float(LA), ",", float(UA), "]  contains 0?", LA<=0<=UA)
print("B(1,1,a4*) in [", float(LB), ",", float(UB), "]")
# a5* = -B/A; A<0? B? signs:
# a5* in interval: -[LB,UB]/[LA,UA]
if LA>0:
    iv5 = (-UB/LA, -LB/UA)
elif UA<0:
    iv5 = (-LB/UA, -UB/LA)
print("a5* in [", float(iv5[0]), ",", float(iv5[1]), "]")
print("certified a4*>0:", lo>0, " a5*>0:", iv5[0]>0)
print("width a5 interval:", float(iv5[1]-iv5[0]))
