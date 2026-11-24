# step_53: decagon pair verification + load OCT9c data for task 3
import mpmath as mp
mp.mp.dps = 120
def Li3(x): return mp.re(mp.polylog(3, mp.mpf(x)))
def Li2(x): return mp.re(mp.polylog(2, mp.mpf(x)))
def L3(x):
    x = mp.mpf(x)
    return Li3(x) - mp.log(abs(x))*Li2(x) - mp.mpf(1)/3*mp.log(abs(x))**2*mp.log(1-x)
z3 = mp.zeta(3)
L5 = (mp.zeta(3, mp.mpf(1)/5) - mp.zeta(3, mp.mpf(2)/5) - mp.zeta(3, mp.mpf(3)/5) + mp.zeta(3, mp.mpf(4)/5))/125
import sympy as sp
s5 = sp.sqrt(5)
phi = (1+mp.sqrt(5))/2
c1 = sp.Rational(-75,32) + sp.Rational(81,64)*s5
c2 = sp.Rational(-35,32) + sp.Rational(53,64)*s5
c3 = sp.Rational(33,64) - sp.Rational(19,80)*s5
c3loc = sp.Rational(5,32) - sp.Rational(1,8)*s5
a = sp.Rational(1,64) + sp.Rational(31,320)*s5
b = sp.Rational(875,768) - sp.Rational(500,768)*s5
X1, Xm1 = L3(1/phi), L3(-1/phi)
def ev(cf, val):
    return mp.mpf(str(sp.N(cf.subs(s5, mp.sqrt(5)), 40)))*val
valN = ev(c1,X1)+ev(c2,Xm1)+ev(c3-c3loc,z3)
valab = ev(a,z3)+ev(b,L5)
print("decagon: val(N3) =", valN)
print("         a*zeta3 + b*L5 =", valab, " diff:", valN-valab)
print("a =", a, " b =", b)

# ---- OCT9c data ----
import pickle, os
from fractions import Fraction as Fr
base = "<project>"
O9B = os.path.join(base,"research/08_theory/OCT9b_octagon_debug")
s123 = pickle.load(open(os.path.join(O9B,"sol_B123.pkl"),"rb"))
lists = pickle.load(open(os.path.join(O9B,"Sunit_lists.pkl"),"rb"))
print("sol_B123 keys:", list(s123.keys()) if isinstance(s123,dict) else type(s123))
print("lists keys:", list(lists.keys()))
