import sympy as sp, pickle, os, sys
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
t = sp.Symbol('t'); r5 = sp.sqrt(5); phi = (1+r5)/2
kap1 = 1/phi; kap4 = -1/phi**2
u1 = 1+kap1*t; u4 = 1+kap4*t
w1 = sp.log(u1); w4 = sp.log(u4)
lf = sp.log(phi)
import dec_engine as E
# --- D9: IBP: d(w1^2 w4 /2) = kap1 w1 w4/u1 + (kap4/2) w1^2/u4
expr = sp.expand(sp.diff(w1**2*w4/2, t) - kap1*(w1*w4/u1) - sp.Rational(1,2)*kap4*(w1**2/u4))
print("D9 IBP identity residual (0 expected):", sp.simplify(expr))
# boundary: [w1^2 w4 /2]_0^1 = w1(1)^2 w4(1)/2 = lf^2*(-lf)/2 = -lf^3/2
bnd = sp.expand((w1**2*w4/2).subs(t,1) - (w1**2*w4/2).subs(t,0))
print("D9 boundary:", sp.simplify(bnd.subs(sp.log(phi), E.lf)))
D9 = sp.expand((-E.lf**3/2 - kap4*E.D7/2)/kap1)
D9 = sp.expand(D9.subs(E.lf, E.lf))
print("D9 diff vs engine:", sp.expand(D9 - E.D9))
# --- D10: d(w1 w4^2 /2) = (kap1/2) w4^2/u1 + kap4 w1 w4/u4
expr2 = sp.expand(sp.diff(w1*w4**2/2, t) - sp.Rational(1,2)*kap1*(w4**2/u1) - kap4*(w1*w4/u4))
print("D10 IBP identity residual (0 expected):", sp.simplify(expr2))
bnd2 = sp.expand((w1*w4**2/2).subs(t,1) - (w1*w4**2/2).subs(t,0))
print("D10 boundary:", sp.simplify(bnd2.subs(sp.log(phi), E.lf)))
D10 = sp.expand((E.lf**3/2 - kap1*E.D8/2)/kap4)
print("D10 diff vs engine:", sp.expand(D10 - E.D10))
