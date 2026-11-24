import os, pickle
import sympy as sp, mpmath as mp
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
with open(f'{WD}/coef_acc_raw.pkl','rb') as f: coef_acc = pickle.load(f)
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
# zeta5 monomial: ((),0,(0,0,1)); z2z3: ((),0,(1,1,0))
mk_z5 = ((),0,(0,0,1)); mk_z2z3 = ((),0,(1,1,0))
print("z5 monomial present:", mk_z5 in coef_acc, " z2z3 present:", mk_z2z3 in coef_acc)
# coefficient of z5 (before 120x)
co_z5 = coef_acc.get(mk_z5, sp.Integer(0))
co_z2z3 = coef_acc.get(mk_z2z3, sp.Integer(0))
# For edge 3, logs don't matter for zeta coeffs? They might. Substitute with logs anyway.
rd3 = {A0:sp.Rational(1,4),A1:sp.Rational(1,4),B:sp.Rational(3,4),kap:sp.Integer(0),
      L2:sp.Integer(2),q0:sp.Integer(-1),V02:sp.Integer(1),
      lA0:sp.log(sp.Rational(1,4)),lA1:sp.log(sp.Rational(1,4)),lB:sp.log(sp.Rational(1,4))}
z5coef = sp.simplify(120*co_z5.subs(rd3))
z23coef = sp.simplify(120*co_z2z3.subs(rd3))
print("z5 coeff (edge3) =", z5coef, "=", sp.N(z5coef,20))
print("Z_e^(5) expected 3195*sqrt(2)/2048 =", sp.N(3195*sp.sqrt(2)/2048,20))
print("z2z3 coeff (edge3) =", z23coef, "=", sp.N(z23coef,20))
print("Z_e^(2,3) expected 1515*sqrt(2)/2048 =", sp.N(1515*sp.sqrt(2)/2048,20))
