import os, pickle
import sympy as sp, mpmath as mp
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
with open(f'{WD}/coef_acc_raw.pkl','rb') as f: coef_acc = pickle.load(f)
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
# edge 1 (r=2), polar logs
rd1 = {A0:sp.Rational(1,2),A1:sp.Rational(3,4),B:sp.Rational(1,4),kap:sp.Integer(1),
      L2:sp.Integer(5),q0:sp.Integer(1),V02:sp.Integer(1),
      lA0:sp.log(sp.Rational(1,2)),lA1:sp.log(sp.Rational(1,4)),lB:sp.log(sp.Rational(1,4))}
for name,mk in [('z5',((),0,(0,0,1))),('z2z3',((),0,(1,1,0)))]:
    co = coef_acc.get(mk, sp.Integer(0))
    v = sp.simplify(120*co.subs(rd1))
    print(f"edge1 {name} coeff (P_arc) = {v} = {sp.N(v,20)}")
# E5e: X = P_arc/sqrt(5); X z5 coord = -20187/2000, z2z3 = -858/125
# => P_arc z5 = sqrt(5)*(-20187/2000), z2z3 = sqrt(5)*(-858/125)
print("\nE5e expectations for P_arc:")
print("  z5  = sqrt(5)*(-20187/2000) =", sp.N(sp.sqrt(5)*sp.Rational(-20187,2000),20))
print("  z2z3= sqrt(5)*(-858/125)   =", sp.N(sp.sqrt(5)*sp.Rational(-858,125),20))
print("\n local Z_e: z5=-2079*sqrt(5)/250 =", sp.N(-2079*sp.sqrt(5)/250,20),
      " z2z3=-477*sqrt(5)/250 =", sp.N(-477*sp.sqrt(5)/250,20))
