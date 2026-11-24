import os, pickle
import sympy as sp
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
with open(f'{WD}/coef_acc_raw.pkl','rb') as f: coef_acc = pickle.load(f)
A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
# generators
l2,pi,z3,z5,Li4h,Li5h = sp.symbols('l2 pi z3 z5 Li4h Li5h')
pi2=pi**2; pi4=pi**4
# atom values at x=-1 in generators (z2z3 -> pi^2*z3/6)
def Z23(cf): return cf*pi**2*z3/6
atomval = {
 'H1': -l2,
 'Li2': -pi**2/12,
 'Li3': -3*z3/4,
 'H011': z3/8,
 'Li4': -7*pi**4/720,
 'H0011': -pi**4/48 + sp.Rational(7,4)*z3*l2 - sp.Rational(1,12)*pi**2*l2**2 + sp.Rational(1,12)*l2**4 + 2*Li4h,
 'H0111': -pi**4/90 + sp.Rational(7,8)*z3*l2 - sp.Rational(1,24)*pi**2*l2**2 + sp.Rational(1,24)*l2**4 + Li4h,
 'Li5': -sp.Rational(15,16)*z5,
 'H00011': -sp.Rational(29,32)*z5 + Z23(sp.Rational(1,2)),
 'H00101': sp.Rational(41,32)*z5 - Z23(sp.Rational(5,8)),
 'H00111': sp.Rational(33,32)*z5 + Z23(sp.Rational(1,2)) - sp.Rational(7,8)*z3*l2**2 + sp.Rational(1,18)*pi**2*l2**3 - sp.Rational(1,15)*l2**5 - 2*Li5h - 2*Li4h*l2,
 'H01011': -sp.Rational(73,64)*z5 - Z23(sp.Rational(3,2)) + sp.Rational(7,4)*z3*l2**2 - sp.Rational(1,9)*pi**2*l2**3 + sp.Rational(2,15)*l2**5 + 4*Li5h + 4*Li4h*l2,
 'H01111': z5 - sp.Rational(7,16)*z3*l2**2 + sp.Rational(1,36)*pi**2*l2**3 - sp.Rational(1,30)*l2**5 - Li5h - Li4h*l2,
}
# edge1 exact: A0=1/2,A1=3/4,B=1/4,kap=1,L2=5,q0=1,V02=1; lA0=-l2,lA1=-2l2,lB=-2l2; x=-1
rd1 = {A0:sp.Rational(1,2),A1:sp.Rational(3,4),B:sp.Rational(1,4),kap:sp.Integer(1),
      L2:sp.Integer(5),q0:sp.Integer(1),V02:sp.Integer(1),lA0:-l2,lA1:-2*l2,lB:-2*l2}
# z-symbols in coef_acc are z2,z3,z5 (sympy) -> z2=pi^2/6, z3, z5
z2s,z3s,z5s = sp.symbols('z2 z3 z5')
zsub = {z2s: pi**2/6, z3s: z3, z5s: z5}
import time; t0=time.time()
folded = sp.Integer(0)
for (at,xe,ze),co in coef_acc.items():
    cov = co.subs(rd1)  # exact in l2, sqrt5
    # monomial value: atom product * (-1)^xe * z^ze
    mv = sp.Integer((-1)**xe)
    for nm,e in at:
        mv *= atomval[nm]**e
    mv *= z2s**ze[0]*z3s**ze[1]*z5s**ze[2]
    mv = mv.subs(zsub)
    folded += cov*mv
folded = sp.expand(120*folded)
print("folding done, %.1f s"%(time.time()-t0))
# extract coeff of z5 and of pi^2*z3
c_z5 = sp.expand(folded).coeff(z5,1)
# ensure no other gens in the z5 coeff
c_z5 = sp.simplify(c_z5.subs({l2:0,pi:0,z3:0,Li4h:0,Li5h:0}) + (c_z5 - c_z5.subs({l2:0,pi:0,z3:0,Li4h:0,Li5h:0})))
# Actually isolate terms that are exactly (const)*z5
cz5_terms = [t for t in sp.Add.make_args(sp.expand(folded)) if t.has(z5) and not any(t.has(g) for g in [l2,pi,z3,Li4h,Li5h])]
c_z5 = sp.simplify(sum(cz5_terms)/z5)
# pi^2 z3 terms
cz23_terms = [t for t in sp.Add.make_args(sp.expand(folded)) if t.has(z3) and t.has(pi) and not any(t.has(g) for g in [l2,z5,Li4h,Li5h])]
c_pi2z3 = sp.simplify(sum(cz23_terms)/(pi**2*z3))
c_z2z3 = sp.simplify(6*c_pi2z3)  # z2z3 = pi^2 z3/6
print("edge1 folded z5 coeff (P_arc)  =", c_z5, "=", sp.N(c_z5,20))
print("edge1 folded z2z3 coeff(P_arc) =", c_z2z3, "=", sp.N(c_z2z3,20))
print("\nE5e: P_arc z5  = sqrt5*(-20187/2000) =", sp.N(sp.sqrt(5)*sp.Rational(-20187,2000),20))
print("E5e: P_arc z2z3= sqrt5*(-858/125)    =", sp.N(sp.sqrt(5)*sp.Rational(-858,125),20))
