
import pickle, sympy as sp
WORK='research/08_theory/K1e_hexagon_closed'
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
A1,B,kap,q0,V02,L2 = sp.symbols('A1 B kap q0 V02 L2')
A0 = A1 - B*kap; r = kap+1
PhiA = 6*B*(q0**2-L2*V02)/L2**sp.Rational(5,2)*(3*A0*kap*q0 - L2*(2*A0+r*A1))
PhiL = 18*B*(q0**2-L2*V02)/L2**sp.Rational(5,2)*(kap*(-A1**2+A1*B*(2*kap-1)+2*B**2*kap)*q0 + L2*(A0*B+r*A1**2))
cL3 = sp.expand(clean['Li3(x)'][0]/L2**sp.Rational(5,2))
cL3L = sp.expand(clean['Li3(x/(x-1))'][0]/L2**sp.Rational(5,2))
d1 = sp.simplify(sp.expand(cL3L - PhiA).subs({A1:1-B}))
d2 = sp.simplify(sp.expand(cL3 - (2*PhiA+PhiL)).subs({A1:1-B}))
print("after A1=1-B: coef[Li3(x/(x-1))]-Phi_A =", d1)
print("after A1=1-B: coef[Li3(x)]-(2Phi_A+Phi_L) =", d2)
# also show factorization of the difference as multiple of (A1+B-1)
diff = sp.expand(cL3L - PhiA)
num,den = sp.fraction(sp.together(diff))
print("difference numerator has factor (A1+B-1):", sp.rem(sp.Poly(sp.expand(num), A1), sp.Poly(A1+B-1, A1)).as_expr()==0)
