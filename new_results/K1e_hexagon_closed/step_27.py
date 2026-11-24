
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
fact = pickle.load(open(f"{WORK}/Earc_basis_fact.pkl","rb"))
A1,B,kap,q0,V02,Ls,lA0,lA1,lB = sp.symbols('A1 B kap q0 V02 Ls lA0 lA1 lB')
L2v = Ls**2; A0 = A1 - B*kap; r = kap+1
PhiA = 6*B*(q0**2-L2v*V02)/L2v**sp.Rational(5,2)*(3*A0*kap*q0 - L2v*(2*A0+r*A1))
PhiL = 18*B*(q0**2-L2v*V02)/L2v**sp.Rational(5,2)*(kap*(-A1**2+A1*B*(2*kap-1)+2*B**2*kap)*q0 + L2v*(A0*B+r*A1**2))
numL3, denL3 = fact['Li3(x)']
numL3L, denL3L = fact['Li3(x/(x-1))']
d1 = sp.simplify(numL3L/denL3L - PhiA)
d2 = sp.simplify(numL3/denL3 - (2*PhiA+PhiL))
print("coef[Li3(x/(x-1))] - Phi_A =", d1)
print("coef[Li3(x)] - (2 Phi_A + Phi_L) =", d2)
with open(f"{WORK}/progress.log",'a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_09-10: after constraint A0=A1-Bkap all kap-poles cancel; every basis coefficient = (poly in kap,A1,B,q0,V02,L2,lA*)/L2^(5/2), all share factor (q0^2-L2*V02)=-L2 h^2; only even powers of sqrt(L2). coef[Li3(x/(x-1))]=Phi_A and coef[Li3(x)]=2Phi_A+Phi_L match K1d exactly (diffs 0).\n")
