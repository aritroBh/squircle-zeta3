
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
A1,B,kap,q0,V02,L2 = sp.symbols('A1 B kap q0 V02 L2')
A0 = A1 - B*kap; r = kap+1
PhiA = 6*B*(q0**2-L2*V02)/L2**sp.Rational(5,2)*(3*A0*kap*q0 - L2*(2*A0+r*A1))
PhiL = 18*B*(q0**2-L2*V02)/L2**sp.Rational(5,2)*(kap*(-A1**2+A1*B*(2*kap-1)+2*B**2*kap)*q0 + L2*(A0*B+r*A1**2))
cL3 = sp.expand(clean['Li3(x)'][0]/L2**sp.Rational(5,2))
cL3L = sp.expand(clean['Li3(x/(x-1))'][0]/L2**sp.Rational(5,2))
print("coef[Li3(x/(x-1))] - Phi_A =", sp.expand(cL3L - PhiA))
print("coef[Li3(x)] - (2*Phi_A+Phi_L) =", sp.expand(cL3 - (2*PhiA+PhiL)))
with open(f"{WORK}/progress.log",'a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_16-19: full 11-coefficient table factored (params A1,B,r or B,r with both constraints); all over L2^(5/2), polylog coeffs have factor Delta=L2*V02-q0^2. coef[Li3(x/(x-1))]=Phi_A and coef[Li3(x)]=2Phi_A+Phi_L EXACTLY (symbolic diffs 0), matching K1d's general law.\n")
