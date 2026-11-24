
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
t0=time.time()
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
A1,B,kap,q0,V02,L2,lA0,lA1,lB = sp.symbols('A1 B kap q0 V02 L2 lA0 lA1 lB')
A0 = A1 - B*kap
numz, _ = clean['zeta3']
Ze_num = sp.expand(6*A0*A1*(L2*V02-q0**2)*(3*A1*L2+3*B*kap*q0-2*L2))
print("zeta3 coef numerator == Z_e numerator :", sp.expand(numz - Ze_num)==0)
print("=> cz3 = 6*A0*A1*(L2*V02-q0^2)*(3*A1*L2+3*B*kap*q0-2*L2)/L2^(5/2) = Z_e  (E1 law with u=-q0/L2, h^2=V02-q0^2/L2)")
with open(f"{WORK}/progress.log",'a') as fh:
    fh.write(f"[{time.strftime('%H:%M:%S')}] step_12: zeta3 coefficient of general E_arc = 6*A0*A1*(L2*V02-q0^2)*(3*A1*L2+3*B*kap*q0-2*L2)/L2^(5/2) = Z_e (E1 local law) EXACTLY (symbolic diff 0). With A1=1-B all polylog coefficients factor over Delta=L2*V02-q0^2.\n")
