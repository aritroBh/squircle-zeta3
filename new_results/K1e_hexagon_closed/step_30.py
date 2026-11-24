
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
t0=time.time()
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
A1,B,kap,q0,V02,L2,lA0,lA1,lB = sp.symbols('A1 B kap q0 V02 L2 lA0 lA1 lB')
A0 = A1 - B*kap
# 1) zeta3 coefficient == Z_e (E1 local law) ?
numz, _ = clean['zeta3']
numz = sp.expand(numz)
h2 = V02 - q0**2/L2
Ze = 6*A0*A1*h2/sp.sqrt(L2)*(3*(A1 + B*kap*q0/L2) - 2)*L2**sp.Rational(5,2)/L2**sp.Rational(5,2)
Ze_num = sp.expand(6*A0*A1*(L2*V02-q0**2)*(3*A1*L2+3*B*kap*q0-2*L2))
print("zeta3 coef * L2^(5/2) - Z_e numerator =", sp.expand(numz - (-Ze_num)))
# i.e. numz should equal -6*A0*A1*(L2*V02-q0**2)*(3*A1*L2+3Bkapq0-2L2) = 6*A0*A1*(q0^2-L2V02)*(...)
# 2) apply A1 = 1-B and refactor the small ones
for name in ['zeta3','Li3(x)','Li3(x/(x-1))','Li2(x)*ln(1-x)','zeta2*ln(1-x)','zeta2','Li2(x)','ln^3(1-x)']:
    num, den = clean[name]
    n2 = sp.factor(sp.expand(num.subs({A1: 1-B})))
    print(f"--- {name} (A1=1-B) ---"); print(n2)
print(f"({time.time()-t0:.0f}s)")
