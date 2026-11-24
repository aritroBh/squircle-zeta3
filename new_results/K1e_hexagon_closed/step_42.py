
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
A1,B,kap,q0,V02,L2,lA0,lA1,lB,rr = sp.symbols('A1 B kap q0 V02 L2 lA0 lA1 lB r')
final={}
for name in ['zeta3','Li3(x)','Li3(x/(x-1))','Li2(x)*ln(1-x)','zeta2*ln(1-x)','zeta2','Li2(x)','ln^3(1-x)']:
    num,den = clean[name]
    n2 = sp.factor(sp.expand(num.subs({kap: rr-1, A1: 1-B})))
    final[name]=n2
    print(f"=== {name}  (A1=1-B, kap=r-1; over L2^(5/2)) ===")
    print(n2); print()
pickle.dump(final, open(f"{WORK}/Earc_basis_final.pkl","wb"))
