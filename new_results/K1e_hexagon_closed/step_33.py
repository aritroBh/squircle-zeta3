
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
t0=time.time()
clean = pickle.load(open(f"{WORK}/Earc_basis_clean.pkl","rb"))
A1,B,kap,q0,V02,L2,lA0,lA1,lB,rr = sp.symbols('A1 B kap q0 V02 L2 lA0 lA1 lB r')
for name in ['ln^3(1-x)']:
    num, den = clean[name]
    n2 = sp.expand(num.subs({kap: rr-1, A1: 1-B}))
    n2f = sp.factor(n2)
    print(f"--- {name} (A1=1-B, kap=r-1), {len(sp.Add.make_args(n2))} terms ---")
    print(n2f)
print()
for name in ['ln^2(1-x)','zeta2']:
    num, den = clean[name]
    n2 = sp.expand(num.subs({kap: rr-1, A1: 1-B}))
    print(f"--- {name}: {len(sp.Add.make_args(n2))} terms")
print(f"({time.time()-t0:.0f}s)")
