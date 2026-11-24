
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
org = pickle.load(open(f"{WORK}/Earc_by_logmon.pkl","rb"))
A1,B,q0,V02,L2,rr = sp.symbols('A1 B q0 V02 L2 r')
for exps in [(3,0,0),(0,3,0),(0,0,3),(2,1,0),(1,2,0)]:
    c = org[exps]
    c2 = sp.factor(c.subs({A1:1-B}))
    print(f"lA0^{exps[0]} lA1^{exps[1]} lB^{exps[2]}: {c2}")
    print()
