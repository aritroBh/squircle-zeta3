
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
t0=time.time()
mons = pickle.load(open(f"{WORK}/Earc_basis_raw.pkl","rb"))
A0,A1,B,kap,Lp2,q0,V02 = sp.symbols('A0 A1 B kap L2 q0 V02')
# apply edge constraint A0 = A1 - B*kap, then factor coefficients of kap powers
simp = {}
for name, c in mons.items():
    c2 = sp.expand(c.subs({A0: A1 - B*kap}))
    # collect by powers of kap (including negative)
    c2 = sp.expand(c2)
    p = sp.Poly(c2, kap) if not c2.has(1/kap) else None
    simp[name] = c2
    print(f"{name:>16}: after A0-constraint #terms {len(sp.Add.make_args(c2))}  ({time.time()-t0:.0f}s)")
pickle.dump(simp, open(f"{WORK}/Earc_basis_A0sub.pkl","wb"))
print("saved", f"({time.time()-t0:.0f}s)")
