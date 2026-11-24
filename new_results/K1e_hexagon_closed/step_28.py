
import pickle, sympy as sp, time
WORK='research/08_theory/K1e_hexagon_closed'
t0=time.time()
fact = pickle.load(open(f"{WORK}/Earc_basis_fact.pkl","rb"))
A1,B,kap,q0,V02,Ls,L2,lA0,lA1,lB = sp.symbols('A1 B kap q0 V02 Ls L2 lA0 lA1 lB')
clean = {}
for name,(num,den) in fact.items():
    e = sp.expand(num/den).subs({sp.sqrt(Ls**2): Ls, (Ls**2)**sp.Rational(5,2): Ls**5, (Ls**2)**sp.Rational(3,2): Ls**3})
    e = sp.together(e)
    n,d = sp.fraction(e)
    n = sp.expand(n)
    # substitute Ls^2 -> L2 (only even powers expected)
    assert sp.expand(n.subs({Ls:-Ls})-n)==0, (name, "odd Ls!")
    n = sp.expand(n.subs({Ls: sp.sqrt(L2)}))
    d = sp.expand(d.subs({sp.sqrt(Ls**2): Ls, Ls**5: L2**sp.Rational(5,2), Ls**4*sp.sqrt(Ls**2): L2**sp.Rational(5,2)}))
    clean[name] = (sp.factor(n), d)
    print(f"=== {name}: den = {d}")
for name in clean:
    print(f"--- {name} ---")
    print(clean[name][0])
pickle.dump(clean, open(f"{WORK}/Earc_basis_clean.pkl","wb"))
print(f"({time.time()-t0:.0f}s)")
